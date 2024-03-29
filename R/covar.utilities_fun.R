## Internal for sauron.plot and covar.utilties

## Redimension
redimension <- function(one_post, dimensions) {
    if(!is.null(one_post$VCV)) {
        one_post$VCV <- one_post$VCV[dimensions, dimensions, drop = FALSE]
    }
    if(!is.null(one_post$Sol)) {
        one_post$Sol <- one_post$Sol[dimensions]
    }
    return(one_post)
}

## Selecting n elements from a covar object
sample.n <- function(covar, n, selected_n, dimensions = NULL) {
    ## Get the posterior size
    n_post <- length(covar[[1]])
    if(!missing(n) && n == n_post) {
        if(is.null(dimensions)) {
            ## Short cut
            return(covar)    
        } else {
            return(lapply(covar, lapply, redimension, dimensions = dimensions))
        }
    } else {
        ## Sample n
        if(missing(selected_n)) {
            selected_n <- sample.int(n_post, n, replace = n > n_post)
        }
        ## Return n for each group
        if(is.null(dimensions)) {
            return(lapply(covar, function(group, n) group[n], n = selected_n))
        } else {
            return(lapply(lapply(covar, function(group, n) group[n], n = selected_n), lapply, redimension, dimensions = dimensions))
        }
    }
}

## Internal: get the coordinates of one axes
get.one.axis <- function(data, axis = 1, level = 0.95, dimensions) {

    # The magic: https://stackoverflow.com/questions/40300217/obtain-vertices-of-the-ellipse-on-an-ellipse-covariance-plot-created-by-care/40316331#40316331

    ## VCVing the matrix
    if(!is(data, "list")) {
        data <- list(VCV = data)
    } else {
        if(is.null(data$VCV)) {
            data$VCV <- data
        }
    }
    ## adding a loc
    if(is.null(data$loc)) {
        data$loc <- rep(0, nrow(data$VCV))
    }

    ## Select the right dimensions
    data$VCV <- data$VCV[dimensions, dimensions, drop = FALSE]

    ## Get the data dimensionality
    dims <- length(diag(data$VCV))

    ## Create the unit hypersphere (a hypersphere of radius 1) for the scaling
    unit_hypersphere1 <- unit_hypersphere2 <- matrix(0, ncol = dims, nrow = dims)
    ## The "front" (e.g. "top", "right") units
    diag(unit_hypersphere1) <- 1
    ## The "back" (e.g. "bottom", "left") units
    diag(unit_hypersphere2) <- -1
    unit_hypersphere <- rbind(unit_hypersphere1, unit_hypersphere2)
    ## Scale the hypersphere (where level is the confidence interval)
    unit_hypersphere <- unit_hypersphere * sqrt(qchisq(level, dims))

    ## Do the eigen decomposition (symmetric - faster)
    eigen_decomp <- eigen(data$VCV, symmetric = TRUE)

    ## Re-scaling the unit hypersphere
    scaled_edges <- unit_hypersphere * rep(sqrt(abs(eigen_decomp$values)), each = dims*2)
    ## Rotating the edges coordinates
    edges <- tcrossprod(scaled_edges, eigen_decomp$vectors)

    ## Move the matrix around
    edges <- edges + rep(data$loc[dimensions, drop = FALSE], each = dims*2)

    ## Get the edges coordinates
    return(edges[c(axis, axis+dims), , drop = FALSE])
}

## Internal: get the coordinates of one axes (second try) from: https://stackoverflow.com/a/40304894/9281298
# get.one.axis2 <- function(data, axis = 1, level = 0.95, dimensions) {

#     # location along the ellipse
#     # linear algebra lifted from the code for ellipse()
#     ellipse.loc <- function(theta, center, shape, radius)
#     {
#         vert <- cbind(cos(theta), sin(theta))
#         Q <- chol(shape, pivot=TRUE)
#         ord <- order(attr(Q, "pivot"))
#         t(center + radius*t(vert %*% Q[, ord]))
#     }

#     # distance from this location on the ellipse to the center 
#     ellipse.rad <- function(theta, center, shape, radius)
#     {
#         loc <- ellipse.loc(theta, center, shape, radius)
#         (loc[,1] - center[1])^2 + (loc[,2] - center[2])^2
#     }

#     # ellipse parameters
#     center <- c(-0.05, 0.09)
#     A <- matrix(c(20.43, -8.59, -8.59, 24.03), nrow=2)
#     radius <- 1.44

#     # solve for the maximum distance in one hemisphere (hemi-ellipse?)
#     t1 <- optimize(ellipse.rad, c(0, pi - 1e-5), center=center, shape=A, radius=radius, maximum=TRUE)$m
#     l1 <- ellipse.loc(t1, center, A, radius)

#     # solve for the minimum distance
#     t2 <- optimize(ellipse.rad, c(0, pi - 1e-5), center=center, shape=A, radius=radius)$m
#     l2 <- ellipse.loc(t2, center, A, radius)

#     # other points obtained by symmetry
#     t3 <- pi + t1
#     l3 <- ellipse.loc(t3, center, A, radius)

#     t4 <- pi + t2
#     l4 <- ellipse.loc(t4, center, A, radius)
#     return(cbind())

#     ## Visualise
#     # plot(ellipse::ellipse(x = A, centre = center), type = "l")
#     # points(rbind(l1, l2, l3, l4), cex=2, col="blue", lwd=2)
#     # lines(rbind(l1, l3))
#     # lines(rbind(l2, l4))
# }



## Internal: summarising distributions
summarise.fun <- function(one_group, fun) {
    output <- list()
    ## Summarise all elements
    if("VCV" %in% names(one_group[[1]])) {
        if(all(dim(one_group[[1]]$VCV) == c(1,1))) {
            output$VCV <- matrix(apply(do.call(cbind, lapply(one_group, `[[`, "VCV")), 1, fun))
        } else {
            output$VCV <- apply(simplify2array(lapply(one_group, `[[`, "VCV")), 1:2, fun)
        }
    }
    if("Sol" %in% names(one_group[[1]])) {
        output$Sol <- apply(do.call(rbind, lapply(one_group, `[[`, "Sol")), 2, fun)
    }                
    return(output)
}