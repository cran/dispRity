# ' Create a curve set out of a list in the right form.
# '
# ' @param curve_set A list containing elements r, obs, sim_m and optionally
# '   the element theo. r must be a vector describing the radius vector. obs
# '   must be a vector containing the summary function values for the
# '   observed pattern. obs must have same length as r. sim_m must be a
# '   matrix containing summary function values for all the simulated
# '   patterns. Each column corresponds to a pattern. The number of rows must
# '   match the length of r. If included, theo corresponds to the theoretical
# '   summary function curve. If present, its length must match the length of
# '   r.
# ' @param ... Do not use. (For internal use only.)
# ' @return The given list adorned with the proper class name.
# ' @export
create_curve_set <- function(curve_set, ...) {
    check_curve_set_content(curve_set, ...)
    class(curve_set) <- 'curve_set'
    curve_set
}


# ' Check the content validity of a potential curve_set object.
# '
# ' @param curve_set A curve_set object to be checked.
# ' @param allow_Inf_values Logical, for internal use. Can be used to allow infinite or nonnumeric
# ' values in an \code{\link[spatstat]{envelope}} object at the first place, if those are cropped
# ' away (in \code{\link{crop_curves}}).
check_curve_set_content <- function(curve_set, allow_Inf_values = FALSE) {
    possible_names <- c('r', 'obs', 'sim_m', 'theo', 'is_residual')

    n <- length(curve_set)
    if (n < 1L) {
        stop('curve_set must have some elements.')
    }
    if (!is.list(curve_set)) {
        stop('curve_set must be a list.')
    }

    name_vec <- names(curve_set)
    if (length(name_vec) != n) {
        stop('Every element of curve_set must be named.')
    }
    if (!all(name_vec %in% possible_names)) {
        stop('curve_set should contain only elements with some of these ',
             ' names: ', paste(possible_names, collapse = ', '))
    }

    r <- curve_set[['r']]
    n_r <- length(r)
    if (n_r < 1L) {
        stop('curve_set[["r"]] must have at least one element.')
    }
    if (!is.vector(r)) {
        stop('curve_set[["r"]] must be a vector.')
    }
    if (!all(is.numeric(r)) || !all(is.finite(r))) {
        stop('curve_set[["r"]] must have only finite numeric values.')
    }

    obs <- curve_set[['obs']]
    if (length(obs) != n_r) {
        stop('curve_set[["obs"]] must have as many values as ',
             'curve_set[["r"]].')
    }
    if (!is.vector(obs)) {
        stop('curve_set[["obs"]] must be a vector.')
    }
    if (!(allow_Inf_values | ( all(is.numeric(obs)) && all(is.finite(obs)) ))) {
        stop('curve_set[["obs"]] must have only finite numeric values.')
    }

    sim_m <- curve_set[['sim_m']]
    dim_sim_m <- dim(sim_m)
    if (!is.matrix(sim_m)) {
        stop('curve_set[["sim_m"]] must be a matrix.')
    }
    if (dim_sim_m[1] != n_r) {
        stop('curve_set[["sim_m"]] must have as many rows as there are ',
             'elements in curve_set[["r"]].')
    }
    if (dim_sim_m[2] < 1L) {
        stop('curve_set[["sim_m"]] must have at least one column.')
    }
    if (!(allow_Inf_values | ( all(is.numeric(sim_m)) && all(is.finite(sim_m)) ))) {
        stop('curve_set[["sim_m"]] must have only finite numeric values.')
    }

    theo <- curve_set[['theo']]
    n_theo <- length(theo)
    if (n_theo > 0L) {
        if (n_theo != n_r) {
            stop('curve_set[["theo"]] must have as many values as ',
                 'curve_set[["r"]].')
        }
        if (!is.vector(theo)) {
            stop('curve_set[["theo"]] must be a vector.')
        }
        if (!(allow_Inf_values | ( all(is.numeric(theo)) && all(is.finite(theo)) ))) {
            stop('curve_set[["theo"]] must have only finite numeric ',
                 'values.')
        }
    }

    is_residual <- curve_set[['is_residual']]
    n_is_residual <- length(is_residual)
    if (n_is_residual > 0L && (n_is_residual != 1L ||
                               !is.logical(is_residual) ||
                               !is.finite(is_residual))) {
        stop('curve_set[["is_residual"]] must be either TRUE or FALSE.')
    }

    if (n_is_residual > 0L && is_residual && n_theo > 0L) {
        stop('A residual curve set must not contain a theoretical curve.')
    }

    curve_set
}

# ' The rank envelope test
# '
# ' The rank envelope test, p-value and global envelope
# '
# '
# ' The rank envelope test is a completely non-parametric test, which provides a p-value
# ' interval given by the most liberal and the most conservative p-value estimate and
# ' the 100(1-alpha)\% global envelope for the chosen test function T(r) on
# ' the chosen interval of distances.
# '
# ' Given a \code{curve_set} (or an \code{\link[spatstat]{envelope}}) object,
# ' the test is carried out as follows.
# '
# ' For each curve in the curve_set, both the data curve and the simulations,
# ' the global rank measure R is determined. If savedevs = TRUE, then the
# ' global rank values R_1, R_2, ..., R_(s+1) are returned in the component 'k',
# ' where k[1] is the value for the data.
# '
# ' Based on R_i, i=1, ..., s+1, the p-interval is calculated. This interval is
# ' by default plotted for the object returned by the rank_envelope function.
# ' Also a single p-value is calculated and returned in component 'p'. By default
# ' this p-value is the mid-rank p-value, but another option can be used by specifying
# ' \code{ties} argument which is passed to \code{\link{estimate_p_value}}. For
# ' options see \code{\link{estimate_p_value}}.
# '
# ' The 100(1-alpha)\% global envelope is given by the 'k_alpha'th lower and
# ' upper envelope. For details see Myllymaki et al. (2013).
# '
# ' The above holds for p-value calculation if \code{lexo == FALSE} and then the test
# ' corresponds to the rank envelope test by Myllymaki et. al (2013). If \code{lexo == TRUE},
# ' then all the pointwise ranks are used to rank the curves by rank count ordering (Myllymaki et al., 2015)
# ' and the single p-value in \code{p} is the p-value based on the rank count ordering.
# '
# ' The rank count ordering test allows in principle a lower number of simulations to be used,
# ' but then the test may no longer be usable as a graphical test.
# '
# ' @references
# ' Myllymaki, M., Mrkvicka, T., Seijo, H., Grabarnik, P. (2013). Global envelope tests for spatial point patterns. arXiv:1307.0239 [stat.ME]
# '
# ' Myllymaki, M., Mrkvicka, T., Grabarnik, P., Seijo, H. and Hahn, U. (2015). Global envelope tests for spatial point patterns. arXiv:1307.0239v4 [stat.ME]
# '
# ' @param curve_set A curve_set (see \code{\link{create_curve_set}}) or an \code{\link[spatstat]{envelope}}
# '  object. If an envelope object is given, it must contain the summary
# '  functions from the simulated patterns which can be achieved by setting
# '  savefuns = TRUE when calling envelope().
# ' @param alpha The significance level. The 100(1-alpha)\% global envelope will be calculated.
# ' @param savedevs Logical. Should the global rank values k_i, i=1,...,nsim+1 be returned? Default: FALSE.
# ' @param alternative A character string specifying the alternative hypothesis. Must be one of the following:
# '         "two.sided" (default), "less" or "greater".
# ' @param lexo Logical, whether or not to use rank count ordering for calculation of the p-value. See details.
# ' @param ties Ties method to be passed to \code{\link{estimate_p_value}}. Used to obtain
# ' a point estimate for the p-value. The default point estimate is the mid-rank p-value.
# '
# ' @return An "envelope_test" object containing the following fields:
# ' \itemize{
# '   \item r = Distances for which the test was made.
# '   \item method = The name of the envelope test.
# '   \item alternative = The alternative specified in the function call.
# '   \item p = A point estimate for the p-value (default is the mid-rank p-value).
# '   \item ties = As the argument \code{ties}.
# '   \item p_interval = The p-value interval [p_liberal, p_conservative].
# '   \item k_alpha = The value of k corresponding to the 100(1-alpha)\% global envelope.
# '   \item k = Global rank values (k[1] is the value for the data pattern). Returned only if savedevs = TRUE.
# '   \item central_curve = If the curve_set (or envelope object) contains a component 'theo',
# '         then this function is used as the central curve and returned in this component.
# '         Otherwise, the central_curve is the mean of the test functions T_i(r), i=2, ..., s+1.
# '         Used for visualization only.
# '   \item data_curve = The test function for the data.
# '   \item lower = The lower envelope.
# '   \item upper = The upper envelope.
# '   \item call = The call of the function.
# ' }
# ' @export
# ' @seealso \code{\link{random_labelling}}, \code{\link{plot.envelope_test}}
# ' @examples
# '
# ' ## Testing complete spatial randomness (CSR)
# ' #-------------------------------------------
# ' require(spatstat)
# ' pp <- unmark(spruces)
# ' # Generate nsim simulations under CSR, calculate L-function for the data and simulations
# ' env <- envelope(pp, fun="Lest", nsim=2499, savefuns=TRUE, correction="translate")
# ' # The rank envelope test
# ' res <- rank_envelope(env)
# ' # Plot the result.
# ' # - The central curve is now obtained from env[['theo']], which is the
# ' # value of the L-function under the null hypothesis (L(r) = r).
# ' plot(res)
# ' # or (requires R library ggplot2)
# ' plot(res, use_ggplot2=TRUE)
# '
# ' ## Advanced use:
# ' # Choose the interval of distances [r_min, r_max] (at the same time create a curve_set from 'env')
# ' curve_set <- crop_curves(env, r_min = 1, r_max = 7)
# ' # For better visualisation, take the L(r)-r function
# ' curve_set <- residual(curve_set, use_theo = TRUE)
# ' # Do the rank envelope test
# ' res <- rank_envelope(curve_set); plot(res, use_ggplot2=TRUE)
# '
# ' ## Random labeling test
# ' #----------------------
# ' # requires library 'marksummary'
# ' mpp <- spruces
# ' # 1) Perform simulations under the random labelling hypothesis and calculate
# ' # the test function T(r) for the data pattern (mpp) and each simulation.
# ' # The command below specifies that the test function is T(r) = \hat{L}_m(r),
# ' # which is an estimator of the mark-weighted L function, L_m(r),
# ' # with translational edge correction (default).
# ' # The random_labelling function returns the centred functions \hat{L}_m(r)-T_0(r),
# ' # where T_0(r) = \hat{L}(r) is the unmarked L function.
# ' curve_set <- random_labelling(mpp, mtf_name = 'm', nsim=2499, r_min=1.5, r_max=9.5)
# ' # 2) Do the rank envelope test
# ' res <- rank_envelope(curve_set)
# ' # 3) Plot the test result
# ' plot(res, use_ggplot2=TRUE, ylab=expression(italic(L[m](r)-L(r))))
# '
# ' # Make the test using instead the test function T(r) = \hat{L}_mm(r);
# ' # which is an estimator of the mark-weighted L function, L_mm(r),
# ' # with translational edge correction (default).
# ' curve_set <- random_labelling(mpp, mtf_name = 'mm', nsim=2499, r_min=1.5, r_max=9.5)
# ' res <- rank_envelope(curve_set)
# ' plot(res, use_ggplot2=TRUE, ylab=expression(italic(L[mm](r)-L(r))))
# '
# ' ## Goodness-of-fit test (typically conservative)
# ' #-----------------------------------------------
# ' pp <- unmark(spruces)
# ' # Minimum distance between points in the pattern
# ' min(nndist(pp))
# ' # Fit a model
# ' fittedmodel <- ppm(pp, interaction=Hardcore(hc=1)) # Hardcore process
# '
# ' \dontrun{
# ' # Simulating Gibbs process by 'envelope' is slow, because it uses the MCMC algorithm
# ' #env <- envelope(fittedmodel, fun="Jest", nsim=999, savefuns=TRUE,
# '                  correction="none", r=seq(0, 4, length=500))
# '
# ' # Using direct algorihm can be faster, because the perfect simulation is used here.
# ' simulations <- NULL
# ' for(j in 1:2499) {
# '    simulations[[j]] <- rHardcore(beta=exp(fittedmodel$coef[1]),
# '                                  R = fittedmodel$interaction$par$hc,
# '                                  W = pp$window);
# '    if(j%%10==0) cat(j, "...", sep="")
# ' }
# ' env <- envelope(pp, simulate=simulations, fun="Jest", nsim=length(simulations),
# '                 savefuns=TRUE, correction="none", r=seq(0, 4, length=500))
# ' curve_set <- crop_curves(env, r_min = 1, r_max = 3.5)
# ' res <- rank_envelope(curve_set); plot(res, use_ggplot2=TRUE)
# ' }
# '
rank_envelope <- function(curve_set, alpha=0.05, savedevs=FALSE,
                          alternative=c("two.sided", "less", "greater"),
                          lexo=FALSE, ties) {
    curve_set <- convert_envelope(curve_set)

    if(alpha < 0 | alpha > 1) stop("Unreasonable value of alpha.")
    if(!is.logical(savedevs)) cat("savedevs should be logical. Using the default FALSE.")
    alternative <- match.arg(alternative)

    # The type of the p-value
    if(missing(ties))
        ties <- p_value_ties_default()
    else if(lexo) cat("The argument ties ignored, because lexo = TRUE. \n")
    if(lexo) ties <- "lexical"

    # data_curve = the vector of test function values for data
    # sim_curves = matrix where each row contains test function values of a simulation under null hypothesis
    data_curve <- curve_set[['obs']]
    sim_curves <- t(curve_set[['sim_m']])

    Nsim <- dim(sim_curves)[1];
    nr <- length(curve_set$r)
    # Define the central curve T_0
    T_0 <- get_T_0(curve_set)

    data_and_sim_curves <- rbind(data_curve, sim_curves)
    loranks <- apply(data_and_sim_curves, MARGIN=2, FUN=rank, ties.method = "average")
    hiranks <- Nsim+2-loranks
    # k:
    switch(alternative,
           "two.sided" = {
               allranks <- pmin(loranks, hiranks)
           },
           "less" = {
               allranks <- loranks
           },
           "greater" = {
               allranks <- hiranks
           })

    distance <- apply(allranks, MARGIN=1, FUN=min)
    u <- -distance
    #-- p-interval
    p_low <- estimate_p_value(obs=u[1], sim_vec=u[-1], ties='liberal')
    p_upp <- estimate_p_value(obs=u[1], sim_vec=u[-1], ties='conservative')

    #-- p-value
    if(!lexo) {
        p <- estimate_p_value(obs=u[1], sim_vec=u[-1], ties=ties)
    }
    else { # rank the curves by lexical ordering
        # order ranks within each curve
        sortranks <- apply(allranks, 1, sort) # curves now represented as columns
        lexo_values <- do.call("order", split(sortranks, row(sortranks)))

        # find ties
        sorted <- sortranks[ ,lexo_values]
        dupp <- duplicated(split(sorted, col(sorted)))
        tied <- dupp | c(dupp[-1], FALSE)

        # replace ranks of tied values by mean ranks
        # (maybe a little bit awkward, but isntitcool)
        tie.rle <- rle(tied)
        tie.end <- cumsum(tie.rle$lengths)
        tie.start <- cumsum(c(1,tie.rle$lengths))
        tie.start <- tie.start[-length(tie.start)]
        rank.rle <- tie.rle
        rank.rle$values <- (tie.start + tie.end)/2
        tieranks <- inverse.rle(rank.rle)
        newranks <- 1:(Nsim+1)
        newranks[tied] <- tieranks[tied]

        distance_lexo <- newranks[order(lexo_values)]
        #-- calculate the p-value
        u_lexo <- -distance_lexo
        p <- estimate_p_value(obs=u_lexo[1], sim_vec=u_lexo[-1])
    }

    #-- calculate the 100(1-alpha)% global envelope
    distancesorted <- sort(distance, decreasing=TRUE)
    kalpha <- distancesorted[floor((1-alpha)*(Nsim+1))]
    LB <- array(0, nr);
    UB <- array(0, nr);

    for(i in 1:nr){
        Hod <- sort(data_and_sim_curves[,i])
        LB[i]<- Hod[kalpha];
        UB[i]<- Hod[Nsim+1-kalpha+1];
    }

    res <- list(r=curve_set[['r']], method="Rank envelope test", alternative = alternative,
                p=p, p_interval=c(p_low,p_upp), ties=ties,
                k_alpha=kalpha,
                central_curve=T_0, data_curve=data_curve, lower=LB, upper=UB,
                call=match.call())
    if(savedevs) res$k <- distance
    class(res) <- "envelope_test"
    res
}

# ' Convert an envelope object to a curve_set object.
# '
# ' If given an envelope object, convert it into a curve_set object. If given
# ' a curve_set object, check its correctness and give it back.
# '
# ' @param curve_set A curve_set or an \code{\link[spatstat]{envelope}}
# '   object. If an envelope object is given, it must contain the summary
# '   functions from the simulated patterns which can be achieved by setting
# '   savefuns = TRUE when calling envelope().
# ' @param ... Allows to pass arguments to \code{\link{check_curve_set_content}}
# ' and \code{\link{envelope_to_curve_set}} (to be passed further through
# ' \code{\link{create_curve_set}} to \code{\link{check_curve_set_content}}).
# ' @return If an \code{\link[spatstat]{envelope}} object was given, return a
# '   corresponding curve_set object. If a curve_set object was given, return
# '   it unharmed.
convert_envelope <- function(curve_set, ...) {
    if (inherits(curve_set, 'envelope')) {
        curve_set <- envelope_to_curve_set(curve_set, ...)
    } else if (!is(curve_set, 'curve_set')) {
        stop('curve_set must either have class "envelope" (from spatstat) ',
             'or class "curve_set".')
    }
    check_curve_set_content(curve_set, ...)
    curve_set
}

# ' Turn an \code{\link[spatstat]{envelope}} object into a curve_set object.
# '
# ' @param env An \code{\link[spatstat]{envelope}} object. The envelope()
# '   functions must have been called with savefuns = TRUE.
# ' @return A corresponding curve_set object.
# ' @param ... Do not use. (For internal use only.)
# ' @export
envelope_to_curve_set <- function(env, ...) {
    if (!inherits(env, 'envelope')) {
        stop('env is not an envelope object.')
    }

    r <- env[['r']]
    n_r <- length(r)
    if (n_r < 1L) {
        stop('env[["r"]] must exist.')
    }

    obs <- env[['obs']]
    if (length(obs) != n_r) {
        stop('env[["obs"]] and env[["r"]] must have the same length.')
    }

    simfuns_fv <- attr(env, 'simfuns', exact = TRUE)
    if (length(simfuns_fv) < 1L) {
        stop('env did not include the simulated curve set. envelope must ',
             'be run with savefuns = TRUE.')
    }
    simulation_df <- as.data.frame(simfuns_fv)
    simulation_df_names <- names(simulation_df)

    if (!('r' %in% simulation_df_names)) {
        stop('The env attribute simfuns did not include r.')
    }
    if (!identical(r, simulation_df[['r']])) {
        stop('env[["r"]] must be equal to ',
             'attributes(env)[["simfuns"]][["r"]].')
    }

    # Dimensions: r_idx, sim_idx.
    sim_m <- as.matrix(simulation_df[, !(simulation_df_names %in% 'r')])
    if (nrow(sim_m) != n_r) {
        stop('simulated curves must have the same length as r.')
    }

    theo <- env[['theo']]
    n_theo <- length(theo)
    # If theo exists it must be correct.
    if (n_theo > 0L && n_theo != n_r) {
        stop('env[["theo"]] and env[["r"]] must have the same length.')
    }

    res <- list(r = r, obs = obs, sim_m = sim_m)
    if (n_theo > 0L) {
        res[['theo']] <- theo
    }
    res[['is_residual']] <- FALSE

    res <- create_curve_set(res, ...)
    res
}

# ' The default ties method for the p-value
# '
# ' The default ties method for the p-value calculated by estimate_p_value
p_value_ties_default <- function() {
    'midrank'
}


# ' Define T_0 from a curve_set object
# '
# ' Define T_0, the expectation of the test function under H_0, from a curve_set object.
# '
# ' @inheritParams convert_envelope
# ' @export
get_T_0 <- function(curve_set) {
    curve_set <- convert_envelope(curve_set)

    if(with(curve_set, exists('is_residual'))) {
        if(!curve_set[['is_residual']]) {
            if(with(curve_set, exists('theo'))) {
                T_0 <- curve_set[['theo']]
            }
            else {
                sim_curves <- t(curve_set[['sim_m']])
                T_0 <- colMeans(sim_curves)
            }
        }
        else {
            T_0 <- rep(0, times=length(curve_set$r))
        }
    }
    else { # Assume curve_set does not contain residuals
        if(with(curve_set, exists('theo'))) {
            T_0 <- curve_set[['theo']]
        }
        else {
            sim_curves <- t(curve_set[['sim_m']])
            T_0 <- colMeans(sim_curves)
        }
    }
    T_0
}

# ' Estimate p-value.
# '
# ' @param x The first argument.
# ' @param ... Additional arguments.
# ' @export
# ' @seealso \code{\link{estimate_p_value.default}}
# estimate_p_value <- function (x, ...) UseMethod('estimate_p_value')


# FIXME: Do we need to consider NA values at some point in life? Our
# methods should not produce them but if we export this function, they
# should be handled.
# FIXME: How about both two-sided and one-sided?
# FIXME: north_note_2002 uses the term anti-conservative for the incorrect p-value.
# ' Estimate p-value.
# '
# ' Estimates the p-value of the given observation for the given set of Monte
# ' Carlo samples. User can choose which method is used to treat possible
# ' tied values.
# '
# ' @usage \method{estimate_p_value}{default}(obs, sim_vec, ties = 'midrank')
# '
# ' @param obs The data sample. A scalar real value. Must not be
# '   NULL.
# ' @param sim_vec The Monte Carlo samples. A vector of real values.
# '   Must not be NULL.
# ' @param ties The method to treat tied values with. If one or more of the
# '   elements of sim_vec are equal to obs, how should the rank of obs be
# '   determined? For 'conservative' the resulting p-value will be the
# '   highest possible. For 'liberal' the p-value will be the lowest
# '   possible. For 'random' the rank of the obs within the tied values is
# '   uniformly sampled so that the resulting p-value is at most the
# '   conservative option and at least the liberal option. For 'midrank'
# '   the mid-rank within the tied values is taken. 'midrank' is the default.
# ' @return The p-value estimate. A scalar real value between 0 and 1.
# '
# ' @references Hajek & Sidak & Sen. Theory of Rank Tests. 1999. ff. 130.
# ' @export
estimate_p_value <- function(obs, sim_vec, ties = 'midrank') {
    if (length(obs) != 1L || !is.finite(obs) || !is.numeric(obs)) {
        stop('obs must be a scalar finite real value.')
    }
    n_sim <- length(sim_vec)
    if (n_sim < 1L || !all(is.finite(sim_vec)) ||
        !all(is.numeric(sim_vec))) {
        stop('sim_vec must have at least one element and the ',
             'elements must be finite and real.')
    }
    possible_ties <- c('midrank', 'random', 'conservative', 'liberal')
    if (length(ties) != 1L || !(ties %in% possible_ties)) {
        stop('ties must be exactly one of the following: ',
             paste(possible_ties, collapse=', '))
    }

    # FIXME: Ugly redundancy.
    if (ties %in% 'conservative') {
        n_less <- sum(sim_vec < obs)
    } else if (ties %in% 'liberal') {
        n_less <- sum(sim_vec <= obs)
    } else if (ties %in% 'midrank') {
        n_smaller <- sum(sim_vec < obs)
        n_equal <- sum(sim_vec == obs)
        n_less <- n_smaller + n_equal / 2L
    } else if (ties %in% 'random') {
        n_smaller <- sum(sim_vec < obs)
        n_equal <- sum(sim_vec == obs)
        n_rand_smaller <- sample(seq(0L, n_equal, by = 1L), size = 1L)
        n_less <- n_smaller + n_rand_smaller
    } else {
        stop('ties argument not recognized: ', ties)
    }
    n_all <- n_sim + 1L
    p_estimate <- 1 - n_less / n_all
    p_estimate
}

# ' Print method for the class 'envelope_test'
# ' @usage \method{print}{envelope_test}(x, ...)
# '
# ' @param x an 'envelope_test' object
# ' @param ... Ignored.
# '
# ' @method print envelope_test
# ' @export
print.envelope_test <- function(x, ...) {
    with(x, cat(method, "\n",
                " p-value of the test: ", p, sep=""))
    with(x, if(exists('ties')) cat(" (ties method: ", ties, ")\n", sep="")
            else cat("\n"))
    with(x, if(exists('p_interval'))
            cat(" p-interval         : (", p_interval[1], ", ", p_interval[2],")\n", sep=""))
}
