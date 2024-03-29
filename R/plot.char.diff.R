#' @title Plots character differences
#'
#' @description Plots a character difference matrix from a discrete character matrix or its character differences density profile.
#'
#' @param x A discrete matrix or an already computed character difference matrix of class \code{char.diff}.
#' @param ... Any additional graphical arguments to be passed to \code{image}.
#' @param type Either \code{"matrix"} (or \code{"m"}) or \code{"density"} (or \code{"d"}) for respectively plotting the matrix of character differences or its character differences density profile.
#' @param col Two colors for forming the gradient if \code{type = "correlation"} or for the density lines colors if \code{type = "density"}.
#' @param legend A logical value stating whether to print the legend or not (default = \code{TRUE}).
#' @param legend.title A \code{character} string to be displayed as the title of the legend (default = \code{Difference}).
#' @param legend.pos The position of the legend. Can be two \code{numeric}. Default is \code{"topleft"}.
#' @param legend.round A \code{numeric} value for digits up legend values. Default is \code{0}.
#' @param axis A logical value stating whether to print the axis or not (default = \code{TRUE}).
#' @param xlim Two \code{numeric} values to determine the x axis limits. If missing (default), the limits are calculated automatically to fit the plot window.
#' @param ylim Two \code{numeric} values to determine the y axis limits. If missing (default), the limits are calculated automatically to fit the plot window.
#' @param xlab A \code{character} string for the the x axis. Can be missing.
#' @param ylab A \code{character} string for the the y axis. Can be missing.
#' @param main An overall title for the plot.
#' 
#' @examples
#' ## Comparing two characters
#' char.diff(list(c(0, 1, 0, 1), c(0, 1, 1, 1)))
#' 
#' ## Pairwise comparisons in a morphological matrix
#' morpho_matrix <- matrix(sample(c(0,1), 100, replace = TRUE), 10)
#' 
#' ## Plotting a matrix
#' plot.char.diff(morpho_matrix)
#' 
#' ## Plotting the density profile of a char.diff object
#' char.diff_matrix <- char.diff(morpho_matrix)
#' plot(char.diff_matrix, type = "density")
#' 
#' @seealso \code{\link{char.diff}}
#' 
#' @author Thomas Guillerme
# @export
#' 

plot.char.diff <- function(x, ..., type = "matrix", legend = TRUE, legend.title = "Difference", legend.pos = "topleft", legend.round = 0, axis = TRUE, xlim, ylim, xlab, ylab, col, main) {

    matrix <- x

    ## Saving the call
    match_call <- match.call()

    ## Checking the input type
    if(!is(matrix, "matrix")) {
        stop.call(match_call$matrix, " must be a matrix.")
    } else {
        ## If the input is just a matrix, calculate the characters differences
        if(!is(matrix, "char.diff")) {
            matrix <- char.diff(matrix)
        }
    }
    class(matrix) <- "matrix"

    ## type must be either "matrix", "m", "density", or "d"
    all_types <- c("matrix", "m", "density", "d")
    ## type must be a character string
    check.class(type, "character")
    type <- tolower(type)
    ## type must have only one element
    check.length(type, 1, paste(" argument must only one of the following:\n", paste(all_types, collapse = ", "), ".", sep = ""))
    check.method(type, all_types, "type argument")
    ## if type is a letter change it to the full word (lazy people...)
    type <- ifelse(type == "m", "matrix", type)
    type <- ifelse(type == "d", "density", type)

    ## legend
    legend_class <- check.class(legend, c("logical", "character"))
    # if(legend_class == "character") {
    #     check.length(legend.pos, 2, " must be a logical value or two character strings.")
    # }
    check.class(legend.title, "character")
    check.length(legend.title, 1, " must be a single character string.")
    legend_pos_class <- check.class(legend.pos, c("character", "numeric"))
    if(legend_pos_class == "character") {
        check.length(legend.pos, 1, " must be a single character string or a pair of coordinates.")
    } else {
        check.length(legend.pos, 2, " must be a single character string or a pair of coordinates.")
    }
    check.class(legend.round, "numeric")
    check.length(legend.round, 1, " must be a single numeric value.")

    ## axis
    check.class(axis, "logical")

    ## plot arguments
    if(type == "density") {
        ## Default options:
        if(missing(main)) {
            main <- "Character differences profile"
        } else {
            check.class(main, "character")
            check.length(main, 1, " must be a single character string.", errorif = FALSE)
        }
        if(is(legend, "logical") && legend == TRUE) {
            legend <- c("Combined", "Individual")
        }
        if(missing(col)) {
            col <- c("black", "grey")
        } else {
            check.class(col, "character")
            check.length(col, 1, " must have at least two colours.", errorif = TRUE)
        }
        if(missing(xlab)) {
            xlab <- "Character differences"
        } else {
            check.class(xlab, "character")
            check.length(xlab, 1, " must be a single character string.", errorif = FALSE) 
        }
        if(missing(ylab)) {
            ylab <- "Density"
        } else {
            check.class(ylab, "character")
            check.length(ylab, 1, " must be a single character string.", errorif = FALSE) 
        }

    } else {
        ## Default options
        if(missing(col)) {
            col <- c("orange", "blue")
        } else {
            check.class(col, "character")
            check.length(col, 1, " must have at least two colours.", errorif = TRUE)
        }
        if(missing(main)) {
            main <- "Character differences matrix"
        }
    }


    ## Plotting the matrix
    if(type == "matrix") {
        ## Remove the upper triangle
        matrix[upper.tri(matrix)] <- NA

        ## Setting the colours
        colfunc <- grDevices::colorRampPalette(col)
        colheat <- rev(colfunc(10))

        ## Plotting the heat map
        image(matrix, col = colheat, axes = FALSE, main = main, ...)
        if(axis) {
            axis(1, at = seq(from = 0, to = 1, length.out = ncol(matrix)), labels = FALSE, tick = TRUE)
            axis(2, at = seq(from = 0, to = 1, length.out = ncol(matrix)), labels = FALSE, tick = TRUE)
        }

        ## Adding the legend
        if(legend) {
            legend("topleft", legend = c(as.character(round(max(matrix, na.rm = TRUE), legend.round)), as.character(round(min(matrix, na.rm = TRUE), legend.round))), title = legend.title, col = col, pch = 19)
        }
    } else {
        ## Plotting the density profile
        do.plot.char.diff.density(matrix, main, legend, col, xlim, ylim, legend.pos, xlab, ylab)
    }
}