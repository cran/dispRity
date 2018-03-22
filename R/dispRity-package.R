#' Measuring Disparity in R
#' 
#' A modular package for measuring disparity from multidimensional matrices. Disparity can be calculated from any matrix defining a multidimensional space. The package provides a set of implemented metrics to measure properties of the space and allows users to provide and test their own metrics. The package also provides functions for looking at disparity in a serial way (e.g. disparity through time) or per groups as well as visualising the results. Finally, this package provides several basic statistical tests for disparity analysis.
#' 
#' @name dispRity-package
#'
#' @docType package
#'
#' @author Thomas Guillerme <guillert@@tcd.ie>
#'
#' @keywords disparity, ordination, phylogeny, cladistic, morphometric, ecology
#'
#' @import ape
#' @import stats
#' @importFrom paleotree timeSliceTree 
#' @importFrom geometry convhulln 
#' @importFrom ade4 randtest as.randtest 
#' @importFrom grDevices colorRampPalette grey 
#' @importFrom caper comparative.data 
#' @importFrom graphics axis boxplot hist image lines mtext par plot points polygon text legend
#' @importFrom utils combn data capture.output
#' @importFrom phyclust gen.seq.HKY 
#' @importFrom phangorn dist.hamming NJ RF.dist CI RI optim.parsimony parsimony
#' @importFrom vegan adonis vegdist
#' @importFrom geiger dtt ratematrix sim.char
#' @importFrom parallel parLapply detectCores makeCluster clusterExport stopCluster
##' @importFrom Claddis MorphDistMatrixFast


NULL

#' Beck and Lee 2014 datasets
#'
#' Example datasets from Beck and Lee 2014.
#'
#' \itemize{
#'   \item \code{BeckLee_tree} A phylogenetic tree with 50 living and fossil taxa
#'   \item \code{BeckLee_mat50} The ordinated matrix based on the 50 taxa cladistic distances
#'   \item \code{BeckLee_mat99} The ordinated matrix based on the 50 taxa + 49 nodes cladistic distances
#'   \item \code{BeckLee_ages} A list of first and last occurrence data for fossil taxa
#' }
#'
#' @format three matrices and one phylogenetic tree.
#' @source \url{http://rspb.royalsocietypublishing.org/content/281/1793/20141278.short}
#' @references Beck RMD & Lee MSY. 2014. Ancient dates or accelerated rates?
#' Morphological clocks and the antiquity of placental mammals.
#' Proc. R. Soc. B 2014 281 20141278; DOI: 10.1098/rspb.2014.1278
#' @name BeckLee
#' @aliases BeckLee_tree BeckLee_mat50 BeckLee_mat99 BeckLee_ages
#' @seealso McClean_data disparity
NULL

# #' McClean dataset
# #'
# #' Example datasets from McClean (unpublished).
# #'
# #' \itemize{
# #'   \item \code{ordination} The ordinated matrix based on 40 sites and 20 dimensions
# #'   \item \code{treatment} A vector of \code{character} treatment types (\code{"a"} or \code{"b"})
# #'   \item \code{depth} A vector of \code{numeric} depth types (\code{1} or \code{2})
# #' }
# #'
# #' @format one matrix and two vectors.
# #' @name McClean_data
# #' @seealso BeckLee_data disparity
# NULL

#' disparity
#'
#' And example of a \code{dispRity} object.
#'
#' This matrix is based on the \code{\link{BeckLee}} dataset and split into seven continuous subsets (\code{\link{chrono.subsets}}).
#' It was bootstrapped 100 times (\code{\link{boot.matrix}}) with four rarefaction levels.
#' Disparity was calculated as the \code{\link[stats]{median}} of the \code{\link{centroids}} (\code{\link{dispRity}}).
#'
#' @format one \code{dispRity} object.
#' @name disparity
#' @seealso McClean_data BeckLee_data
NULL