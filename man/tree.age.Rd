% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tree.age.R
\name{tree.age}
\alias{tree.age}
\title{Extracting the age of nodes and tips in a tree.}
\usage{
tree.age(tree, age, order = "past", fossil = TRUE, digits = 3)
}
\arguments{
\item{tree}{A \code{phylo} object.}

\item{age}{The age of the tree. If missing the age is set to be the tree height.}

\item{order}{Either "past" if the units express time since the present (e.g. million years ago), or "present" if the unit is expressed in time since the root.}

\item{fossil}{\code{logical}, whether to always consider the tree as containing at least one living taxa (\code{TRUE}) or allowing only fossil taxa (\code{FALSE} - default), see details.}

\item{digits}{A \code{numeric} value or \code{integer} for the precision of the output.}
}
\description{
Extracting the age of each node and tip in a tree give the height of the tree or some specified age.
}
\details{
When \code{fossil = TRUE}, if the \code{tree} contains a \code{tree$root.time} element (for tree's root age), and that \code{order} is set to \code{"past"}, the output ages are adjusted to be starting from the root.time. Else, if no \code{tree$root.time} exists or \code{fossil = FALSE}, tips and nodes age is relative from the tip furthest away from the root. \emph{THIS FUNCTION DOES NOT ESTIMATE TREE AGES}, it just extracts branch length information and converts it into time units. For basic dating functions in R, check \code{\link[ape]{chronos}}, \code{\link[ape]{chronopl}},  \code{\link[ape]{chronoMPL}} or use more specialised dating software (e.g. MrBayes, BEAST, RAxML, etc.).
}
\examples{
## A dated random phylogeny with a root 50 units of time old.
tree.age(rtree(10), age = 50)
## A random tree with the distance since the root.
tree.age(rtree(10), order = 'present')

}
\seealso{
\code{\link{slice.tree}}, \code{\link{chrono.subsets}}.
}
\author{
Thomas Guillerme
}
