\name{plot.fsa_stored}
\alias{plot.fsa_stored}
\title{plot form fsa files stored with storing.inds}
\description{\code{plot} method for class \code{"fsa_stored"}.
}
\usage{
\method{plot}{fsa_stored}(x, lay=c(2,1), channel=NULL, cex.legend=.5, ncol.legend=4,lims=NULL, color=NULL, \dots)
}
\arguments{
 \item{x}{an object of class \code{"fsa_stored"}}
 \item{lay}{layout for the numberof plots to visualize.}
 \item{channel}{if preferred a single channel can be specified.}
 \item{cex.legend}{value of the size of the legend.}
 \item{ncol.legend}{number of columns to divide the legend.}
 \item{lims}{equivalen to ylim argument in plot, xlim still can be used.}
 \item{color}{specific colors to use in case the user don't want the default palete.}
 \item{\dots}{Further arguments to be passed}
}
\value{
  vector of plot
 }
\author{Giovanny Covarrubias \email{covarrubiasp@wisc.edu}}
\seealso{\code{\link{Fragman}}}
\keyword{models}

