\name{lapply_pb}
\alias{lapply_pb}
\title{complementary tools for Fragma}
\description{
This function is a wrapper of lapply function that allows the drawing of a progress bar to assess the speed of the process.
}
\usage{
lapply_pb (X, FUN, ...)
}
\arguments{
\item{X}{a vector (atomic or list) or an expression object. See see \code{\link{lapply}}.}
\item{FUN}{the function to be applied to each element of X: see \code{\link{lapply}}.}
\item{...}{passes another arguments to the typical lapply function.}

}
\details{
No major details
}
\value{
Performs lapply drawing a progress bar
\describe{
\item{res}{the same result than using lapply}
}
}

\references{
See see \code{\link{lapply}}
}

\examples{
l <- sapply(1:200, function(x) list(rnorm(1000)))
lapply_pb(l, mean)
}
