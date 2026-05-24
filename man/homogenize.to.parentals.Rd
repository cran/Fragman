\name{homogenize.to.parentals}
\alias{homogenize.to.parentals}
\title{complementary tools}

\description{
This functions takes a list of positions, heights and weights for ssr calls of a certain plant and uses parental information and a window to homogenize to the parent calls.
}
\usage{
homogenize.to.parentals(x, parents, window)
}

\arguments{
  \item{x}{list of positions, heights and weights}

\item{parents}{different parental calls}

\item{window}{window in base pairs indicating when should be considered the same allele}
}
\details{

No major details.
}
\value{
If arguments are correct the function returns a list containing
\describe{
\item{$pos}{the index positions for the intensities}
\item{$hei}{the intensities for the fragments found}
\item{$wei}{the putative weights in base pairs based on the ladder provided}
}
}

\references{
Robert J. Henry. 2013. Molecular Markers in Plants. Wiley-Blackwell. ISBN 978-0-470-95951-0.

Ben Hui Liu. 1998. Statistical Genomics. CRC Press LLC. ISBN 0-8493-3166-8.
}

\examples{
#No example provided, internally working to round to the closest parent ssr call.
x <- 1:10

}
