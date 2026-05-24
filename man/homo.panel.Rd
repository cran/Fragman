\name{homo.panel}
\alias{homo.panel}
\title{complementary tools}

\description{
This functions takes a list of positions, heights and weights for ssr calls of a certain plant and uses panel information and a window to homogenize the weights of base pairs to the panel calls provided.
}
\usage{
homo.panel(x, panel, window)
}

\arguments{
  \item{x}{list of positions, heights and weights}

\item{panel}{different dna sizes usually obtained by using overview and locator functions}

\item{window}{window in base pairs indicating when should be considered the same panel allele}
}
\details{

No major details.
}
\value{
If arguments are correct the function returns a list containing
\describe{
\item{$pos}{the index positions for the intensities}
\item{$hei}{the intensities for the fragments found}
\item{$wei}{the putative weights in base pairs based on the panel provided}
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
