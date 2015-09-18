\name{saturate}
\alias{saturate}
\title{Checking and correcting saturated peaks}

\description{
This function takes a vector of intensities and looks for peaks above 8000 RFUs and correct for possible splits at the top of the peaks by inverting the vally between splitted peaks and correcting the peak. 
}
\usage{
saturate(y)
}

\arguments{
  \item{y}{ a vector containing the DNA intensities for the capillary electrophoresis.}

}
\details{
No major details.
}
\value{
If arguments are correctly specified the function returns:
\describe{
\item{$y}{ A new vector of DNA intensities adjusted for saturated peaks over 8000 RFUs.}
}
}

\references{
Robert J. Henry. 2013. Molecular Markers in Plants. Wiley-Blackwell. ISBN 978-0-470-95951-0.

Ben Hui Liu. 1998. Statistical Genomics. CRC Press LLC. ISBN 0-8493-3166-8.
}

\examples{
data(my.plants)
y <- my.plants[[1]][,3]
layout(matrix(1:2,2,1))
plot(y, type="l", xlim=c(2750,2850))
y2 <- saturate(y=y)
plot(y2, type="l", xlim=c(2750,2850))
}
