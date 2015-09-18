\name{big.peaks.col}
\alias{big.peaks.col}
\title{Peak search by first derivatives}

\description{
This function find all peaks by taking the first derivative based on 'rle' function. Is used in different Fragma functions.
}
\usage{
big.peaks.col(x,tre)
}

\arguments{
 \item{x}{A vector of heights or intensities}

\item{tre}{A scalar value deciding when peaks will be ignored}
}
\details{
No major details.
}
\value{
Retuns the biggest peaks for a vector of intensities.
\describe{
\item{out}{ a vector of positions where the derivative is zero and therefore a peak was found}
}
}

\references{
Robert J. Henry. 2013. Molecular Markers in Plants. Wiley-Blackwell. ISBN 978-0-470-95951-0.

Ben Hui Liu. 1998. Statistical Genomics. CRC Press LLC. ISBN 0-8493-3166-8.
}

\examples{
data(my.plants)
big.peaks.col(my.plants[[1]][,1],100)#for any color
}