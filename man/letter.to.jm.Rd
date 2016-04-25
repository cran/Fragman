\name{letter.to.jm}
\alias{letter.to.jm}
\title{
Letter to JoinMap code converter
}
\description{
This function converts a vector of ssr calls in letter format to joinmap code using information from mother and father provided in first and second row respectively
}
\usage{
letter.to.jm(x)
}

\arguments{
  \item{x}{A vector of ssr calls in letter format or snp types, mother and father of the population should be in 1st and 2nd position respectively}
}
\details{

If numeric data exists first needs to be converted to letter code in order to use this function.
}
\value{
If arguments are correct the function returns a list containing
\describe{
\item{$y}{ A vector with ssr calls in joinmap format}
}
}

\references{
Covarrubias-Pazaran G, Diaz-Garcia L, Schlautman B, Salazar W, Zalapa J. Fragman: An R package for fragment analysis. 2016. BMC Genetics 17(62):1-8.

Robert J. Henry. 2013. Molecular Markers in Plants. Wiley-Blackwell. ISBN 978-0-470-95951-0.

Ben Hui Liu. 1998. Statistical Genomics. CRC Press LLC. ISBN 0-8493-3166-8.
}

\examples{
xx <- data.frame(cbind(a=rep(150, 96), b=c(rep(100,48), rep(150,48))))
xx1 <- num.to.lett(xx)
letter.to.jm(unlist(xx1))
# try using apply to a dataframe
}
