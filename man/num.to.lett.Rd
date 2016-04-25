\name{num.to.lett}
\alias{num.to.lett}
\title{
Number to Letter code converter
}
\description{
This function converts dataframes with rounded calls (numeric format in 2 cells), to letter format based on the GBS pipeline developed by Elshire et al. (2011) which can be used as intermediate step to transform to joinmap and Onemap formats, it requires mother and father in first and second row
}
\usage{
num.to.lett(xx)
}

\arguments{
  \item{xx}{matrix with numbers, every 2 columns is a marker and each row is an individual, parents are located in the first 2 rows}
}
\details{
No major details.
}
\value{
If arguments are correct the function returns a list containing
\describe{
\item{xx2}{ matrix coded in letter format where each column is a marker and each row is an individual}
}
}

\references{
Covarrubias-Pazaran G, Diaz-Garcia L, Schlautman B, Salazar W, Zalapa J. Fragman: An R package for fragment analysis. 2016. BMC Genetics 17(62):1-8.

Robert J. Henry. 2013. Molecular Markers in Plants. Wiley-Blackwell. ISBN 978-0-470-95951-0.

Ben Hui Liu. 1998. Statistical Genomics. CRC Press LLC. ISBN 0-8493-3166-8.
}

\examples{
xx <- data.frame(cbind(a=rep(150, 96), b=c(rep(100,48), rep(150,48))))
num.to.lett(xx)
}
