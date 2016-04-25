\name{transfft}
\alias{transfft}
\title{Applying the fourier transformation to a data frame}

\description{
This function takes a vector and applies a fourier transformation in order to smooth the peaks usinf the fft function in the base package. Use only top 40 percent of the lowest frequencies. 
}
\usage{
transfft(sn, top=0.3)
}

\arguments{
  \item{sn}{ a numeric vector containing the DNA intensities for the capillary electrophoresis.}

  \item{top}{ percent of lowest frequencies that should be used for the fourier transformation.}
}
\details{
No major details.
}
\value{
If arguments are correctly specified the function returns:
\describe{
\item{$y}{ A new vector of DNA intensities smoothed to avoid extra noisy peaks.}
}
}

\references{
Covarrubias-Pazaran G, Diaz-Garcia L, Schlautman B, Salazar W, Zalapa J. Fragman: An R package for fragment analysis. 2016. BMC Genetics 17(62):1-8.

Robert J. Henry. 2013. Molecular Markers in Plants. Wiley-Blackwell. ISBN 978-0-470-95951-0.

Ben Hui Liu. 1998. Statistical Genomics. CRC Press LLC. ISBN 0-8493-3166-8.
}

\examples{
data(my.plants)
g1 <- transfft(my.plants[[1]][,4], top=0.8)
g2 <- transfft(my.plants[[1]][,4], top=0.4)
g3 <- transfft(my.plants[[1]][,4], top=0.1)
layout(matrix(1:3,3,1))
plot(g1, type="l")
lines(g2, col="red")
lines(g3, col="blue")
par1 <- c("top=0.8", "top=0.4", "top=0.1")
par2 <- c("black", "red", "blue")
par3 <- c(1,1,1)
legend("topright", legend=par1, col=par2, bty = "n", lty=par3, lwd=par3, cex=0.75)
}
