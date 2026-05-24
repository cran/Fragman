\name{pullup}
\alias{pullup}
\title{Applying pullup to channels/colors}

\description{
This function takes a matrix of DNA intensities and merge all the channels (columns) to identify overall peaks and then creates a window moving from peak to peak looking for the channel where this peak is real and adjust the intensities in the other channels. 
}
\usage{
pullup(mati, plotting=FALSE, channel=4)
}

\arguments{
  \item{mati}{ matrix of intensities where each column is a channel/color for a given sample.}

\item{plotting}{ a TRUE/FALSE value indicating if the results from adjusting the intensities should be drawn or not.}
  
\item{channel}{ a numeric value indicating which of the channles/color (column) allocates the ladder intensties.}
}
\details{
No major details.
}
\value{
If arguments are correctly specified the function returns:
\describe{
\item{mati}{ A new matrix of DNA intensities corrected for overlapping of wavelenth readings in different channels.}
}
}

\references{
Covarrubias-Pazaran G, Diaz-Garcia L, Schlautman B, Salazar W, Zalapa J. Fragman: An R package for fragment analysis. 2016. BMC Genetics 17(62):1-8.

Robert J. Henry. 2013. Molecular Markers in Plants. Wiley-Blackwell. ISBN 978-0-470-95951-0.

Ben Hui Liu. 1998. Statistical Genomics. CRC Press LLC. ISBN 0-8493-3166-8.
}

\examples{
data(my.plants)
layout(matrix(1:2,2,1))
# without pull up adjustment
plot(my.plants[[1]][,1], type="l", col="blue", xlim=c(2750,2850))
lines(my.plants[[1]][,2], col="green")
lines(my.plants[[1]][,3], col="gold")
## adjusted
yy <- pullup(my.plants[[1]])
plot(yy[,1], type="l", col="blue", xlim=c(2750,2850))
lines(yy[,2], col="green")
lines(yy[,3], col="gold")
# general view
yy1 <- pullup(my.plants[[1]], plotting=TRUE)
}
