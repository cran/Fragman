\name{separate}
\alias{separate}
\title{Separating peaks by a shift window}

\description{
This function takes a list with positions, heights and weights called "g" and using a shift in base pairs determines when 2 neighboring peaks should be considered only one by getting the tallest peak. For example two peaks found at 173 and 173.5 base pairs are unlikely to be 2 different peaks, therefore only the tallest peak will pe chosen.
}
\usage{
separate(g, shift=1, type="bp")
}

\arguments{
  \item{g}{List with 3 elements; the information of positions, heights and weights for an individual in certain channel}

\item{shift}{The number of base pairs to be used for discarding neighboring peaks to the tallest peaks, i.e. if 2 peaks are 0.3 bp together the smalles will be discarded}

\item{type}{A word indicating if the shift to be used should be used in base pairs or in index. The use is "bp" or "ind"}
}
\details{

No major details.
}
\value{
If arguments are correct the function returns a new list containing
\describe{
\item{$pos}{the index positions for the intensities}
\item{$hei}{the intensities for the fragments found}
\item{$wei}{the putative weights in base pairs based on the ladder provided}
}
}

\references{
Covarrubias-Pazaran G, Diaz-Garcia L, Schlautman B, Salazar W, Zalapa J. Fragman: An R package for fragment analysis. 2016. BMC Genetics 17(62):1-8.

Robert J. Henry. 2013. Molecular Markers in Plants. Wiley-Blackwell. ISBN 978-0-470-95951-0.

Ben Hui Liu. 1998. Statistical Genomics. CRC Press LLC. ISBN 0-8493-3166-8.
}

\examples{
data(my.plants)
x <- big.peaks.col(my.plants[[1]][,1],100)#for any color
#separate(x, shift=1, type="bp") #still needs weight information
}
