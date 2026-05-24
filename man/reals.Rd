\name{reals}
\alias{reals}
\title{Finding the real peaks}

\description{
This function takes a list with the information of positions, heights and weights for an individual and using the panel information finds the real peaks by using the separate function and getting the tallest peaks in the confidence interval constructed for the heights in the inteval of interest.
}
\usage{
reals(x, panel=c(100:400), shi=1, ploidy=2, left.cond=c(0.4,3), 
      right.cond=0.2, window=0.5)
}

\arguments{
  \item{x}{List with 3 elements; the information of positions, heights and weights for an individual in certain channel}

\item{panel}{A vector containing the base pair interval where the peaks should be searched for}

\item{shi}{The number of base pairs to be used for discarding neighboring peaks to the tallest peaks, i.e. if 2 peaks are 0.3 bp together the smalles will be discarded}

\item{ploidy}{A scalar value indicating the ploidy of the organism to be scored}

\item{left.cond}{A percentage value indicating when peaks to the leaft of the tallest peaks should be considered real based on the height, i.e. a very close peak right before the tallest peak if smaller than the tallest (half the size of the tallest one will be real or not)}

\item{right.cond}{A percentage value indicating when peaks to the right of the tallest peaks should be considered real based on the height, i.e. a very close peak right after the tallest peak if smaller than the tallest (half the size of the tallest one will be real or not)}

\item{window}{A value in base pairs indicating how much is the error for detecting a peak in a sample that was provided in the panelas a real peak}
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
data(my.plants)
x <- big.peaks.col(my.plants[[1]][,1],100)#for any color
#reals(x, panel=c(260,280), shi=1, ploidy=2) 
#still needs weight information in order to find the reals, 
#works internally of score.easy function
}
