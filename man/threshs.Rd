\name{threshs}
\alias{threshs}
\title{Customizing thresholds}

\description{
This function takes data contained in a list with 2 elements, the first containing the position in base pairs and the heights of the positions and given a panel finds the best minimum threshold by creating a confidence interval in order to get only the real peaks. Implemented internally of auto,score and score.easy function.
}
\usage{
threshs(my.plant, min.thre=200, panel, ci=1.9)
}

\arguments{
  \item{my.plant}{List with 2 elements; the position in base pairs and the heights associated with them}

\item{min.thre}{Starting threshold to find the peaks and refine the search}

\item{panel}{Two scalar values indicatin the interval where the peaks should be found}

\item{ci}{Scalar value indicating the number of standard errors to be used when creating the confidence interval for ssr calling}
}
\details{
No major details.
}
\value{
If arguments are correctly specified the function returns a list containing
\describe{
\item{$newthre}{ A scalar value indicating the new threshold to be used when calling the real peaks}
}
}

\references{
Robert J. Henry. 2013. Molecular Markers in Plants. Wiley-Blackwell. ISBN 978-0-470-95951-0.

Ben Hui Liu. 1998. Statistical Genomics. CRC Press LLC. ISBN 0-8493-3166-8.
}

\examples{
data(my.plants)
# implemented internally of auto,score and score.easy function
# threshs(my.plant, min.thre=200, panel=(260,290), ci=1.9)
}
