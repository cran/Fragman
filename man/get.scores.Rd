\name{get.scores}
\alias{get.scores}
\title{complementary tools}

\description{
This function extracts the information from auto.score and score.easy functions and fits everything in a dataframe.
}
\usage{
get.scores(my.scores, mark = "mark")
}

\arguments{
  \item{my.scores}{List of individuals which contains at the same time a list of positions, heights and weights for the fragment analysis}

\item{mark}{A vector of names of the markers scored in the panel or the parents}
}
\details{
Nomajor details.
}
\value{
If arguments are correct the function returns a data frame containing
\describe{
\item{$da}{A dataframe containing the ssr calls}
}
}

\references{
Covarrubias-Pazaran G, Diaz-Garcia L, Schlautman B, Salazar W, Zalapa J. Fragman: An R package for fragment analysis. 2016. BMC Genetics 17(62):1-8.

Robert J. Henry. 2013. Molecular Markers in Plants. Wiley-Blackwell. ISBN 978-0-470-95951-0.

Ben Hui Liu. 1998. Statistical Genomics. CRC Press LLC. ISBN 0-8493-3166-8.
}

\examples{

## here "a" is the ouptput of `score.easy` function
par1 <- list(pos=c(3100, 3240), hei=c(22917,20563), wei=c(202,212))
par2 <- list(pos=c(3100, 3240), hei=c(22917,20563), wei=c(202,214))
a <- list(i1=par1, i2=par2)
get.scores(a)
}
