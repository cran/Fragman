\name{jm.conv}
\alias{jm.conv}
\title{
Scores to JoinMap converter
}
\description{
This function converts a data frame containing the scores from score.easy to a new data frame with joinmap calls. The parents need to be provided in the first and second row. Currently only works for CP type of crosses. From CP type to F2 and BC2 is straight forward.
}
\usage{
jm.conv(a)
}

\arguments{
  \item{a}{A data frame containing the scores from score.easy function extracted from get.scores function containing the parental calls in the first 2 rows.
  }
}
\details{

No major details.
}
\value{
If arguments are correct the function returns a new data frame
\describe{
\item{res}{ A new data frame with markers in joinmap format}
}
}

\references{
Covarrubias-Pazaran G, Diaz-Garcia L, Schlautman B, Salazar W, Zalapa J. Fragma: An R package for fragment analysis. http://horticulture.wisc.edu/cggl/ZalapaLab/People.html. 2015.

Robert J. Henry. 2013. Molecular Markers in Plants. Wiley-Blackwell. ISBN 978-0-470-95951-0.

Ben Hui Liu. 1998. Statistical Genomics. CRC Press LLC. ISBN 0-8493-3166-8.
}

\examples{
xx <- data.frame(cbind(a=rep(150, 96), b=c(rep(100,48), rep(150,48))))
jm.conv(xx)
# try using apply to a dataframe
}
