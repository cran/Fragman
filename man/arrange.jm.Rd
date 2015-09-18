\name{arrange.jm}
\alias{arrange.jm}
\title{
Arrange data converted to joinmap code into a joinmap readable file 
}
\description{
This function converts a data frame containing the joinmap code into the readable file for joinmap. This format still needs some extra information, specifically the header indicating the population type, no.loci and no. of individuals, please check file examples included in JoinMap software. 
}
\usage{
arrange.jm(x, par=FALSE) 
}

\arguments{
  \item{x}{A data frame containing the scores in jm coding.
  }
  
  \item{par}{A TRUE/FALSE value indicating if the data returned should include the parents or not, the dafault value is FALSE, indicating that the first 2 individuals will not be included.
  }
}
\details{

No major details.
}
\value{
If arguments are correct the function returns a new data frame
\describe{
\item{joinmap}{ A new data frame with markers in joinmap readable format}
}
}

\references{
Robert J. Henry. 2013. Molecular Markers in Plants. Wiley-Blackwell. ISBN 978-0-470-95951-0.

Ben Hui Liu. 1998. Statistical Genomics. CRC Press LLC. ISBN 0-8493-3166-8.
}

\examples{
xx <- data.frame(cbind(a=rep(150, 96), b=c(rep(100,48), rep(150,48)))); xx[1,] <- c(150,150)
xx2 <- cbind(jm.conv(xx), jm.conv(xx), jm.conv(xx))
xx3 <- arrange.jm(xx2, par=FALSE) 
xx3[,1:10]
}
