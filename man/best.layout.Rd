\name{best.layout}
\alias{best.layout}
\title{complementary tools for layout}
\description{
This function just find the best layout fit for a number of plots desired.
}
\usage{
best.layout(x)
}
\arguments{
\item{x}{A scalar value indicating the number of plots desired}
}
\details{
No major details
}
\value{
Returns the best layout
\describe{
\item{res}{the number of rows and columns giving the best fit}
}
}

\references{

Covarrubias-Pazaran G, Diaz-Garcia L, Schlautman B, Salazar W, Zalapa J. Fragman: An R package for fragment analysis. 2016. BMC Genetics 17(62):1-8.

Robert J. Henry. 2013. Molecular Markers in Plants. Wiley-Blackwell. ISBN 978-0-470-95951-0.

Ben Hui Liu. 1998. Statistical Genomics. CRC Press LLC. ISBN 0-8493-3166-8.
}

\examples{
best.layout(9)
}
