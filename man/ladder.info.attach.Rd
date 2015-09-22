\name{ladder.info.attach}
\alias{ladder.info.attach}
\title{Ladder detection and attachment to R environment}

\description{
This function uses the information stored by the \code{\link{storing.inds}} function and a vector specifying the ladder/standard and finds the real peaks corresponding to the expected weights. The user may use this function to be able to load the ladder information in the global environment of R, so when using the \code{\link{overview}} or \code{\link{score.easy}} functions calculations will be performed faster, if the function is not used the program will calculate the ladder information each time \code{\link{overview}} or \code{\link{score.easy}} functions are used. 

Please! if using the confidence interval method ("ci"), which is NOT the default, once you have found the best parameters for the arguments to match your ladder using the \code{\link{detect.ladder}} function, please pass those values to all the posterior functions, such as the 'dev' argument.
}
\usage{
ladder.info.attach(stored, ladder, channel.ladder=NULL, 
                  ci.upp=1.96, ci.low=1.96, dev=50, warn= FALSE, method="red", 
                  ladd.init.thresh=200, env = parent.frame(), prog=TRUE, draw=TRUE)
}

\arguments{
 \item{stored}{List with the channels information from the individuals specified, usually coming from the \code{\link{storing.inds}} function output.}

\item{ladder}{Vector containing the expected weights of the dna fragments of the ladder in use.}

\item{channel.ladder}{A scalar value indicating in which channel or color the ladder was read}

\item{ci.upp}{A scalar value indicating how many standar errors will be used to detect peaks when checking the height of the ladder peaks(upper bound). To be used in the \code{\link{find.ladder}} function.}

\item{ci.low}{A scalar value indicating how many standar errors will be used to detect peaks when checking the height of the ladder peaks(lower bound). To be used in the 'find.ladder' function}

\item{dev}{A scalar value indicating the number of indexes to be used as peak separation when deciding the ladder peaks. Some ladders contain dna fragments of very closed weights and modifying this parameter helps to detect them correctly}

\item{warn}{A TRUE/FALSE value indicating if warnings should be provided when detecting the ladder}

\item{method}{An argument indicating one of the 2 methods available; "cor" makes all possible combination of peaks and searches exhaustive correlations to find the right peaks corresponsding to the expected DNA weights, or "ci" constructing confidence intervals to look for peaks meeting the conditions specified in the previous arguments}

\item{ladd.init.thresh}{ A value of intensity to detect peaks in the internal use of the \code{\link{find.ladder}} function. We recommend not to deal to much with it unless you identified special situations with your ladder}

\item{env}{ this is used to detect the environment of the user and load the result in the same environment.}

\item{prog}{ A TRUE/FALSE value indicating if a progress bar should be drawn while processing the samples in order to assess the time it takes to find the ladder. The dafault value is TRUE but usually this makes the process slower. Please feel free to set it equal FALSE if the number of samples is quite large and speed is a concern.}

\item{draw}{ A TRUE/FALSE value indicating if a plot showing the peaks matched with your ladder should be drawn. The dafault value is FALSE to avoid extra delay. Please feel free to set it equal TRUE if you prefer to assess the sizing process.}
}
\details{
The peaks are detected by default using a correlation method bu the user can use confidence intervals if desired.

}
\value{
If parameters are indicated correctly the function returns:
\describe{
\item{$pos}{the index positions for the intensities}
\item{$hei}{the intensities for the fragments found}
\item{$wei}{the putative weights in base pairs based on the ladder provided}
}
}

\references{

We have spent valuable time developing this package, please cite it in your publication:

Covarrubias-Pazaran G, Diaz-Garcia L, Schlautman B, Salazar W, Zalapa J. (2015) Fragma: An R package for fragment analysis. R package version 1.0. URL https://cran.r-project.org/web/packages/Fragman/.

Robert J. Henry. 2013. Molecular Markers in Plants. Wiley-Blackwell. ISBN 978-0-470-95951-0.

Ben Hui Liu. 1998. Statistical Genomics. CRC Press LLC. ISBN 0-8493-3166-8.
}

\examples{
data(my.plants)
my.plants <- my.plants[1:2]
my.ladder <- c(120, 125, 129, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375)
ladder.info.attach(stored=my.plants, ladder=my.ladder)
}
