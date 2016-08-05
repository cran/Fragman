\name{ladder.info.attach}
\alias{ladder.info.attach}
\title{Ladder detection and attachment to R environment}

\description{
This function uses the information stored by the \code{\link{storing.inds}} function and a vector specifying the ladder/standard and finds the real peaks corresponding to the expected weights. The user may use this function to be able to load the ladder information in the global environment of R, so when using the \code{\link{overview}} or \code{\link{score.easy}} functions calculations will be performed faster, if the function is not used the program will calculate the ladder information each time \code{\link{overview}} or \code{\link{score.easy}} functions are used. 

NOTE: THE STEP OF MATCHING THE LADDER WITH YOUR SAMPLES USING THE `ladder.info.attach` FUNCTION IS CRITICAL. IF YOU HAVE ANY PROBLEM TRY MODIFYING THE ARGUMENT 'method', WITH THE 2 MOST EFFECTIVE METHODS method="iter" OR method="iter2"

}
\usage{
ladder.info.attach(stored, ladder, channel.ladder=NULL, 
                  method=NULL, ladd.init.thresh=NULL, 
                  env = parent.frame(), prog=TRUE, 
                  draw=TRUE, attempt=10)
}

\arguments{
 \item{stored}{List with the channels information from the individuals specified, usually coming from the \code{\link{storing.inds}} function output.}

\item{ladder}{Vector containing the expected weights of the dna fragments of the ladder in use.}

\item{channel.ladder}{A scalar value indicating in which channel or color the ladder was read}

\item{method}{An argument indicating one of the methods available; "iter" makes an iterative procedure to find the ladder, "iter2" the same but backwards, "cor" makes all possible combination of peaks and searches exhaustive correlations to find the right peaks corresponsding to the expected DNA weights, or "ci" constructing confidence intervals to look for peaks meeting the conditions specified in the previous arguments. The default allows the program to pick among "iter" and "iter2". Older methods have been depreciated.}

\item{ladd.init.thresh}{ A value of intensity to detect peaks in the internal use of the \code{\link{find.ladder}} function. We recommend not to deal to much with it unless you identified special situations with your ladder}

\item{env}{ this is used to detect the environment of the user and load the result in the same environment.}

\item{prog}{ A TRUE/FALSE value indicating if a progress bar should be drawn while processing the samples in order to assess the time it takes to find the ladder. The dafault value is TRUE but usually this makes the process slower. Please feel free to set it equal FALSE if the number of samples is quite large and speed is a concern.}

\item{draw}{ A TRUE/FALSE value indicating if a plot showing the peaks matched with your ladder should be drawn. The dafault value is FALSE to avoid extra delay. Please feel free to set it equal TRUE if you prefer to assess the sizing process.}

\item{attempt}{A scalar value indicating how many attempts should be made to find the real ladder peaks when using the "iter" and "iter2" methods. By default is 7 attempts, which means that will try to build the model assuming that the first peak found in the ladder is the corresponding first peak of the expected ladder, then moves to the 2nd peak until the 7th and the seven models are compared picking the most likely model based on the R2 value for each of the models.}

}
\details{

Method "ci" has been depreciated, currently the method "iter2" is the default and uses the ladder provided and observed peaks to match them using an iterative procedure based on least squares.

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

Covarrubias-Pazaran G, Diaz-Garcia L, Schlautman B, Salazar W, Zalapa J. Fragman: An R package for fragment analysis. 2016. BMC Genetics 17(62):1-8.

Robert J. Henry. 2013. Molecular Markers in Plants. Wiley-Blackwell. ISBN 978-0-470-95951-0.

Ben Hui Liu. 1998. Statistical Genomics. CRC Press LLC. ISBN 0-8493-3166-8.
}

\examples{
data(my.plants)
my.plants <- my.plants[1:2]
my.ladder <- c(50, 75, 100, 125, 129, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375)
ladder.info.attach(stored=my.plants, ladder=my.ladder)
}
