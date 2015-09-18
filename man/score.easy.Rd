\name{score.easy}
\alias{score.easy}
\title{Fragment analysis scoring}

\description{
This function uses information from the fsa files read from \code{\link{storing.inds}} function and does the ssr calling in the channel specified and returns the index position, height and base pair position.

Please! if using the confidence interval method ("ci"), which is NOT the default, once you have found the best parameters for the arguments to match your ladder using the \code{\link{detect.ladder}} function, please pass those values to this function since it will use the same function internally and your dna sizes will depend on that, so if that's the case make sure the 'dev' argument is passed to the new functions.
}
\usage{
score.easy(my.inds, cols = 1, n.inds = NULL, panel=NULL, 
          thresh=NULL, shift=0.8, ladder, channel.ladder=NULL, 
          ploidy=2, ci.upp=1.96, ci.low=1.96, dev=50, left.cond=c(0.6,3), 
          right.cond=0.35, warn=FALSE, window=0.5, init.thresh=200, 
          ladd.init.thresh=200, method="cor", env = parent.frame(), 
          plotting=TRUE, electro=TRUE)
}

\arguments{
  \item{my.inds}{List with the channels information from the individuals specified, usually coming from the \code{\link{storing.inds}} function output}

\item{cols}{The channel you wish to analyze, usually 1 is blue, 2 is green, 3 is yellow, 4 is red and so on}

\item{n.inds}{Vector specifying the plants to be scored}

\item{panel}{A vector containing the base pair interval where the peaks should be searched for}

\item{thresh}{These are the values(number of standard errors) to be used to calculate the heights confidendence interval for ssr calling and be able to customize the thresholds when using the \code{\link{threshs}} function}

\item{shift}{The number of base pairs to be used for discarding neighboring peaks to the tallest peaks, i.e. if 2 peaks are 0.3 bp together the smalles will be discarded}

\item{ladder}{A vector containing the expected weights for the ladder peaks that will be found the using the \code{\link{find.ladder}} function}

\item{channel.ladder}{A scalar value indicating in which channel or color the ladder was read}

\item{ploidy}{A scalar value indicating the ploidy of the organism to be scored}

\item{ci.upp}{A scalar value indicating how many standar errors will be used to detect peaks when checking the height of the ladder peaks(upper bound). To be used in the \code{\link{find.ladder}} function}

\item{ci.low}{A scalar value indicating how many standar errors will be used to detect peaks when checking the height of the ladder peaks(lower bound). To be used in the \code{\link{find.ladder}} function}

\item{dev}{A scalar value indicating the number of indexes to be used as peak separation when deciding the ladder peaks, for more details check \code{\link{find.ladder}} function}

\item{left.cond}{A percentage value indicating when peaks to the leaft of the tallest peaks should be considered real based on the height, i.e. a very close peak right before the tallest peak if smaller than the tallest (half the size of the tallest one will be real or not)}

\item{right.cond}{A percentage value indicating when peaks to the right of the tallest peaks should be considered real based on the height, i.e. a very close peak right after the tallest peak if smaller than the tallest (half the size of the tallest one will be real or not)}

\item{warn}{A TRUE/FALSE value indicating if warnings should be provided when detecting the ladder}

\item{window}{A value in base pairs indicating how much is the error for detecting a peak in a sample that was provided in the panelas a real peak}

\item{init.thresh}{ An initial value of intensity to detect peaks. We recommend not to deal to much with unless you have highly controlled dna concentrations in your experiment}

\item{ladd.init.thresh}{ A value of intensity to detect peaks in the internal use of the \code{\link{find.ladder}} function. We recommend not to deal to much with it unless you identified special situations with your ladder}

\item{method}{An argument indicating one of the 2 methods available; "cor" makes all possible combination of peaks and searches exhaustive correlations to find the right peaks corresponsding to the expected DNA weights, or "ci" constructing confidence intervals to look for peaks meeting the conditions specified in the previous arguments. Passed to the find.ladder function working internally.}

\item{env}{ this is used to detect the environment of the user and load the result in the same environment.}

\item{plotting}{ a TRUE/FALSE value indicating if the plots should be drawn or not. The default value is TRUE.}

\item{electro}{A TRUE/FALSE value indicating if the electrogram/gel should be drawn or not. The default value is TRUE.}
}
\details{
 Remember!! once you have found the best parameters for the arguments to match your ladder using the \code{\link{detect.ladder}} function, please pass those values to this function since it will use the same function internally and your dna sizes will depend on that, please make sure the 'dev' argument is passed to the new functions since could be very specific for each ladder.
}
\value{
If arguments are correct the function returns a plot and a list containing
\describe{
\item{$pos}{the index positions for the intensities}
\item{$hei}{the intensities for the fragments found}
\item{$wei}{the putative weights in base pairs based on the ladder provided}
}
}

\references{

We have spent valuable time developing this package, please cite it in your publication:

Covarrubias-Pazaran G, Diaz-Garcia L, Schlautman B, Salazar W, Zalapa J. Fragma: An R package for fragment analysis. http://horticulture.wisc.edu/cggl/ZalapaLab/People.html. 2015.

Robert J. Henry. 2013. Molecular Markers in Plants. Wiley-Blackwell. ISBN 978-0-470-95951-0.

Ben Hui Liu. 1998. Statistical Genomics. CRC Press LLC. ISBN 0-8493-3166-8.
}

\examples{
data(my.plants)
my.plants <- my.plants[1]
my.ladder <- c(120, 125, 129, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375)
my.panel <- c(190:230)
# make sure to have sufficient plot area or will get an error message
score.easy (my.inds=my.plants, cols = 1, panel=my.panel, ladder=my.ladder)


}
