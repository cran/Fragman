\name{overview2}
\alias{overview2}
\title{Assesing several plants with an overview}

\description{
This function uses information from the FSA files read from \code{\link{storing.inds}} function and creates an overlapping plot to assess graphically the peaks of several plants in certain channel in order to create a panel for the next functions \code{\link{score.easy}}. The function contains several defaults in most of the arguments, please check arguments but in general.

Please! once you have found the best parameters for the arguments to match your ladder using the \code{\link{detect.ladder}} function, please pass those values to this function since it will use the same function internally and your dna sizes will depend on that, please make sure the 'dev' argument is passed to the new functions.
}
\usage{
overview2(my.inds, cols = 1, n.inds = NULL, xlim = NULL, 
        ylim = NULL, ladder, channel.ladder = NULL, ploidy = 2, 
        ci.upp = 1.96, ci.low = 1.96, dev = 50, method="cor", 
        init.thresh=200, ladd.init.thresh=200, lwd=.1, warn=TRUE, 
        min.panel=100, env = parent.frame()) 

}

\arguments{
  \item{my.inds}{List with the channels information from the individuals specified, usually coming from the \code{\link{storing.inds}} function output}

\item{cols}{The channel/color you wish to analyze, usually 1 is blue, 2 is green, 3 is yellow, 4 is red and so on}

\item{n.inds}{Vector specifying the plants to be scored}

\item{xlim}{A vector containing the base pair interval where the plot should be drawn}

\item{ylim}{A vector containing the intensity interval where the plot should be drawn}

\item{ladder}{A vector containing the expected weights for the ladder peaks that will be found the using the \code{\link{find.ladder}} function}

\item{channel.ladder}{A scalar value indicating in which channel or color the ladder was read}

\item{ploidy}{A scalar value indicating the ploidy of the organism to be scored}

\item{ci.upp}{A scalar value indicating how many standar errors will be used to detect peaks when checking the height of the ladder peaks(upper bound). To be used in the \code{\link{find.ladder}} function}

\item{ci.low}{A scalar value indicating how many standar errors will be used to detect peaks when checking the height of the ladder peaks(lower bound). To be used in the \code{\link{find.ladder}} function}

\item{dev}{A scalar value indicating the number of indexes to be used as peak separation when deciding the ladder peaks, for more details check \code{\link{find.ladder}} function}

\item{method}{An argument indicating one of the 2 methods available; "cor" makes all possible combination of peaks and searches exhaustive correlations to find the right peaks corresponsding to the expected DNA weights, or "ci" constructing confidence intervals to look for peaks meeting the conditions specified in the previous arguments}

\item{init.thresh}{ An initial value of intensity to detect peaks. We recommend not to deal to much with unless you have highly controlled dna concentrations in your experiment}

\item{ladd.init.thresh}{ A value of intensity to detect peaks in the internal use of the \code{\link{find.ladder}} function. We recommend not to deal to much with it unless you identified special situations with your ladder}

\item{lwd}{ The width of the line}

\item{warn}{A TRUE/FALSE value indicating if warnings should be provided when detecting the ladder}

\item{min.panel}{A scalar value indicating which peak values should be ignored when creating a panel. If 'xlim' values are specified the 'min.panel' value is ignored and instead the panel peaks provided by the program are based in the region where you are zooming in.}

\item{env}{ this is used to detect the environment of the user and load the result in the same environment.}
}
\details{
No major details.
}
\value{
If rarguments are correct the function returns a list containing
\describe{
\item{$plot}{ Returns a plot joining the channel for the plants specified for the color desired and the peaks found by the function using the parameters specified}
  \item{$nana}{ Returns a vector with the names of the plants specified in the function}
}
}

\references{

We have spent valuable time developing this package, please cite it in your publication:

Covarrubias-Pazaran G, Diaz-Garcia L, Schlautman B, Salazar W, Zalapa J. Fragman: An R package for fragment analysis. http://horticulture.wisc.edu/cggl/ZalapaLab/People.html. 2015.

Robert J. Henry. 2013. Molecular Markers in Plants. Wiley-Blackwell. ISBN 978-0-470-95951-0.

Ben Hui Liu. 1998. Statistical Genomics. CRC Press LLC. ISBN 0-8493-3166-8.
}

\examples{
data(my.plants)
my.plants <- my.plants[1]
my.ladder <- c(120, 125, 129, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375)
overview2(my.inds=my.plants, cols = 1, ladder=my.ladder, lwd=1)
# now use:
# my.panel <- locator(type="p", pch=20, col="red")$x
# to click over the peaks and get the sizes in base pairs
# when you are done make sure you press the "Esc" key, do not push the stop button

}
