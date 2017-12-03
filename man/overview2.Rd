\name{overview2}
\alias{overview2}
\title{Assesing several plants with an overview}

\description{
This function uses information from the FSA files read from \code{\link{storing.inds}} function and creates an overlapping plot to assess graphically the peaks of several plants in certain channel in order to create a panel for the scoring functions \code{\link{score.easy}}. The function contains several defaults in most of the arguments, please check arguments but in general you only need the first 4 arguments to create a panel.

}
\usage{
overview2(my.inds, cols = 1, ladder, xlim = NULL, ylim = NULL, 
        n.inds = NULL, channel.ladder = NULL, ploidy = 2, 
        method="iter2", init.thresh=NULL, ladd.init.thresh=200, 
        lwd=.25, warn=TRUE, min.panel=100, suggested=TRUE, 
        env = parent.frame(), my.palette=NULL) 

}

\arguments{
  \item{my.inds}{List with the channels information from the individuals specified, usually coming from the \code{\link{storing.inds}} function output}

\item{cols}{The channel/color you wish to analyze, usually 1 is blue, 2 is green, 3 is yellow, 4 is red and so on}

\item{ladder}{A vector containing the expected weights for the ladder peaks that will be found the using the \code{\link{find.ladder}} function}

\item{xlim}{A vector containing the base pair interval where the plot should be drawn}

\item{ylim}{A vector containing the intensity interval where the plot should be drawn}

\item{n.inds}{Vector specifying the plants to be scored}

\item{channel.ladder}{A scalar value indicating in which channel or color the ladder was read}

\item{ploidy}{A scalar value indicating the ploidy of the organism to be scored}

\item{method}{An argument indicating one of the 2 methods available; "cor" makes all possible combination of peaks and searches exhaustive correlations to find the right peaks corresponsding to the expected DNA weights, or "ci" constructing confidence intervals to look for peaks meeting the conditions specified in the previous arguments}

\item{init.thresh}{ An initial value of intensity to detect peaks. We recommend not to deal to much with unless you have highly controlled dna concentrations in your experiment}

\item{ladd.init.thresh}{ A value of intensity to detect peaks in the internal use of the \code{\link{find.ladder}} function. We recommend not to deal to much with it unless you identified special situations with your ladder}

\item{lwd}{ The width of the line}

\item{warn}{A TRUE/FALSE value indicating if warnings should be provided when detecting the ladder}

\item{min.panel}{A scalar value indicating which peak values should be ignored when creating a panel. If 'xlim' values are specified the 'min.panel' value is ignored and instead the panel peaks provided by the program are based in the region where you are zooming in.}

\item{suggested}{ a TRUE/FALSE value statement declaring if you want the program to return suggested peaks for your panel. The default is TRUE but can be anoying if the program draws too many peaks.}

\item{env}{ this is used to detect the environment of the user and load the result in the same environment.}

\item{my.palette}{A character vector with the colors to be used when drawing the RFU plots. If NULL it will use the programmed palette.}

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

Covarrubias-Pazaran G, Diaz-Garcia L, Schlautman B, Salazar W, Zalapa J. Fragman: An R package for fragment analysis. 2016. BMC Genetics 17(62):1-8.

Robert J. Henry. 2013. Molecular Markers in Plants. Wiley-Blackwell. ISBN 978-0-470-95951-0.

Ben Hui Liu. 1998. Statistical Genomics. CRC Press LLC. ISBN 0-8493-3166-8.
}

\examples{
data(my.plants)
my.plants <- my.plants[1]
my.ladder <- c(50, 75, 100, 125, 129, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375)
overview2(my.inds=my.plants, cols = 1, ladder=my.ladder, lwd=1)
# now use:
# my.panel <- locator(type="p", pch=20, col="red")$x
# to click over the peaks and get the sizes in base pairs
# when you are done make sure you press the "Esc" key, do not push the stop button

## to look at many channels at the same time you 
## can use the par(new=TRUE) and a for loop

for(u in 1:4){
  overview2(my.inds=my.plants, cols = u, ladder=my.ladder, lwd=1,
            xlim=c(240,350), ylim=c(0,30000))
  par(new=TRUE)
}

}
