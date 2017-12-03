\name{overview}
\alias{overview}
\title{Assesing several plants with an overview}

\description{
This function uses information from the FSA files read from \code{\link{storing.inds}} function and creates a plot to assess graphically the peaks of several plants in certain channel in order to score manually or assess the parental fragments in the case of biparentla ppulations. If you desire to create a panel you may want to take a look at \code{\link{overview2}}. The function contains several defaults in most of the arguments, please check arguments but in general.

}
\usage{
overview(my.inds, cols = 1, n.inds = c(1:length(my.inds)), 
        xlimi=c(min(ladder),max(ladder)), ladder, channel.ladder=dim(my.inds[[1]])[2], 
        ploidy=2, ci.upp=1.96, ci.low=1.96, dev=50, method="iter", 
        init.thresh=200, ladd.init.thresh=200, warn=TRUE, my.palette=NULL,
        env = parent.frame())
}

\arguments{
  \item{my.inds}{List with the channels information from the individuals specified, usually coming from the \code{\link{storing.inds}} function output}

\item{cols}{The channel you wish to analyze, usually 1 is blue, 2 is green, 3 is yellow, 4 is red and so on}

\item{n.inds}{Vector specifying the plants to be scored}

\item{xlimi}{A vector containing the base pair interval where the plot should be drawn}

\item{ladder}{A vector containing the expected weights for the ladder peaks that will be found the using the \code{\link{find.ladder}} function}

\item{channel.ladder}{A scalar value indicating in which channel or color the ladder was read}

\item{ploidy}{A scalar value indicating the ploidy of the organism to be scored}

\item{ci.upp}{A scalar value indicating how many standar errors will be used to detect peaks when checking the height of the ladder peaks(upper bound). To be used in the \code{\link{find.ladder}} function}

\item{ci.low}{A scalar value indicating how many standar errors will be used to detect peaks when checking the height of the ladder peaks(lower bound). To be used in the \code{\link{find.ladder}} function}

\item{dev}{A scalar value indicating the number of indexes to be used as peak separation when deciding the ladder peaks, for more details check \code{\link{find.ladder}} function}

\item{method}{An argument indicating one of the 2 methods available; "cor" makes all possible combination of peaks and searches exhaustive correlations to find the right peaks corresponsding to the expected DNA weights, or "ci" constructing confidence intervals to look for peaks meeting the conditions specified in the previous arguments}

\item{init.thresh}{ An initial value of intensity to detect peaks. We recommend not to deal to much with unless you have highly controlled dna concentrations in your experiment}

\item{ladd.init.thresh}{ A value of intensity to detect peaks in the internal use of the \code{\link{find.ladder}} function. We recommend not to deal to much with it unless you identified special situations with your ladder}

\item{warn}{A TRUE/FALSE value indicating if warnings should be provided when detecting the ladder}

\item{my.palette}{A character vector with the colors to be used when drawing the RFU plots. If NULL it will use the programmed palette.}

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
Covarrubias-Pazaran G, Diaz-Garcia L, Schlautman B, Salazar W, Zalapa J. Fragman: An R package for fragment analysis. 2016. BMC Genetics 17(62):1-8.

Robert J. Henry. 2013. Molecular Markers in Plants. Wiley-Blackwell. ISBN 978-0-470-95951-0.

Ben Hui Liu. 1998. Statistical Genomics. CRC Press LLC. ISBN 0-8493-3166-8.
}

\examples{
data(my.plants)
my.plants <- my.plants[1:10]
my.ladder <- c(50, 75, 100, 125, 129, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375)
overview(my.inds=my.plants, cols = 1, n.inds = c(1:5), ladder=my.ladder, xlim=c(200,220))
# now use:
# locator(type="p", pch=20, col="red")$x
# to click over the peaks and get the sizes in base pairs
# when you are done make sure you press the "Esc" key, 
# do not push the stop button, some versions of R usually crash 
# by stopping instead of pressing 'Esc'.
}
