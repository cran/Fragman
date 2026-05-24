\name{detect.ladder}
\alias{detect.ladder}
\title{Ladder detection by correlation or confidence intervals}

\description{
This function takes a vector of color heights/intensities from the fragment analysis containing the ladder/standard channel, and detects the biggest peaks where the derivative is equal zero and uses the information from the expected weights for the ladder to construct confidence intervals in order to detect the ladder peaks.

Please! if using the confidence interval method ("ci"), which is NOT the default,  once you have found the best parameters for the arguments to match your ladder using this function, please pass those values to all the posterior functions, making sure the 'dev' argument is passed to the new functions. If using the correlation method ("cor"), don't worry about it.
}
\usage{

detect.ladder(stored, ind=1, ladder, channel.ladder=dim(stored[[1]])[2],
              ci.upp=1.96, ci.low=1.96, draw=TRUE, dev=50, warn=TRUE, 
              init.thresh=250,sep.index=8, method="cor", avoid=1500, who="sample")

}

\arguments{
 \item{stored}{Lis of dataframes obtained by using the \code{\link{storing.inds}} function.}
 
\item{ind}{The individual that you wish to analyze for assessing that the ladder was correctly detected.}

\item{ladder}{Vector containing the expected weights of the dna fragments of the ladder in use}

\item{channel.ladder}{A scalar value indicating in which channel or color the ladder was read}

\item{ci.upp}{A scalar value indicating how many standar errors will be used to detect peaks when checking the height of the ladder peaks(upper bound). To be used in the \code{\link{find.ladder}} function}

\item{ci.low}{A scalar value indicating how many standar errors will be used to detect peaks when checking the height of the ladder peaks(lower bound). To be used in the 'find.ladder' function}

\item{draw}{A TRUE/FALSE value indicating if the plot for the ladder found should be printed or not}

\item{dev}{A scalar value indicating the number of indexes to be used as peak separation when deciding the ladder peaks. Some ladders contain dna fragments of very closed weights and modifying this parameter helps to detect them correctly}

\item{warn}{A TRUE/FALSE value indicating if warnings should be provided when detecting the ladder}

\item{init.thresh}{An initial value of color intensity to be used when detecting the ladder, could be really important for the correlation method}

\item{sep.index}{A scalar value indicating how many indexes should be allowed to considered a true peak from noisy peaks}

\item{method}{An argument indicating one of the 2 methods available; "cor" makes all possible combination of peaks and searches exhaustive correlations to find the right peaks corresponsding to the expected DNA weights, or "ci" constructing confidence intervals to look for peaks meeting the conditions specified in the previous arguments}

\item{who}{A name to indicate which sample is being analyzed}

\item{avoid}{A scalar value indicating how many indexes should be avoided when the method of correlation fails to find peaks and a random sample will be drawn from the existing peaks. The default is 1500 indexes which will samples peaks avoiding the first 1500 indexes which is usually related to noisy area in some ladders.}
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

Covarrubias-Pazaran G, Diaz-Garcia L, Schlautman B, Salazar W, Zalapa J. Fragman: An R package for fragment analysis. 2016. BMC Genetics 17(62):1-8.

Robert J. Henry. 2013. Molecular Markers in Plants. Wiley-Blackwell. ISBN 978-0-470-95951-0.

Ben Hui Liu. 1998. Statistical Genomics. CRC Press LLC. ISBN 0-8493-3166-8.
}

\examples{
data(my.plants)
my.ladder <- c(120, 125, 129, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375)
# looking at the first individual
detect.ladder(stored=my.plants, ind=1, ladder=my.ladder)
}
