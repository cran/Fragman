\name{find.ladder}
\alias{find.ladder}
\title{Ladder detection by correlation or confidence intervals}

\description{
This function takes a vector of color heights/intensities from the fragment analysis containing the ladder/standard channel, and detects the biggest peaks where the derivative is equal zero and uses the information from the expected weights for the ladder to construct confidence intervals in order to detect the ladder peaks.

Please! if using the confidence interval method ("ci"), which is not the default, once you have found the best parameters for the arguments to match your ladder using this function, please pass those values to all the posterior functions, please make sure the 'dev' argument is passed to the new functions.
}
\usage{
find.ladder(x, ladder, draw=TRUE, dev=50, warn=TRUE, init.thresh=NULL, 
            sep.index=8, method=NULL, reducing=NULL, who="sample", 
            attempt=10, cex.title=0.8)
}

\arguments{
 \item{x}{Vector of heights from the ladder channel. See example to see how to access to it.}

\item{ladder}{Vector containing the expected weights of the dna fragments of the ladder in use}

\item{draw}{A TRUE/FALSE value indicating if the plot for the ladder found should be printed or not}

\item{dev}{A scalar value indicating the number of indexes to be used as peak separation when deciding the ladder peaks. Some ladders contain dna fragments of very closed weights and modifying this parameter helps to detect them correctly}

\item{warn}{A TRUE/FALSE value indicating if warnings should be provided when detecting the ladder}

\item{init.thresh}{An initial value of color intensity to be used when detecting the ladder}

\item{sep.index}{A scalar value indicating how many indexes should be allowed to considered a true peak from noisy peaks}

\item{method}{An argument indicating one of the 2 methods available; "cor" makes all possible combination of peaks and searches exhaustive correlations to find the right peaks corresponsding to the expected DNA weights, or "ci" constructing confidence intervals to look for peaks meeting the conditions specified in the previous arguments}

\item{who}{A name to indicate which sample is being analyzed}

\item{attempt}{A scalar value indicating how many attempts should be made to find the real ladder peaks. By default is 7 attempts, which means that will try to build the model assuming that the first peak found in the ladder is the corresponding first peak of the expected ladder, then moves to the 2nd peak until the 7th and the seven models are compared picking the most likely model based on the R2 value for each of the models.}

\item{reducing}{A vector of values to reduce the search of peaks to certain indexes in the x axis. Default is NULL so it looks for all peaks for matching the ladder.}

\item{cex.title}{A scalar value indicating how big the title (name of the sample) in the plot should be.}

}
\details{
We have implemented 3 methods for sizing the ladder, each with their advantages and disadvantages. The default method named "red" which stands for "reduction" detect the region where peaks exist (in indexes) in the ladder channel and assumes that your ladder should have some equivalence in indexes and creates an 'expected ladder', then the putative ladder moves along the peak region and correlations and squared distances to the closest peaks are calculated. We have define the coefficient of similarity (CS) as cor(x,y)/var(z), where:

cor(x,y) are the correlations between expected and observed peaks, and
var(z)  is the sum of squares between the differences of expected and observed peaks.

This value usually let us identify the most likely peaks and then all possible combinations for those peaks are computed followed by exhaustive correlations of those combinations with the actual ladder. The highest correlation usually points to the right peaks, which is selected.

In addition the method "cor" is the previous version to "red" which doesn't reduce the search of peaks and computes all possible combinations of peaks from the beggining, with the drawback that slows down the detection process especially when the ladder intensities are low and noisy peaks exist in abundance.

The last method that has been superseded by the previous 2 is the "ci" method based on confidence intervals, which assumes that real ladder peaks have more or less the same intensity and a they can be found by finding the median intensity and computing a 90 percent confidence interval to find the rest of the peaks. This method has been proved to fail when the first condition is broken and ladder have real peaks with intensities greater than the expected.

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
my.ladder <- c(50, 75, 100, 125, 129, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375)
find.ladder(my.plants[[1]][,4], ladder=my.ladder)
}
