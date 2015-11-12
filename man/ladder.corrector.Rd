\name{ladder.corrector}
\alias{ladder.corrector}
\title{Ladder corrector attached to R environment}

\description{
This function was designed to correct manually samples that could not be correctly detected by the \code{\link{ladder.info.attach}} function and allows the user to select manually the peaks he knows are the correct peaks. This function uses the output of the \code{\link{ladder.info.attach}} function which is basically a vector with the names of the samples that were too dificult for the algorithm to find. The console will draw a plot and will ask the user to click over the peaks expected for a ladder provided, once the user is done should press the 'esc' key and continue to the next sample. This process is repeated until the all samples with the names provided are adjusted. 

}
\usage{
ladder.corrector(stored, to.correct, ladder, 
      thresh=200, env = parent.frame())
}

\arguments{
 \item{stored}{List with the channels information from the individuals specified, usually coming from the \code{\link{storing.inds}} function output.}
 
 \item{to.correct}{Vector containing the names ofthe samples to be corrected and usually the output of the \code{\link{ladder.info.attach}} function.}

\item{ladder}{Vector containing the expected weights of the dna fragments of the ladder in use.}

\item{thresh}{A scalar value indicating the minimum value in RFUs to look for peaks to match the user-selected peaks with the program recognized peaks.}

\item{env}{ this is used to detect the environment of the user and load the result in the same environment.}

}
\details{
Once the user has selected the right peaks the function will fix the ladder and attach such information to the R environment.

}
\value{

This function does not produce any output.
\describe{
\item{$data}{the program will attach the information to the R environment}
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
my.ladder <- c(50, 75, 100, 125, 129, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375)
ladder.info.attach(stored=my.plants, ladder=my.ladder, ladd.init.thresh=300)
## now if something goes wrong use the corrector:
#ladder.corrector(stored=my.plants, 
#to.correct="FHN152-CPN01_01A_GH1x35_152-148-209_717-704-793_367-382-381.fsa", 
#ladder=my.ladder)
}
