\name{Fragman-package}
\alias{Fragman}
\docType{package}
\author{
Giovanny Covarrubias-Pazaran, Luis Diaz-Garcia, Brandon Schlautman, Walter Salazar, Juan Zalapa.
}
\title{Fragment analysis and automatic scoring}
\description{Fragman is a package designed for Fragment analysis and automatic scoring of biparental populations (such as F1, F2, BC types) and populations for diversity studies. The program is designed to read files with FSA extension (which stands for FASTA-type file and contains lectures for DNA fragments) and extract the DNA intensities from the channels/colors where they are located, based on ABi machine plattforms to perform sizing and allele scoring. The core of the package relays in 4 functions; 1) \code{\link{storing.inds}} is the function in charge of reading the FSA files and storing them with a list structure, 2) \code{\link{ladder.info.attach}} uses the information read from the FSA files and a vector containing the ladder information (DNA size of the fragments) and matches the peaks from the channel where the ladder was run with the DNA sizes for all samples. Then loads such information in the R environment for the use of posterior functions, 3) \code{\link{overview}} & \code{\link{overview2}} create friendly plots for any number of individuals specified  and can be used to design panels (overview2) for posterior automatic scoring, or make manual scoring (overview) of individuals such as parents of biparental populations or diversity panels, 4) The \code{\link{score.easy}} function score the alleles by finding all regions where the first derivative of the intensity vector iz zero and reduces the search of peaks using a panel (if provided) otherwise returns all peaks present. This function can be automatized if several markers are located in the same channel by creating lists of panels taking advantage of R capabilities and data structures (see vignettes). 

Once the calls have been obtained we can extract a data frame with the \code{\link{get.scores}} function. In addition if a mapping population is being analyzed we can transform those calls to the joinmap format using the \code{\link{jm.conv}} function.

When automatic scoring is not desired the function \code{\link{overview}} can be used for getting an interactive session and click over the peaks (using the \code{\link{locator}} function) in order to get the allele sizes.
	
Vignettes illustrating some of the features of this package can be found at `http://cggl.horticulture.wisc.edu/home-page/`.

We have spent valuable time developing this package, please cite it in your publication:

Covarrubias-Pazaran G, Diaz-Garcia L, Schlautman B, Salazar W, Zalapa J. Fragman: An R package for fragment analysis. http://horticulture.wisc.edu/cggl/ZalapaLab/People.html. 2015.
}

\references{
Covarrubias-Pazaran G, Diaz-Garcia L, Schlautman B, Salazar W, Zalapa J. (2015) Fragma: An R package for fragment analysis. R package version 1.0. URL https://cran.r-project.org/web/packages/Fragman/.

Robert J. Henry. 2013. Molecular Markers in Plants. Wiley-Blackwell. ISBN 978-0-470-95951-0.

Ben Hui Liu. 1998. Statistical Genomics. CRC Press LLC. ISBN 0-8493-3166-8.
}
\keyword{ package }
\seealso{
http://cggl.horticulture.wisc.edu/home-page/
}
\examples{
#####################
## LOAD YOUR DATA ###
#####################

### you would use:
# my.plants <- storing.inds(folder)
### where folder is the path where your samples are, i.e. "~/Documents"
### here we just load our example data
?my.plants
data(my.plants)
my.plants <- my.plants[1:2]

#######################
## MATCH YOU LADDER ###
#######################

### create a vector indicating the sizes of your ladder
my.ladder <- c(120, 125, 129, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375)
### match your ladder to the peaks and attach the information 
### to the R environment using the function:
ladder.info.attach(stored=my.plants, ladder=my.ladder)

#######################
## CREATE A PANEL   ###
#######################

### you may use overview2 or overview to create you customized panel using:
### here we select the channel 3 (yellow) by setting 'cols=3' 
### and providing the samples and ladder
overview2(my.inds=my.plants, cols = 3, ladder=my.ladder, init.thresh=5000)
### you could also click on the peaks you think are real 
### by using the 'locator' function and press 'Esc' when you're done:
# my.panel <- locator(type="p", pch=20, col="red")$x
### so you can click over the peaks and get the sizes
### in base pairs stored in a vector named my.panel
### Instead of doing that I will use the suggested peaks by 
### the program using overview2, which provides a vector with 
### expected DNA sizes to be used in the next step for scoring
### we'll do it in the 160-190 bp region
my.panel <- overview2(my.inds=my.plants, cols = 3, 
 ladder=my.ladder, init.thresh=7000, xlim=c(160,190)); my.panel

##########################
## SCORE YOUR SAMPLES  ###
##########################
a <- score.easy (my.inds=my.plants, cols = 3, panel=my.panel, ladder=my.ladder, electro=FALSE)
### extract your peaks in a data.frame
final.results <- get.scores(a)
final.results 
}