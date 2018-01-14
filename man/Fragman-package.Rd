\name{Fragman-package}
\alias{Fragman}
\docType{package}
\author{
Giovanny Covarrubias-Pazaran, Luis Diaz-Garcia, Brandon Schlautman, Walter Salazar, Juan Zalapa.
}
\title{Fragment analysis and automatic scoring}
\description{Fragman is a package designed for Fragment analysis and automatic scoring of biparental populations (such as F1, F2, BC types) and populations for diversity studies. The program is designed to read files with FSA extension (which stands for FASTA-type file and contains lectures for DNA fragments), and .txt files from Beckman CEQ 8000 system, and extract the DNA intensities from the channels/colors where they are located, based on ABi machine plattforms to perform sizing and allele scoring. 

The core of the package and the workflow of the fragment analysis rely in the following 4 functions; 

1) \code{\link{storing.inds}}(function in charge of reading the FSA or txt(CQS) files and storing them with a list structure)

2) \code{\link{ladder.info.attach}} (uses the information read from the FSA files and a vector containing the ladder information (DNA size of the fragments) and matches the peaks from the channel where the ladder was run with the DNA sizes for all samples. Then loads such information in the R environment for the use of posterior functions)

3) \code{\link{overview2}} (create friendly plots for any number of individuals specified  and can be used to design panels (\code{\link{overview2}}) for posterior automatic scoring (like licensed software does), or make manual scoring (\code{\link{overview}}) of individuals such as parents of biparental populations or diversity populations)

4) The \code{\link{score.markers}} (function score the alleles by finding the peaks provided in the panel (if provided), otherwise returns all peaks present in the channel). Thisfinal function can be automatized if several markers are located in the same channel by creating lists of panels taking advantage of R capabilities and data structures.

** Sometimes during the ladder sizing process some samples can go wrong for several reasons related to the sample quality (low intensity in ladder channel, extreme number of noisy peaks, etc.), because of that we have introduced \code{\link{ladder.corrector}} function which allows the user to correct the bad samples by clicking over the real peaks, by default the \code{\link{ladder.info.attach}} function returns the names of the samples that had a low correlation with the expected peaks.

When automatic scoring is not desired the function \code{\link{overview}} can be used for getting an interactive session and click over the peaks (using the \code{\link{locator}} function) in order to get the allele sizes.
}

\section{Contact}{

Feel free to contact us with questions and improvement suggestions at:

covarrubiasp@wis.edu
	
Just send a sample file with your question to recreate the issue or bug reported along with vector for your ladder.
}

\section{Citation}{

We have spent valuable time developing this package, please cite it in your publication:

Covarrubias-Pazaran G, Diaz-Garcia L, Schlautman B, Salazar W, Zalapa J. Fragman: An R package for fragment analysis. 2016. BMC Genetics 17(62):1-8.
}

\references{
Covarrubias-Pazaran G, Diaz-Garcia L, Schlautman B, Salazar W, Zalapa J. Fragman: An R package for fragment analysis. 2016. BMC Genetics 17(62):1-8.

Robert J. Henry. 2013. Molecular Markers in Plants. Wiley-Blackwell. ISBN 978-0-470-95951-0.

Ben Hui Liu. 1998. Statistical Genomics. CRC Press LLC. ISBN 0-8493-3166-8.
}
\keyword{ package }
\seealso{
http://cggl.horticulture.wisc.edu/home-page/
}
\examples{
## ================================= ##
## ================================= ##
## Fragment analysis requires 
## 1) loading your data
## 2) matching your ladder
## 3) define a panel for scoring
## 4) score the samples
## ================================= ##
## ================================= ##

#####################
## 1) Load your data
#####################

### you would use something like:
# folder <- "~/myfolder"
# my.plants <- storing.inds(folder)
### here we just load our sample data and use the first 2 plants

?my.plants
data(my.plants)
my.plants <- my.plants[1:2]
class(my.plants) <- "fsa_stored"
# plot(my.plants) # to visualize the raw data

#######################
## 2) Match your ladder
#######################

### create a vector indicating the sizes of your ladder and do the match

my.ladder <- c(50, 75, 100, 125, 129, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375)
ladder.info.attach(stored=my.plants, ladder=my.ladder)

### matching your ladder is a critical step and should only happen once per batch of 
### samples read

###****************************************************************************************###
### OPTIONAL:
### If the ladder.info attach function detects some bad samples 
### that you can correct them manually using
### the ladder.corrector() function
### For example to correct one sample in the previous data
### ladder.corrector(stored=my.plants, 
#to.correct="FHN152-CPN01_01A_GH1x35_152-148-209_717-704-793_367-382-381.fsa", 
#ladder=my.ladder)
###****************************************************************************************###

#######################
## 3) Define a panel
#######################

### In fragment analysis you usually design a panel where you indicate
### which peaks are real. You may use the overview2 function which plots all the
### plants in the channel you want in the base pair range you want

overview2(my.inds=my.plants, channel = 2:3, ladder=my.ladder, init.thresh=5000)

### You can click on the peaks you think are real, given that the ones
### suggested by the program may not be correct. This can be done by using the 
### 'locator' function and press 'Esc' when you're done, i.e.:
# my.panel <- locator(type="p", pch=20, col="red")$x
### That way you can click over the peaks and get the sizes
### in base pairs stored in a vector named my.panel

### Just for demonstration purposes I will use the suggested peaks by 
### the program using overview2, which will return a vector with 
### expected DNA sizes to be used in the next step for scoring
### we'll do it in the 160-190 bp region

my.panel <- overview2(my.inds=my.plants, channel = 3, 
                    ladder=my.ladder, init.thresh=7000, 
                    xlim=c(160,190)); my.panel

##########################
## 4) Score the samples
##########################

### When a panel is created is time to score the samples by providing the initial
### data we read, the ladder vector, the panel vector, and our specifications
### of channel to score (other arguments are available)

### Here we will score our samples for channel 3 with our panel created previously

res <- score.markers(my.inds=my.plants, channel = 3, panel=my.panel$channel_3,
                ladder=my.ladder, electro=FALSE)

### Check the plots and make sure they were scored correctly. In case some samples 
### are wrong you might want to use the locator function again and figure out 
### the size of your peaks. To extract your peaks in a data.frame do the following:

final.results <- get.scores(res)
final.results 
}