\name{score.easy}
\alias{score.easy}
\title{Fragment analysis scoring}

\description{
This function uses information from the fsa files read from \code{\link{storing.inds}} function and does the ssr calling in the channel specified and returns the index position, height and base pair position.

}
\usage{
score.easy(my.inds, cols = 1, n.inds = NULL, panel=NULL, shift=0.8,
          ladder, channel.ladder=NULL, 
          ploidy=2, left.cond=c(0.6,3), right.cond=0.35, warn=FALSE, 
          window=0.5, init.thresh=200, ladd.init.thresh=200, 
          method="iter2", env = parent.frame(), 
          plotting=TRUE, electro=FALSE, pref=3)
}

\arguments{
  \item{my.inds}{List with the channels information from the individuals specified, usually coming from the \code{\link{storing.inds}} function output}

\item{cols}{The channel you wish to analyze, usually 1 is blue, 2 is green, 3 is yellow, 4 is red and so on}

\item{n.inds}{Vector specifying the plants to be scored}

\item{panel}{A vector containing the base pair interval where the peaks should be searched for}

\item{shift}{The number of base pairs to be used for discarding neighboring peaks to the tallest peaks, i.e. if 2 peaks are 0.3 bp together the smalles will be discarded}

\item{ladder}{A vector containing the expected weights for the ladder peaks that will be found the using the \code{\link{find.ladder}} function}

\item{channel.ladder}{A scalar value indicating in which channel or color the ladder was read}

\item{ploidy}{A scalar value indicating the ploidy of the organism to be scored to decide the maximum number of peaks the program should look for. TO BE IMPLEMENTED SOON. STILL NOT FUNCTIONAL.}

\item{left.cond}{A percentage value (0-1) indicating when peaks to the left of the tallest peaks should be considered real based on the height, i.e. a value of 0.5 would mean that a close peak (to the left of the tallest peak) will be picked only if such peak is at least 50 percent as tall with respect to the tallest peak. The second argument is the number of base pair indicating when peaks to the left of the tallest peaks should be considered real based on the distance, i.e. a value of 3 would mean that a close peak (to the left of the tallest peak) will be picked only if such peak is at least 3 base pairs far away from the tallest peak}

\item{right.cond}{A percentage value (0-1) indicating when peaks to the right of the tallest peaks should be considered real based on the height, i.e. a value of 0.5 would mean that a close peak (to the right of the tallest peak) will be picked only if such peak is at least 50 percent as tall with respect to the tallest peak.}

\item{warn}{A TRUE/FALSE value indicating if warnings should be provided when detecting the ladder}

\item{window}{A value in base pairs indicating how much is the error for detecting a peak in a sample when providing a panel with expected peaks.}

\item{init.thresh}{ An initial value of intensity to detect peaks. We recommend not to deal to much with it unless you have highly controlled dna concentrations in your experiment.}

\item{ladd.init.thresh}{ If samples were not sized using the info.ladder.attach function this value will be used to detect ladder peaks. Internally the program will use the \code{\link{find.ladder}} function. We recommend not to deal to much with it unless you identified special situations with your ladder}

\item{method}{If samples were not sized using the info.ladder.attach function this method will be used to detect ladder peaks. An argument indicating one of the 3 methods available; "cor" makes all possible combination of peaks and searches exhaustive correlations to find the right peaks corresponsding to the expected DNA weights, or "ci" constructing confidence intervals to look for peaks meeting the conditions specified in the previous arguments, "iter2" an iterative procedure looking for the most likely peaks meeting your ladder expectation. Default is "iter2".}

\item{env}{ this is used to detect the environment of the user and load the result in the same environment. Don't mess with it please.}

\item{plotting}{ a TRUE/FALSE value indicating if the plots should be drawn or not. The default value is TRUE.}

\item{electro}{A TRUE/FALSE value indicating if the electrogram/gel should be drawn or not. The default value is FALSE.}

\item{pref}{A scalar value indicating how many plots should be drawn in the output plotting. The dafault is 3.}
}
\details{

Method "ci" has been depreciated, currently the method "iter2" is the default and uses the ladder provided and observed peaks to match them using an iterative procedure based on least squares.
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

Covarrubias-Pazaran G, Diaz-Garcia L, Schlautman B, Salazar W, Zalapa J. Fragman: An R package for fragment analysis. 2016. BMC Genetics 17(62):1-8.

Robert J. Henry. 2013. Molecular Markers in Plants. Wiley-Blackwell. ISBN 978-0-470-95951-0.

Ben Hui Liu. 1998. Statistical Genomics. CRC Press LLC. ISBN 0-8493-3166-8.
}

\examples{
## ================================= ##
## ================================= ##
##    FIRST PART OF THE ANALYSIS
## LOAD DATA, SET LADDER, MATCH LADDER 
## ================================= ##
## ================================= ##

#####################
## LOAD YOUR DATA ###
#####################

### you would use:
# my.plants <- storing.inds(folder)
### where folder is the path where your samples are, i.e. "~/Documents"
### here we just load our example data and use the first 2 plants

?my.plants
data(my.plants)
my.plants <- my.plants[1:2]

#######################
## MATCH YOU LADDER ###
#######################

### create a vector indicating the sizes of your ladder

my.ladder <- c(50, 75, 100, 125, 129, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375)

### match your ladder to the peaks and attach the information 
### to the R environment using the function: (DO ONLY ONCE PER BATCH)

ladder.info.attach(stored=my.plants, ladder=my.ladder)

###****************************************************************************************###
### OPTIONAL:
### If the function detects some bad samples you can correct them manually using
### the ladder.corrector() function, i.e.:
### ladder.corrector(stored=my.plants, 
#to.correct="FHN152-CPN01_01A_GH1x35_152-148-209_717-704-793_367-382-381.fsa", 
#ladder=my.ladder)
###****************************************************************************************###

## ================================= ##
## ================================= ##
##    SECOND PART OF THE ANALYSIS
## CREATE PANEL, SCORE SAMPLES 
## ================================= ##
## ================================= ##

#######################
## CREATE A PANEL   ###
#######################

### In fragment analysis you usually design a panel where you indicate
### which peaks are real. You may use the overview2 function which plots all the
### plants in the channel you want in the base pair range you want

### Just to show the uptput. Here we select the channel 3 (yellow) by setting 'cols=3' 
### and providing the samples (my.plants) and ladder (my.ladder)

overview2(my.inds=my.plants, cols = 3, ladder=my.ladder, init.thresh=5000)

### You could also click on the peaks you think are real if the ones
### selected by the program are not correct. This can be done by using the 
### 'locator' function and press 'Esc' when you're done, i.e.:

# my.panel <- locator(type="p", pch=20, col="red")$x

### That way you can click over the peaks and get the sizes
### in base pairs stored in a vector named my.panel
### Just for demonstration purposes I will use the suggested peaks by 
### the program using overview2, which will return a vector with 
### expected DNA sizes to be used in the next step for scoring
### we'll do it in the 160-190 bp region
### KEEP IN MIND THIS IS NOT THE BEST WAY TO DO IT, BETTER
### USE "my.panel <- locator(type="p", pch=20, col="red")$x" AND SELECT MANUALLY

my.panel <- overview2(my.inds=my.plants, cols = 3, 
                    ladder=my.ladder, init.thresh=7000, 
                    xlim=c(160,190)); my.panel

##########################
## SCORE YOUR SAMPLES  ###
##########################

### When a panel is created is time to score the samples by providing the initial
### data we read, the ladder vector, the panel vector, and our specifications
### of channel to score (other arguments are available)

### Here we will score our samples for channel 3 with our panel created previously

a <- score.easy (my.inds=my.plants, cols = 3, panel=my.panel,
                ladder=my.ladder, electro=FALSE)

### Check the plots and make sure they were scored correctly. In case some samples 
### are wrong you might want to use the locator function again and figure out 
### the size of your peaks. To extract your peaks in a data.frame do the following:

final.results <- get.scores(a)
final.results 
}
