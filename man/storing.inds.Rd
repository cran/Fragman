\name{storing.inds}
\alias{storing.inds}
\title{Extracting channel information}

\description{
This function reads the FSA files using a function named 'read.abif' from another R package called seqinr. This will extract the information of the DNA intensities of the capillary electrophoresis and will store it in a data structure know in R as a list. The usage of the function and the arguments it takes are as follows:
}
\usage{
storing.inds(folder, channels=NULL, fourier=TRUE, 
            saturated=TRUE, lets.pullup=TRUE, 
            plotting=FALSE, rawPlot=FALSE,
            llength=3000, ulength=80000)
}

\arguments{
  \item{folder}{A path/directory where the FSA files are located. We recommend to use the Seession tab -> Set working directory and provide that folder as an argument}

\item{channels}{A scalar value indicating how many channels/colors should be found in the FSA files, usually people using the rox375 ladder in the red channel have 4 channels, people using the LIZ ladder have 5 channels. The default is the last channel.}

\item{fourier}{A FALSE/TRUE value indicating if data should be smooth aplying a Fourier transormation using 40 percent of the lowest frequencies. The dafault is TRUE}

\item{saturated}{A FALSE/TRUE value indicating if data should be checked and treated for saturated peaks above 8000 RFU which usually split at the top in 2 different peaks. The dafault is TRUE}

\item{lets.pullup}{A FALSE/TRUE value indicating if data should be treated for noise from channel to channel known as pull up or pull down peaks since wavelengths where the dyes are read usually overlap (blue->green->yellow->red->orange. The dafault is TRUE}

\item{plotting}{A FALSE/TRUE value indicating if results after data cleaning steps should be plotted to asses graphically how data was handled. The dafault is FALSE}

\item{rawPlot}{A FALSE/TRUE value indicating if a plot drawing all vectors read should be plotted. The dafault is FALSE since this consumes a lot of memory.}

\item{llength}{A numeric value indicating how small can be the number of indexes in each channel. Default is 3000 indexes for small runs.}

\item{ulength}{A numeric value indicating how big can be the number of indexes in each channel. Default is 80000 indexes for long runs.}

}
\details{
No major details.
}
\value{
If arguments are correct the function returns a list containing
\describe{
\item{all.inds.mats}{A list where each element is a data frame containing the "n" channels of an individual}
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
### the correct way to do it for a population of inds with 
### 4 colors + ladder= 5 would be:
# my.plants <- storing.inds(folder)
}
