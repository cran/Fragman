\name{read.abif}
\alias{read.abif}
\title{Read ABIF formatted files}
\description{
ABIF stands for Applied Biosystem Inc. Format, a binary fromat modeled after TIFF format.
Corresponding files usually have an .ab1 or .fsa extension.
}
\usage{
read.abif(filename, max.bytes.in.file = file.info(filename)$size, 
 pied.de.pilote = 1.2, verbose = FALSE)
}
\arguments{
  \item{filename}{The name of the file.}
  \item{max.bytes.in.file}{The size in bytes of the file, defaulting to what is returned by file.info}	
  \item{pied.de.pilote}{Safety factor: the argument readBin is set as pied.de.pilote*max.bytes.in.file.}	
  \item{verbose}{logical [FALSE]. If TRUE verbose mode is on.}	
}
\details{
All data are imported into memory, there is no attempt to read items on the fly.
}
\value{
A list with three components: Header which is a list that contains various low-level information,
among which numelements is the number of elements in the directory and dataoffset
the offset to find the location of the directory. Directory is a data.frame for the directory
of the file with the number of row being the number of elements in the directory and the 7
columns describing various low-level information about the elements.  
}
\references{ 

Anonymous (2006) Applied Biosystem Genetic Analysis Data File Format.
Available at \url{http://www.appliedbiosystems.com/support/software_community/ABIF_File_Format.pdf}.
Last visited on 03-NOV-2008.

The figure in the example section is an attempt to reproduce figure 1A from:

Krawczyk, J.,  Goesmann, A., Nolte, R., Werber, M., Weisshaar, B. (2009)
Trace2PS and FSA2PS: two software toolkits for converting trace and fsa files to PostScript format.
\emph{Source Code for Biology and Medicine}, \bold{4}:4.

}
\author{J.R. Lobry}


\examples{
#No examples provided, please download seqinr package for going deeper in this function.

}
