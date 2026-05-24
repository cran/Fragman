##################################################################################################
#Startup function
#this function is executed once the library is loaded
.onAttach = function(library, pkg)
{
  Rv = R.Version()
  if(!exists("getRversion", baseenv()) || (getRversion() < "2.1"))
    stop("This package requires R 2.1 or later")
  if(interactive())
  {
    desc <- utils::packageDescription(pkg)
    packageStartupMessage(paste("## ============================================================= ## "),appendLF=TRUE)
    packageStartupMessage(paste("## ============================================================= ## "),appendLF=TRUE)
    packageStartupMessage(paste("# Fragman: An R package for Fragment Analysis ", desc$Version, ". ",sep=""),appendLF=TRUE)
    packageStartupMessage("# Author: Covarrubias-Pazaran et al. (2016)",appendLF=TRUE)
    packageStartupMessage("# Published: BMC Genetics 17(62):1-8",appendLF=TRUE)
    packageStartupMessage("# Supported and partially funded by:", appendLF=TRUE)
    packageStartupMessage("#    + Council of Science and Technology (CONACYT)", appendLF=TRUE)
    packageStartupMessage("#    + US Department of Agriculture (USDA-ARS)", appendLF=TRUE)
    packageStartupMessage("# Type 'help(Fragman)' for summary information",appendLF=TRUE)
    packageStartupMessage("# Type 'citation(Fragman)' to know how to cite this package",appendLF=TRUE)
    packageStartupMessage("# Update Fragman from time to time using 'install.packages('Fragman')'",appendLF=TRUE)
    packageStartupMessage(paste("## ============================================================= ## "),appendLF=TRUE)
    packageStartupMessage(paste("## ============================================================= ## "),appendLF=TRUE)
  }
  invisible()
}

