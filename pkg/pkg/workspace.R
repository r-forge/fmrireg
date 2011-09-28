library(gtools)
library(numDeriv)
library(yaml)
library(splines)
library(gregmisc)
library(gsubfn)

wd <- getwd()

if (basename(wd) == "R") {
	prefix = ""
} else if (basename(wd) == "pkg") {
	prefix = "R"	
} else {
	stop(paste("non-standard directory ", wd))
}


fnames <- c("AllGeneric.R",
  "AllClass.R",
  "DataSource-methods.R",
  "EventVector-methods.R",
  "EventTable-methods.R",
  "EventModel-methods.R",
  "FMRIDesign-methods.R",
  "Regressor-methods.R",
  "HRF-methods.R",
  "waver.R",
  "AFNIStim-methods.R",
  "AFNI_GLM.R",
  "ConvolvedTerm-methods.R",
  "Contrast-methods.R")


 

fnames <- paste(prefix, "/", fnames, sep="")
lapply(fnames, source)
