

afni.glm <- function(fmridesign, workingDir=".", mask, options=list()) {

  defopts <- list(fout=TRUE, rout=TRUE, tout=TRUE, float=TRUE, cbucket="coefout", bucket="statout", jobs=2, polort=5)

  for (optname in names(options)) {
    defopts[[optname]] <- options[[optname]]
  }

  
  new("AFNICommand", design=fmridesign, workingDir=workingDir, mask=mask, options=defopts)
}

.makeCommandStr <- function(cmdlines) {
  cmdstr <- lapply(names(cmdlines), function(optname) {
    entry <- cmdlines[[optname]]
    if (is.list(entry)) {
      switchnames <- rep(paste("-", optname, sep=""), length(entry))
      paste(switchnames, entry, collapse=" ")
    } else if (is.null(entry) || length(entry) == 0) {
      ""
    } else if (entry[[1]] == TRUE) {
      paste("-", optname, sep="")
    } else {
      paste(paste("-", optname, sep=""), paste(entry, collapse=" "))
    }
  })
  cmdstr <- Filter(function(x) !is.null(x) & x != "", cmdstr)
  cmdstr <- paste(cmdstr, collapse=" ")
  cmdstr <- paste("3dDeconvolve", cmdstr)
  cmdstr

}
  

setMethod("build", signature(x="AFNICommand"),
          function(x) {

            opts <- getOptions(x)

            desmat <- designMatrix(x@design)
           
            
            #stimlabels <- gsub("\\]", "_", gsub("\\[", "_", colnames(desmat)))
            #stimlabels <- gsub("\\(", "+", gsub("\\)", "+", stimlabels))
            #stimlabels <- gsub("\\,", "@", stimlabels)
            #stimlabels <- gsub(" ", "", stimlabels)

           
            stimlabels <- shortnames(x@design)
            stopifnot(length(unique(stimlabels)) == length(stimlabels))
            
             
            stimfiles <- paste(stimlabels, "_reg.1D", sep="")

            
            cmdlines <- list(input=filelist(x@design, full.names=T),
                             mask=x@mask,
                             polort=opts[["polort"]],
                             num_stimts=length(colnames(x@design)),
                             stim_file=lapply(1:NCOL(desmat), function(i) paste(i, stimfiles[i], collapse=" ")),
                             stim_label=lapply(1:NCOL(desmat), function(i) paste(i, stimlabels[i], collapse=" ")),
                             fout=opts[["fout"]],
                             rout=opts[["rout"]],
                             tout=opts[["tout"]],
                             bout=opts[["bout"]],
                             cbucket=opts[["cbucket"]],
                             bucket=opts[["bucket"]],
                             jobs=opts[["jobs"]],
                             float=TRUE)

         
            cmdstr <- .makeCommandStr(cmdlines)
             
            ret <- list()
            wd <- workingDir(x)
            
            ret$run <- function() {
              
              print(paste("setting directory:", wd))
              setwd(wd)
              
              sapply(1:NCOL(desmat), function(i) {
                write(desmat[,i], file=stimfiles[i], ncolumns=1)
              })

              write(cmdstr, "3ddeconvolve.sh")
              
              
              system(cmdstr)
            }

            ret$command <- cmdstr
            ret
                            
            
            
          })

setMethod("workingDir", signature(x="AFNICommand"),
          function(x) {
            x@workingDir
            
          })

setMethod("getOptions", signature(x="AFNICommand"),
          function(x) {
            x@options           
          })



getfun <- function() {
  cmd <- "ls -lrt"
  x <- list()
  
  x$run <- function() {
    system(cmd)
  }

  x
}
  
  
  
  

  
