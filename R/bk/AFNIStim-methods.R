AFNIStim <- function(events, basis="GAM", basis.params=list()) {
  new("AFNIStim", events=events, basis=basis, basis.params=basis.params)
}


.buildBLOCK <- function(stimlabel, stimfile, stimnum, duration) {
  paste(paste("-stim_times", stimnum, stimfile,
        paste("\'BLOCK(", duration, ")\'"),
        paste("-stim_label", stimnum, stimlabel)))
}

.buildGAM <- function(stimlabel, stimfile, stimnum) {
  paste(paste("-stim_times", stimnum, stimfile,
        paste("\"GAM\""),
        paste("-stim_label", stimnum, stimlabel)))
}

.buildSPMG <- function(stimlabel, stimfile, stimnum) {
  paste(paste("-stim_times", stimnum, stimfile,
        paste("\"SPMG\""),
        paste("-stim_label", stimnum, stimlabel)))
}
.buildSPMG3 <- function(stimlabel, stimfile, stimnum) {
  paste(paste("-stim_times", stimnum, stimfile,
        paste("\"SPMG3\""),
        paste("-stim_label", stimnum, stimlabel)))
}

.buildBasis <- function(basisname, stimlabel, stimfile, stimnum, start, end, nparams) {
  paste("-stim_times", stimnum, stimfile,
    paste("\'", basisname, "(", start, ",", end, ",", nparams, ")\'", sep=""),
    paste("-stim_label", stimnum, stimlabel),
    paste("-iresp", stimnum, paste(stimlabel, "_iresp", sep="")))

}

setMethod("stimlabels", signature(x="AFNIStim"),
          function(x) {
            return(names(cells(x@events)))
          })

setMethod("stimfiles", signature(x="AFNIStim"),
          function(x) {
            return(paste(names(cells(x@events)), ".reg.1D", sep=""))
          })


setMethod("writeStims", signature(x="AFNIStim", dir="character"),
          function(x, dir) {
            
            .writeOnsets <- function(outname, eventList) {
              hfile = file(outname, "w")
              lapply(eventList, function(line) {
                if (length(line) == 0) {
                  cat("*", "\n", file=hfile)
                } else {
                  cat(line, "\n", file=hfile)
                }
              })

              close(hfile)
            }
            
            
            fnames <- stimfiles(x)
            slabels <- stimlabels(x)

            cell.list <- cells(x@events)

            
            lapply(1:length(fnames), function(i) {
              cell <- cell.list[[i]]
              fname <- fnames[i]
              eventList <- split(cell$onsets, cell$blocknum)
              .writeOnsets(paste(dir, "/", fname, sep=""), eventList)
              
            })
          })
                   
            
         
           
setMethod("basis", signature(x="AFNIStim"), function(x) x@basis)

setMethod("buildArgs", signature(x="AFNIStim"),
          function(x) {

            slabels <- stimlabels(x)
            snames <- stimfiles(x)
            switch(basis(x),
                   BLOCK=lapply(1:length(slabels), function(i) .buildBlock(slabels[i], snames[i], i, 2)),
                   
                   GAM=lapply(1:length(slabels), function(i) .buildGAM(slabels[i], snames[i], i)),
                   
                   SPMG=lapply(1:length(slabels), function(i) .buildSPMG(slabels[i], snames[i], i)),
                   
                   SPMG3=lapply(1:length(slabels), function(i) .buildSPMG(slabels[i], snames[i], i)),
                   
                   TENT=lapply(1:length(slabels), function(i) .buildBasis(x@basis, slabels[i], snames[i], i,
                     x@basis.params$start, x@basis.params$end, x@basis.params$nparams)),
                   
                   CSPLIN=lapply(1:length(slabels), function(i) .buildBasis(basis(x), slabels[i], snames[i], i,
                     x@basis.params$start, x@basis.params$end, x@basis.params$nparams)),
                   
                   POLY=lapply(1:length(slabels), function(i) .buildBasis(basis(x), slabels[i], snames[i], i,
                     x@basis.params$start, x@basis.params$end, x@basis.params$nparams)))
            
          })
          
