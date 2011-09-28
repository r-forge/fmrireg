setMethod("buildAFNIStims", signature(x="AFNITerm"),
	function(x, iresp=FALSE, tr.times=1) {
		
		#stimlabels <- shortnames(x, exclude.basis=TRUE)
		split.ons <- splitOnsets(x, global=TRUE)
		stimlabels <- names(split.ons)
		stimfiles <- paste(stimlabels, "_reg.1D", sep="")		
		
		
		lapply(1:length(stimlabels), function(i) {
			AFNIStimTimes(stimlabels[i], stimfiles[i], x@hrf, split.ons[[stimlabels[[i]]]], iresp, tr.times)
		})
	
		#labels <- unlist(lapply(1:length(stimlabels), function(i) paste(i+k-1, stimlabels[i], collapse=" ")))
		#stimes <- unlist(lapply(1:length(stimlabels), function(i) paste(i+k-1, stimfiles[i], as.character(x@hrf), collapse=" ")))
		#if (iresp) {
		#	iresp <- unlist(lapply(1:length(stimlabels), function(i) paste(i+k-1, paste(paste(stimlabels[i], "_iresp", sep=""), collapse=" "))))
		#	list(stim_label=labels, stim_times=stimes, stim_file=NULL, iresp=iresp)
		#} else {
		#	list(stim_label=labels, stim_times=stimes, stim_file=NULL, iresp=NULL)
		#}
				
	})


setMethod("buildGLT", signature(con="Contrast"),
	function(con) {
	 		
		gltlabels <- row.names(con) #if (nbasis(term) > 1) { 
			## this may be necessary for AFNI terms, but screws things up otherwise
			#gsub("(^.*):basis([0-9]+)$","\\1[\\2]", shortnames(term, exclude.basis=FALSE))
		#} else {
			#shortnames(term, exclude.basis=FALSE)
		#}
		
		glts <- lapply(1:NCOL(con), function(i) {
			weights <-zapsmall(con[,i], digits=4)	
				
			paste(paste(round(con[weights != 0,i],digits=3), "*", gltlabels[weights!=0], sep=""), collapse=" ")
		})
		
		names(glts) <- colnames(con)
		glts
	})

setMethod("buildAFNIStims", signature(x="RegressionTerm"),
	function(x,...) {
		stimlabels <- longnames(x)
		stimfiles <- paste(stimlabels, "_reg.1D", sep="")
		desmat <- convolve(x)
		
		lapply(1:length(stimlabels), function(i) {
			AFNIStimFile(stimlabels[i], stimfiles[i], desmat[,i])
		})	
		
	})
	
setMethod("buildCommandSwitch", signature(x="AFNIStimTimes", k="numeric", type="character"),
	function(x, k, type) {
		switch(type,
			label=paste(k, x@label, collapse=" "),
			times=paste(k, x@fileName, as.character(x@hrf), collapse=" "),
			file=NULL,
			iresp=if (x@iresp) paste(k, paste(paste(x@label, "_iresp", sep=""), collapse=" ")) else NULL
		)
	})

setMethod("buildCommandSwitch", signature(x="AFNIStimFile", k="numeric", type="character"),
	function(x, k, type) {
		switch(type,
			label=paste(k, x@label, collapse=" "),
			file=paste(k, x@fileName, collapse=" "),
			times=NULL,
			iresp=NULL
		)
	})


	

AFNIStimFile <- function(label, fileName, values) {
  new("AFNIStimFile", label=label, fileName=fileName, values=values)
}

AFNIStimTimes <- function(label, fileName, hrf, onsets, iresp=FALSE, tr.times=1) {
  new("AFNIStimTimes", label=label, fileName=fileName, hrf=hrf, onsets=onsets, iresp=iresp, tr.times=tr.times)
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

setMethod("writeAFNIStim", signature(x="AFNIStimFile", dir="character"),
          function(x, dir) {
            
            .writeOnsets <- function(outname, vals) {
              hfile = file(outname, "w")
              write(x@values, file=hfile, ncolumns=1)
              close(hfile)
            }

            .writeOnsets(paste(dir, "/", x@fileName, sep=""), x@values)
          })


setMethod("writeAFNIStim", signature(x="AFNIStimTimes", dir="character"),
          function(x, dir) {
            
            .writeOnsets <- function(outname, onsets) {
              hfile = file(outname, "w")
              write(onsets, file=hfile, ncolumns=1)
              close(hfile)
            }

            .writeOnsets(paste(dir, "/", x@fileName, sep=""), x@onsets)
          })
                   
            
         
           

