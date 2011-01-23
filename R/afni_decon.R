

afni.glm <- function(fmridesign, workingDir=".", mask, options=list(), gltlist=list()) {
  
  defopts <- list(noFDR=FALSE, fout=TRUE, rout=TRUE, tout=TRUE, float=TRUE, cbucket="coefout", bucket="statout", jobs=2, polort=5, iresp=FALSE, TR_times=1)
  
  for (optname in names(options)) {
    defopts[[optname]] <- options[[optname]]
  }

  new("AFNICommand", design=fmridesign, workingDir=workingDir, mask=mask, options=defopts, glts=gltlist)
}

.makeCommandStr <- function(cmdlines) {
  cmdstr <- lapply(names(cmdlines), function(optname) {
    entry <- cmdlines[[optname]]
    if (is.list(entry)) {
      switchnames <- rep(paste("-", optname, sep=""), length(entry))
      paste(switchnames, entry, collapse=" ")
    } else if (is.null(entry) || length(entry) == 0) {
      ""
    } else if (is.numeric(entry[[1]])) {
		paste(paste("-", optname, sep=""), entry[[1]])
	} else if (entry[[1]] == TRUE) {
      paste("-", optname, sep="")
    } else if (entry[[1]] == FALSE) {
      paste("")
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
            stimlabels <- longnames(x@design)
            
            stopifnot(length(unique(stimlabels)) == length(stimlabels))
            stopifnot(length(stimlabels) == length(conditions(x@design)))
             
            #stimfiles <- paste(stimlabels, "_reg.1D", sep="")
			
            gltlist <- glts(x)
            gltnames <- if (length(gltlist) > 0) names(gltlist) else NULL
            gltfiles <- if (length(gltnames) > 0) paste("glt_", gltnames, ".txt", sep="") else NULL
			
			funcTerms <- terms(x@design)
			#browser()
			afni.stims <- unlist(lapply(funcTerms, function(term) { buildAFNIStims(term, opts$iresp, opts$TR_times ) }))
			
			purgeNulls <- function(A) {
				A[!sapply(A, is.null)]				
			}
			
			
			opt_stim_labels <-  purgeNulls(lapply(seq_along(afni.stims), function(i) buildCommandSwitch(afni.stims[[i]], i, "label")))
			opt_stim_files  <-  purgeNulls(lapply(seq_along(afni.stims), function(i) buildCommandSwitch(afni.stims[[i]], i, "file")))
			opt_stim_times  <-  purgeNulls(lapply(seq_along(afni.stims), function(i) buildCommandSwitch(afni.stims[[i]], i, "times")))			
			opt_stim_iresp  <-  purgeNulls(lapply(seq_along(afni.stims), function(i) buildCommandSwitch(afni.stims[[i]], i, "iresp")))
			
			
            cmdlines <- list(input=filelist(x@design, full.names=T),
                             mask=x@mask,
                             polort=opts[["polort"]],
                             num_stimts=length(afni.stims),
							 num_glt=length(gltlist),
                             stim_file=opt_stim_files,
                             stim_label=opt_stim_labels,
							 stim_times=opt_stim_times,
							 TR_times=opts[["TR_times"]],
							 iresp=opt_stim_iresp,
                             gltsym=lapply(seq_along(gltfiles), function(i) paste(gltfiles[i], collapse=" ")),
                             glt_label=lapply(seq_along(gltnames), function(i) paste(i, gltnames[i], collapse=" ")),						 
                             nofullf_atall=opts[["nofullf_atall"]],
                             fout=opts[["fout"]],
                             rout=opts[["rout"]],
                             tout=opts[["tout"]],
                             bout=opts[["bout"]],
							 noFDR=opts[["noFDR"]],
                             cbucket=opts[["cbucket"]],
                             bucket=opts[["bucket"]],
                             jobs=opts[["jobs"]],
                             float=TRUE)

         
            cmdstr <- .makeCommandStr(cmdlines)
             
            ret <- list()
            wd <- workingDir(x)

            nextDirName <- function(wd) {
              nd <- paste(wd, "+", sep="")
              if (!file.exists(nd)) {
                nd
              } else {
                Recall(nd)
              }
            }

            writeStimFiles <- function() {
              sapply(afni.stims, function(stim) {
                  writeAFNIStim(stim, ".")
                })
            }

            writeGLTs <- function() {
              lapply(seq_along(gltlist), function(i) {
                fout <- file(gltfiles[i], "w")
                .glt <- gltlist[[i]]
				
				write(unlist(.glt), file=fout, sep="\n")
                
                close(fout)
              })
            }
              
                          
            ret$run <- function() {
              startDir <- getwd() 

              res <- try({

                if (!file.exists(wd)) {
                  dir.create(wd)
                } else {
                  warning(paste("glm output directory: ", wd, " already exists"))
                  wd <- nextDirName(wd)
                  dir.create(wd)
                  warning(paste("outputting to: ", wd))
                  
                }
                            
                print(paste("setting directory:", wd))
                setwd(wd)
				
                writeStimFiles()
                writeGLTs()
              
                

                write(cmdstr, "3ddeconvolve.sh")             
                system(cmdstr)
              })

              setwd(startDir)

              
            }

            ret$command <- cmdstr
            #ret$designMat <- desmat
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

setMethod("glts", signature(x="AFNICommand"),
          function(x) {
            x@glts
          })




  
  
  

  
