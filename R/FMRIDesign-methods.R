



ConvolveModel <- function(eventModel, hrf, TR, blocklens, granularity=.33, normalize.height=FALSE, drop.unused.levels=TRUE, use.multicore=TRUE) {
  stopifnot(numBlocks(eventModel) == length(blocklens))
  if (!inherits(hrf, "HRF") && is.function(hrf)) {
     hrf <- HRF(hrf)
   }

  termlist <- lapply(eventModel@eventTerms, ConvolvedTerm, hrf, blocklens, TR, start=TR/2, 
	granularity=granularity, normalize.height=normalize.height, drop.unused.levels=drop.unused.levels, use.multicore=use.multicore)
	
  new("ConvolvedModel", eventModel=eventModel, functionalTerms=termlist)
   
}
  
  #FMRIDesign(eventModel, dataSource,termlist,hrf)


AFNIModel <- function(eventModel, afni_hrf, dataSource) {
	stopifnot(length(filelist(dataSource)) == length(unique(blockids(eventModel))))
    stopifnot(is(afni_hrf, "AFNI_HRF"))

	termlist <- lapply(eventModel@eventTerms, AFNITerm, afni_hrf, blocklens(dataSource), TR(dataSource))
    afmod <- new("AFNIModel", eventModel=eventModel, functionalTerms=termlist)
    afmod
    	
}


FMRIDesign <- function(fmriModel, dataSource) {
  stopifnot(length(filelist(dataSource)) == length(unique(blockids(fmriModel))))
  stopifnot(length(conditions(fmriModel)) == length(unique(conditions(fmriModel))))
  
  new("FMRIDesign", fmriModel=fmriModel, dataSource=dataSource)
 
}

.convolveRegressors <- function(eventModel, dataSource, hrf, granularity=.1, normalize.height=FALSE) {

   bl <- blocklens(dataSource)
   ons <- onsets(eventModel)
   ids <- blockids(eventModel)

   ids <- rep(1:length(unique(ids)), table(ids))
   tr <- TR(dataSource)
   
   globons <- sapply(1:length(ons),function(i) {
     ons[i] + ((ids[i]-1) * bl[ids[i]]) * tr
   })

   durs <- durations(eventModel)
   
   cell.list <- cells(eventModel)
   dmat <- designMatrix(eventModel)

   reglist <- lapply(1:NCOL(dmat), function(i) {
     amp <- zapsmall(dmat[,i])
     Regressor(hrf, globons, cell.list[[i]], amp, durs, granularity=granularity, normalize.height=normalize.height)
   })

   reglist
 }
   
 
setMethod("filelist", signature(x = "FMRIDesign", full.names="logical"),
    function(x, full.names) {
      filelist(x@dataSource, full.names)
    })

setMethod("filelist", signature(x = "FMRIDesign", full.names="missing"),
    function(x) {
      filelist(x@dataSource)
    })



setMethod("TR", signature(x = "FMRIDesign"),
    function(x) {
      TR(x@dataSource)
    })

setMethod("path", signature(x = "FMRIDesign"),
    function(x) {
      path(x@dataSource)
    })

setMethod("formula",  signature(x = "FMRIModel"),
          function(x) {
            formula(x@eventModel)
         })

setMethod("formula",  signature(x = "FMRIDesign"), function(x) formula(x@fmriModel) )

setMethod("terms", signature(x="FMRIModel"), function(x) x@functionalTerms )
setMethod("terms", signature(x="FMRIDesign"), function(x) terms(x@fmriModel) )

setMethod("eventVars", signature(x="FMRIModel"), function(x) names(eventTable(x@eventModel)) )
setMethod("eventVars", signature(x="FMRIDesign"), function(x) eventVars(x@fmriModel))

setMethod("onsets", signature(x="FMRIModel"), function(x) onsets(x@eventModel) )
setMethod("onsets", signature(x="FMRIDesign"), function(x) onsets(x@fmriModel) )

setMethod("blocklens", signature(x = "FMRIDesign"), function(x) blocklens(x@dataSource))

setMethod("blockids", signature(x = "FMRIModel"), function(x) blockids(x@eventModel))
setMethod("blockids", signature(x = "FMRIDesign"), function(x) blockids(x@fmriModel))

setMethod("conditions", signature(x = "FMRIModel"),
          function(x) {
            unlist(lapply(terms(x), function(cterm) {
              conditions(cterm)
            }))           
          })

setMethod("conditions", signature(x = "FMRIDesign"), function(x) conditions(x@fmriModel))

setMethod("designMatrix",  signature(x = "FMRIModel"),
          function(x) {
            mat <- do.call(cbind, lapply(x@functionalTerms, designMatrix))            
          })

setMethod("designMatrix",  signature(x = "FMRIDesign"), function(x) designMatrix(x@fmriModel))


setMethod("colnames", signature(x="FMRIModel", do.NULL="missing", prefix="missing"),
          function(x) {  
	        conditions(x)
          })

setMethod("colnames", signature(x="FMRIDesign", do.NULL="missing", prefix="missing"), function(x) conditions(x@fmriModel))


## need to distinguish betweens "continuous" and "factor" terms
setMethod("cells", signature(x="FMRIModel"),
          function(x) {
            cell.list <- lapply(terms(x), cells)
            names(cell.list) <- terms(x)
            cell.list
            
          })

setMethod("cells", signature(x="FMRIDesign"), function(x) cells(x@fmriModel))


setMethod("shortnames", signature(x="FMRIModel"),
          function(x) {
            unlist(lapply(terms(x), function(cterm) {
              row.names(cells(cterm))
			}))                                      
          })
  
setMethod("longnames", signature(x="FMRIModel"),
		  function(x) {
			  unlist(lapply(terms(x), longnames))                               
		  })

setMethod("shortnames", signature(x="FMRIDesign"),function(x) shortnames(x@fmriModel))
setMethod("longnames", signature(x="FMRIDesign"),function(x) longnames(x@fmriModel))

      
setMethod("globalOnsets", signature(x= "FMRIModel",blockDurations="missing"),
          function(x) {
            stop("not implemented")          
 
         })

setMethod("globalOnsets", signature(x= "FMRIDesign",blockDurations="missing"), function(x) stop("not implemented"))

setMethod("[[", signature(x="FMRIModel", i = "character", j = "missing"),
          function(x, i, j) {
            x@eventModel[[i]]
          })


setMethod("[[", signature(x="FMRIDesign", i = "character", j = "missing"), function(x,i,j)  { x@fmriModel[[i]] })

setMethod("[[", signature(x="FMRIModel", i = "numeric", j = "missing"),
          function(x, i, j) {
            x@eventModel[[i]]
          })

setMethod("[[", signature(x="FMRIDesign", i = "numeric", j = "missing"), function(x,i,j)  { x@fmriModel[[i]] })

setMethod("show", signature(object="ConvolvedModel"),
          function(object) {
            cat("formula: ", as.character(formula(object@eventModel)), "\n")
            cat("columns: ", colnames(object@eventModel), "\n")
                             
          })


setMethod("show", signature(object="FMRIModel"),
          function(object) {
            cat("formula: ", as.character(formula(object@eventModel)), "\n")
            cat("columns: ", colnames(object@eventModel), "\n")                           
          })

setMethod("show", signature(object="FMRIDesign"),
          function(object) {
            cat("formula: ", as.character(formula(object@fmriModel@eventModel)), "\n")
            cat("columns: ", colnames(object@fmriModel), "\n")
            cat("num blocks: ", length(blocklens(object)), "\n")
            cat("block lengths: ", blocklens(object), "\n")                    
          })
