
FMRIDesign <- function(eventModel, dataSource, hrf) {
  stopifnot(length(filelist(dataSource)) == length(unique(blockids(eventModel))))
  new("FMRIDesign", eventModel=eventModel, dataSource=dataSource, hrf=hrf)
}

setMethod("initialize", "FMRIDesign", function(.Object, eventModel, dataSource, hrf) {
  .Object <- callNextMethod()
  .Object@eventModel <- eventModel
  .Object@dataSource <- dataSource

  cnames <- colnames(eventModel)
  dmat <- designMatrix(eventModel)

 

  bl <- blocklens(dataSource)
  ons <- onsets(eventModel)
  ids <- blockids(eventModel)

  
  ## should block ids truly denote the run number?
  ids <- rep(1:length(unique(ids)), table(ids))

  tr <- TR(dataSource)
  globons <- sapply(1:length(ons),function(i) {
    ons[i] + ((ids[i]-1) * bl[ids[i]]) * tr
  })

  durs <- durations(eventModel)


  if (!inherits(hrf, "HRF") && is.function(hrf)) {
    hrf <- HRF(hrf)
  }

   
  cell.list <- cells(eventModel)

  reglist <- lapply(1:NCOL(dmat), function(i) {
    amp <- zapsmall(dmat[,i])
    Regressor(cell.list[[i]], hrf, globons, amp, durs)
  })

  .Object@HRF <- hrf
  .Object@regressors <- reglist

  .Object

})

  


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

setMethod("formula",  signature(x = "FMRIDesign"),
          function(x) {
            formula(x@eventModel)
          })

setMethod("terms", signature(x="FMRIDesign"),
          function(x) {
            ret <- terms(formula(x))
            attr(ret, "term.labels")
          })

setMethod("eventVars", signature(x="FMRIDesign"),
          function(x) {
            names(eventTable(x@eventModel))
          })
          

setMethod("onsets", signature(x="FMRIDesign"),
          function(x) {
            onsets(x@eventModel)
          })

setMethod("blocklens", signature(x = "FMRIDesign"),
    function(x) blocklens(x@dataSource))

setMethod("blockids", signature(x = "FMRIDesign"),
          function(x) blockids(x@eventModel))

setMethod("conditions", signature(x = "FMRIDesign"),
          function(x) conditions(x@eventModel))

setMethod("designMatrix",  signature(x = "FMRIDesign"),
          function(x) {

             samps <- samples(x@dataSource)
            res <- lapply(x@regressors, function(reg) {
              reg(samps)
            })

            regnames <- unlist(lapply(x@regressors, names))
            mat <- do.call("cbind", res)
            colnames(mat) <- regnames

            mat
            
          })

setMethod("colnames", signature(x="FMRIDesign", do.NULL="missing", prefix="missing"),
          function(x) {
            unlist(lapply(x@regressors, names))
                   
          })

setMethod("shortnames", signature(x="FMRIDesign"),
          function(x) {
            unlist(lapply(x@regressors, shortnames))                    
          })
      
setMethod("globalOnsets", signature(x= "FMRIDesign",blocklens="missing"),
          function(x) {
            browser()
            bl <- blocklens(x)
            ons <- onsets(x)
            ids <- blockids(x)

            ## ids must be strictly ordinal             
            globons <- sapply(1:length(ons),function(i) {
              ons[i] + ((ids[i]-1) * bl[ids[i]])
            })

            
          })


setMethod("[[", signature(x="FMRIDesign", i = "character", j = "missing"),
          function(x, i, j) {
            x@eventModel[[i]]
          })

setMethod("[[", signature(x="FMRIDesign", i = "numeric", j = "missing"),
          function(x, i, j) {
            x@eventModel[[i]]
          })
