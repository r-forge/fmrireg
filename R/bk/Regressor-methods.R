


  



makeEventHRF <- function(eventOnset, HRF, amp=1) {

  localHRF <- HRF
  localOnset <- eventOnset
  localAmp <- amp

  if (is.na(amp)) {
    ## this is wrong for basis functions
    function(t) { rep(0, length(t)) }
  } else {
    localHRF(t-localOnset)*amp   
  }
}




makeBlockHRF <- function(eventOnset, HRF, amp=1, duration, granularity=.1) {
  if (duration == 0) {
    return(makeEventHRF(eventOnset, HRF, amp))
  }
           
  localHRF <- HRF
  localOnset <- eventOnset
  onsets <- seq(eventOnset, eventOnset+duration, granularity)
  funlist <- lapply(onsets, function(onset) makeEventHRF(onset, localHRF, amp))
  
  f <- function(t) {
    ret <- lapply(funlist, function(fun) fun(t))   
    Reduce("+", ret)
  }

  
}

EventHRF <- function(hrf, onset, amplitude=1, duration=0, granularity=.1, ...) {
  stopifnot(is.function(hrf))
  stopifnot(duration >= 0)

 
  
  if (!inherits(hrf, "HRF")) {
    hrf <- HRF(hrf, ...)
  }

  erf <- NULL
  if (duration < granularity) {
    ## single impulese event
    ehrf <- function(t) {
      hrf(t-onset)*amplitude
    }
  } else {

    ### block event
    .onsets <- seq(onset, onset+duration, granularity)
     funlist <- lapply(.onsets, function(ons) {
      ons <- ons
      function(t) {
         hrf(t-ons)*amplitude
      }
    })

    ehrf <- function(t) {
      ret <- lapply(funlist, function(fun) fun(t))
      Reduce("+", ret)
    }
  }
    
    
  new("EventHRF", ehrf, name=hrf@name, nbasis=nbasis(hrf), onset=onset, amplitude=amplitude, duration=duration)
     
    
}
  

Regressor <- function(cells, hrf, onsets, amp=1, durations=0, granularity=.1) {
                    
  if (length(amp) == 1) {
    .amp <- rep(amp, length(onsets))
  } else {
    .amp <- amp
  }

 

  
  stopifnot(length(.amp) == length(onsets))
  stopifnot(all(durations >= 0))
  
  keep <- .amp != 0
  nzonsets <- onsets[keep]
  nzamp <- .amp[keep]

  funlist <- NULL

  eventMaker <- function(onsets, hrf, amp, durations) {
    lapply(1:length(onsets), function(i) {
      EventHRF(hrf, onsets[i], amp[i], durations[i])
    })
  }
    

  funlist <- if (length(durations == 1) && durations==0) {
    nzdurations <- rep(0, length(nzonsets))
    eventMaker(nzonsets, hrf, nzamp, nzdurations)
  } else if (length(durations) == 1 && durations != 0) {
    nzdurations <- rep(durations, length(nzonsets))
    eventMaker(nzonsets, hrf, nzamp, nzdurations)
  } else {  
    stopifnot(length(durations) == length(onsets))
    nzdurations <- durations[keep]
    eventMaker(nzonsets, hrf, nzamp, nzdurations)
  }
             
                 
  fun <- function(t) {
    ret <- lapply(funlist, function(fun) fun(t))
    zapsmall(Reduce("+", ret))
  }

  if (nbasis(hrf) > 1) {
    cells$basis <- seq(1, nbasis(hrf))
  }

 
  new("Regressor", fun, cells=cells, onsets=nzonsets, hrf=hrf, amplitudes=nzamp, durations=nzdurations)
}

#setMethod("initialize", "Regressor", function(.Object, fun, labels, onsets, hrf, amplitudes, durations=0) {
#  .Object@fun <- fun
#  .Object@labels <- labels
#  .Object@onsets <- onsets
#  .Object@hrf <- hrf
#  .Object@amplitudes <- amplitudes
#  .Object@durations <- durations
#  .Object
#}


setMethod("onsets", signature(x="Regressor"),
          function(x) {
            x@onsets
          })

setMethod("durations", signature(x="Regressor"),
          function(x) {
            x@durations
          })

setMethod("amplitudes", signature(x="Regressor"),
          function(x) {
            x@amplitudes
          })

setMethod("cells", signature(x="Regressor"),
          function(x) {
            x@cells
          })

setMethod("names", signature(x="Regressor"),
          function(x) {
            fnames <- names(x@cells)

                       
            egrid <- expand.grid(x@cells)

            ret <- apply(egrid, 1, function(row) {
              levs <- ifelse(row == "", "", sapply(row, function(el) paste("[", el, "]", sep="")))           
              paste(fnames, levs, sep="")
            })

            if (is.matrix(ret)) {
              apply(ret, 2, paste, collapse=":")
            } else {
              ret
            }
             
          })

setMethod("shortnames", signature(x="Regressor"),
          function(x) {
            fnames <- names(x@cells)
            egrid <- expand.grid(x@cells)
            
            ret <- apply(egrid, 1, function(row) {
              levs <- ifelse(row == "", names(row), sapply(row, function(el) el))           
            })

            if (is.matrix(ret)) {
              apply(ret, 2, paste, collapse=":")
            } else {
              ret
            }
                        
          })
            

setMethod("nbasis", signature(x="Regressor"),
          function(x) {
            nbasis(x@hrf)
          })

#setMethod("names", signature(x="RegressorSet"),
#          function(x) {
#            x@labels
#          })
    
setMethod("show", signature(object="Regressor"),
          function(object) {
            cat("Class: ", class(object))
            cat("\n")
            cat("name: ", names(object))
            cat("\n")

           
            N <- min(c(6, length(onsets(object))))
            cat(paste("onsets: ", paste(onsets(object)[1:N], collapse=" "), "..."))
            cat("\n")

            

            durs <- durations(object)
            amps <- amplitudes(object)
           
            if (all(durs == durs[1])) {
              cat(paste("durations: ", durs[1], "for all events"))
            } else {
              cat(paste("durations: ", paste(durs[1:N], collapse=" "), "..."))
            }

            cat("\n")

            

            if (all(amps == amps[1])) {
              cat(paste("amplitudes: ", amps[1], "for all events"))
            } else {     
              cat(paste("amplitudes: ", paste(amps[1:N], collapse=" "), "..."))
             
            }

            cat("\n")

            
                        
              
          })
            
