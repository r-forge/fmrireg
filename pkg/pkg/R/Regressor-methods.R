


  
EventHRF <- function(hrf, onset, amplitude=1, duration=0, span=20, granularity=.1, height=NULL, ...) {
  stopifnot(is.function(hrf))
  stopifnot(duration >= 0)

        
  if (!inherits(hrf, "HRF")) {
    hrf <- HRF(hrf, ...)
  }

  erf <- NULL
  if (duration < granularity) {
	  
	### FAILS TO SCALE TO MAX_HEIGHT
	  
    ## single impulse event
    ehrf <- function(t) {
      hrf(t-onset)*amplitude
    }
	
	### FAILS TO SCALE TO MAX_HEIGHT
      
  } else {

    ### block event
    .onsets <- seq(onset, onset+duration, granularity)
    funlist <- lapply(.onsets, function(ons) {
      ons <- ons
      function(t) {
        hrf(t-ons)*amplitude
      }
    })

	if (!is.null(height)) {
		ret <- lapply(funlist, function(fun) fun(seq(.onsets[1], .onsets[1]+span, by=granularity)))
		vals <- if (nbasis(hrf) > 1) {
					Reduce("+", ret)
				} else {
					zapsmall(rowSums(do.call(cbind, ret)))
				}
		MAX_HEIGHT <- max(vals)
		#print(paste("MAX_HEIGHT", MAX_HEIGHT))
				
	}

    ehrf <- function(t) {
		ret <- lapply(funlist, function(fun) fun(t))		
	  
      vals <- if (nbasis(hrf) > 1) {
        Reduce("+", ret)
      } else {
        zapsmall(rowSums(do.call(cbind, ret)))
      }
	  
	  if (!is.null(height)) {  
	
		  sf <- if (MAX_HEIGHT > 0) {
		  	height/MAX_HEIGHT
		  } else {
			  1
		  }
  		  sf*vals
	  } else {
		  vals
	  }
      
      #if (normalize.height) {
      #  sf <- abs(diff(range(vals)))
      #  if (sf == 0) {
      #    vals
      #  } else {
      #    vals/max(vals)
      #  }
      #} else {
      #  vals
      #}
    }
  }
    
   
  new("EventHRF", ehrf, name=hrf@name, nbasis=nbasis(hrf), onset=onset, amplitude=amplitude, duration=duration, span=span+duration)
         
}

 
  
makeEventList <- function(onsets, hrf, amp, durations, granularity=.1,span=20, height=NULL) {

  sapply(1:length(onsets), function(i) {
      EventHRF(hrf, onsets[i], amp[i], durations[i], span, granularity=granularity, height=height)
    })
}

convolveBlockEvent <- function(hrf, onset, amp, duration, granularity) {
	stopifnot(duration > granularity)
	.onsets <- seq(onset, onset+duration, granularity)
    funlist <- lapply(.onsets, function(ons) {
      ons <- ons
      function(t) {
        hrf(t-ons)*amp
      }
    })

	function(t) {
     ret <- lapply(funlist, function(fun) fun(t))
      if (nbasis(hrf) > 1) {
        vals <- Reduce("+", ret)
      } else {
        mat <- do.call("cbind", ret)
        vals <- zapsmall(rowSums(mat))
      }
	}	
}

convolveBlockEvent2 <- function(hrf, onset, amp, duration, granularity) {
	stopifnot(duration > granularity)
	.onsets <- seq(onset, onset+duration, granularity)
    funlist <- lapply(.onsets, function(ons) {
	  ons <- ons
      function(t) {
        ifelse((ons > t) | (t > (ons + duration + 20)), 0, hrf(t-ons)*amp)
      }
    })

	function(t) {
	
		ret <- lapply(funlist, function(fun) fun(t))
	
      	if (nbasis(hrf) > 1) {
        	Reduce("+", ret)
      	} else {
        	mat <- do.call(cbind, ret)
        	zapsmall(rowSums(mat))
      	}
			
	}
}

convolveSimpleEvent <- function(hrf, onset, amp) {
	function(t) {
        hrf(t-onset)*amp
    }
}

createFun <- function(funlist) {
	function(t) {	
		if (require(multicore)) {
			ret <- mclapply(funlist, function(fun) fun(t))
		} else {
			ret <- lapply(funlist, function(fun) fun(t))
		}
		zapsmall(rowSums(do.call(cbind, ret)))
    
	}
}

createFun1 <- function(funlist) {
	
	function(t) {
		ons <- sapply(funlist, function(f) f@onset)	
		offset <- sapply(funlist, function(f) f@span) + ons	
		
		
				
		ret <- mclapply(seq_along(t), function(i) {
			idx <- which(t[i] >= ons & t[i] <= offset)
			vals <- numeric(length(t))
			if (length(idx) > 0) {
				vals[i] <- sum(sapply(funlist[idx], function(f) f(t[i])))
			}
			vals
			
		})
		
		rowSums(do.call(cbind, ret))	
		
	}
}

createFun2 <- function(funlist) {
	function(t) {
		sapply(t, function(x) {
			sum(sapply(funlist, function(f) {
				if (x >= f@onset && x < (f@onset + f@span)) {
					f(x)
				} else {
					0
				}
			}
		))
	})
}}
  
Regressor <- function(hrf, onsets, amp=1, durations=0, granularity=.1, span=20, height=NULL, splitBy=NULL) {
  if (!inherits(hrf, "HRF")) {
    hrf <- HRF(hrf)
  }
 
  .amp <- if (length(amp) == 1) {
    rep(amp, length(onsets))
  } else {
    amp
  }

  
  stopifnot(length(.amp) == length(onsets))
  stopifnot(all(durations >= 0))
  
  keep <- !is.na(.amp)
  nzonsets <- onsets[keep]
  nzamp <- .amp[keep]

  nzdurations <- if (length(durations) == 1) {
    rep(durations, length(nzonsets))
  } else {
    stopifnot(length(durations) == length(onsets))
    durations[keep]
  }
  
  funlist <- makeEventList(nzonsets, hrf, nzamp, nzdurations, granularity, span, height)
  

  
  fun <- function(t) {	
		ret <- lapply(funlist, function(fun) {				
				fun(t)
		})
				
		zapsmall(Reduce("+", ret))
		
	}
		
  new("Regressor", fun, onsets=nzonsets, hrf=hrf, amplitudes=nzamp, durations=nzdurations)
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


#setMethod("names", signature(x="Regressor"),
#          function(x) {
#            fnames <- names(x@cells)                   
#            egrid <- expand.grid(x@cells)
#
#            ret <- apply(egrid, 1, function(row) {
#              levs <- ifelse(row == "", "", sapply(row, function(el) paste("[", el, "]", sep="")))           
#              paste(fnames, levs, sep="")
#            })
#
#            if (is.matrix(ret)) {
#              apply(ret, 2, paste, collapse=":")
#            } else {
#              ret
#            }
#             
#          })

#setMethod("shortnames", signature(x="Regressor"),
#          function(x) {
#            fnames <- names(x@cells)
#            egrid <- expand.grid(x@cells)
#            
#            ret <- apply(egrid, 1, function(row) {
#              levs <- ifelse(row == "", names(row), sapply(row, function(el) el))           
#            })
#
#            ret <- if (is.matrix(ret)) {
#               apply(ret, 2, function(sname) {
#                oname <- paste(sname, collapse=":")
#                gsub("[\\(\\)]", "_", oname)
#              })
#            } else {
#              gsub("[\\(\\)]", "_", ret)            
#            }
#
#            ret
#                        
#          })
            

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
            
