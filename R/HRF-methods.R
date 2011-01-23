



AFNI_SPMG1 <- function(d=1) { new("AFNI_HRF", name="SPMG1", nbasis=as.integer(1), params=list(d=d)) }
AFNI_SPMG2 <- function(d=1) { new("AFNI_HRF", name="SPMG2", nbasis=as.integer(2), params=list(d=d)) }
AFNI_SPMG3 <- function(d=1) { new("AFNI_HRF", name="SPMG3", nbasis=as.integer(3), params=list(d=d)) }
AFNI_BLOCK <- function(d=1,p=1) { new("AFNI_HRF", name="BLOCK", nbasis=as.integer(1), params=list(d=d,p=p)) } 
AFNI_TENT <- function(b=0,c=18, n=10) { new("AFNI_HRF", name="TENT", nbasis=as.integer(n), params=list(b=b,c=c,n=n)) } 
AFNI_CSPLIN <- function(b=0,c=18, n=6) { new("AFNI_HRF", name="CSPLIN", nbasis=as.integer(n), params=list(b=b,c=c,n=n)) }
AFNI_POLY <- function(b=0,c=18, n=10) { new("AFNI_HRF", name="POLY", nbasis=as.integer(n), params=list(b=b,c=c,n=n)) }
AFNI_SIN <- function(b=0,c=18, n=10) { new("AFNI_HRF", name="SIN", nbasis=as.integer(n), params=list(b=b,c=c,n=n)) }
AFNI_GAM <- function(p=8.6,q=.547) { new("AFNI_HRF", name="GAM", nbasis=as.integer(1), params=list(p=p,q=q)) }
AFNI_WAV <- function(d=1) { new("AFNI_HRF", name="WAV", nbasis=as.integer(1), params=list(d=1)) }


get_AFNI_HRF <- function(name) {
	hrf <- switch(name,
			gamma=AFNI_GAM,
			spmg1=AFNI_SPMG1,
			spmg2=AFNI_SPMG2,
			spmg3=AFNI_SPMG3,
			csplin=AFNI_CSPLIN,
			poly=AFNI_POLY,
			sine=AFNI_SIN,
			wav=AFNI_WAV,
			block=AFNI_BLOCK)
	
	if (is.null(hrf)) {
		stop("could not find hrf named: ", name)
	}
	
	hrf
	
}


makeDeriv <- function(HRF, n=1) {
  if (n == 1) {
    function(t) grad(HRF, t)
  } else {
    Recall(function(t) grad(HRF,t), n-1)
  }
}



HRF <- function(hrf, name=NULL, ...) {
  .orig <- list(...)

  if (is.list(hrf)) {
    stopifnot(all(sapply(hrf, class) == "function"))
    if (length(.orig) != 0) { warning(paste("extra args ", .orig, "ignored")) }
    
    fun <- function(t) {
      do.call("cbind", lapply(hrf, function(fun) fun(t)))
    }
  } else if (length(.orig) > 0) {
    fun <- function(t) do.call(hrf, c(list(t), .orig))
  } else {
    fun <- function(t) hrf(t)
  }

  test <- fun(1)
  
  if (is.matrix(test)) {
    nbasis <- ncol(test)
  } else {
    nbasis <- 1
  }

  name <- attr(hrf, "name")
  if (is.null(name)) {
    name = "unknown_hrf"
  }
  
  new("HRF", fun, name=name, nbasis=as.integer(nbasis))
}

setMethod("nbasis", signature(x="HRF"),
	function(x) {
		x@nbasis
    })

setMethod("nbasis", signature(x="AFNI_HRF"),
	function(x) {
		x@nbasis
	})

setMethod("params", signature(x="AFNI_HRF"),
	function(x) {
		x@params
	})
	
setAs("AFNI_HRF", "character", function(from) {
	paste(from@name, "\\(", paste(from@params, collapse=","), "\\)", sep="")
})

setMethod("as.character", signature(x="AFNI_HRF"),
	function(x) {
		as(x, "character")
	})


setMethod("names", signature(x="HRF"),
	function(x) {
		x@name
	})

setMethod("names", signature(x="AFNI_HRF"),
	function(x) {
		x@name
	})


setMethod("show", signature(object="AFNI_HRF"),
	function(object) {
			cat(object@name, "\n")
			cat("nbasis ", object@nbasis, "\n")
			cat("params: ", "\n")
			print(object@params)
	})


HRFDerivBasis <- function(hrf, nderiv=2) {
  if (nderiv == 0) {
    return(hrf)
  }

  hrflist <- list()
  hrflist[[1]] <- hrf
  for (i in 1:nderiv) {
    hrflist[[i+1]] <- makeDeriv(hrf, i)
  }

  HRF(hrflist)
}
    
  
