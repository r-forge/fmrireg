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
    
  
