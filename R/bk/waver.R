
library(gtools)
library(numDeriv)


create.HRF <- function(HRF, ...) {
  .orig <- list(...)

  if (length(.orig) > 0) {
    ret <- function(t) {
      do.call(HRF, c(list(t), .orig))
    }
    attr(ret, "params") <- .orig
    ret
  } else {
     function(t) {     
      HRF(t)
    }
  }
}


HRF.TIME <- function(t, maxt) {
  ifelse(t > 0 & t < maxt, t, 0)
}

attr(HRF.TIME, "name") <- "time.hrf"

.HRF.GAMMA <- function(t, shape, rate) {
  dgamma(t, shape=shape, rate=rate)
}

HRF.GAMMA <- create.HRF(.HRF.GAMMA, shape=6, rate=1)


 
.HRF.GAUSSIAN <- function(t, mean=6, sd=2) {
  dnorm(t, mean=mean, sd=sd)
}

HRF.GAUSSIAN <- create.HRF(.HRF.GAUSSIAN, mean=6, sd=2)




HRF.BSPLINE <- function(t, width=20, N=5, degree=3) {
  
  ord <- 1 + degree
  nIknots <- N - ord + 1
  if (nIknots < 0) {
    nIknots <- 0
    #warning("'df' was too small; have used  ", ord - (1 - intercept))
  }
  

  
  knots <- if (nIknots > 0) {
    knots <- seq.int(from = 0, to = 1, length.out = nIknots + 2)[-c(1, nIknots + 2)]
    stats::quantile(0:width, knots)
  } else {
    0
  }

 if (any(t < 0)) {
    t[t < 0] <- 0
  }

  if(any(t > width)) {
    t[t > width] <- 0
  }
 
  bs(t, df=N, knots=knots, degree=degree, Boundary.knots=c(0,width))
}

attr(HRF.GAUSSIAN, "name") <- "gaussian.hrf"
attr(HRF.GAMMA, "name") <- "gamma.hrf"
attr(HRF.TIME, "name") <- "time.hrf"
attr(HRF.BSPLINE, "name") <- "bspline.hrf"

HRF.LOGIT <- function(t, a1=1, T1, T2, T3, D1, D2, D3) {

  a2 <- a1 * (((inv.logit(-T3)/D3) - (inv.logit(-T1)/D1)) / ((inv.logit(-T3)/D3) + (inv.logit(-T2/D2))))
 
  a3 <- abs(a2) - abs(a1)
  
  a1*inv.logit((t-T1)/D1) + a2*inv.logit((t-T2)/D2) + a3*inv.logit((t-T3)/D3)
}

HRF.LOGIT2 <- function(t, a1, a2, a3, T1, T2, T3, D1, D2, D3) {
  
  a1*inv.logit((t-T1)/D1) + a2*inv.logit((t-T2)/D2) + a3*inv.logit((t-T3)/D3)
  
}


getfun <- function(time) {

  TIME <- time
  ret <- function(par) {
    HRF.LOGIT2(TIME, par[1], par[2], par[3], par[4], par[5], par[6], par[7], par[8], par[9])
  }

  return(ret)
}

minimize.fun <- function(yvals, fun, par) {
  ypred <- fun(par)
  return(sum((yvals-ypred)^2))
}

get.minimizer <- function(yfun, yvals) {
  YFUN <- yfun
  YVALS <- yvals
  ret <- function(par) {
    #browser()
    ypred <- YFUN(par)
    return(sum((YVALS-ypred)^2))
  }

  return(ret)
}

    

create.HRF.SET <- function(hrflist) {
  function(t) {
    do.call("cbind", lapply(hrflist, function(fun) fun(t)))
  }
}

shift.HRF <- function(HRF, shift) {
  localShift <- shift
  function(t) {
    HRF(t+localShift)
  }
}

makeDeriv <- function(HRF, n=1) {
  if (n == 1) {
    function(t) grad(HRF, t)
  } else {
    Recall(function(t) grad(HRF,t), n-1)
  }
}

makeBlock <- function(HRF, duration) {
  d1 <- duration
  if (duration < 2) {
    stop("duration must be greater than 1")
  }

  funlist <- c(HRF, lapply(seq(-1, -(duration-1)), function(i) shift.HRF(HRF, i)))
  function(t) {
    ret <- numeric(length(t))
    for (i in 1:length(t)) {
      ret[i] <- sum(unlist(lapply(funlist, function(fun) fun(t[i]))))
    }
    ret
    
  }
}


makeBlockHRF <- function(eventOnset, duration, HRF) {
  
  localHRF <- HRF
  localOnset <- eventOnset
  onsets <- seq(eventOnset, eventOnset+duration, 1)
  funlist <- lapply(onsets, function(onset) makeEventHRF(onset, localHRF))
  
  function(t) {
    ret <- lapply(funlist, function(fun) fun(t)) 
    Reduce("+", ret)
  }
}    
  
  
makeEventHRF <- function(eventOnset, HRF, amp=1) {

  localHRF <- HRF
  localOnset <- eventOnset
  localAmp <- amp
  function(t) {
    #browser()
    
    localHRF(t-localOnset)*amp
    #for (i in 1:length(t)) {
    #  ret[i] <- localHRF(t[i]-localOnset)*amp
    #}
 
  }
}

.makeEventHRF2 <- function(eventOnset, HRF, amp=1) {

  localHRF <- HRF
  localOnset <- eventOnset
  localAmp <- amp
  function(t) {
    ret <- numeric(length(t))
    for (i in 1:length(t)) {
      if (t[i] < localOnset) {
        ret[i] <- 0
      } else {
        ret[i] <- localHRF(t[i]-localOnset)*amp
      }
    }
    ret
  }
}

.makeEventHRF3 <- function(eventOnset, HRF, amp=1) {

  localHRF <- HRF
  localOnset <- eventOnset
  localAmp <- amp
  function(t) {
    if (t < localOnset) {
      0
    } else {
      localHRF(t-localOnset)*amp
    }

  }
}


    
      



    






    
      
      
