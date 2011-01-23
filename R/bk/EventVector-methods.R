

EV <- function(vals, name, onsets, blockids=1, durations=NULL) {
  if (is.matrix(vals)) {
    return(EventVariableSet(vals, name, onsets, blockids, durations))
  }

  if (is.factor(vals)) {
    return(EventFactor(vals, name, onsets, blockids, durations))
  }

  if (is.vector(vals)) {
    return(EventVariable(vals, name, onsets, blockids, durations))
  }
}

.checkEVArgs <- function(vals, onsets, blockids, durations) {

  ret <- list(vals=vals, onsets=onsets, blockids=blockids, durations=durations)

  stopifnot(length(onsets) == length(vals))
   
  
  if (is.null(durations) || length(durations) == 1) {
    durations <- rep(0, length(onsets))
  }

  
  stopifnot(length(durations) == length(vals))
   
  if (length(blockids) != length(onsets)) {
    blockids <- rep(blockids[1], length(onsets))
  }

  list(vals=vals, onsets=onsets, blockids=blockids, durations=durations)

}

  
EventFactor <- function(facvar, name, onsets, blockids=1, durations=NULL) {
   if (!is.factor(facvar)) {
    facvar <- as.factor(facvar)
  }

  ret <- .checkEVArgs(facvar, onsets, blockids, durations)
  
  new("EventFactor", value=ret$vals, varname=name, onsets=ret$onsets, durations=ret$durations, blockids=as.integer(ret$blockids))
   
}          

EventVariable <- function(vals, name, onsets, blockids=1, durations=NULL) {

  ret <- .checkEVArgs(vals, onsets, blockids, durations)
  
  new("EventVariable", value=ret$vals, varname=name, onsets=ret$onsets, durations=ret$durations,blockids=as.integer(ret$blockids))
}


EventVariableSet <- function(mat, name, onsets, blockids=1, durations=NULL) {
  stopifnot(is.matrix(mat))

  ret <- .checkEVArgs(mat[,1], onsets, blockids, durations)

  class(mat) <- "matrix"

  if (is.null(colnames(mat))) {
    colnames(mat) <- 1:NCOL(mat)
  }
  
  new("EventVariableSet", value=mat, varname=name, onsets=ret$onsets, durations=ret$durations, blockids=as.integer(ret$blockids))
}






setMethod("[", signature(x="EventVector", i="logical", j="missing", drop="missing"),
          function(x,i,j) {
            EventFactor(elements(x)[i], name=name(x), onsets=onsets(x)[i], durations=durations(x)[i],
                        blockids=blockids(x)[i])
          })

setMethod("[", signature(x="EventVector", i="numeric", j="missing", drop="missing"),
          function(x,i,j) {
            EventFactor(elements(x)[i], name=name(x), onsets=onsets(x)[i], durations=durations(x)[i],
                        blockids=blockids(x)[i])
          })

setMethod("cells", signature(x="EventFactor"),
          function(x) {
            f1 <- elements(x)

            onsplit <- split(onsets(x), f1)
            dursplit <- split(durations(x), f1)
            blocksplit <- split(blockids(x), f1)
            
            cells <- lapply(1:length(onsplit), function(i) {
              dat <- data.frame(onsets=onsplit[[i]], durations=dursplit[[i]], blockids=blocksplit[[i]])
              dat <- dat[order(dat$onsets),]
            })

            cellnames <- names(onsplit)
            names(cells) <- cellnames
            return(cells)
          })

  

.checkCompat <- function(e1, e2) {

  
   
   v1 <- elements(e1)
   v2 <- elements(e2)

   if (is.matrix(v1)) {
     v1 <- v1[,1]
   }

   if (is.matrix(v2)) {
     v2 <- v2[,1]
   }
   
   if (length(v1) != length(v2)) {
     stop("error: event vectors must have same length to be crossed")
   }
   
   
   on1 <- onsets(e1)
   on2 <- onsets(e2)

   if (!all(on1==on2)) {
     stop("error: event vectors must have identical onsets to be crossed")
   }
 }
  

setMethod("*", signature(e1 = "EventFactor", e2="EventFactor"),
  function(e1,e2) {
    .checkCompat(e1,e2)

    if (identical(e1, e2)) {
      stop(paste("cannot cross identical factors", varname(e1), " and ", varname(e2)))
    }
    

    f3 <- interaction(elements(e1),elements(e2), sep=":", drop=T) 
    name <- paste(varname(e1), ":", varname(e2), sep="")

   
    ev <- EventFactor(f3, name, onsets(e1), blockids=blockids(e1), durations=durations(e1))
    
    attr(ev, "parents") <- list(varname(e1), varname(e2))
    ev
  })

setMethod("*", signature(e1 = "EventVariable", e2="EventFactor"),
          function(e1,e2) {
            callGeneric(e2,e1)
          })

setMethod("*", signature(e1 = "EventFactor", e2="EventVariable"),
          function(e1,e2) {
            .checkCompat(e1,e2)

            on1 <- onsets(e1)
            fac <- elements(e1)
            var <- elements(e2)

            omat <- matrix(0, length(on1), length(levels(fac)))
            for (i in 1:length(levels(fac))) {
              lev <- levels(fac)[i]
              idx <- which(fac == lev)
              omat[idx, i] <- var[idx]
            }

            colnames(omat) <- paste(levels(fac), ":", varname(e2), sep="")
                       
            ev <- EV(omat, paste(varname(e1), ":", varname(e2), sep=""),
                     on1, blockids=as.integer(blockids(e1)), durations=as.integer(durations(e1)))
            
            attr(ev, "parents") <- list(varname(e1), varname(e2))
            ev
            
          })


setMethod("*", signature(e1="EventVariableSet", e2="EventVariableSet"),
  function(e1,e2) {
    .checkCompat(e1,e2)

    on1 <- onsets(e1)
    mat1 <- elements(e1)
    mat2 <- elements(e2)

    
    rmat <- model.matrix( ~ mat1 * mat2)
    offset <- NCOL(mat1) + NCOL(mat2) + 2
    rmat <- rmat[,offset:NCOL(rmat)]

    cnames <- gsub("mat1", varname(e1), colnames(rmat))
    cnames <- gsub("mat2", varname(e2), cnames)
   

    colnames(rmat) <- cnames
    oname <- paste(varname(e1), ":", varname(e2), sep="")
    ev <- EV(rmat, oname, on1, blockids=as.integer(blockids(e1)),
             durations=as.integer(durations(e1)))
         
    attr(ev, "parents") <- list(varname(e1), varname(e2))
    ev
          
  })


setMethod("*", signature(e1="EventFactor", e2="EventVariableSet"),
  function(e1,e2) {
    callGeneric(e2,e1)
  })
    

setMethod("*", signature(e1="EventVariableSet", e2="EventFactor"),
  function(e1,e2) {
    .checkCompat(e1,e2)
    on1 <- onsets(e1)
    fac <- elements(e2)
    mat <- elements(e1)
    ocols <- length(levels(fac)) * NCOL(mat)
    omat <- matrix(0, length(on1), ocols)
    nlevs <- length(levels(fac))
 
    cnames <- paste(varname(e1), ":", levels(e1), sep="")
    cnames <- unlist(lapply(levels(fac), function(lev) paste(lev, ":", cnames, sep="")))
    
    for (i in 1:nlevs) {
      cstart <- (i-1)*NCOL(mat) + 1
      cend <- i*NCOL(mat)
      
      lev <- levels(fac)[i]
      idx <- which(fac == lev)
      omat[idx, cstart:cend] <- mat[idx,]
    }

    colnames(omat) <- cnames
    oname <- paste(varname(e1), ":", varname(e2), sep="")
    
    ev <- EV(omat, oname, on1, blockids=as.integer(blockids(e1)), durations=as.integer(durations(e1)))

    attr(ev, "parents") <- list(varname(e1), varname(e2))
    ev

    
  })
  

  

setMethod("*", signature(e1 = "EventVariable", e2="EventVariable"),
  function(e1,e2) {
    .checkCompat(e1,e2)
    

    v3 <- elements(e1)*elements(e2)

    
    name <- paste(varname(e1), ":", varname(e2), sep="")
    
    EventVariable(v3, name, onsets(e1),
                  blockids=as.integer(blockids(e1)),
                  durations=durations(e1))
    

  })
         


                         
setMethod("merge", signature(x = "EventVariable", y="EventVariable"),
  function(x,y, ...) {
    namex <- name(x)
    namey <- name(y)
    if (namex != namey) {
      warning("append: arguments have different \"name\" attributes")
    }

    .name <- namex
    
    .blockids <- c(blockids(x), blockids(y))
    .onsets <- c(onsets(x), onsets(y))
    .durations <- c(durations(x), durations(y))
    .value <- c(elements(x), elements(y))

    rest <- list(...)
    ret <- EventVariable(.value, .name, .onsets, durations=.durations, blockids=.blockids)

    if (length(rest) >= 1) {
      ret <- Reduce("merge", c(ret, rest))
    }

    return(ret)
                    
  })
     


setMethod("merge", signature(x = "EventFactor", y="EventFactor"),
  function(x,y, ...) {
    namex <- name(x)
    namey <- name(y)
    if (namex != namey) {
      warning("append: arguments have different \"name\" attributes")
    }

    .name <- namex
    
    .blockids <- c(blockids(x), blockids(y))
    .onsets <- c(onsets(x), onsets(y))
    .durations <- c(durations(x), durations(y))
    .value <- factor(c(elements(x), elements(y)))

    rest <- list(...)
    ret <- EventFactor(.value, .name, .onsets, durations=.durations, blockids=.blockids)

    if (length(rest) >= 1) {
      ret <- Reduce("merge", c(ret, rest))
    }

    return(ret)
                    
  })
     

setMethod("convolve", signature(x = "EventFactor", y="HRF", blocklens="numeric"),
          function(x, y, blocklens) {
            levs <- levels(x)          
            ons <- split(globalOnsets(x, blocklens), elements(x))
            RegressorSet(y, ons)
          })
           
           
setMethod("convolve", signature(x = "EventVariableSet", y="HRF", blocklens="numeric"),
          function(x, y, blocklens) {
            ons <- globalOnsets(x, blocklens)
            mat <- elements(x)
                     
            ampList <- lapply(1:NCOL(mat), function(i) mat[,i])
            onsetList <- lapply(1:NCOL(mat), function(i) ons)
            names(onsetList) <- paste(varname(x), ".", colnames(mat), sep="")
            
            RegressorSet(y, onsetList, ampList)
          })
           
           
                   
            
setMethod("blockids", signature(x = "EventVector"),
          function(x) x@blockids)
  
setMethod("varname", signature(x = "EventVector"),
          function(x) x@varname)

setMethod("levels", signature(x="EventVector"),
          function(x) {
            levels(as.factor(x@value))
          })

setMethod("conditions", signature(x = "EventFactor"),
          function(x) {
            olevs <- paste("[", levels(x), "]", sep="")
            paste(varname(x), olevs, sep="")
          })

setMethod("conditions", signature(x = "EventVariable"),
          function(x) {
            varname(x)
          })

setMethod("conditions", signature(x = "EventVariableSet"),
          function(x) {
            olevs <- paste("[", levels(x), "]", sep="")
            paste(varname(x), olevs, sep="")
          })
                
            

setMethod("levels", signature(x="EventVariableSet"),
          function(x) {
            colnames(elements(x))
          })


setMethod("elements", signature(x="EventVector"),
          function(x) {
            x@value
          })


setMethod("onsets", signature(x = "EventVector"),
    function(x) x@onsets)

setMethod("globalOnsets", signature(x= "EventVector", blocklens="numeric"),
          function(x, blocklens) {
            ons <- onsets(x)
            ids <- blockids(x)
            if (length(blocklens) != length(ids)) {
              stop("length of blocklens (run lengths)  must equal length of blokids (sessions)")
            }
            
            globons <- sapply(1:length(ons),function(i) {
              ons[i] + ((ids[i]-1) * blocklens[ids[i]])
            })
          })

setMethod("durations", signature(x = "EventVector"),
    function(x) x@durations)

