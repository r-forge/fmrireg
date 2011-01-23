

.parseTypes <- function(varnames) {
  sapply(varnames, function(vn) {
    class(parse(text=vn)[[1]])
  })
}
       
      
  


EventModel <- function(formula, data, block, durations=0, center=T, scale=F) {

  if (!inherits(formula, "terms")) {
    formula <- terms(formula, data = data)
  }

 
  env <- environment(formula)

  #rownames <- .row_names_info(data, 0L)
  vars <- attr(formula, "variables")
  
  predvars <- attr(formula, "predvars")
  if (is.null(predvars)) {
    predvars <- vars
  }


  block <- eval(substitute(block), envir=data)
    
  durations <- if (is.name(quote(durations))) {
    eval(substitute(durations), envir=data)
  } else {
    as.numeric(durations)
  }

  #browser()

 
  ### may need way to include "intercept" i.e. column of all ones for main effects in the case where no factors are in model
  
  varnames <- sapply(predvars, deparse, width.cutoff = 500)[-1]
  resp <- attr(formula, "response")
  
  if (center || scale) {
    for (i in 1:length(varnames)) {
      vn <- varnames[i]
      if (!(i == resp || is.factor(data[[vn]]))) {
        data[[vn]] <- scale(data[[vn]], center=center, scale=scale)[,,drop=TRUE]     
      }
      
    }
  }

  
  variables <- eval(predvars, data, env)  
  names(variables) <- varnames  
  # check if resp > 0
  lhs <- variables[[resp]]

              
  vartypes <- .parseTypes(varnames)
 
  
  #figure out matrix of interactions
  expmat <- attr(formula, "factors")

  etab <- .createEventTable(expmat, varnames, variables, lhs, block, durations)

 
  facnames <- names(factors(etab))   
  varnames <- names(covariates(etab))

  ## find the highest order interaction term and add
  facind <- row.names(expmat) %in% facnames
 
  fac.terms <- paste(facnames, collapse=":")
  var.terms <- .extractVarTerms(varnames, facnames, expmat)
  

  formula.terms <- paste(unique(c(fac.terms, var.terms)), collapse=" + ")  
  eventformula <- as.formula(paste("~", formula.terms, " -1"), env=environment(formula))
  
  new("EventModel", eventTable=etab, formula=eventformula, context=data) 
}

setMethod("initialize", "EventModel", function(.Object, eventTable, formula, context) {
  .Object <- callNextMethod()

  .Object@eventTable <- eventTable
  .Object@formula <- formula
  .Object@context <- context
  .Object@designMatrix <- .createDesignMatrix(eventTable, formula, context)
  .Object
})
  




.extractVarTerms <- function(varnames, facnames, expmat) {
  varind <- row.names(expmat) %in% varnames
  varmat <- expmat[varind,, drop=F]
  var.terms <- unlist(lapply(1:NROW(varmat), function(i) {
    idx <- which(varmat[i,] > 0)
    if (length(idx) > 0) {
      vnames <- unique(unlist((apply(expmat[,idx, drop=F], 2, function(vals) names(vals[vals>0])))))
      if (all(vnames %in% varnames)) {
        ## all terms are varibles, therefore need to keep main effects
        colnames(varmat)[idx]
      } else {
        #browser()
        ## variable is crossed with factor, need only to keep higher order terms       
        term.count <- t(apply(expmat, 2, function(c1) {
          parts <- names(c1[c1>0])
          nvars <- sum(parts %in% varnames)
          nfacs <- sum(parts %in% facnames)
          c(nvars, nfacs)
        }))

        maxfac <- max(term.count[,2])
        idx <- which(term.count[,2] == maxfac & term.count[,1] > 0)
        term.names <- rownames(term.count[idx,])
        term.names                      
     
      }
    }
  }))

}
  


.createEventTable <- function(expmat, varnames, variables, lhs, block, durations) {
  nterms <- 0
  evs <- list()

  #nfactors <- sum(sapply(variables, is.factor))
  #if (nfactors == 0) {
  #  nterms <- 1
  #  evs[[nterms]] <- EV(factor(rep("ON", length(lhs))), "__STIMULUS__", lhs, block, durations)
  #}
  
  for (i in 1:NCOL(expmat)) {
    idx <- which(expmat[,i] == 1)
    if (length(idx) == 1) {
      nterms <- nterms + 1
      vname <- varnames[idx]      
      evs[[nterms]] <- EV(variables[[vname]], vname, lhs, block, durations)      
    }
  }

 
  EventTable(evs)
}
  



setMethod("formula",  signature(x = "EventModel"),
          function(x) {
            x@formula
          })

setMethod("eventTable", signature(x="EventModel"),
          function(x) {
            x@eventTable
          })


setMethod("onsets", signature(x="EventModel"),
          function(x) {
            onsets(eventTable(x))
          })

setMethod("durations", signature(x="EventModel"),
          function(x) {
            durations(eventTable(x))
          })

setMethod("blockids", signature(x="EventModel"),
          function(x) {
            blockids(eventTable(x))
          })

setMethod("[[", signature(x="EventModel", i = "character", j = "missing"),
          function(x, i, j) {
            eventTable(x)[[i]]
          })

setMethod("[[", signature(x="EventModel", i = "numeric", j = "missing"),
          function(x, i, j) {
            eventTable(x)[[i]]
          })

setMethod("conditions", signature(x="EventModel"),
          function(x) {
            conditions(eventTable(x))
          })

setMethod("cells", signature(x="EventModel"),
          function(x) {
            condlist <- .levelTables(x@eventTable, x@formula, x@context)
            ret <- unlist(lapply(condlist, function(condmat) {
              if (prod(dim(condmat)) == 1) {
                ret <- list()
                ret[[names(condmat)]] = ""
                list(ret)
              } else {
                nlevs <- apply(condmat, 2, function(x) length(unique(x)))             
                ret <- apply(as.matrix(condmat), 1, as.list)
                if (any(nlevs == 1)) {
                  ret <- lapply(ret, function(el) {
                    tmp <- el; tmp[which(nlevs==1)] <- ""; tmp
                  })
                }

                ret
                      
                
                
              }
            }), recursive=FALSE)
                          
            
          })
                                 
setMethod("colnames", signature(x="EventModel", do.NULL="missing", prefix="missing"),
          function(x) {
            colnames(designMatrix(x))         
          })

setMethod("shortnames", signature(x="EventModel"),
          function(x) {
            condlist <- .levelTables(x@eventTable, x@formula, x@context)
            cnames <- lapply(condlist, function(condmat) {
             
              apply(as.matrix(condmat), 1, function(row) {
                nlevs <- apply(condmat, 2, function(v) length(unique(v)))
                olevs <- sapply(1:length(row), function(i) {
                  if (nlevs[i] > 1) {
                    row[i]
                  } else {
                    ""
                  }
                })

                paste(olevs, sep="", collapse=":")
              })
            })

            unlist(cnames)

          })


         
setMethod("terms", signature(x="EventModel"),
          function(x) {
            ret <- terms(formula(x))
            attr(ret, "term.labels")
          })

setMethod("factors", signature(x="EventModel"),
          function(x) {
            names(factors(eventTable(x)))
          })

.levelTables <- function(etab, formula, context) {
  mat <- model.matrix(formula, data=context)
  term.idx <- attr(mat, "assign")
  term.names <- colnames(attr(terms(formula), "factors"))
  runs <- rle(term.idx)
  ctabs <- sapply(runs$values, function(val) {
    runlen <- runs$lengths[val]
    facname <- term.names[val]
    if (runlen == 1) {
      names(facname) <- facname
      facname
    } else {

       #must disallow factor names with colons ...
      facparts <- strsplit(facname, ":")[[1]]
      levlist <- lapply(facparts, function(fpart) {
        els <- elements(etab[[fpart]])
        if (is.factor(els) || NCOL(els) > 1) {
          levels(etab[[fpart]])
        } else {
          fpart
        }
      })

      condmat <- expand.grid(levlist)
      colnames(condmat) <- facparts
      condmat
    }
  })
}
  
  


.createDesignMatrix <- function(etab, formula, context) {
 
  mat <- model.matrix(formula, data=context)
  condlist <- .levelTables(etab, formula, context)
  cnames <- lapply(condlist, function(condmat) {
    condmat <- as.matrix(condmat)
    if (prod(dim(condmat) == 1)) {
      return(rownames(condmat))
    }

    nlevs <- apply(condmat, 2, function(v) length(unique(v)))
    apply(as.matrix(condmat), 1, function(row) {
      olevs <- sapply(1:length(row), function(i) {
        if (nlevs[i] > 1) {
          paste("[", row[i], "]", sep="")
        } else {
            ""
          }
        })
      paste(names(row), olevs, sep="", collapse=":")
    })
  })
                   
  colnames(mat) <- unlist(cnames)
  mat
}
  

                   
          
setMethod("designMatrix",  signature(x = "EventModel"),
          function(x) {
            x@designMatrix            
          })


          
