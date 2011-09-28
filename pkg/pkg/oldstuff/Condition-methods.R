

ConditionModel <- function(conlist) {
  if (all(sapply(conlist, inherits, "data.frame"))) {
    conlist <- lapply(conlist, function(clist) ConditionList(clist))
  }
  
  stopifnot(all(sapply(conlist, inherits, "ConditionList")))
  
  new("ConditionModel", conditionList=conlist)
}
    

ConditionList <- function(cell.frame) {
  cell.frame <- as.data.frame(cell.frame)
  
  varnames <- names(cell.frame)
  levs <- apply(cell.frame, 2, unique)
  
  conlist <- list()
  
  for (i in 1:NROW(cell.frame)) {
    ret <- lapply(cell.frame[i,], as.character)
    conlist[[i]] <- Condition(names(ret), as.character(ret))
  }

  new("ConditionList", conditions=conlist, table=cell.frame)
}
  


Condition <- function(names, levels, types=NULL) {
  if (is.null(types)) {
    types <- rep("factor", length(levels))
  }  

  new("Condition", names=names, levels=levels, types=types)
}

setMethod("names", signature(x="ConditionModel"),
          function(x) {
           unique(unlist(sapply(x@conditionList, names)))
          })

setMethod("names", signature(x="ConditionList"),
          function(x) {
            names(x@table)
          })

setMethod("levels", signature(x="ConditionList"),
          function(x) {
			ret <- lapply(x@conditions, levels)
	        unlist(lapply(ret, paste, collapse=":"))
          })

setMethod("levels", signature(x="ConditionModel"),
	function(x) {
		sapply(x@conditionList, levels, simplify=FALSE)
	})
	
setMethod("[[", signature(x="ConditionList", i = "character", j = "missing"),
	function(x, i, j) {
		x@table[[i]]
	})
	
setMethod("[[", signature(x="ConditionModel", i = "character", j = "missing"),
	function(x, i, j) {
		ret <- lapply(x@conditionList, "[[", i)
		ret <- ret[!sapply(ret, is.null)]
		ret[[1]]
	})

	

setMethod("cells", signature(x="ConditionList"),
          function(x) {
            x@table
          })

setMethod("cells", signature(x="ConditionModel"),
          function(x) {
             lapply(x@conditionList, cells)
          })

setMethod("conditions", signature(x="ConditionList"),
          function(x) {
            sapply(x@conditions, longform)
          })

setMethod("conditions", signature(x="ConditionModel"),
          function(x) {
            ret <- unlist(lapply(x@conditionList, longform))
            names(ret) <- NULL
            ret           
          })

.eval.seq <- function(string.seq, eval.seq=NULL) {
	if (is.null(eval.seq)) {  
		stopifnot(length(string.seq) >= 3)
		first.seq <- .eval.args(string.seq[1], string.seq[2], string.seq[3])
		Recall(string.seq[4:length(string.seq)], first.seq)
	} else {    
		if (length(string.seq) == 1) {
			return(eval.seq)
		} else {
			op <- string.seq[1]
		    arg2 <- string.seq[2]     
		    Recall(string.seq[3:length(string.seq)], .eval.args(eval.seq, op, arg2))
		 }
	}
}



expandPattern <- function(components) {
	
	exp.comps <- lapply(components, function(comp) {
		.expandComponent(comp)
	})
	
	comp.grid <- expand.grid(exp.comps)
	apply(comp.grid, 1, paste, collapse=":")
}

#.makeExpression(comp) {
#	f1 <- strsplit(comp, "\\|")[[1]]
#	f2 <- paste("x==", f1, sep="", collapse=" | ")
#}


.expandComponent <- function(comp) {

  len <- nchar(comp)
  res <- regexpr("(^\\{)([\\w|\\.\\.|\\|]+)(\\}$)", comp, perl=TRUE)
  
  if (res > 0 && attr(res, "match.length") != len) {
    stop("illegal subexpression:", comp)
  }
  if (res > 0 && attr(res, "match.length") == len) {
    ret <- strapply(comp, "(\\{)(.+)(\\})", c)

    stopifnot(length(ret) == 1)
    stopifnot(length(ret[[1]]) == 3)

    meat <- ret[[1]][[2]]

    string.seq <- strapply(meat, "(\\w+\\s*)(\\.\\.|\\|)*", FUN <- function(x,...) {
      args <- list(...)

      if (args[[1]] != "") {
        c(x,args[[1]])
      } else {
        x
      }
    })

    string.seq <- c(string.seq[[1]], "")   
    .eval.seq(string.seq)

  } else {
    comp
  }
}

# should allow longform pattern selection, e.g.:      
# SEQ1:basis[{1..2}]  
# FACTOR[{pattern}]
# FACTOR[{!(expression)}]

# :level: notation means match any levels at any position, e.g. :left: would match left:basis1, left:basis2 etc.
# :factor[level]: should works as well
# :{level1 | level2}: 
# :factor[{level1 | level2}]:

.isPatternExpr <- function(patstr) {
	ret <- regexpr("(^{)([\\w&|!\\.]+)(}$)", patstr, perl=T)	
	attr(ret, "match.length") == nchar(patstr)
}

.isLongForm <- function(patstr) {
	pat <- "(^\\w+)(\\[)([^[\\]]+)(\\])"
	ret <- regexpr(pat, patstr, perl=T)
	attr(ret, "match.length") == nchar(patstr)
}

.extractLongForm <- function(patstr) {
	ret <- strapply(patstr, "(^\\w+)(\\[)([^[\\]]*)(\\])", c, perl=T)
	if (is.null(ret[[1]])) {
		stop(paste("failed to parse expression: ", patstr))
	}
	
	list(factor=ret[[1]][1], level=ret[[1]][3])
}

.lookupFactor <- function(conmod, match.levs) {
	varnames <- names(conmod)
	ret <- lapply(varnames, function(fac) {
		flevs <- levels(conmod[[fac]])
		if (all(match.levs %in% flevs)) {
			fac
		} else {
			NULL
		}
	})
	
	unique(unlist(ret))
}

resolvePattern <- function(conmod, constr) {
   
	if (substr(constr, 1,1) == ":") {
		comps <- strsplit(substr(constr, 2, nchar(constr)), ":")[[1]]		
	} else {
		comps <- strsplit(constr, ":")[[1]]
	}
		
	conlist <- list()
	for (i in 1:length(comps)) {
		if (.isLongForm(comps[i])) {
			parts <- .extractLongForm(comps[i])
					
		    if (!any(names(conmod) %in% parts$factor)) {	
				stop(paste("variable ", facname, " not in model with variables: ", names(conmod), "\n"))						
			} else {
				idx <- which(names(conmod) %in% parts$factor)
				conlist[[parts$factor]] <- .expandComponent(gsub(" ", "", parts$level))				
			}
							
		} else {
		
			excomp <- .expandComponent(gsub(" ", "", comps[i]))
			fac <- .lookupFactor(conmod, excomp)
			if (is.null(fac)) {
				stop(paste("cannot find matching variable for level pattern: ", comps[i]))
			} else {
				conlist[[fac]] <- excomp		
			}	
		}
	}
	
	#ord.idx <- sapply(names(conlist), function(nam) which(names(conmod) %in% nam))
	#conlist[ord.idx]
	conlist
}

setMethod("patternMatch", signature(x="ConditionModel", y="character"),
          function(x,y, ...) {
			# remove whitespace
	        constr <- gsub(" ", "", y)
			resolved <- resolvePattern(x, constr)
			browser()
            patterns <- expandPattern(resolved)
            print(paste("patterns are: ", patterns))
            ## need to check that all expanded levels are valid

            res <- lapply(patterns, function(pattern) {
              unlist(lapply(x@conditionList, patternMatch, pattern))
            })

            ret <- as.logical(apply(do.call("rbind", res), 2, sum))
            names(ret) <- conditions(x)
            ret
            
          })

setMethod("patternMatch", signature(x="ConditionList", y="character"),
          function(x,y, ...) {
            sapply(x@conditions, isMatch, y)
          })

setMethod("names", signature(x="Condition"),
          function(x) {
            x@names
          })

setMethod("levels", signature(x="Condition"),
          function(x) {
            x@levels
          })

setMethod("shortform", signature(x="Condition"),
          function(x) {
            paste(levels(x), collapse=":")
          })

setMethod("longform", signature(x="Condition"),
          function(x) {          
            levs <- paste("[", levels(x), "]", sep="")
            paste(names(x), levs, sep="", collapse=":")
          })

setMethod("longform", signature(x="ConditionList"),
          function(x) {
            sapply(x@conditions, longform)
            
          })

setMethod("longform", signature(x="ConditionModel"),
          function(x) {
            sapply(x@conditionList, longform)
            
          })

setMethod("shortform", signature(x="ConditionList"),
          function(x) {
            sapply(x@conditions, shortform)
            
          })

setMethod("length", signature(x="Condition"),
          function(x) {
            length(names(x))
          })

.transformConStr <- function(condition, constr) {
	#browser()
  comps <- strsplit(constr, ":")[[1]]
  conspec <- list()

  for (i in 1:length(comps)) {
    if (regexpr("(^\\w+)(\\[)(\\w+)(\\])", comps[i])[1] >0) {
      
      ret <- strapply(comps[i], "(^\\w+)(\\[)(\\w+)(\\])", c)
      
      if (length(ret) != 1 || length(ret[[1]]) != 4) {
        stop(paste("invalid constrast representation:", constr))
      }
      
      facname <- ret[[1]][1]
      levname <- ret[[1]][3]

      if ((facname %in% names(condition))) {
        conspec[[facname]] <- levname
        #stop(paste("invalid contrast spec: ", constr, " -- factor ", facname, " invalid for condition: ", longform(condition)))
      }     
        
    } else {
      ### shorthand form, must be in order
      if (comps[i] == "_") {
        conspec[[names(condition)[i]]] <- levels(condition)[i]
      } else if (comps[i] == levels(condition)[i]){
        
        ## valid level
        facname <- names(condition)[i]
        if (!is.null(conspec[[facname]])) {
          stop(paste("invalid constrast spec: ", constr, " -- two components refer to same factor: ", facname))
        }
        conspec[[facname]] <- comps[i]
      } #else {
        #stop(paste("invalid contrast spec: ", constr, "--level ", comps[i], " invalid for condition: ", longform(condition)))
      #}
    }
  }

  conspec <- conspec[match(names(condition), names(conspec))]
  conspec <- conspec[!sapply(conspec, is.null)]

}
      
        
    
setMethod("isMatch", signature(x="Condition", y="character"),
          function(x,y, ...) {
			ncomps <- length(strsplit(y, ":")[[1]])
	        if (ncomps != length(levels(x))) {
				# number of pattern components does not equal number of levels for condition
				return(FALSE)
			}
            
            comps <- .transformConStr(x, y)
            levs <- levels(x)
            
            if (length(comps) != length(x)) {
              return(FALSE)
            }
            
            for (i in 1:length(comps)) {
              if (comps[i] == "_") {
                next
              }
              if (levs[i] != comps[i]) {
                return(FALSE)
              }
            }

            TRUE
          })


setMethod("show", signature(object="Condition"),
          function(object) {
            cat("condition: ", longform(object), "\n")
            
          })

setMethod("show", signature(object="ConditionList"),
          function(object) {
            cat("conditions: ", sapply(object@conditions, longform), "\n")
            
          })

setMethod("show", signature(object="ConditionModel"),
          function(object) {
            cat("conditions: ", unlist(lapply(object@conditionList, longform)), "\n")
            
          })



		.eval.args <- function(arg1, op, arg2) {
		  if (op == "|") {
		    c(arg1,arg2)
		  } else if (op == "..") {

		    a1 <- as.numeric(arg1[length(arg1)])
		    a2 <- as.numeric(arg2)
		    if (is.na(a1) || is.na(a2)) {
		      stop(paste("cannot evaluate sequence arguments: ", paste(arg1[length(arg1)], "and", arg2, collapse=" ")))
		    }
		    piece <- seq(as.numeric(arg1[length(arg1)]), as.numeric(arg2))
		    if (length(piece) < 2) {
		      arg1
		    } else {
		      c(arg1, piece[2:length(piece)])
		    }

		  } else {
		    stop("illegal operator: ", op)
		  }
		}

		


		

		
          
          
            
            
            
