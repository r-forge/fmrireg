

EventTerm <- function(..., subset=NULL) {
	events <- list(...)
	stopifnot(all(sapply(events, function(ev) inherits(ev, "EventVector"))))
	
	cnames <- sapply(events, varname)
	names(events) <- cnames
	
	if (is.null(subset)) {
		## assumes all events have same length
		subset=rep(TRUE, length(events[[1]]))
	}
	new("EventTerm", events=events, subset=subset)
}




EV <- function(vals, name, onsets, blockids = 1, durations = NULL) {
	
	if (inherits(vals, "ParametricBasis")) {
		### omit name
		return(EventBasis(vals, onsets, blockids, durations))		
	}
	
	if (is.matrix(vals) && NCOL(vals) == 1) {
		vals <- vals[, 1, drop=TRUE]
	}
	
	if (is.vector(vals)) {
		return(EventVariable(vals, name, onsets, blockids, durations))
	}
	
	if (is.matrix(vals)) {
		return(EventVariableSet(vals, name, onsets, blockids, durations))
	}
	if (is.factor(vals)) {
		return(EventFactor(vals, name, onsets, blockids, durations))
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

EventBasis <- function(basis, onsets, blockids=1, durations=NULL) {
	stopifnot(inherits(basis, "ParametricBasis"))
	ret <- .checkEVArgs(basis@x, onsets, blockids, durations)
	
	#name <- paste(basis@fun, "_", basis@argname, "_", postfix, sep="")
	new("EventBasis", value=basis, varname=basis@name, onsets=ret$onsets, durations=ret$durations, blockids=as.integer(ret$blockids))
}



EventFactor <- function(facvar, name, onsets, blockids=1, durations=NULL) {
	if (!is.factor(facvar)) {
		facvar <- as.factor(facvar)
	}
	
	
	ret <- .checkEVArgs(facvar, onsets, blockids, durations)
	
	new("EventFactor", value=ret$vals, varname=name, onsets=ret$onsets, durations=ret$durations, blockids=as.integer(ret$blockids))
	
}          

EventVariable <- function(vals, name, onsets, blockids=1, durations=NULL) {
	stopifnot(is.vector(vals))
	
	if (is.factor(vals)) {
		stop("cannot create an EventVariable from a factor")
	}
	
	ret <- .checkEVArgs(vals, onsets, blockids, durations)
	
	new("EventVariable", value=ret$vals, varname=name, onsets=ret$onsets, durations=ret$durations,blockids=as.integer(ret$blockids))
}


EventVariableSet <- function(mat, name, onsets, durations=NULL, blockids=1 ) {
	stopifnot(is.matrix(mat))
	
	ret <- .checkEVArgs(mat[,1], onsets, blockids, durations)
	
	
	if (is.null(colnames(mat))) {
		colnames(mat) <- 1:NCOL(mat)
	}
	
	new("EventVariableSet", value=mat, varname=name, onsets=ret$onsets, durations=ret$durations, blockids=as.integer(ret$blockids))
}

setValidity("EventFactor", function(object) {
			(length(object@onsets) == length(object@value)) && (length(object@durations) == length(object@onsets))
		})


setMethod("[[", signature(x="EventTerm", i="character", j="missing"),
		function(x,i,j) {
			x@events[[i]]
		})

setMethod("[[", signature(x="EventTerm", i="numeric", j="missing"),
		function(x,i,j) {
			x@events[[i]]
		})

setMethod("[", signature(x="EventFactor", i="logical", j="missing", drop="missing"),
		function(x,i,j) {
			EventFactor(elements(x)[i, drop=TRUE], name=name(x), onsets=onsets(x)[i], durations=durations(x)[i],
					blockids=blockids(x)[i])
		})

setMethod("[", signature(x="EventFactor", i="numeric", j="missing", drop="missing"),
		function(x,i,j) {
			EventFactor(elements(x)[i, drop=TRUE], name=name(x), onsets=onsets(x)[i], durations=durations(x)[i],
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


setMethod("eventTable",	signature(x="EventTerm"),
		function(x) {
			
			pterms <- parentTerms(x)
			ret <- lapply(pterms, function(nam) {
						if (isContinuous(x[[nam]])) {
							rep(.sanitizeName(nam), length(onsets(x)))
						} else {
							x[[nam]]@value
						}			
					})
			
			ret <- as.data.frame(ret)
			names(ret) <- sapply(pterms, .sanitizeName)
			ret
		})


setMethod("cells", signature(x="EventTerm"),
		function(x, drop.unused=TRUE, expand=FALSE) {    
			
			
			evtab <- eventTable(x)
			evset <- expand.grid(lapply(x@events, levels))
			which.cat <- which(!sapply(x@events, isContinuous))
			
			### counts should only be based on factors
			counts <- apply(evset[,which.cat,drop=F], 1, function(row1) {
						sum(apply(evtab[x@subset,which.cat,drop=F], 1, function(row2) {										
											all(row1 == row2)
										}))
					})
			
			
			if (drop.unused) {				
				evset <- evset[counts > 0,,drop=F]
				attr(evset, "count") <- counts[counts > 0]
				
			} else {
				attr(evset, "count") <- counts
			}
			
			evset
		})

setMethod("covariates", signature(x="EventTerm"),
		function(x) {
			keep <- sapply(x@events, isContinuous)
			if (!any(keep)) {
				warning("this term has no covariates, returning empty list")
				list()
			} else {
				ret <- x@events[keep]
				names(ret) <- parentTerms(x)[keep] 
				ret
			}
			
		})

setMethod("factors", signature(x="EventTerm"),
		function(x) {
			keep <- sapply(x@events, function(el) !isContinuous(el))
			ret <- if (!any(keep)) {
						warning("this term has no factors, returning empty list")
						list()
					} else {
						ret <- x@events[keep]
						names(ret) <- parentTerms(x)[keep] 
						ret
					}
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




.sanitizeName <- function(name) {
	name <- gsub(":", ".", name)
	name <- gsub(" ", "", name)
	name <- gsub("[\\(\\)]", ".", name, perl=TRUE)
	name <- gsub(",", ".", name)
	name <- gsub("\\.$", "", name)
	name
}


setMethod("predict",  signature(object="EventVector"), 
		function(object, y, onsets=0, durations=0) {
			stopifnot(length(onsets) == length(y))
			EV(y, varname(object), onsets, durations)
		})

setMethod("predict",  signature(object="EventVariableSet"), 
		function(object, y, onsets=0, durations=0) {
			stopifnot(length(onsets) == length(y))
			els <- elements(object)
			ret <- try(predict(els,  y))
			if (inherits(ret,  "try-error")) {
				ret <- y
			}
			EventVariableSet(ret, varname(object), onsets, durations)
		})

setMethod("predict",  signature(object="EventFactor"), 
		function(object, y, onsets=0, durations=0) {
			stopifnot(length(onsets) == length(y))
			stopifnot(all(levels(y) %in% levels(object)))
			stopifnot(length(levels(y)) <= length(levels(object)))
			
			vy <- as.character(y)
			fac <- factor(vy, levels=levels(object@value))
			#if (length(levels(y)) < length(levels(object))) {
			#	levels(y) <- c(levels(y),  setdiff(levels(object),  levels(y)))			
			#}
			
			#.modelMat(varname(object),  y)
			EventFactor(fac, varname(object), onsets, durations)
		})

setMethod("predict", signature(object="EventBasis"), 
		function(object, y, onsets=0, durations=0) {
			stopifnot(length(onsets) == length(y))
			mat <- predict(object@value, y)
			EV(mat, varname(object), onsets, durations)
		})



setMethod("predict",  signature(object="scaled"), 
		function(object,  y) {
			center <- attr(object, "scaled:center")
			sc <- attr(object, "scaled:scale")
			if (is.null(center) && is.null(sc)) {
				y
			} else if (is.null(center)) {
				y/sc
			} else if (is.null(sc)) {
				(y-center)		
			} else {
				(y-center)/sc
			}
		})


.modelMat <- function(vname,  vals) {
	locenv <- new.env()
	vname <- gsub(":", "_", vname)
	assign(vname, vals, envir = locenv)
	form <- as.formula(paste("~", vname, "-1"))
	model.matrix(form, data = locenv)
}	

setMethod("formula",  signature(x = "EventTerm"),
		function(x) {
			as.formula(paste("~ ", "(", paste(parentTerms(x), collapse=":"), "-1", ")"))
		})

setMethod("designMatrix", signature(x = "EventVector"), function(x) {
			mat <- .modelMat(varname(x),  elements(x))
			colnames(mat) <- conditions(x)
			mat
		})


setMethod("designMatrix", signature(x="EventTerm"),
		function(x, drop.unused.levels=TRUE) {
			
			locenv <- new.env()
			pterms <- sapply(parentTerms(x), .sanitizeName)					
			for (ev in x@events) {
				vname <- .sanitizeName(varname(ev))
				els <- elements(ev)
				lapply(names(els), function(n) assign(n, els[[n]],envir=locenv))			
			}
			
			
			els <- as.data.frame(elements(x))		
			nas <- try(apply(els,1, function(vals) any(is.na(vals))))			
			
			counts <- attr(cells(x, drop=FALSE), "count")
			
			if (ncol(els) == 1 && is.factor(els[,1]) && length(levels(els[,1])) == 1) {
				### this is a term with only an intercept
				### column of 1s
				mat <- cbind(rep(1, NROW(els)))
				#mat <- els[,1,drop=F]
			} else { 		
				form <- as.formula(paste("~", paste(pterms, collapse=":"), "-1"))
				mat <- model.matrix(form, data=locenv)
			}
			
			### multiply design matrix by subset
			rmat <- mat * x@subset
			
			#remove rows with NAS
			if (any(nas)) {
				rmat <- matrix(0, length(x), length(conditions(x)))
				rmat[!nas,] <- mat
				rmat[nas,] <- NA				
			} 
			
			#browser()
			
			# remove columns with no events (postpone this to later stage) 
			if (any(counts == 0) && (length(conditions(x, drop=F)) == length(counts)) && drop.unused.levels) {
				rmat <- rmat[, !(counts==0), drop=FALSE]
				colnames(rmat) <- conditions(x, drop=T)
			} else {
				colnames(rmat) <- conditions(x, drop=F)			
			}
			
			
			
			rmat
			
		})




setMethod("*", signature(e1 = "EventFactor", e2="EventFactor"),
		function(e1,e2) {
			.checkCompat(e1,e2)
			
			if (identical(e1, e2)) {
				stop(paste("cannot cross identical factors", varname(e1), " and ", varname(e2)))
			}
			
			
			f3 <- interaction(elements(e1), elements(e2), sep=":", drop=T) 
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

setMethod("convolve", signature(x= "EventTerm"),
		function(x,hrf,blocklens, TR) {
			ids <- rep(1:length(unique(blockids(x))), table(blockids(x)))
			onsets <- onsets(x)
			globons <- sapply(1:length(onsets),function(i) {
						onsets[i] + ((ids[i]-1) * blocklens[ids[i]]) * TR
					})
			
			nimages <- sum(blocklens)
			start <- TR/2
			samples <- seq(start, length.out=nimages, by=TR)
			
			dmat <- designMatrix(x)
			
			convolvedMat <- apply(dmat, 2, function(vals) {
						Regressor(hrf, onsets, amp=vals*x@subset, durations=durations(x), granularity=.1)(samples)
					})
			
		})


setMethod("blockids", signature(x = "EventVector"),
		function(x) x@blockids)


setMethod("blockids", signature(x = "EventTerm"),
		function(x) blockids(x@events[[1]]))

setMethod("blockids", signature(x = "EventFunction"),
		function(x) blockids(x@value))

setMethod("varname", signature(x = "EventVector"),
		function(x) x@varname)

setMethod("varname", signature(x = "EventTerm"),
		function(x) {
			paste(sapply(x@events, varname), collapse=":")
		})

setMethod("argname", signature(x="EventVector"),
		function(x) {
			varname(x)
		})
setMethod("argname", signature(x="EventBasis"),
		function(x) {
			x@value@argname
		})

setMethod("parentTerms", signature(x="EventTerm"),
		function(x) {
			lapply(x@events, varname)
		})

setMethod("names", signature(x="EventTerm"),
		function(x) {
			as.vector(unlist(lapply(x@events, names)))		
		})

setMethod("names", signature(x="EventVector"),
		function(x) {
			varname(x)
		})

setMethod("names", signature(x="EventVariableSet"),
		function(x) {
			paste(.sanitizeName(varname(x)), ".", levels(x), sep="")
		})


setMethod("isContinuous", signature(x="EventBasis"), 
		function(x) {
			TRUE
		})

setMethod("isContinuous", signature(x="EventFactor"),
		function(x) {
			FALSE
		})

setMethod("isContinuous", signature(x="EventVariable"),
		function(x) {
			TRUE
		})

setMethod("isContinuous", signature(x="EventVariableSet"),
		function(x) {
			TRUE
		})

setMethod("isContinuous", signature(x="EventFunction"),
		function(x) {
			TRUE
		})

setMethod("isContinuous", signature(x="EventTerm"),
		function(x) {
			any(sapply(x@events, isContinuous))
		})


setMethod("levels", signature(x="EventVector"),
		function(x) {
			varname(x)
		})

setMethod("levels", signature(x="EventVariableSet"),
		function(x) {
			colnames(x@value)
		})

setMethod("levels", signature(x="EventTerm"),
		function(x) {
			facs <- x@events[!sapply(x@events, isContinuous)]
			if (length(facs) == 1) {
				levels(facs[[1]])
			} else {
				facs <- lapply(facs, elements)
				f <- function(...) {
					interaction(..., drop=TRUE, sep=":")
				}
				levels(do.call(f, facs))			
			}
						
		})

setMethod("levels", signature(x="EventFactor"),
		function(x) {
			levels(x@value)
		})

setMethod("levels", signature(x="EventBasis"),  
		function(x) {
			colnames(x@value@y)
		})

setMethod("nlevels", signature(x="EventTerm"),
		function(x) {
			length(levels(x))
		})

setMethod("nlevels", signature(x="EventFactor"),
		function(x) {
			length(levels(unlist(elements(x))))
		})

setMethod("nlevels", signature(x="EventVariable"),
		function(x) {
			1
		})

setMethod("nlevels", signature(x="EventVariableSet"),
		function(x) {
			NCOL(elements(x)[[1]])
		})

setMethod("nlevels", signature(x="EventBasis"),
		function(x) {
			NCOL(x@value@y)
		})


setMethod("conditions", signature(x="EventTerm"),
		function(x, drop=TRUE) {
			
			.cells <- cells(x, drop=drop)
			pterms <- parentTerms(x)
			levs <- apply(.cells, 1, paste, collapse=":")
			
			splitlevs <- strsplit(levs, ":")
			ret <- lapply(1:length(pterms), function(i) {
						lev <- sapply(splitlevs, "[[", i)
						if (nlevels(x[[pterms[[i]]]]) > 1) {
							paste(.sanitizeName(pterms[i]), "[", lev, "]", sep="")
						} else {
							.sanitizeName(pterms[i])
						}
					})
			
			
			do.call(function(...) paste(..., sep=":"), ret)
			
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
			paste(.sanitizeName(varname(x)), olevs, sep="")
		})

setMethod("conditions", signature(x = "EventFunction"),
		function(x) {
			paste(x@funcname,  "_",  conditions(x@value),  sep="")
		})

setMethod("elements", signature(x="EventVariable"),
		function(x, values=TRUE,...) {
			if (values) {
				ret <- list(x@value)
				names(ret) <- varname(x)
				ret
			} else {
				ret <- list(rep(varname(x), length(x)))
				names(ret) <- varname(x)
				ret		
			}
		})


.vset.names <- function(x) {
	paste(.sanitizeName(varname(x)), ".", levels(x), sep="")	
}



setMethod("elements", signature(x="EventVariableSet"),
		function(x, values=TRUE,...) {
			if (values) {
				ret <- x@value
				colnames(ret) <- .vset.names(x)
				n <- .sanitizeName(varname(x))
				ret <- list(ret)
				names(ret) <- n
				ret
			} else {
				N <- length(x)
				vnames <- .vset.names(x)
				res <- lapply(vnames, function(el) rep(el, N))
				mat <- do.call(cbind, res)
				colnames(mat) <- vnames			
				ret <- list(mat)
				names(ret) <- .sanitizeName(varname(x))
				ret			
			}
		})

setMethod("elements", signature(x="EventVector"),
		function(x, values=TRUE,...) {		
			ret <- list(x@value)
			names(ret) <- varname(x)
			ret
		})

### code duplication with EventVariableSet -- perhaps EventBasis should wrap an EventVariableSet?

setMethod("elements", signature(x="EventBasis"), 
		function(x, values=TRUE, transformed=TRUE) { 
			if (values && !transformed) {
				x@value@x				
			}
			else if (values) {
				ret <- x@value@y
				colnames(ret) <- .vset.names(x)
				n <- .sanitizeName(varname(x))
				ret <- list(ret)
				names(ret) <- n
				ret
			} else {
				N <- length(x)
				vnames <- .vset.names(x)
				res <- lapply(vnames, function(el) rep(el, N))
				mat <- do.call(cbind, res)
				colnames(mat) <- vnames			
				ret <- list(mat)
				names(ret) <- .sanitizeName(varname(x))
				ret			
			}
			
		})

setMethod("elements", signature(x="EventTerm"), 
		function(x, values=TRUE) {
			els <- lapply(x@events, elements, values=values)
			n <- sapply(names(els), function(nam) .sanitizeName(nam))
			names(els) <- as.vector(n)
			els
			
		})

setMethod("splitOnsets", signature(x="EventTerm"), 
		function(x, global=FALSE, blockDurations=NULL) {
			
			### need to check for 0 factors
			facs <- x@events[!sapply(x@events, isContinuous)]
			facs <- lapply(facs, function(fac) unlist(elements(fac)))
			
			f <- function(...) {
				interaction(..., drop=TRUE, sep=":")
			}
			
			ret <- try(do.call(f, facs))
			
			if (inherits(ret, "try-error")) {
				browser()
			}
			
			if (global) {
				if (is.null(blockDurations)) {
					stop("must supply blockDurations argument to compute global onsets")
				}
				split(globalOnsets(x@events[[1]], blockDurations=blockDurations), ret)
			} else {
				split(onsets(x), ret)
			}
		})


setMethod("length", signature(x="EventVector"), function(x) {
			length(onsets(x))
		})

setMethod("length", signature(x="EventTerm"), function(x) {
			length(parentTerms(x))
		})

setMethod("onsets", signature(x = "EventVector"),
		function(x) x@onsets)


setMethod("onsets", signature(x = "EventTerm"),
		function(x) onsets(x@events[[1]]))

setMethod("globalOnsets", signature(x= "EventVector", blockDurations="numeric"),
		function(x, blockDurations) {
			onsets <- onsets(x)			
			ids <- rep(1:length(unique(blockids(x))), table(blockids(x)))
			
			sapply(1:length(onsets),function(i) {
						onsets[i] + (((ids[i]-1) * blockDurations[ids[i]]))
					})
			
			
		})



setMethod("durations", signature(x = "EventVector"),
		function(x) x@durations)

setMethod("durations", signature(x = "EventTerm"),
		function(x) durations(x@events[[1]]))


setMethod("show", signature(object="EventTerm"),
		function(object) {
			cat("Event Term", "\n")
			cat(" ", "Term Name: ", names(object), "\n")
		})



