

.parseTypes <- function(varnames) {
	sapply(varnames, function(vn) {
				class(parse(text=vn)[[1]])
			})
}



center <- function(vals, modefun=mean) {
	vals - modefun(vals)
}

.extractTerms <- function(formula, data) {
	if (!inherits(formula, "terms")) {
		terms(formula, data = data)
	} else {
		formula
	}	
}


extractCovariates <- function(.terms, variables, resp, etab) {
	vars <- attr(.terms, "variables") 
	varnames <- sapply(vars, deparse, width.cutoff = 500)[-1]
	ind.vars <- varnames[-resp] 
	orig.covar.names <- ind.vars[which(sapply(variables[-resp], function(obj) is.numeric(obj) || .is.parametric.basis(obj)))]
	new.covar.names <- names(events(etab))[sapply(events(etab), isContinuous)]
	covar.names <- as.list(orig.covar.names)
	names(covar.names) <- new.covar.names
	covar.names
}

.is.parametric.basis <- function(obj) { inherits(obj, "ParametricBasis") }

extractVariables <- function(.terms, data) {
	
	env <- environment(.terms)
	varnames <- sapply(attr(.terms, "variables") , deparse, width.cutoff = 500)[-1]
	variables <- eval(attr(.terms, "variables"), data, env)  
	names(variables) <- varnames
	variables
	
}

createEventTerms <- function(.terms, variables, resp, etab, facnames, expmat) {
	covar.names <- extractCovariates(.terms, variables, resp, etab)
	facterm <- do.call(EventTerm, lapply(facnames, function(fac) etab[[fac]]))
	var.terms <- if(length(covar.names) > 0) .extractVarTerms(covar.names, facnames, expmat, etab) else NULL
	c(facterm,var.terms)
}


afni_hrf <- function(..., conv="gamma", onsets=NULL, durations=NULL, drop.unused.levels=TRUE) {
	varlist <- list(...)
	
	anames <- parse(text=match.call())
	anames <- as.list(anames)[2:length(anames)]
	
	.onsets=onsets
	.durations=durations
	
	f <- function(onsets, TR, blocklens, blockids, durations, data) {
		
		if (!is.null(.onsets)) {
			onsets <- .onsets
		}		
		
		if (!is.null(.durations)) {
			durations <- .durations
		}
		
		evs <- lapply(1:length(varlist), function(i) {
					EV(eval(varlist[[i]], data, enclos=parent.frame()), as.character(anames[[i]]), onsets, blockids, durations)
				})
		
		eterm <- do.call(EventTerm, evs)
		
		if (is.character(conv)) {
			conv <- get_AFNI_HRF(conv)() 
		}
		
		AFNITerm(eterm, conv, blocklens, TR)					   		
	}
	
	class(f) <- c("function", "functerm")
	f
	
}
	




hrf <- function(..., conv="gamma", onsets=NULL, durations=NULL, granularity=.1, height=NULL, drop.unused.levels=TRUE, subset=NULL, labelPrefix=NULL) {
	varlist <- list(...)
	
	anames <- parse(text=match.call())
	anames <- as.list(anames)[2:length(anames)]
	
	if (!is.null(labelPrefix)) {
		anames <- paste(labelPrefix, "_", anames, sep="")
	}
		
	.onsets <- onsets
	.durations <- durations
	.subset <- subset
	
	f <- function(onsets, TR, blocklens, blockids, durations, data) {
		
		if (!is.null(.onsets)) {
			onsets <- .onsets
		}		
		
		if (!is.null(.durations)) {
			durations <- .durations
		}
		
		if (is.null(.subset)) {
			.subset <- rep(TRUE, length(onsets))				
		} else {		
			
			if (sum(.subset) < 1) {
				stop(paste("Error: provided subset contains no cases, aborting"))
			}
			#stopifnot(sum(.subset) > 1)
		}
		
		
		
		
		evs <- lapply(seq_along(varlist), function(i) {							
				EV(base:::eval(varlist[[i]], envir=data, enclos=parent.frame()), as.character(anames[[i]]), onsets, blockids, durations)
		})

		
		eterm <- do.call(EventTerm, evs)
		## hack until I can figure out how to add extra arg to do.call
		eterm@subset <- .subset
		## hack until I can figure out how to add extra arg to do.call
		
		.hrf <- if (is.character(conv)) {
			HRF(getHRF(conv))		
		} else if (is(conv, "HRF")){
			conv			
		} else if (is(conv, "function")) {
			HRF(conv)
		} else {
			stop(paste("illegal argument: conv ", conv))
		}

		ConvolvedTerm(eterm, .hrf, blocklens, TR, start=TR/2, granularity=granularity, height=height, drop.unused.levels=drop.unused.levels)					   		
	}
	
	class(f) <- c("function", "functerm")
	f
	
}

fromNuisanceList <- function(matlist, keep.names = FALSE, labelPrefix="Stim", blockPrefix="Run") {
	
	if (length(matlist) == 1) {
		return(matlist[[1]])
	}
	
	stopifnot(all(sapply(matlist, ncol) == ncol(matlist[[1]])))
		
	blocklens <- sapply(matlist, nrow)
	NC <- ncol(matlist[[1]])
	
	mat <- matrix(0, sum(blocklens), NC * length(matlist))
	
	
	gen.names <- function() {
		unlist(lapply(1:length(matlist), function(i) { 
					P1 <- paste(blockPrefix, i, "_", sep="")
					P2 <- paste(labelPrefix, "#", 1:NC, sep="")
					paste(P1, P2, sep="")
				}))
		
	}
	cnames <- if(keep.names) {
		unlist(lapply(matlist, colnames))
	} else {
		gen.names()
		
	}

	if (is.null(cnames)) {
		stop("matrix columns must have names if keep.names = TRUE ")
	}
	
	if (!length(unique(cnames)) == length(cnames)) {
		stop(paste("column names must be unique: found duplicate names in: ", paste(cnames, colapse=" ")))
	}
	
				
	
	for (i in seq_along(matlist)) {		
		rowstart <- sum(blocklens[1:i]) - blocklens[i] + 1
		rowend <- rowstart + blocklens[i] -1
		colstart <- (NC*(i-1)) + 1
		colend <- colstart + NC -1
		mat[rowstart:rowend, colstart:colend] <- as.matrix(matlist[[i]])	
	}
	
	colnames(mat) <- cnames
	mat
}

nuisance <- function(matlist, keep.names = FALSE, labelPrefix="Stim", blockPrefix="Run") {
	anames <- parse(text=match.call())
	varname <-  as.character(as.list(anames)[[2]])
	f <- function(onsets, TR, blocklens, blockids, durations=NULL, data) {
		#browser()
		if (is.matrix(matlist) || is.data.frame(matlist)) {
			mat <- as.matrix(matlist)
		} else {
			#mat <- fromNuisanceList(list(matlist), keep.names, labelPrefix, blockPrefix)
			mat <- fromNuisanceList(matlist, keep.names, labelPrefix, blockPrefix)
		}
		
		MatrixTerm(varname, mat, blocklens, TR)
	}
	
	class(f) <- c("function", "functerm")
	f
}

block <- function(blockterm, blocklens, TR) {
	.blockterm <- blockterm
	.blocklens <- blocklens
	.TR <- TR
	
	ret <- list(blockids=.blockterm, blocklens=.blocklens, TR=.TR)
	class(ret) <- c("list", "blockstruc")
	ret
}

baseline <- function(N=5, basis=c("bs", "poly", "ns")[1], name=paste("Baseline_", basis, "_", N, sep="")) {
	bfun <- switch(basis,
			bs=bs,
			ns=ns,
			poly=poly)
	f <- function(onsets=.onsets, TR, blocklens, blockids, durations=NULL, data) {
		ret <- lapply(blocklens, function(bl) cbind(rep(1, bl), bfun(seq(1, bl), N)))
		mat <- matrix(0, sum(blocklens), (N+1)*length(blocklens))
		for (i in seq_along(ret)) {
			rowstart <- sum(blocklens[1:i]) - blocklens[1] + 1
			rowend <- rowstart + blocklens[i] -1
			colstart <- (N+1)*(i-1) + 1
			colend <- colstart + N
			mat[rowstart:rowend, colstart:colend] <- ret[[i]]
		}
				
		cnames <- apply(expand.grid(paste("Basis",1:length(blocklens),sep=""), paste("Run", 1:(N+1), sep="")), 1, paste, collapse="_")
		colnames(mat) <- cnames
		MatrixTerm(name, mat, blocklens, TR)			
	}
	
	class(f) <- c("function", "functerm")
	f
		
}


fmrireg <- function(formula, data, hrf.fun=HRF.GAMMA, block=NULL, durations=0, blocklens=NULL, TR=NULL, drop.unused.levels=TRUE) {
	stopifnot(inherits(formula, "formula"))
	.terms <- .extractTerms(formula, data)
	resp <- attr(.terms, "response")
	stopifnot(resp > 0)
	
	if (is.null(resp)) {
		stop("need to provide onset vector on left side of formula, e.g. Onsets ~  a + b")
	}
	
	durations <- base:::eval(substitute(durations), envir=data, enclos=parent.frame())
	
	
	variables <- extractVariables(.terms, data)
	lhs <- variables[[resp]]
	
	rhs <- variables[(resp+1):length(variables)]
	vclass <- sapply(rhs, class)
	
	
	block <- Filter(function(x) is(x, "blockstruc"), rhs) 
	
	if (length(block) == 0) {
		block <- list(blockids=eval(substitute(block), envir=data), blocklens=blocklens, TR=TR)
		class(block) <- c("list", "blockstruc")
	} else {
		block <- block$block
	} 
	
	cfuncs <- Filter(function(x) is(x, "functerm"), rhs)
	
	
	cterms <- lapply(cfuncs, function(func) func(lhs, block$TR, block$blocklens, block$blockids, durations, data))
	
	eterms <- sapply(cterms, function(term) {
				if (inherits(term, "EventRegressionTerm")) term@eventTerm else NULL 
			})
	
	eterms <- eterms[!sapply(eterms, is.null)]
	
	## temporary restriction
	stopifnot(length(eterms) >= 1)
	
	
	evlist <- unique(unlist(lapply(eterms, function(et) et@events)))
		
    ##any(sapply(set, function(s) identical(elements(s, transformed=FALSE), elements(item, transformed=FALSE))))

	
	## temporary restriction
	#browser()
		
	eventModel <- new("EventModel", eventTable=EventTable(evlist), eventTerms=eterms, context=data, call.formula=formula, drop.unused.levels=drop.unused.levels) 
	new("ConvolvedModel", eventModel=eventModel, functionalTerms=cterms)
	
}



EventModel <- function(formula, data, block, durations=0, drop.unused.levels=TRUE) {
	if (missing(block)) {
		stop("error: need to provide \"block\" variable that divides the design in to sessions or runs")
	}
	
	
	
	.terms <- .extractTerms(formula, data)
	
	resp <- attr(.terms, "response")
	stopifnot(resp > 0)
	
	if (is.null(resp)) {
		stop("need to provide onset vector on left side of formula, e.g. Onsets ~  a + b")
	}
	
	block <- eval(substitute(block), envir=data)
	
	
	### may need way to include "intercept" i.e. column of all ones for main effects in the case where no factors are in model
	
	durations <- if (is.name(quote(durations))) {
				eval(substitute(durations), envir=data)
			} else {
				as.numeric(durations)
			}
	
	#figure out matrix of interactions
	variables <- extractVariables(.terms, data)
	lhs <- variables[[resp]]
	if (any(is.na(lhs))) {
		stop("onset vector cannot contain NA's")	
	}
	
	expmat <- attr(.terms, "factors")
	if (length(expmat) == 0 && length(variables) == 1) {
		#browser()
		### intercept only model
		Intercept <- EventFactor(rep(1, length(lhs)), "Intercept", lhs, blockids=block, durations=durations)
		iterm <- EventTerm(Intercept)
		etab <- EventTable(Intercept)
		facnames <- names(factors(etab))
		new("EventModel", eventTable=etab, eventTerms=list(iterm), context=data, call.formula=formula)
		
	} else {	
		etab <- .createEventTable(expmat, names(variables), variables, lhs, block, durations)
		facnames <- names(factors(etab))
	
		if (length(facnames) == 0) {
			stop("must provide at least one factor variable at present")
		}
	
		eventTerms <- createEventTerms(.terms, variables, resp, etab, facnames, expmat)
	
		new("EventModel", eventTable=etab, eventTerms=eventTerms, context=data, call.formula=formula, drop.unused.levels=drop.unused.levels) 
	}
}


setMethod("predict",  signature(object="EventModel"), 
		function(object,  newdata) {
			evs <- events(object@eventTable)
			evlist <- lapply(evs, function(ev) {
						predict(ev, newdata[[argname(ev)]])			
					})
			
			
			etab <- EventTable(evlist)
			facnames <- names(factors(etab))
			.terms <- .extractTerms(object@call.formula, object@context)
			orig.variables <- extractVariables(.terms, object@context)
			expmat <- attr(.terms, "factors")			
			resp <- attr(.terms, "response")
			eventTerms <- createEventTerms(.terms, orig.variables, resp, etab, facnames, expmat)
			
			### there is a difference between non-existent levels and unused levels??
			new("EventModel", eventTable=etab, eventTerms=eventTerms, context=newdata, call.formula=object@call.formula, drop.unused.levels=FALSE) 
			
		})



setMethod("initialize", "EventModel", function(.Object, eventTable, eventTerms, context,  call.formula, drop.unused.levels) {
			.Object <- callNextMethod()
			.Object@eventTable <- eventTable
			.Object@eventTerms <- eventTerms
			.Object@context <- context
			.Object@call.formula <- call.formula
			.Object@designMatrix <- as.matrix(do.call(cbind, lapply(eventTerms, designMatrix, drop.unused.levels=drop.unused.levels)))
			.Object
		})





.extractVarTerms <- function(varnames, facnames, expmat, eventTable) {
	
	varind <- row.names(expmat) %in% varnames
	varmat <- expmat[varind,, drop=F]
	
	term.count <- t(apply(expmat, 2, function(c1) {
						parts <- names(c1[c1>0])
						nvars <- sum(parts %in% varnames)
						nfacs <- sum(parts %in% facnames)
						c(nvars, nfacs)
					}))
	
	
	
	var.terms <- lapply(1:NROW(varmat), function(i) {
				idx <- which(varmat[i,] > 0)
				if (length(idx) > 0) {
					vnames <- unique(unlist((apply(expmat[,idx, drop=F], 2, function(vals) names(vals[vals>0])))))
					
					vterms <- varmat[i,idx,drop=FALSE]
					vterm.count <- term.count[idx,,drop=FALSE]
					
					
					if (all(vnames %in% varnames)) {
						## all terms are continuous variables, therefore need to keep main effects
						keep.idx <- which(vterm.count[,1] > 0)
						term.sets <- sapply(idx[keep.idx], function(id) {
									rownames(expmat[which(expmat[,id] == 1),,drop=FALSE])
								})
						
						term.sets      
						
					} else {
						
						## variable is crossed with factor, need only to keep higher order terms.      
						maxfac <- max(vterm.count[,2])
						keep.idx <- which(vterm.count[,2] == maxfac & vterm.count[,1] > 0)
						
						
						term.sets <- lapply(idx[keep.idx], function(id) {
									rownames(expmat[which(expmat[,id] == 1),,drop=FALSE])
								})
						
						
						term.sets                
						
					}
				}
			})
	
	# duplicate terms can occur in higher order interactions, need to find the unique terms
	var.terms <- unique(unlist(var.terms, recursive=F))
	
	# here we sort the terms according to their position in the contingency table (expmat)
	term.order <- order(sapply(var.terms, function(vt) {
						sum(sapply(vt, function(v) which(rownames(expmat) %in% v)))
					}))
	
	# a hack to handle the fact that we transform (name of) function call poly(time, 3) to an event name of poly_time
	name.map <- c(facnames, names(varnames))
	names(name.map) <- c(facnames, varnames)
	
	# here we create the higher order terms as "EventTerm" instances
	event.terms <- lapply(var.terms[term.order], function(terms) {
				evs <- eventTable[name.map[terms]]
				do.call(EventTerm, evs)
			})
	
	#should be the unique terms in the design
	event.terms
	
}



.createEventTable <- function(expmat, varnames, variables, lhs, block, durations) {
	
	nterms <- 0
	evs <- list()
	
	for (i in 1:NCOL(expmat)) {
		idx <- which(expmat[,i] == 1)
		if (length(idx) == 1) {
			nterms <- nterms + 1
			vname <- varnames[idx]      
			ev <- EV(variables[[varnames[idx]]], vname, lhs, block, durations) 
			evs[[nterms]] <- ev   
		}
	}
	
	
	EventTable(evs)
}



setMethod("formula",  signature(x = "EventModel"),
		function(x) {
			as.formula(paste("~", paste(sapply(x@eventTerms, varname), collapse=" + "), " -1"))
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

setMethod("blockids", signature(x="FMRIModel"),
		function(x) {
			blockids(x@eventModel)
		})


setMethod("[[", signature(x="EventModel", i = "character", j = "missing"),
		function(x, i, j) {
			eventTable(x)[[i]]
		})

setMethod("[[", signature(x="EventModel", i = "numeric", j = "missing"),
		function(x, i, j) {
			eventTable(x)[[i]]
		})

setMethod("numBlocks", signature(x="EventModel"),
		function(x) {
			numBlocks(x@eventTable)
		})

setMethod("conditions", signature(x="EventModel"),
		function(x) {
			unlist(lapply(model1@eventTerms, conditions))
		})

setMethod("levels", signature(x="EventModel"),
		function(x) {
			unlist(lapply(x@eventTerms, levels))
		})



setMethod("cells", signature(x="EventModel"),
		function(x) {
			cell.list <- lapply(x@eventTerms, cells)
			names(cell.list) <- terms(x)      
			cell.list 
		})


setMethod("colnames", signature(x="EventModel", do.NULL="missing", prefix="missing"),
		function(x) {
			colnames(designMatrix(x))         
		})


setMethod("terms", signature(x="EventModel"),
		function(x) {
			#sapply(x@eventTerms, varname)
			x@eventTerms
		})

setMethod("factors", signature(x="EventModel"),
		function(x) {
			names(factors(eventTable(x)))
		})

setMethod("designMatrix",  signature(x = "EventModel"),
		function(x) {
			x@designMatrix            
		})

setMethod("show", signature(object="EventModel"),
		function(object) {
			cat("formula: ", as.character(formula(object)), "\n")
			#print("conditions: \n"); print(conditions(object))
			#print("cells: ", cells(object), "\n")
			#cat("block lengths: ", blocklens(object), "\n")        
			
		})


