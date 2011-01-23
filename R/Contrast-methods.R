


#match.list <- lapply(1:length(levs), function(i) {
#    lapply(cell.list, function(cell) {
#      ifelse(cell[[i]] == levs[i], TRUE, FALSE)
#    })
#  })

.conFromString <- function(cell.list, constr) {
  facs <- names(cell.list)

  levs <- strsplit(constr, ":")[[1]] 
  
} 

.computeWeights <- function(cmod, pat) {
	patmat <- patternMatch(cmod, pat)
	if (sum(patmat) <= 0) {
		stop(paste("failed to match any conditions with pattern: ", pat))
	}
	
	cvec <- ifelse(patmat, 1, 0)
	cvec/sum(cvec)	
}



setMethod("-", signature(e1="SumToOneContrast", e2="SumToOneContrast"),
	function(e1,e2) {
		stopifnot(all(rownames(e1) == rownames(e2)))
		res <- e1@.Data -(e2@.Data)
		new("SumToZeroContrast", res)
	})
	
setMethod("-", signature(e1="SumToZeroContrast", e2="SumToZeroContrast"),
	function(e1,e2) {
		stopifnot(all(rownames(e1) == rownames(e2)))
		res <- e1@.Data -(e2@.Data)
		new("SumToZeroContrast", res)
	})
	
setMethod("*", signature(e1="SumToZeroContrast", e2="SumToZeroContrast"),
	function(e1,e2) {
		stopifnot(all(rownames(e1) == rownames(e2)))
		res <- e1@.Data * (e2@.Data)
		new("SumToZeroContrast", res)
	})
	

setMethod("FContrast",signature(x="MatrixTerm"),
		function(x) {        
			els <- cells(x) 
			cmat <- diag(NROW(els))
			row.names(cmat) <- rownames(els)
			new("SumToOneContrast", cmat)
		})
	
.FContrast <- function(term) {
	celltab <- cells(term)
	
	M <- unlist(apply(replicate(ncol(celltab), c(0,1)), 2, function(col) list(col)), recursive=F)
	M <- do.call(expand.grid, M)
	
	###### was in the middle of fixing the factorial contrasts algorithm, see Mukerjee and Wu
	
	effect.vec <- apply(celltab, 2, function(vals) {
		rep(1, length(levels(as.factor(vals))))
	})
	
	diff.vec <- apply(celltab, 2, function(vals) {
		contr.sum(length(levels(as.factor(vals))))
	})
	
	main.effects <- lapply(seq_along(diff.vec), function(i) {
		kmat <- diff.vec[[i]]
		if (length(effect.vec) > 1) {
			for (j in 1:length(effect.vec)) {
				if (j == i) next
				kmat <- kronecker(kmat, effect.vec[[j]])
			}
		}
		kmat
	})
	
	names(main.effects) <- colnames(celltab)
	
	ieffects <- if (length(diff.vec) > 1) {		
		X <- unlist(lapply(2:length(diff.vec), function(i) {
			cmat <- combinations(length(diff.vec), i)
			ret <- lapply(1:NROW(cmat), function(j) {
				row.idx <- cmat[j,]
				dvecs <- lapply(row.idx, function(k) diff.vec[[k]])
				kmat <- Reduce(kronecker, dvecs)
				attr(kmat, "label") <- paste(names(diff.vec)[row.idx], collapse=":")
				kmat
			})
		}), recursive=F)
			
	    names(X) <- unlist(lapply(X, attr, "label"))	
		X
	}
	
	ret <- lapply(c(main.effects, ieffects), function(mat) {
		new("SumToZeroContrast", mat)
		
	})
}
	
.factorMap <- function(fac) {
	as.numeric(as.character(fac))
}	

makeWeights <- function(keepA, keepB=NULL) {
	weights <- matrix(0, length(keepA), 1)
	numA <- sum(keepA)
	
	weights[keepA,1] <- rep(1/numA, numA)
	if (!is.null(keepB)) {
		numB <- sum(keepB)
		weights[keepB,1] <- -rep(1/numB, numB)
	}
	
	weights
}


baselineContrast <- function(term, A, where=TRUE, splitBy=NULL) {
	stopifnot(inherits(term, "EventRegressionTerm"))
	term.cells <- cells(term)
	
	row.names(term.cells) <- longnames(term)
	
	count <- attr(term.cells, "count")		
	term.cells <- subset(term.cells, count > 0)
	
	keep <- eval(substitute(where), envir=term.cells, enclos=parent.frame())	
	keepA <- eval(substitute(A), envir=term.cells, enclos=parent.frame())
	
	split.fac <- eval(substitute(splitBy), envir=term.cells, enclos=parent.frame())
	weights <- if (!is.null(split.fac)) {
				split.levs <- levels(split.fac)
				do.call(cbind,lapply(split.levs, function(lev) {					
									keepC <- split.fac == lev 
									makeWeights(keep & keepA & keepC)
								}))
			} else { 
				makeWeights(keep & keepA) 
			}			
	
	
				
	row.names(weights) <- row.names(term.cells)
	new("SumToOneContrast", weights)
	
}



pairContrast <- function(term, A, B, where=TRUE, splitBy=NULL) {
	stopifnot(inherits(term, "EventRegressionTerm"))
	
	term.cells <- cells(term)
	row.names(term.cells) <- longnames(term)
	
	count <- attr(term.cells, "count")		
	term.cells <- subset(term.cells, count > 0)
	
	keep <- eval(substitute(where), envir=term.cells, enclos=parent.frame())	
	keepA <- eval(substitute(A), envir=term.cells, enclos=parent.frame())	
	keepB <- eval(substitute(B), envir=term.cells, enclos=parent.frame())
	
	split.fac <- eval(substitute(splitBy), envir=term.cells, enclos=parent.frame())
	weights <- if (!is.null(split.fac)) {
					split.levs <- levels(split.fac)
					do.call(cbind,lapply(split.levs, function(lev) {					
						keepC <- split.fac == lev 
						makeWeights(keep & keepA & keepC, keep & keepB & keepC)
					}))
				} else { 
					makeWeights(keep & keepA, keep & keepB) 
				}			
				

	row.names(weights) <- row.names(term.cells)
	new("SumToZeroContrast", weights)
		
}

polyContrast <- function(term, over, where=TRUE, degree=1, valueMap=NULL) {
	stopifnot(inherits(term, "EventRegressionTerm"))


	term.cells <- cells(term)
	row.names(term.cells) <- longnames(term)
	#count <- attr(term.cells, "count")		
	#term.cells <- subset(term.cells, count > 0)
		
	keep <- eval(substitute(where), envir=term.cells)	
	reduced.term.cells <- subset(term.cells, keep)
		
	vals <- if (is.null(valueMap)) {
		as.numeric(as.character(reduced.term.cells[[over]]))
	} else {
		valueMap[as.character(reduced.term.cells[[over]])]
	}
	
	weights <- matrix(0, NROW(term.cells), degree)	
	pvals <- poly(vals, degree=degree)
	row.names(weights) <- row.names(term.cells)
	colnames(weights) <- paste("poly", 1:degree, sep="")
	
	weights[keep, ] <- pvals
	new("SumToZeroContrast", weights)
	
}
	
		
Contrast <- function(conditionModel, a, b=NULL) {		
	Acon <- if (is.character(a)) {
		Amat <- .computeWeights(conditionModel, a)
		mat <- matrix(Amat, 1, length(Amat))
		colnames(mat) <- names(Amat)
		new("SumToOneContrast", mat)
	} else if (inherits(a, "Contrast")) {
		a
	} else {
		stop("argument \'a\' must be a of type character or Contrast.")
	}
		
	Bcon <- if (!is.null(b) && is.character(b)) {
		Bmat <- .computeWeights(conditionModel, b)	
		mat <- matrix(Bmat, 1, length(Bmat))	
		colnames(mat) <- names(Bmat)
		new("SumToOneContrast", mat)
	} else if (inherits(b, "Contrast")) {
		b
	} else {
		NULL
	}
	
	if (!is.null(Bcon)) {
		Acon - Bcon
	} else {
		Acon
	}
		
	
}
	

  


  
