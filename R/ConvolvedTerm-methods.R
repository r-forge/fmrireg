



MatrixTerm <- function(varname, mat, blocklens, TR) {
	stopifnot(is.matrix(mat))
	stopifnot(all(TR < blocklens))
	
	if (length(TR) == 1) {
		TR <- rep(TR, length(blocklens))
	}
	
	new("MatrixTerm", varname=varname, blocklens=blocklens, TR=TR, designMatrix=mat)
}


AFNITerm <- function(eventTerm, hrf, blocklens, TR) {
	globons <- .globalOnsets(eventTerm, blocklens, TR)
	stopifnot(is(hrf, "AFNI_HRF"))
	
	new("AFNITerm", varname=names(eventTerm), blocklens=blocklens,TR=TR, eventTerm=eventTerm, onsets=globons, hrf=hrf)
	
}

.globalOnsets <- function(eventTerm, blocklens, TR) {
	
	ids <- rep(1:length(unique(blockids(eventTerm))), table(blockids(eventTerm)))
	onsets <- onsets(eventTerm)
	globons <- sapply(1:length(onsets),function(i) {
				onsets[i] + (((ids[i]-1) * blocklens[ids[i]]) * TR)
			})
	
}

ConvolvableTerm <- function(eventTerm, hrf, blocklens, TR) {
	globons <- .globalOnsets(eventTerm, blocklens, TR)
	nimages <- sum(blocklens)
	if (!is(hrf, "HRF") && is.function(hrf)) {
		hrf <- HRF(hrf)
	} else {
		stop("hrf must be a function or of class HRF")
	}
	
	new("ConvolvableTerm", varname=names(eventTerm), blocklens=blocklens,TR=TR, eventTerm=eventTerm, onsets=globons)
}


convolveRegressor <- function(hrf, dmat, globons, durations, granularity, samples, height=NULL) {
	
	cond.names <- colnames(dmat)
	if (any(is.na(dmat)) || any(is.na(globons))) {
		cat("starting convolution of ", NCOL(dmat), "regressors", "\n")
		keep <- apply(dmat, 1, function(vals) all(!is.na(vals)))
		keep[is.na(globons)] <- FALSE
		ret <- lapply(1:NCOL(dmat), function(i) {
					cat("convolving regressor: ", cond.names[i], "\n")
					Regressor(hrf, globons[keep], amp=dmat[keep,i], durations=durations[keep], granularity=granularity, height=height)(samples)
					
				})	
	} else {
		ret <- lapply(1:NCOL(dmat), function(i) {
					cat("convolving regressor: ", cond.names[i], "\n")
					Regressor(hrf, globons, amp=dmat[,i], durations=durations, granularity=granularity, height=height)(samples)
				})
		
	}
	
	ret
	
}

ConvolvedTerm <- function(eventTerm, hrf, blocklens, TR, start=NULL, granularity=.1, height=NULL, drop.unused.levels=TRUE, use.multicore=TRUE) {
	#### onsets should not exceed block lengths -- this would cause onsets to run into next block
	
	ids <- rep(1:length(unique(blockids(eventTerm))), table(blockids(eventTerm)))
	onsets <- onsets(eventTerm)
	
	if (max(ids) > length(blocklens)) {
		stop("there are more block ids than block lengths, cannot compute global onsets")
	}
	
	globons <- sapply(1:length(onsets),function(i) {
				blocknum <- ids[i]
				offset <- (sum(blocklens[1:blocknum]) - blocklens[blocknum])*TR
				if (onsets[i] > blocklens[blocknum]*TR) {
					NA
				} else {
					onsets[i] + offset
				}
				
			})
	
	
	nimages <- sum(blocklens)
	
	if (is.null(start)) {
		start <- TR/2
	}
	
	samples <- seq(start, length.out=nimages, by=TR)
	
	dmat <- designMatrix(eventTerm, drop.unused.levels)
	
	split.dmat <- split(as.data.frame(dmat), factor(blockids(eventTerm)))
	split.ons <- split(globons, factor(blockids(eventTerm)))
	split.durations <- split(durations(eventTerm), factor(blockids(eventTerm)))
	split.samples <- split(samples, rep(1:length(blocklens), blocklens))
	
	do_norm <- function(vals) {	  
		if (normalize.height) {
			sf <- diff(range(vals))
			if (sf == 0) vals else vals/sf		 
		} else {
			vals
		}
	}
	
	
	if (use.multicore && require(multicore) && length(globons) > 100) {
		print("using mclapply to convolve regressors")
		.lapply <- mclapply
	} else {
		print("using lapply to convolve regressors")
		.lapply <- lapply
	}
	
	reglist <- .lapply(1:length(split.dmat), function(i) {
				do.call(cbind, convolveRegressor(hrf, split.dmat[[i]], split.ons[[i]], split.durations[[i]], granularity, split.samples[[i]], height))
			})
	
	
	
	convolvedMat <- do.call(rbind, reglist)
			
	new("ConvolvedTerm", varname=names(eventTerm), blocklens=blocklens,TR=TR, hrf=hrf, eventTerm=eventTerm, onsets=globons, samples=samples, designMatrix=convolvedMat)
	
}

setMethod("convolve", signature(x="AFNITerm"),
		function(x) {
			stop("cannot convolve an AFNITerm.  An AFNITerm is a \"dummy\" term that serves only as a specification for the AFNI 3dDeconvolve program")
		})


setMethod("convolve", signature(x="MatrixTerm"),
		function(x) {
			designMatrix(x)
		})


setMethod("convolve", signature(x="ConvolvedTerm"),
		function(x) {
			x@designMatrix
		})

setMethod("convolve", signature(x="ConvolvableTerm"),
		function(x, start=TR(x)/2, granularity=.1, normalize.height=FALSE) {
			nimages <- sum(blocklens(x))
			samples <- seq(start, length.out=nimages, by=TR(x))
			dmat <- designMatrix(x@eventTerm)
			globons <- onsets(x)
			durs <- durations(x)
			if (any(is.na(dmat))) {
				keep <- apply(dmat, 1, function(vals) all(!is.na(vals)))
				dmat <- dmat[keep,]
				globons <- globons[keep]
				durs <- durs[keep]
			}
			
			reglist <- lapply(1:NCOL(dmat), function(i) {
						Regressor(x@hrf, globons, amp=dmat[,i], durations=durs, granularity=granularity, normalize.height=normalize.height)(samples)
					})
			
			convolvedMat <- do.call(cbind, reglist)
			new("ConvolvedTerm", varname=names(eventTerm), blocklens=blocklens(x),TR=TR(x), eventTerm=x@eventTerm, onsets=globons, hrf=x@hrf, samples=samples, designMatrix=convolvedMat)		
			
		})


setMethod("isContinuous", signature=(x="MatrixTerm"),
	function(x) {
		TRUE
	})

setMethod("isContinuous", signature(x="EventRegressionTerm"),
		function(x) {
			isContinuous(x@eventTerm)
		})

setMethod("cells", signature(x="MatrixTerm"),
		function(x) {
			cnames <- colnames(x@designMatrix)
			vname <- varname(x)
			df1 <- data.frame(cnames)
			colnames(df1) <- vname
			rownames(df1) <- cnames
			df1
			
		})

setMethod("cells", signature(x="EventRegressionTerm"),
		function(x) {
			
			evtab <- eventTable(x@eventTerm)
			
			evset <- if (nbasis(x) > 1) {
						evlist <- c(list(factor(paste("basis", 1:nbasis(x), sep=""))), cells(x@eventTerm))
						names(evlist) <- c("basis", parentTerms(x@eventTerm))
						evlist <- lapply(evlist,levels)
						ret <- expand.grid(evlist, stringsAsFactors=TRUE)
						ret[c(2:length(ret), 1)]
					} else {
						cells(x@eventTerm)
					}
			
			
			
			strels <- apply(apply(evtab,2,trim), 1, paste, collapse=":")
			strlevs <- apply(apply(evset,2,trim), 1, paste, collapse=":")
			row.names(evset) <- strlevs
			counts <- rep(attr(cells(x@eventTerm), "count"), each=nbasis(x))
			
			ret <- evset[counts > 0,,drop=F]			
			attr(ret, "count") <- counts[counts > 0]
			
			ret
			
		})

setMethod("elements", signature(x="EventRegressionTerm"), 
		function(x) {
			eterm <- x@eventTerm
			els <- lapply(eterm@events, elements)
			
			if (nbasis(x) > 1) {
				els.basis <- as.data.frame(lapply(els, rep, each=nbasis(x)))
				els.basis$basis <- paste("basis", rep(1:nbasis(x), length.out=nbasis(x)*length(els[[1]])), sep="")
				names(els.basis) <- c(unlist(parentTerms(eterm)), "basis")	
				as.data.frame(els.basis)	
			} else {
				names(els) <- unlist(parentTerms(eterm))
				as.data.frame(els)
			}	
		})

setMethod("conditions", signature(x="MatrixTerm"),
		function(x, drop=TRUE) {
			colnames(designMatrix(x))
		})

setMethod("conditions", signature(x="EventRegressionTerm"),
		function(x, drop=TRUE) {
			ret <- conditions(x@eventTerm, drop=drop)
			if (nbasis(x) > 1) {
				blevs <- paste("[", 1:nbasis(x), "]", sep="")
				unlist(lapply(ret, function(prefix) paste(prefix, ":basis", blevs, sep="")))
			} else {
				ret
			}           
		})

setMethod("levels", signature(x="EventRegressionTerm"),
		function(x) {
			ret <- levels(x@eventTerm)
			if (nbasis(x) > 1) {
				blevs <- 1:nbasis(x)
				unlist(lapply(ret, function(prefix) paste(prefix, ":basis", blevs, sep="")))
			} else {
				ret
			}
			
		})

setMethod("shortnames", signature(x="RegressionTerm"),
		function(x,exclude.basis=FALSE) {
			# ignores exclude.basis
			row.names(cells(x))	
			# ignores exclude.basis
		})

setMethod("longnames", signature(x="RegressionTerm"),
		function(x,exclude.basis=FALSE) {
			# ignores exclude.basis
			term.cells <- cells(x)
			# ignores exclude.basis
			apply(sapply(1:ncol(term.cells), 
							function(i) {
								paste(names(term.cells)[i], "#", term.cells[,i], sep="")
							}), 1, paste, collapse=":")
		})



setMethod("shortnames", signature(x="AFNITerm"),
		function(x,exclude.basis=FALSE) {
			ret <- row.names(cells(x@eventTerm))
			if (nbasis(x) > 1 && !exclude.basis) {
				blevs <- 1:nbasis(x)
				unlist(lapply(ret, function(prefix) paste(prefix, ":basis", blevs, sep="")))
			} else {
				ret
			}
		})

#setMethod("[[", signature(x="EventRegressionTerm", i="character", j="missing"),
#	function(x,i,j) {
#		x@eventTerm[[i]]
#	})

#setMethod("[[", signature(x="EventRegressionTerm", i="numeric", j="missing"),
#	function(x,i,j) {
#		x@eventTerm[[i]]
#	})


setMethod("nlevels", signature(x="EventRegressionTerm"),
		function(x) {
			nlevels(x@eventTerm)*nbasis(x)
			
		})

setMethod("nbasis", signature(x="ConvolvedTerm"),
		function(x) {
			nbasis(x@hrf)
		})

setMethod("nbasis", signature(x="AFNITerm"),
		function(x) {
			nbasis(x@hrf)
		})

setMethod("isConvolved", signature(x="AFNITerm"),
		function(x) {
			FALSE
		})

setMethod("isConvolved", signature(x="AFNITerm"),
		function(x) {
			TRUE
		})

setMethod("TR", signature(x="EventRegressionTerm"),
		function(x) {
			x@TR
		})

setMethod("blocklens", signature(x="EventRegressionTerm"),
		function(x) {
			x@blocklens
		})

setMethod("onsets", signature(x="EventRegressionTerm"),
		function(x) {
			x@onsets
		})

setMethod("splitOnsets", signature(x="EventRegressionTerm"), 
		function(x, global=FALSE, blockDurations=x@blocklens*x@TR) {
			### what about basis functions? currently ignores them ....
			splitOnsets(x@eventTerm, global, blockDurations)
		})

setMethod("durations", signature(x="EventRegressionTerm"),
		function(x) {
			durations(x@eventTerm)
		})

setMethod("designMatrix", signature(x="ConvolvedTerm"),
		function(x) {     
			## hack to add colnames to design matrix
			ret <- x@designMatrix    
			
			####### TEMPORARY HACK !!!
			n1 <- conditions(x, drop=FALSE)
			n2 <- conditions(x, drop=TRUE)
			####### TEMPORARY HACK !!!
			if (ncol(ret) == length(n1))
				colnames(ret) <- n1
			else {
				colnames(ret) <- n2
			}
			ret
		})

setMethod("designMatrix", signature(x="MatrixTerm"),
		function(x) {        
			x@designMatrix                      
		})


setMethod("designMatrix", signature(x="AFNITerm"),
		function(x) {        
			designMatrix(x@eventTerm)                        
		})

setMethod("varname", signature(x="RegressionTerm"),
		function(x) {
			x@varname
		})

setMethod("blocklens", signature(x="RegressionTerm"),
		function(x) {
			x@blocklens
		})

setMethod("plot", signature(x="RegressionTerm"),
		function(x, splitByBlock=TRUE) {
			if (splitByBlock) {
				blockFactor <- factor(rep(1:length(blocklens(x)), blocklens(x)))
				mats <- split(as.data.frame(designMatrix(x)), blockFactor)
				
				NC <- 3
				NR <- max(length(mats)/NC,1)
				par(mfrow=c(NR, NC))
				lapply(mats, function(mat) ts.plot(ts(mat), col=1:NCOL(mat)))
			} else {
				ts.plot(ts(designMatrix(x)), col=1:NCOL(mat))		
			}
		})
		
setMethod("show", signature(object="MatrixTerm"),
		function(object) {
			cat("MatrixTerm", "\n")
			cat(" ", "Term Name: ", varname(object), "\n")
			cat(" ", "Num Columns", ncol(designMatrix(object)), "\n")
			
		})

setMethod("show", signature(object="AFNITerm"),
		function(object) {
			cat("AFNITerm", "\n")
			cat(" ", "Term Name: ", names(object@eventTerm), "\n")
			cat(" ", "HRF: ", object@hrf@name, "\n")
			cat(" ", "nbasis: ", object@hrf@nbasis, "\n")
			#print(object@hrf)
		})

setMethod("show", signature(object="EventRegressionTerm"),
		function(object) {
			cat("EventRegressionTerm", "\n")
			cat(" ", "Term Name: ", names(object@eventTerm), "\n")
			cat(" ", "HRF: ", object@hrf@name, "\n")
			cat(" ", "nbasis: ", object@hrf@nbasis, "\n")
			
		})


