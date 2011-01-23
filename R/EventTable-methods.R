



EventTable <- function(...) {
	
	## is there a better way? how does data.frame do it?
	vlist <- list(...)
	
	if (is.list(unlist(vlist))) {
		vlist <- unlist(vlist)
	}
	
	if (!all(sapply(vlist, function(obj) is(obj, "EventVector")))) {
		stop("arguments to EventTable must inherit from class EventVector")
	}
	
	
	ret <- sapply(vlist, varname)
	#browser()
	
	if (length(ret) > 1) {
		if (length(unique(ret)) != length(ret)) {
			stop("event vectors names are not unique: ", print(ret))
		}
	}
	
	names(vlist) <- ret
	
	new("EventTable", eventVars=vlist)
}

setMethod("initialize", "EventTable", function(.Object, eventVars) {
			#eframe <- as.data.frame(do.call("cbind", lapply(eventVars, elements)))
			.Object@eventVars <- eventVars
			
			enames <- sapply(eventVars, varname)
			
			### this doesn't properly account for eventVectors that are matrices, e.g. poly(motion,2)
			### need way to derive sensible names from "call" like poly(motion,2) -- could be transformed to motion:poly1
			### motion:poly2, ...
			
			
			ons  <- onsets(eventVars[[1]])
			durs <- durations(eventVars[[1]])
			
			## onsets should(?) be the sorted union of the onsets of all evs to allow for non-factorial designs
			
			.Object@onsets <- onsets(eventVars[[1]])
			.Object@durations <- durations(eventVars[[1]])
			
			.Object@blockids <- blockids(eventVars[[1]])
			.Object
			
		})


setMethod("names", signature(x="EventTable"),
		function(x) {
			sapply(x@eventVars, varname)
		})

setMethod("show", signature(object="EventTable"),
		function(object) {
			cat(paste("variables:", paste(names(object), collapse=", ")))
			cat("\n")
			for (name in names(object)) {
				cat(paste(name, ": ", paste(levels(object[[name]]), collapse=" "), sep=""))
				cat("\n")
				
			}
			
		})      

setMethod("[[", signature(x="EventTable", i = "character", j = "missing"),
		function(x, i, j) {
			x@eventVars[[i]]
		})

setMethod("[[", signature(x="EventTable", i = "numeric", j = "missing"),
		function(x, i, j) {
			x@eventVars[[i]]
		})

setMethod("[", signature(x="EventTable", i = "numeric", j = "missing", drop="missing"),
		function(x, i, j) {
			x@eventVars[i]
		})

setMethod("[", signature(x="EventTable", i = "character", j = "missing", drop="missing"),
		function(x, i, j) {
			x@eventVars[i]
		})

setMethod("events", signature(x="EventTable"),
		function(x) {
			x@eventVars
			#evList <- lapply(x@eventVars, elements)
			#names(evList) <- enames
			#evList
		})

setMethod("onsets", signature(x="EventTable"),
		function(x) {
			onsets(x@eventVars[[1]])
		})

setMethod("durations", signature(x="EventTable"),
		function(x) {
			durations(x@eventVars[[1]])
		})

setMethod("blockids", signature(x="EventTable"),
		function(x) {
			blockids(x@eventVars[[1]])
		})

setMethod("numBlocks", signature(x="EventTable"),
		function(x) {
			length(unique(blockids(x@eventVars[[1]])))
		})

setMethod("conditions", signature(x="EventTable"),
		function(x) {
			facparts <- names(factors(x))
			
			levlist <- lapply(facparts, function(fpart) {
						els <- elements(x[[fpart]])
						levels(x[[fpart]])
					})
			
			condmat <- expand.grid(levlist)
			colnames(condmat) <- facparts
			
			cnames <- apply(as.matrix(condmat), 1, function(row) {
						nlevs <- apply(condmat, 2, function(v) length(unique(v)))
						olevs <- sapply(1:length(row), function(i) {
									if (nlevs[i] > 1) {
										return(paste("[", row[i], "]", sep=""))
									} else {
										return("")
									}
								})
						
						paste(names(row), olevs, sep="", collapse=":")
					})
			
			rownames(condmat) <- cnames
			condmat
		})



setMethod("factors", signature(x="EventTable"),
		function(x) {
			ret <- Filter(function(ev) is(ev,"EventFactor"), x@eventVars)
			lapply(ret, function(o) elements(o))
		})

setMethod("covariates", signature(x="EventTable"),
		function(x) {
			ret <- Filter(function(ev) is(ev,"EventVariable") || is(ev, "EventVariableSet"), x@eventVars)
			lapply(ret, function(o) elements(o))
		})

setMethod("subset", signature(x="EventTable"),
		function(x, f, rowsub=NULL) {
			if (length(f) != 2) {
				stop("subset: formula must be of the form : ~ fac1 + fac2 ...")
			}
			
			
			
			facstr <- deparse(f[[2]])
			facnames <- strsplit(facstr, ".\\+.")[[1]]
			faclist <- lapply(facnames, function(fn) x[[fn]])
			
			if (!is.null(rowsub)) {
				keep <- eval(substitute(rowsub), events(x))
				
			}
			
			EventTable(faclist)
			
		})


setMethod("as.data.frame", signature(x = "EventTable", row.names="missing", optional="missing"),
		function(x) as(x, "data.frame"))

setAs("EventTable", "data.frame",
		function(from) {
			### todo unravel any EventRegressorSets before converting to data.frame
			as.data.frame(events(from))
		})



setAs("data.frame", "EventTable",
		function(from) {  
			onsets <- from$onsets
			blocknum <- from$blocknum
			
			if ("durations" %in% names(from)) {
				durations <- from$durations
			} else {
				durations <- rep(length(onsets), 0)
			}
			
			durations <- from$duration
			varnames <- names(from)
			olist <- list()
			for (i in 1:length(varnames)) {
				vname <- varnames[i]
				if (toupper(vname) %in% c("ONSETS", "DURATIONS", "BLOCKNUM")) {
					next
				}
				
				vals <- from[[vname]]
				if (is.character(vals) || is.factor(vals)) {
					vals <- as.factor(vals)
					olist[[vname]] <- EventFactor(vals, vname, onsets, blocknum, durations)
				} else {
					vals <- as.numeric(vals)
					olist[[vname]] <- EventVariable(vals, vname, onsets, blocknum, durations)
				}
			}
			
			EventTable(olist)
			
			
		})








