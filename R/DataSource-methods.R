
DataSource <- function(path, filelist, blocklens, TR) {
  if (length(blocklens) != length(filelist)) {
    blocklens <- rep(blocklens, length.out=length(filelist))
  }

  if (length(path) == 1) {
	  path <- rep(path, length(filelist))
  } 
  
  stopifnot(length(path) == length(filelist))
  
  new("DataSource", path=path, filelist=basename(filelist), blocklens=blocklens, TR=TR)
}

  

setMethod("path", signature(x = "DataSource"),
    function(x) x@path)


setMethod("filelist", signature(x = "DataSource", full.names="logical"),
    function(x, full.names) {
      if (full.names) {
        paste(path(x), .Platform$file.sep, x@filelist, sep="")
      } else {
        x@filelist
      }
    })

setMethod("filelist", signature(x = "DataSource", full.names="missing"),
    function(x) {
      callGeneric(x, FALSE)
      
    })




setMethod("length", signature(x="DataSource"),
          function(x) {
            sum(blocklens(x))
          })


setMethod("TR", signature(x = "DataSource"),
    function(x) x@TR)

setMethod("samples", signature(x="DataSource", global="missing", start="missing"),
          function(x) {
            nimages <- length(x)
            seq(TR(x)/2, length.out=nimages, by=TR(x))
          })

setMethod("samples", signature(x="DataSource", global="missing", start="numeric"),
          function(x, start) {
            callGeneric(x, global=TRUE, start=TR(x)/2)
          })


setMethod("samples", signature(x="DataSource", global="logical", start="missing"),
          function(x, global) {
            callGeneric(x, global, start=TR(x))
          })

setMethod("samples", signature(x="DataSource", global="logical", start="numeric"),
          function(x, global, start) {
            if (global) {
              nimages <- length(x)
              seq(start, length.out=nimages, by=TR(x))
            } else {
              lapply(blocklens(x), function(bl) {
                seq(start, length.out=bl, by=TR(x))
              })
            }
          })




            
setMethod("blocklens", signature(x = "DataSource"),
    function(x) x@blocklens)

setMethod("c", signature(x="DataSource"),
          function(x, ..., recursive=F) {
            varargs <- list(...)
            
            bind <- function(x, y) {           
              if (!inherits(y, "DataSource")) {
                stop(paste("cannot append ", class(y), " to a DataSource"))
              }
              
              DataSource(c(path(x), path(y)), c(filelist(x), filelist(y)), c(blocklens(x), blocklens(y)), c(TR(x), TR(y)))
            }

            ret <- x

            for (i in 1:length(varargs)) {
              ret <- bind(ret, varargs[[i]])
            }

            ret          
              
          })



