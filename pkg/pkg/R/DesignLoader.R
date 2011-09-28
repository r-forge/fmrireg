library(yaml)

source("rlist.SOURCE")



.blockPath <- function(basepath, inode) {
  if (!is.null(inode$path)) {
    if (substr(inode$path, 1,1) == "/") {
      return(path.expand(inode$path))
    }
  }
      
  blockPath <- ifelse(!is.null(inode$path),
                        paste(basepath, "/", inode$path, sep=""),
                        basepath)

  
}



.parseImageData <- function(basepath, imnode) {
  blockPath <- .blockPath(basepath, imnode)

  blocksizes <- imnode[["block-sizes"]]
  blocks <- lapply(seq_along(imnode$blocks), function(i) {   
    new("DataBlock", path=blockPath, filelist=imnode$blocks[[i]], blocknum=i, blocksize=blocksizes[i])
  })

  datalist <- DataBlockList(blocks)
}

.parseTableData <- function(basepath, tabnode) {
  #browser()
  tabfilenode <- tabnode[["table-files"]]
  tabpath <- .blockPath(basepath, tabfilenode)

  tabnames <- paste(tabpath, "/", tabfilenode$files, sep="")

  fulltable <- do.call("rbind", lapply(tabnames, read.table, header=T))

  colnode <- tabnode[["column-header"]]
  blockfactor <- colnode[["block-factor"]]
  onsetvar <- colnode[["onsets"]]
  durvar <- colnode[["durations"]]
  factors <<- colnode[["factors"]]
  covars <<- colnode[["covariates"]]
  # will make ordinal even if there is a "skip" (i.e. 1, 2, 4, 5, 6) in block factor
  # this will work provided there is a one-to-one mapping between the ordered set of blocks and the list of data blocks.
  # if there is more behavioral data than image data, the we need an "exclude-blocks" tag which will filter out certain runs.
  
  blocks <- sort(unique(fulltable[[blockfactor]]))
  
  if (is.null(onsetvar)) {
    onsetvar = "onsets"
    
  }

  if (is.null(durvar)) {
    durvar = "durations"
  } else if (is.numeric(durvar)) {
    # durvar is a constant value
    fulltable$durations <- rep(durvar, NROW(fulltable))
    durvar <- "durations"
  }

      

  .buildVariable <- function(varnames, blocks, builder) {
    lapply(varnames, function(var) {
      v <<- var
      do.call("merge", lapply(blocks, function(block) {
        idx <- which(fulltable[[blockfactor]] == block)
        vals <- fulltable[[v]][idx]
        onsets <- fulltable[[onsetvar]][idx]
        durations <- fulltable[[durvar]][idx]
        builder(vals, v, onsets, blocknum=block, durations)
      }))
    })
  }

  

  evfacs <- .buildVariable(factors, blocks, function(vals, name, onsets, blocknum, durations) {
    EventFactor(vals, name, onsets, blocknum, durations)
  })

  evvars <- .buildVariable(covars, blocks, function(vals, name, onsets, blocknum, durations) {
    EventVariable(vals, name, onsets, blocknum, durations)
  })           

                                        
  etable <- do.call("EventTable", c(evfacs, evvars))

  
}
  

loadDesign <- function(filename) {
  
  des <- yaml.load(paste(readLines(filename), collapse="\n"))

  Design <- des$Design
  
  basepath <- path.expand(ifelse(!is.null(Design[["base-path"]]), Design[["base-path"]], getwd()))


  imnode <- Design[["image-data"]]
  datalist <- .parseImageData(basepath, imnode)

  tabnode <- Design[["table-data"]]
  etable <- .parseTableData(basepath, tabnode)

  TR <- imnode[["TR"]]
  if (is.null(TR)) {
    warning("no TR set, defaulting to 2 s")
    TR = 2
  }

  FMRIDesign(etable, datalist, TR)
  

}





  
