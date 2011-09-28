
EventVectorList <- function(vlist) {
  #browser()
  if (!all(sapply(vlist, function(obj) is(obj, "EventVector")))) {
    stop("EventVectorList can only contains EventVector instances")
  }

 
  ret <- sapply(vlist, name)
  if (length(ret) > 1) {
    if (!all(ret[1] == ret[2:length(ret)])) {
      stop("each EventVector component must have the same \"name\" attribute")
    }
  }

  name <-name(vlist[[1]])


  new("EventVectorList", name=name, eventlist=as.list(vlist))

  
}


setMethod("filelist", signature(x = "DataBlockList"),
    function(x) {
      sapply(x@list, function(obj) obj@filelist)
    })


setMethod("path", signature(x = "DataBlockList"),
    function(x) {
      sapply(x@list, function(obj) obj@path)
    })



  
          


          
