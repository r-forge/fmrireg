
#### generics related to DataBlock
setGeneric("path",     function(x) standardGeneric("path"))
#setGeneric("fullpath",     function(x) standardGeneric("fullpath"))
setGeneric("filelist",     function(x, full.names) standardGeneric("filelist"))
setGeneric("samples", function(x, global, start) standardGeneric("samples"))

setGeneric("params", function(x) standardGeneric("params"))

setGeneric("onsets",  function(x) standardGeneric("onsets"))
setGeneric("globalOnsets",  function(x, blockDurations) standardGeneric("globalOnsets"))
setGeneric("blockids",  function(x) standardGeneric("blockids"))
setGeneric("blocklens", function(x) standardGeneric("blocklens"))
setGeneric("numBlocks", function(x) standardGeneric("numBlocks"))

#setGeneric("blocksize",  function(x) standardGeneric("blocksize"))
setGeneric("amplitudes", function(x) standardGeneric("amplitudes"))
setGeneric("durations",  function(x) standardGeneric("durations"))
setGeneric("varname",  function(x) standardGeneric("varname"))
setGeneric("argname", function(x) standardGeneric("argname"))
setGeneric("levels",  function(x) standardGeneric("levels"))

setGeneric("nlevels",  function(x) standardGeneric("nlevels"))
setGeneric("elements",  function(x, ...) standardGeneric("elements"))
setGeneric("cells", function(x, ...) standardGeneric("cells"))
setGeneric("conditions", function(x, ...) standardGeneric("conditions"))
setGeneric("events", function(x) standardGeneric("events"))
setGeneric("parentTerms", function(x) standardGeneric("parentTerms"))
           
#setGeneric("concatenate",  function(x, y) standardGeneric("concatenate"))

setGeneric("convolve", function(x, ...) standardGeneric("convolve"))
           
setGeneric("buildArgs", function(x) standardGeneric("buildArgs"))
setGeneric("stimlabels", function(x) standardGeneric("stimlabels"))
setGeneric("stimfiles", function(x) standardGeneric("stimfiles"))
setGeneric("writeStims", function(x, dir) standardGeneric("writeStims"))
setGeneric("basis", function(x) standardGeneric("basis"))
setGeneric("nbasis", function(x) standardGeneric("nbasis"))

setGeneric("isContinuous", function(x) standardGeneric("isContinuous"))
setGeneric("isCategorical", function(x) standardGeneric("isCategorical"))
setGeneric("isConvolved", function(x) standardGeneric("isConvolved"))

setGeneric("factors", function(x) standardGeneric("factors"))
setGeneric("splitOnsets", function(x, ...) standardGeneric("splitOnsets"))
setGeneric("covariates", function(x) standardGeneric("covariates"))
setGeneric("subset", function(x,...) standardGeneric("subset"))
setGeneric("TR", function(x,...) standardGeneric("TR"))

setGeneric("formula", function(x, ...) standardGeneric("formula"))
setGeneric("terms", function(x, ...) standardGeneric("terms"))
setGeneric("designMatrix", function(x, ...) standardGeneric("designMatrix"))
setGeneric("eventTable", function(x) standardGeneric("eventTable"))
setGeneric("eventVars", function(x) standardGeneric("eventVars"))
setGeneric("shortnames", function(x, ...) standardGeneric("shortnames"))
setGeneric("longnames", function(x, ...) standardGeneric("longnames"))

setGeneric("addOption<-", function(x, optstr) standardGeneric("addOption<-"))

setGeneric("getOptions", function(x) standardGeneric("getOptions"))
setGeneric("workingDir", function(x) standardGeneric("workingDir"))
setGeneric("buildAFNIStims", function(x,...) standardGeneric("buildAFNIStims"))
setGeneric("buildCommandSwitch", function(x,k,type) standardGeneric("buildCommandSwitch"))
setGeneric("buildGLT", function(con) standardGeneric("buildGLT"))

setGeneric("writeAFNIStim", function(x, dir, ...) standardGeneric("writeAFNIStim"))

setGeneric("build", function(x) standardGeneric("build"))
setGeneric("glts", function(x) standardGeneric("glts"))

setGeneric("isMatch", function(x, y, ...) standardGeneric("isMatch"))
setGeneric("patternMatch", function(x,y,...) standardGeneric("patternMatch"))
setGeneric("shortform", function(x) standardGeneric("shortform"))
setGeneric("longform", function(x) standardGeneric("longform"))


setGeneric("FContrast", function(x, ...) standardGeneric("FContrast"))




           
