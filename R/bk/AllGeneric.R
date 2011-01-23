
#### generics related to DataBlock
setGeneric("path",     function(x) standardGeneric("path"))
#setGeneric("fullpath",     function(x) standardGeneric("fullpath"))
setGeneric("filelist",     function(x, full.names) standardGeneric("filelist"))
setGeneric("samples", function(x, global, start) standardGeneric("samples"))

setGeneric("onsets",  function(x) standardGeneric("onsets"))
setGeneric("globalOnsets",  function(x, blocklens) standardGeneric("globalOnsets"))
setGeneric("blockids",  function(x) standardGeneric("blockids"))
setGeneric("blocklens", function(x) standardGeneric("blocklens"))

#setGeneric("blocksize",  function(x) standardGeneric("blocksize"))
setGeneric("amplitudes", function(x) standardGeneric("amplitudes"))
setGeneric("durations",  function(x) standardGeneric("durations"))
setGeneric("varname",  function(x) standardGeneric("varname"))
setGeneric("levels",  function(x) standardGeneric("levels"))
setGeneric("elements",  function(x) standardGeneric("elements"))
setGeneric("cells", function(x) standardGeneric("cells"))
setGeneric("conditions", function(x) standardGeneric("conditions"))
setGeneric("events", function(x) standardGeneric("events"))
           
#setGeneric("concatenate",  function(x, y) standardGeneric("concatenate"))

setGeneric("convolve", function(x, y, blocklens) standardGeneric("convolve"))
           
setGeneric("buildArgs", function(x) standardGeneric("buildArgs"))
setGeneric("stimlabels", function(x) standardGeneric("stimlabels"))
setGeneric("stimfiles", function(x) standardGeneric("stimfiles"))
setGeneric("writeStims", function(x, dir) standardGeneric("writeStims"))
setGeneric("basis", function(x) standardGeneric("basis"))
setGeneric("nbasis", function(x) standardGeneric("nbasis"))


setGeneric("factors", function(x) standardGeneric("factors"))
setGeneric("covariates", function(x) standardGeneric("covariates"))
setGeneric("subset", function(x,...) standardGeneric("subset"))
setGeneric("TR", function(x,...) standardGeneric("TR"))

setGeneric("formula", function(x, ...) standardGeneric("formula"))
setGeneric("terms", function(x, ...) standardGeneric("terms"))
setGeneric("designMatrix", function(x, ...) standardGeneric("designMatrix"))
setGeneric("eventTable", function(x) standardGeneric("eventTable"))
setGeneric("eventVars", function(x) standardGeneric("eventVars"))
setGeneric("shortnames", function(x) standardGeneric("shortnames"))

setGeneric("addOption<-", function(x, optstr) standardGeneric("addOption<-"))

setGeneric("getOptions", function(x) standardGeneric("getOptions"))
setGeneric("workingDir", function(x) standardGeneric("workingDir"))

setGeneric("build", function(x) standardGeneric("build"))
           
