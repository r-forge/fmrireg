
setOldClass("bs")
setOldClass("basis")
setOldClass("poly")
setOldClass("ns")
setOldClass("scaled")

setClassUnion("predictable",  c("basis",  "bs",  "poly",  "ns", "scaled"))
setClassUnion("matrixOrbasis", c("matrix",  "basis",  "bs",  "poly",  "ns", "scaled"))


setClass("BlockStructure",
		representation=representation(blocklens="numeric",
				TR="numeric"))

setClass("DataSource",
         representation(path="character",
                        filelist="character",
                        blocklens="numeric",
                        TR="numeric")
         )


setClass("EventVector",
         representation(varname="character",
                        blockids="integer",
                        onsets="numeric",        
                        durations="numeric"),
         contains="VIRTUAL")


setClass("ParametricBasis", representation(x="numeric", y="matrixOrbasis", fun="character", argname="character", name="character"), contains="VIRTUAL")

setClass("Poly", representation(degree="integer"), contains="ParametricBasis")

Poly <- function(x, degree) {
	mc <- match.call()
	pres <- poly(x, degree)
	n <- paste("poly", "_", as.character(mc[["x"]]), "_", degree, sep="")
	new("Poly", x=x, y=pres, fun="poly", argname=as.character(mc[["x"]]), name=n, degree=as.integer(degree))
}

setMethod("predict", signature(object="Poly"),
          function(object, newvals) {
			predict(object@y, newvals)
          
          })


                        
setClass("HRF", representation(name="character", nbasis="integer"), contains="function")

setClass("AFNI_HRF", representation(name="character", nbasis="integer", params="list"))


setClass("EventHRF", representation(onset="numeric", amplitude="numeric", duration="numeric", span="numeric"), contains="HRF")

setClass("Regressor",
         representation(onsets="numeric", hrf="HRF", amplitudes="numeric", durations="numeric"),
         contains="function")

              
setClass("EventFactor",
         representation(value="factor"), contains="EventVector")

setClass("EventBasis", 
		representation(value="ParametricBasis"),  contains="EventVector")
		
setClass("EventVariable",
         representation(value="numeric"), contains="EventVector")
         

setClass("EventVariableSet",
         representation(value="matrixOrbasis"), contains="EventVector")
         

setClass("EventTerm",
         representation(events="list"),
         validity=function(object) {
           if (!all(sapply(object@events, inherits, "EventVector"))) {
             "all events must inherits from EventVector"
           } #else if (!Reduce(function(x,y) all.equal(onsets(x), onsets(y)), object@events)) {
             #"all events must have same onsets"
             #}           
         })


									  
setClass("RegressionTerm", representation(varname="character", blocklens="numeric", TR="numeric"), contains="VIRTUAL")

setClass("MatrixTerm", representation(designMatrix="matrix"), contains="RegressionTerm")

setClass("EventRegressionTerm", representation(eventTerm="EventTerm", onsets="numeric"), contains="RegressionTerm")

setClass("AFNITerm", representation(hrf="AFNI_HRF"), contains="EventRegressionTerm")


setClass("ConvolvableTerm",
	representation(hrf="HRF"),contains="EventRegressionTerm")

setClass("ConvolvedTerm",
         representation(samples="numeric", designMatrix="matrix"),contains="ConvolvableTerm",
         validity=function(object) {
           #if (!(nbasis(object@hrf)*length(conditions(object@eventTerm)) == NCOL(object@designMatrix))) {
           #  "number of basis function != number of columns in design matrix"
           if (NROW(object@designMatrix) != length(object@samples)) {
	          paste("number of samples != rows in design matrix: ", "SAMPLES: ", length(object@samples), "ROWS:", NROW(object@designMatrix))				
		   }
         })



setClass("EventVectorList",
         representation(name="character",
                        eventlist="list"))

setClass("EventTable",
         representation(eventVars="list", onsets="numeric", durations="numeric", blockids="integer"))

setClass("EventModel",
         representation(eventTable="EventTable", eventTerms="list", context="list", call.formula="formula", designMatrix="matrix"))

setClass("FMRIModel", representation(eventModel="EventModel", functionalTerms="list"), contains="VIRTUAL")

setClass("ConvolvedModel", contains="FMRIModel",
	validity=function(object) {
		#ret <- sapply(object@functionalTerms, function(term) {
		#	is(term, "ConvolvedTerm")
		#})
		
		#if (!all(ret)) {
		#	"all functional terms in ConvolvedModel must be of class ConvolvedTerm"
		#}
	})

setClass("AFNIModel", contains="FMRIModel",
	validity=function(object) {
		
		ret <- sapply(object@functionalTerms, function(term) { is(term, "AFNITerm") })
		if (sum(ret) < 1) {
			"at least one term in an AFNIModel must be an AFNITerm"
		}
	})
		
		

setClass("FMRIDesign",
         representation(fmriModel="FMRIModel",
                        dataSource="DataSource"))
         

setClass("AFNIStim",
         representation(label="character", fileName="character"), contains="VIRTUAL")

setClass("AFNIStimFile",  representation(values="numeric"), contains="AFNIStim")

setClass("AFNIStimTimes", representation(hrf="AFNI_HRF", onsets="numeric", iresp="logical", tr.times="numeric"), contains="AFNIStim")

	                     
setClass("AFNICommand", representation(design="FMRIDesign", workingDir="character", mask="character", options="list", glts="list"))


setClass("Contrast", contains=c("matrix"))

setClass("SumToZeroContrast", contains="Contrast",
	validity=function(object) {
		if (!all(round(apply(object,2, sum),digits=4) == 0)) {
			"all columns of a SumToZeroContrast must sum to 0."
		} else {		
			TRUE
		}
	})
	
setClass("SumToOneContrast", contains="Contrast",
		validity=function(object) {
			if (!all(round(apply(object,2, sum),digits=4) == 1)) {
				"all columns of a SumToOneContrast must sum to 1."
			} else {		
				TRUE
			}
		})
	
setMethod("initialize", "Contrast", function(.Object, ...) { 
	value <- callNextMethod() 
	validObject(value) 
	value 
}) 
	

setClass("Condition", representation(names="character", levels="character", types="character"),
         validity=function(object) {
           if (length(object@names) != length(object@levels)) {
             "length of names must equal length of levels"
           } else if (length(object@names) != length(object@types)) {
             "length of names must equal length of types"
           } 
         })


setClass("ConditionList", representation(conditions="list", table="data.frame"),
         validity=function(object) {
           if (!all(sapply(object@conditions, inherits, "Condition"))) {
             "all items in conditions list must inherit from Condition"
           }
          
         })
           
setClass("ConditionModel", representation(conditionList="list"),
         validity=function(object) {
           if (!all(sapply(object@conditionList, inherits, "ConditionList"))) {
             "all items in list argument must inherit from ConditionList"
           }
         })
           
           
