

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

                        
setClass("HRF", representation(name="character", nbasis="integer"), contains="function")

setClass("EventHRF", representation(onset="numeric", amplitude="numeric", duration="numeric"), contains="HRF")

setClass("Regressor",
         representation(cells="list", onsets="numeric", hrf="HRF", amplitudes="numeric", durations="numeric"),
         contains="function")

setClass("RegressorSet",
         representation(labels="character", onsets="list", hrf="HRF"),
         contains="function")
         
         
         
setClass("EventFactor",
         representation(value="factor"), contains="EventVector"
         )



setClass("EventVariable",
         representation(value="numeric"), contains="EventVector"
         )

setClass("EventVariableSet",
         representation(value="matrix"), contains="EventVector"
         )


setClass("EventVectorList",
         representation(name="character",
                        eventlist="list"))

setClass("EventTable",
         representation(eventVars="list", onsets="numeric", durations="numeric", blockids="integer"))

setClass("EventModel",
         representation(eventTable="EventTable", formula="formula", context="list", designMatrix="matrix"))



setClass("FMRIDesign",
         representation(eventModel="EventModel",
                        dataSource="DataSource",
                        HRF="HRF",
                        regressors="list"))

setClass("AFNIStim",
         representation(events="EventVector", basis="character", basis.params="list"))

                       
                       
setClass("AFNICommand", representation(design="FMRIDesign", workingDir="character", mask="character", options="list"))


