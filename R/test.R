
source("rlist.SOURCE")

fpath <- "~/data1/Dots/KAR/Data_RS/epi"
flist <- list.files(fpath, "^t107.*nii.gz")



destable <- read.table("~/data1/Dots/KAR/Data_RS/behavior/alldata.txt", header=T)
destable$image <- rep(flist, each=table(destable$run))
destable$main <- rep(1, NROW(destable))

trainruns <- as.integer(seq(1,48, length.out=25))
traintable <- drop.levels(subset(destable, run %in% trainruns))

emodel <- EventModel(Onset ~ direction*motion*time, data=traintable, block=run)
#emodel <- EventModel(Onset ~ direction*poly(motion,2), data=traintable, block=run)
dsource <- DataSource(fpath, unique(traintable$image), 302, 1.37)
fmridesign <- FMRIDesign(emodel, dsource,  HRF.GAMMA)
aglm <- afni.glm(fmridesign, "junk2", mask="/Users/brad/data1/Dots/KAR/Data_RS/epi/global_mask.nii",options=list(polort=2, jobs=6))
ret <- build(aglm)



### FMRIModel convolves events with HRF, and takes contrasts, F-test specifications

# what to do with basis functions?
# EV --> Regressor or RegressorSet
emodel <- EventModel(Onset ~ direction*day*poly(motion,2), data=destable, block=run, durations=time)
fmridesign <- FMRIDesign(emodel, dsource, HRF.GAMMA)

emodel2 <- EventModel(Onset ~ direction*poly(motion,2), data=subset(destable, run<9), block=run)
#fmridesign2 <- FMRIDesign(emodel2, dsource, HRF.GAMMA)
fmridesign3 <- FMRIDesign(emodel2, dsource, create.HRF(HRF.BSPLINE))

aglm <- afni.glm(fmridesign3, "junk", mask="/Users/brad/data1/Dots/KAR/Data_RS/epi/global_mask.nii")


tmp <- function(expr) {
  browser()
}
