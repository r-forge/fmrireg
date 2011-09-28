#! /usr/bin/env Rscript

library(methods)
library(neuroim)
library(gregmisc)


cat(paste("sourcing", file.path(getwd(), "config.R")))
source(file.path(getwd(), "config.R"))


attach(config)




dsource <- DataSource(EPI_PATH, unique(DESIGN$image), NIMAGES, REPTIME)

hrf.spmg <- HRF.SPMG1


DESIGN$repnum <- factor(ifelse(DESIGN$repnum <= 1, 0, DESIGN$repnum-1))
DESIGN$lag <- factor(DESIGN$lag)

reg1 <- fmrireg(onset ~ hrf(repnum,lag, word_status) +
                nuisance(PCNuisance, label="Nuisance") +
                nuisance(MotNuisance, label="Motion") + 
			    block(run, blocklens(dsource), TR(dsource)), 
			    hrf.fun="spmg1", data=DESIGN)

fmrimod <- FMRIDesign(reg1, dsource)
term1 <- terms(reg1)[[1]]

lag_pseudo <- polyContrast(term1, "lag", where=lag != 0 & word_status == "pseudo", degree=1)
lag_real <- polyContrast(term1, "lag", where=lag != 0 & word_status == "real", degree=1)

lag_by_wtype <- lag_real - lag_pseudo

glt.set <- list(
	main_word = buildGLT(baselineContrast(term1, A=repnum %in% levels(repnum))),
        main_rep = buildGLT(pairContrast(term1, A=repnum != 0, B=repnum == 0)),
        main_wtype = buildGLT(pairContrast(term1, A=word_status == "real", B=word_status == "pseudo")),
	lin_lag = buildGLT(polyContrast(term1, "lag", where=lag != 0, degree=2)),
        lin_repnum = buildGLT(polyContrast(term1, "repnum", where=repnum != 0, degree=2)),
        lin_lag_pseudo = buildGLT(lag_pseudo),
        lin_lag_real = buildGLT(lag_real),
        lin_lag_real_pseudo = buildGLT(lag_by_wtype))
        
                

runmodel <- function(model, outname, glts) {
  glm <- afni.glm(model, outname, mask=MASK, options=list(polort=5, jobs=7, noFDR=TRUE, rout=FALSE, fout=FALSE), gltlist=glt.set)
  ret <- build(glm)
  ret
}

 

afnimod <- runmodel(fmrimod, "glm_basic", glt.set)
afnimod$run()


                
