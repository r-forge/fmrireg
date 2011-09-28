source("../R/rlist.SOURCE")


design <- read.table("all_data2.txt", header=TRUE)
design$repnum <- ifelse(design$repnum <= 1, 0, design$repnum-1)
design <- subset(design, run==1)

reg1 <- fmrireg(onset ~ hrf(word_status) + hrf(word_status, repnum, subset=repnum!=0) + block(run, 100, 1.5), hrf.fun="spmg1", data=design)
