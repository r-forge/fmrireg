.SID=scan("SID")
.BASE_PATH=file.path("/euterpe/data/brad/projects/words_study", .SID)
.MASK = file.path(.BASE_PATH, "epi", "global_mask.nii")
.EPI_PATH = file.path(.BASE_PATH, "epi")

.SCANS = file.path(.EPI_PATH, scan("scans.txt", what="character"))
.RUNS = scan("runs.txt")
.NRUNS = length(.RUNS)
.DESIGN_FILE = file.path(.BASE_PATH, "behavior/all_data2.txt")
.DESIGN = subset(read.table(.DESIGN_FILE, header=T), run %in% .RUNS)
.TR = 1.37
.NIMAGES=302
.PCNuisance=lapply(scan("pc_nuisance.txt", what=character()), read.table, header=TRUE)
.MotNuisance=lapply(scan("mot_nuisance.txt",what=character()), read.table, header=TRUE)


.DESIGN$image <- rep(.SCANS, each=table(.DESIGN$run))
.DESIGN$trial <- factor(seq(1, NROW(.DESIGN)))
.DESIGN$srepnum <- scale(.DESIGN$repnum)

source("~brad/Rcode/FMRIDesign/R/rlist.SOURCE")

library(gregmisc)

config = list(
	SID=.SID,
	BASE_PATH=.BASE_PATH,
	MASK=.MASK,
	EPI_PATH=.EPI_PATH,
	SCANS=.SCANS,
	RUNS=.RUNS,
	NRUNS=.NRUNS,
	DESIGN_FILE=.DESIGN_FILE,
	DESIGN=.DESIGN,
	REPTIME=.TR,
	NIMAGES=.NIMAGES,
	PCNuisance=.PCNuisance,
	MotNuisance=.MotNuisance	
	)
