# TODO: Add comment
# 
# Author: brad
###############################################################################


BlockStructure <- function(blocklens, TR) {
	if (length(TR) == 1) {
		TR <- rep(TR, length(blocklens))
	}
			
	new("BlockStructure", blocklens=blocklens, TR=TR)
	
}
