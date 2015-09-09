#author -- Fredy Siegrist, PhD
#date -- 09.09.2015
#title -- Analysis functions -for- statistical analyses of tef scaffold matches from CoGe's SynMap on Sorghum chromosomes


# compute fold changes from log fold changes -- for testing purpose
calcFC=function(lfc){
	fc=unlist(lapply(lfc, function(x){
					if (is.na(x)) { NA}
					else if(x>=0) {2^x}
					else{-(1/2^x)}}))
	return(fc)
 }