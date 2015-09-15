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


chrDensPlot <- function(chr="1", data=DAG6){
    densChr.10 <- apply(data[data[,1]==chr,2:3], 1, function(cord) seq(as.numeric(cord[1])/10, as.numeric(cord[2])/10) )
    print(data[data[,1]==chr,2:3])
    print(length(floor(unlist(densChr.10))))
    hist(floor(unlist(densChr.10)), freq=TRUE, breaks=diff(range(densChr.10))/10, main=paste("chr.",chr), xlab="10-base occ.")

}
