#author -- Fredy Siegrist, PhD
#date -- 09.09.2015
#title -- Analysis functions -for- statistical analyses of tef scaffold matches from CoGe's SynMap on Sorghum chromosomes

# Some basic ploting functions for tef genome analysis.

# compute fold changes from log fold changes -- for testing purpose

chrDensPlot <- function(chr="1", data=DAG6){
    densChr.10 <- apply(data[data[,1]==chr,2:3], 1, function(cord) seq(floor(as.numeric(cord[1])/10), floor(as.numeric(cord[2])/10)))
    print(data[data[,1]==chr,2:3])
    print(length(floor(unlist(densChr.10))))
    hist(unlist(densChr.10), freq=TRUE, breaks=diff(range(densChr.10))/10, main=paste("chr.",chr), xlab="10-base occ.")
}

chrDensPlot2 <- function(chr="1", data=DAG6){
    densChr <- apply(data[data[,1]==chr,2:3], 1, function(cord) seq(as.numeric(cord[1]), as.numeric(cord[2])))
    hist(unlist(densChr), freq=TRUE, breaks=diff(range(densChr)), main=paste("chr.",chr), xlab="base occ.")
}

densityMap <- function(data, ink, chromosomes = 1:10, overview=TRUE, chrlen = c(73840612, 77932577, 74440842, 68034345, 62352331, 62208772, 64342021, 55460251, 59635592, 60981625)) {
if (overview) {
    par(mfrow=c(2,ceiling(length(chromosomes)/2)))
    subprefix="chr #"; xlabl="nt"; ylabl="scaffolds"; mainl="E. tef on Sorghum"
    }
else {
    subprefix="Sorghum chromosome #"; xlabl="nt position on "; ylabl="runif distributed E. tef scaffolds"; mainl="Mapping of E. tef on Sorghum"
    }
for (chrno in chromosomes) {
    plot(1:chrlen[chrno], (1:(chrlen[chrno])), type='n', sub=paste(subprefix,chrno), xlab=xlabl, ylab=ylabl, main=mainl)
    apply(data[data$chr==chrno, c(2:3,15)], 1, function(z) {n<-sample(1:chrlen[chrno],1); x<- z[1:2]; y<-c(n,n); colr <- ink[z[3]]; lines(x, y, col=colr, lwd=2)})
    }
}

# plot the differences of positions to the position
plotdiff <- function(betterlist, best, slimmaster, connec, penalty) {
    betterlist <- cond[best,]
    itermaster <- betterlist[,c(1, 5, 12)]
    iterconnec <- match(slimmaster$CHR2, itermaster$scaffold, nomatch=NULL)
    diffs <- (1:length(iterconnec)) - iterconnec
    plot(1:length(connec), (1:length(connec))-iterconnec, pch=16, xlab="ordinal # in better master", ylab="position difference", cex=1, main=paste("ordering in master and ", penalty, "master file"), col=colors()[as.numeric(itermaster$chr)+98])
    }

# Calculation of how many nucleotides are attributed to A and B chromosome,
# allowing counting for each reference chromosome separately.

nttable <- function(data) {
 m<-cbind('A'=sum(data[data[,10]==1,6])-sum(data[data[,10]==1,5]),
'B'=sum(data[data[,10]==2,6])-sum(data[data[,10]==2,5]),
'unmapped'=sum(data[data[,10]==3,6])-sum(data[data[,10]==3,5]) )
rownames(m) <- as.character(i)
 return(m) }