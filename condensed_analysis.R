#author -- Fredy Siegrist, PhD
#date -- 09.09.2015
#title -- statistical analyses of tef scaffold matches from CoGe's SynMap on Sorghum chromosomes

source("tef_analysis_functions_2015.R")

# set working directory and prepare output directory
setwd("../../i1sz/")
statdir <- file.path("output")

# Loading the file
# cond <- read.delim(file="output.scoreminus1/stretches_condensed_12mio_scoreFirst.csv", header=FALSE, skip=0, comment.char = "")
cond <- read.delim(file="stretches_condensed_ac.csv", header=FALSE, skip=0, comment.char = "")

colnames(cond) <- c("chr", "start", "end", "ori", "scaffold", "sstart", "send", "sori", "stretch", "number", "score", "reversion", "gir", "block", "position")
cond$scaffold <- as.character(cond$scaffold)
if (max(cond$score) > 1) cond[,11] <- 1/cond[,11]
chrlen <- c(73840612, 77932577, 74440842, 68034345, 62352331, 62208772, 64342021, 55460251, 59635592, 60981625)

# set colors for plotting
require('colorspace')
seqcol <- heat_hcl(16, c=c(80,30), l = c(30,90), power = c(1/5, 1.5))
scorecol <- terrain_hcl(20, c=c(65,0), l=c(45, 95), power=c(1/3, 1.5))
seqcol[1] <- "#3C3CEC"

# This is a very energy consuming thing, don't run it if not necessary
pdf(file=paste(getwd(),"/output/density_map_chromosomes_cond.pdf", sep=""), paper="a4r", width = (2967/150)/2.54, height = (2099/150)/2.54)
densityMap(cond, seqcol, overview=FALSE)
densityMap(cond, seqcol)
dev.off()

# Find the stretches that cover more than 1 billion nucleotides:
cond[(cond[,3]-cond[,2])>2000000,]

pdf(file=paste(getwd(),"/output/bettermaster.pdf", sep=""), paper="a4", width = (2099/100)/2.54, height = (2967/100)/2.54)
par(mfrow=c(3,1))
# Plot the s the match length on Sorghum chromosome against scaffold length of E.tef for found stretches
plot(abs(cond[,3]-cond[,2]), abs(cond[,6]-cond[,7]), pch=16, cex=1, col=seqcol[cond[,15]], xlab="Sorghum chromosome legnth", ylab="E. tef scaffold length", main="scaffold length vs. match length", log="xy", xlim=c(5e+3, 1e+7), ylim=c(1e+3, 1e+6))
abline(0, glm(log(abs(cond[,6]-cond[,7])) ~ log(abs(cond[,3]-cond[,2]))-1)$coefficients,  lwd=2, col="yellow4")
abline(0, 1,  lwd=2, col="yellow2", lty=3)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE)
length(levels(factor(.bincode(cond[,11], seq(75, 3500, by=100)))))
points(abs(cond[,3]-cond[,2]), abs(cond[,6]-cond[,7]), pch=16, cex=0.4, col=scorecol[as.numeric(levels(factor(.bincode(cond[,11], seq(75, 3500, by=100)))))])
legend(1e+6, 1e+4, legend=c("1st stretch","2nd stretch","3rd stretch","score   75 -   175","score 476 -   575","score 976 - 1075"), col=c(seqcol[c(1,2,6)], scorecol[c(1, 5, 10)]), pch=c(rep(c(16, 20), each=3)), ncol=2)
(75+100*(as.numeric(levels(factor(.bincode(cond[,11], seq(75, 3500, by=100)))))-1))[c(1, 5, 10)]

# reconstruct the master file
bestlist <- cond[cond[,15]==1,]
bestlist <- bestlist[order((bestlist[,1]*1e10)+(bestlist[,2]+bestlist[,3])/2),]

# scaffold7112 appears on chromosome 1 while it is attributed to chromosome 2 in the master list.
cond[cond$scaffold=="scaffold7112",]

# RECONSTRUCTION of master file from CoGe output.
remaster <- bestlist[,c(1, 5, 12)]

# Read in the master file for comparision and refiltering.
master <- read.delim(file="master_22790_24796.CDS-CDS.last.tdd10.cs0.filtered.dag.all.go_D20_g10_A3.aligncoords.gcoords_ct0.w1000.spa1.sr.cs1.csoS.log.nsd.spa_info.txt", header=TRUE)
# The most stupid and most often used way to solve the problem
# - easy solution welcome!
levels(master$X.CHR1) <- c("01", "10", "02", "03", "04", "05", "06", "07", "08", "09", "unmapped")
master$X.CHR1 <- factor(master$X.CHR1, levels=c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10","unmapped"))
master$CHR2 <- as.character(master$CHR2)
master <- master[order(master$X.CHR1),]

# Reducing the CoGe master file to have only entries matching my
slimmaster <- master[which(master$CHR2 %in% cond$scaffold),]

# Matching the positions on the two master files
connec <- match(slimmaster$CHR2,remaster$scaffold, nomatch=NULL)

# ploting the positions on the two master files
plot(1:length(connec), connec, pch=16, xlab="ordinal # in new master", ylab="position on SynMap master file", cex=1, main="ordering in master and new master file", col=colors()[as.numeric(remaster$chr)+98])
legend(0, 2000, unique(remaster$chr), col=colors()[99:108], pch=16, ncol=5, title="E. tef Chromosome" )

# Highlight the outliers and set the 2nd ranked Stretch
outliers <- (1:length(connec))[abs((1:length(connec)) - connec) > sd((1:length(connec))- connec)]
points(outliers, connec[outliers], pch="°", col="red", cex=2)

pdf(file=paste(getwd(),"/output/positiondiffs.pdf", sep=""), width = 14, height = 10)
betterlist <- cond[cond[,15]==1,]
outlierpos <- cond[cond[,5] %in% slimmaster[outliers,2],]
outlier <- outlierpos[outlierpos[,15]==2,]
# write.csv(cbind(outlierpos[outlierpos[,15]==2,], outlierpos[outlierpos[,15]==1,]), file="outlier12.csv")
betterlist[match(outlier$scaffold, betterlist$scaffold),] <- outlier
betterlist <- betterlist[order((betterlist[,1]*1e10)+betterlist[,2]+betterlist[,3]),]
betterremaster <- betterlist[,c(1, 5, 12)]
newcon <- match(slimmaster$CHR2, betterremaster$scaffold, nomatch=NULL)
updated <- (1:3018)[abs(newcon - connec) > 100] #match(outlier$scaffold, cond$scaffold)
plot(1:3018, (1:3018)-newcon, pch=16, xlab="placement in better master", ylab="position difference", cex=1, main="ordering in master and 'better' master file", col=colors()[as.numeric(betterremaster$chr)+98])
legend(0, 3, unique(betterremaster$chr), col=colors()[98:107], pch=16, ncol=5)
points(updated, updated-newcon[updated], pch="°", col="green", cex=2)
dev.off()

# Iterate over the Stretches for outliers to find best combination of positions
# overall penalty of the last 500 steps set to 500:1
# poor man algorithm

betterlist <- cond[cond[,15]==1,]
betterlist <- betterlist[order((betterlist[,1]*1e10)+betterlist[,2]+betterlist[,3]),]
best <- match(paste(betterlist$scaffold, betterlist$position), paste(cond$scaffold, cond$position))
itermaster <- betterlist[,c(1, 5, 12)]
iterconnec <- match(slimmaster$CHR2, itermaster$scaffold, nomatch=NULL)
diffs <- (1:length(iterconnec)) - iterconnec
penalty <- sum(diffs^2)
pos_diff <- 500:1
counter <- 0
i <- 0
bestsquare <- 1E10
betterscore <- sum(betterlist$score)
mincounter <- NULL
while (counter < 2501 && length(unique(pos_diff))>1) {
    oldpenalty <- penalty
    # sample from an array of at least 2 in length!
    worst <- sample(c(rep(which(abs(diffs) %in% max(abs(diffs))), max(c(1,var(diffs)))), sample(which(abs(diffs) %in% 1:max(abs(diffs))), max(c(1, sd(diffs)))), sample(dim(slimmaster)[1],3)),1)
    alternativestretches <- cond[cond[,5] %in% slimmaster[worst,2],]
    pos <- match(alternativestretches$scaffold[1],  betterlist$scaffold)
    # find the best position
    ifelse (betterlist[pos,15]==max(alternativestretches[,15]), nextStretch <- 1, nextStretch <- (betterlist[pos,15]+1) )
    betterlist[pos,] <- alternativestretches[alternativestretches[,15]==nextStretch,]
    betterlist <- betterlist[order((betterlist[,1]*1e10)+betterlist[,2]+betterlist[,3]+rnorm(1, 10000, 10000)),]
    itermaster <- betterlist[,c(1, 5, 12)]
    iterconnec <- match(slimmaster$CHR2, itermaster$scaffold, nomatch=NULL)
    diffs <- (1:length(iterconnec)) - iterconnec
    penalty <- sum(diffs^2)
    # Feed the new difference to the end of the 500 last differences.
    if((oldpenalty - penalty)>0) {
        pos_diff <- c(pos_diff[2:500], oldpenalty - penalty)
        }
    # Save best positions for reloading if the algorithm get lost.
    if ( penalty < bestsquare | (penalty <= bestsquare & sum(betterlist$score) < betterscore)) { # & pos_diff[500] == min(pos_diff)
        ifelse(best == match(paste(betterlist$scaffold, betterlist$position), paste(cond$scaffold, cond$position)),  i <- i+1, {betterscore <- sum(betterlist$score); best <- match(paste(betterlist$scaffold, betterlist$position), paste(cond$scaffold, cond$position)); print(paste('saved', penalty, betterscore, counter, i)); i <- 0; bestsquare <- penalty})
        plotdiff(betterlist, best, slimmaster, connec, bestsquare)
        }
    # Check if the best improvement was 500 cycles ago and reload randomly.
    else {
        if (pos_diff[1] == min(pos_diff)) {
            counter <- counter + 1
            if (sample(10, 1)==1) {
                betterlist <- cond[best,]
                print(paste('loaded', bestsquare, counter))
                betterlist <- cond[best,]
                print(quantile(mincounter))
                mincounter <- NULL
                }
            mincounter <- cbind(mincounter, penalty)
            }
        }
}

# reorder them to calculate 'true' sum 1/DAGscore
betterlist <- cond[best,]
betterlist <- betterlist[order((betterlist[,1]*1e10)+betterlist[,2]+betterlist[,3]),]
itermaster <- betterlist[,c(1, 5, 12)]
iterconnec <- match(slimmaster$CHR2, itermaster$scaffold, nomatch=NULL)
diffs <- (1:length(iterconnec)) - iterconnec
penalty <- sum(diffs^2)
( betterscore <- sum(betterlist$score) )

plotdiff(betterlist, best, slimmaster, connec, bestsquare)
dev.off()

# Generate a better master list and write
bettermaster <- cbind(betterlist[,c(1:3, 5:7, 12, 13)], no=as.numeric(rownames(betterlist)), chrt=factor(rep("unmapped", dim(betterlist)[1]), levels=c("A", "B", "unmapped") ))
write.table(bettermaster, file="output/bettermaster_secondatonce.delim", quote=FALSE, sep="\t")

# Definition of colors to display DAGchainer score
max(as.numeric(rownames(bettermaster)))
length(unique(bettermaster[,8]))
gcol <- terrain_hcl(36, c=c(65,0), l=c(45, 95), power=c(1/3, 1.5))
gircol <- rep("#FFFFFF", 68)
girs <- unique(bettermaster[,8])
girorder <- girs[order(unique(bettermaster[,8]))]
transmatrix <- cbind(unique(bettermaster[,8]), order(girorder, decreasing = TRUE))
gircol[transmatrix[,1]] <- gcol[transmatrix[,2]]

pdf(file=paste(getwd(),"/output/bettermaster_secondatonce.pdf", sep=""), paper="a4r", width = (2967/100)/2.54, height = (2099/100)/2.54)

par(mfrow=c(2,5))
for (chrno in 1:10) {
    n<-1
    plot(seq(1, chrlen[chrno], length.out=max(as.numeric(rownames(bettermaster)))), (1:max(as.numeric(rownames(bettermaster)))), type='n', sub=paste("chr #",chrno), xlab="nt", ylab="scaffolds", main="E. tef on Sorghum")
    apply(bettermaster[bettermaster$chr==chrno, c(2:3,8:9)], 1, function(z) {n<-z[4]; x<- z[1:2]; y<-c(n,n); colr <- gircol[z[3]]; lines(x, y, col=colr, lwd=3)})
}
dev.off()

# Split to A and B chromosome
# 'Random' assignement to A and B chromosome starting on A, then B the rest
# remains unmapped.
for (entry in 1:dim(bettermaster)[1]) {
    actchr <- bettermaster[entry, 1]
    if (sum(bettermaster$chrt=="A" & bettermaster$chr==actchr)==0) { bettermaster[entry, 10] <- "A" }
    else { if (bettermaster[entry, 2] > range(bettermaster[bettermaster$chrt=="A" & bettermaster$chr==actchr, 2:3])[2]) { bettermaster[entry, 10] <- "A" }
           else { if (sum(bettermaster$chrt=="B" & bettermaster$chr==actchr)==0)  { bettermaster[entry, 10] <- "B" }
                  else { if (bettermaster[entry, 2] > range(bettermaster[bettermaster$chrt=="B" & bettermaster$chr==actchr, 2:3])[2]) { bettermaster[entry, 10] <- "B" }
                  }
           }
    }
}

bettermaster[,10] <- as.numeric(bettermaster[,10])

pdf(file=paste(getwd(),"/output/firstAB_maxsorted.pdf", sep=""), paper="a4r", width = (2967/100)/2.54, height = (2099/100)/2.54)
for (chrno in 1:10) {
    plot(seq(1, chrlen[chrno], length.out=5), (0:4), type='n', sub=paste("chr #",chrno), xlab="nt", ylab="scaffolds", main="Randomly attributed E. tef scaffolds to A/B/umapped on Sorghum chromosomes")
    apply(bettermaster[bettermaster$chr==chrno, c(2:3,8,10)], 1, function(z) {n<-jitter(as.numeric(z[4]), 5); x<- z[1:2]; y<-c(n,n); colr <- gircol[z[3]]; lines(x, y, col=colr, lwd=5)})
}
par(mfrow=c(2,5))
for (chrno in 1:10) {
    plot(seq(1, chrlen[chrno], length.out=5), (0:4), type='n', sub=paste("chr #",chrno), xlab="nt", ylab="scaffolds", main="E. tef on Sorghum")
    apply(bettermaster[bettermaster$chr==chrno, c(2:3,8,10)], 1, function(z) {n<-jitter(as.numeric(z[4]), 3); x<- z[1:2]; y<-c(n,n); colr <- gircol[z[3]]; lines(x, y, col=colr, lwd=3)})
}

dev.off()

# Bettersplit using number of genes in rows

require(IRanges)
bettermaster <- cbind(betterlist[,c(1:3, 5:7, 12, 13)], no=as.numeric(rownames(betterlist)), chrt=factor(rep("unmapped", dim(betterlist)[1]), levels=c("A", "B", "unmapped") ))
ABmaster <- bettermaster[order(bettermaster$gir, decreasing = TRUE),]

for (entry in 1:dim(ABmaster)[1]) {
    actchr <- ABmaster[entry, 1]
    ranges <- IRanges(ABmaster[entry,2], ABmaster[entry,3])
    if (sum(ABmaster$chrt=="A" & ABmaster$chr==actchr)==0) { ABmaster[entry, 10] <- "A" }
    else {
        rangesA <- IRanges(ABmaster[ABmaster$chrt=="A" & ABmaster$chr==actchr, 2], ABmaster[ABmaster$chrt=="A" & ABmaster$chr==actchr, 3])
        if ( countOverlaps(ranges, rangesA, type="any")==0 ) { ABmaster[entry, 10] <- "A" }
           else {
                rangesB <- IRanges(ABmaster[ABmaster$chrt=="B" & ABmaster$chr==actchr, 2], ABmaster[ABmaster$chrt=="B" & ABmaster$chr==actchr, 3])
                if (sum(ABmaster$chrt=="B" & ABmaster$chr==actchr)==0)  { ABmaster[entry, 10] <- "B" }
                  else {
                        rangesB <- IRanges(ABmaster[ABmaster$chrt=="B" & ABmaster$chr==actchr, 2], ABmaster[ABmaster$chrt=="B" & ABmaster$chr==actchr, 3])
                        if (countOverlaps(ranges, rangesB, type="any")==0 )  { ABmaster[entry, 10] <- "B" }
                  }
           }
    }
}

ABmasterSorted <- ABmaster[order(ABmaster[,1]*1e10+(ABmaster[,2]+ABmaster[,3])/2),]
write.table(ABmasterSorted, file="output/ABmaster_2atonce.delim", quote=FALSE, sep="\t")
ABmaster[,10] <- as.numeric(ABmaster[,10])

# Calculation of how many scaffold are attributed to A and B chromosome
table(ABmasterSorted$chrt)
print(nttable(ABmaster, 'All'))
for (i in 1:10) print(nttable(ABmaster[ABmaster[,1]==i,], i))

pdf(file=paste(getwd(),"/output/ABmaster2_maxsorted_2nd.pdf", sep=""), paper="a4r", width = (2967/100)/2.54, height = (2099/100)/2.54)
for (chrno in 1:10) {
    plot(seq(1, chrlen[chrno], length.out=5), (0:4), type='n', sub=paste("chr #",chrno), xlab="nt", ylab="scaffolds", main="Randomly attributed E. tef scaffolds to A/B/umapped on Sorghum chromosomes")
    apply(ABmasteier[ABmaster$chr==chrno, c(2:3,8,10)], 1, function(z) {n<-jitter(as.numeric(z[4]), 5); x<- z[1:2]; y<-c(n,n); colr <- gircol[z[3]]; lines(x, y, col=colr, lwd=5)})
}

par(mfrow=c(2,5))
for (chrno in 1:10) {
    plot(seq(1, chrlen[chrno], length.out=5), (0:4), type='n', sub=paste("chr #",chrno), xlab="nt", ylab="scaffolds", main="E. tef on Sorghum")
    apply(ABmaster[ABmaster$chr==chrno, c(2:3,8,10)], 1, function(z) {n<-jitter(as.numeric(z[4]), 3); x<- z[1:2]; y<-c(n,n); colr <- gircol[z[3]]; lines(x, y, col=colr, lwd=3)})
}

dev.off()

# Tutorial for usage of IRanges
library(IRanges)
ranges <- IRanges(ABmaster[entry,2], ABmaster[entry,3])
rangesA <- IRanges(ABmaster[ABmaster$chrt=="A" & ABmaster$chr==actchr, 2], ABmaster[ABmaster$chrt=="A" & ABmaster$chr==actchr, 3])
rangesB <- IRanges(ABmaster[ABmaster$chrt=="B" & ABmaster$chr==actchr, 2], ABmaster[ABmaster$chrt=="B" & ABmaster$chr==actchr, 3])
rangesX <- IRanges(ABmaster[ABmaster$chrt=="unmapped" & ABmaster$chr==actchr, 2], ABmaster[ABmaster$chrt=="unmapped" & ABmaster$chr==actchr, 3])
# Which rangesX contain at least one ranges?
countOverlaps(ranges, rangesX, type="any")>0

# Easy plot the occupance on the genome
pdf(file=paste(getwd(),"/output/plotRanges1.pdf", sep=""), paper="a4r", width = (2967/100)/2.54, height = (2099/100)/2.54)

par(mfrow=c(2,5))
for (chrno in 1:10) {
    rangesC <- IRanges(cond[cond$chr==chrno, 2], cond[cond$chr==chrno, 3])
    plotRanges(rangesC, main=paste("All on Chr", chrno))
}

for (chrno in 1:10) {
    rangesC <- IRanges(betterlist[betterlist$chr==chrno, 2], betterlist[betterlist$chr==chrno, 3])
    plotRanges(rangesC, main=paste("Best on Chr", chrno))
}

dev.off()


