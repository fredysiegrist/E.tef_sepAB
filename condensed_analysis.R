#author -- Fredy Siegrist, PhD
#date -- 09.09.2015
#title -- statistical analyses of tef scaffold matches from CoGe's SynMap on Sorghum chromosomes

require(limma)
source("tef_analysis_functions_2015.R")


##########
#Load Data
##########


# set working directory and prepare output directory
setwd("../../i1sz/")
statdir <- file.path("output")


# Adressing the files
cond <- read.delim(file="stretches_condensed3.csv", header=FALSE, skip=0, comment.char = "")
colnames(cond) <- c("chr", "start", "end", "ori", "scaffold", "sstart", "send", "sori", "stretch", "number", "score", "reversion", "gir", "block", "position")
chrlen <- c(73840612, 77932577, 74440842, 68034345, 62352331, 62208772, 64342021, 55460251, 59635592, 60981625)

require('colorspace')
seqcol <- heat_hcl(16, c=c(80,30), l = c(30,90), power = c(1/5, 1.5))
scorecol <- terrain_hcl(20, c=c(65,0), l=c(45, 95), power=c(1/3, 1.5))
seqcol[1] <- "#3C3CEC"

pdf(file=paste(getwd(),"/output/density_map_chromosomes_cond4.pdf", sep=""), paper="a4r", width = (2967/150)/2.54, height = (2099/150)/2.54)
for (chrno in 1:10) {
    n<-1
    plot(1:chrlen[chrno], (n:(chrlen[chrno])), type='n', sub=paste("Sorghum chromosome #",chrno), xlab="nt position on ", ylab="runif distributed E. tef scaffolds", main="Mapping of E. tef on Sorghum")
    apply(cond[cond$chr==chrno, c(2:3,15)], 1, function(z) {n<-sample(1:chrlen[chrno],1); x<- z[1:2]; y<-c(n,n); colr <- seqcol[z[3]]; lines(x, y, col=colr, lwd=3)})
}

par(mfrow=c(2,5))
for (chrno in 1:10) {
    n<-1
    plot(1:chrlen[chrno], (n:(chrlen[chrno])), type='n', sub=paste("chr #",chrno), xlab="nt", ylab="scaffolds", main="E. tef on Sorghum")
    apply(cond[cond$chr==chrno, c(2:3,15)], 1, function(z) {n<-sample(1:chrlen[chrno],1); x<- z[1:2]; y<-c(n,n); colr <- seqcol[z[3]]; lines(x, y, col=colr, lwd=2)})
}
dev.off()

# Find the stretches that cover more than 1 billion nucleotides:
cond[(cond[,3]-cond[,2])>2000000,]

# Plot the s the match length on Sorghum chromosome against scaffold length of E.tef for found stretches
plot(abs(cond[,3]-cond[,2]), abs(cond[,6]-cond[,7]), pch="•", cex=1, col=seqcol[cond[,15]], xlab="Sorghum chromosome legnth", ylab="E. tef scaffold length", main="scaffold length vs. match length", log="xy", xlim=c(5e+3, 1e+7), ylim=c(1e+3, 1e+6))
abline(0, glm(log(abs(cond[,6]-cond[,7])) ~ log(abs(cond[,3]-cond[,2]))-1)$coefficients,  lwd=2, col="yellow4")
abline(0, 1,  lwd=2, col="yellow2", lty=3)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE)
length(levels(factor(.bincode(cond[,11], seq(75, 3500, by=100)))))
points(abs(cond[,3]-cond[,2]), abs(cond[,6]-cond[,7]), pch="•", cex=0.6, col=scorecol[as.numeric(levels(factor(.bincode(cond[,11], seq(75, 3500, by=100)))))])
legend(2e+6, 5e+3, legend=c("први   stretch","други stretch","шести stretch","score   75 -   175","score 476 -   575","score 976 - 1075"), col=c(seqcol[c(1,2,6)], scorecol[c(1, 5, 10)]), pch=c(rep(c(16, 20), each=3)))
(75+100*(as.numeric(levels(factor(.bincode(cond[,11], seq(75, 3500, by=100)))))-1))[c(1, 5, 10)]


# reconstruct the master file
bestlist <- cond[cond[,15]==1,]
bestlist <- bestlist[order((bestlist[,1]*1e10)+bestlist[,2]+bestlist[,3]),]

    # scaffold7112 appears on chromosome 1 while it is attributed to chromosome 2 in the master list.
cond[cond$scaffold=="scaffold7112",]

# RECONSTRUCTION

remaster <- bestlist[,c(1, 5, 12)]

# read in the master file for comparision and refiltering

master <- read.delim(file="master_22790_24796.CDS-CDS.last.tdd10.cs0.filtered.dag.all.go_D20_g10_A3.aligncoords.gcoords_ct0.w1000.spa1.sr.cs1.csoS.log.nsd.spa_info.txt", header=TRUE)
# The most stupid and most often used way to solve the problem - easy solution welcome!
levels(master$X.CHR1) <- c("01", "10", "02", "03", "04", "05", "06", "07", "08", "09", "unmapped")
master$X.CHR1 <- factor(master$X.CHR1, levels=c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10","unmapped"))
master <- master[order(master$X.CHR1),]

connec <- match(master$CHR2,remaster$scaffold, nomatch=NULL)
connec <- connec[!is.na(connec)]

rownames(remaster[order(master[connec,2]),])
plot(1:3014, connec, pch=16, xlab="ordinary # in new master", ylab="position on SynMap master file", cex=1, main="ordering in master and new master file", col=colors()[as.numeric(remaster$chr)+98])
legend(2000, 500, unique(remaster$chr), col=colors()[98:107], pch=16, ncol=4, title="E. tef Chromosome" )

