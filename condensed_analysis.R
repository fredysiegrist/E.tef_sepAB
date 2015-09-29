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
cond <- read.delim(file="stretches_condensed2.csv", header=FALSE, skip=0, comment.char = "")
colnames(cond) <- c("chr", "start", "end", "ori", "scaffold", "sstart", "send", "sori", "stretch", "number", "score", "reversion", "gir", "block", "position")
chrlen <- c(73840612, 77932577, 74440842, 68034345, 62352331, 62208772, 64342021, 55460251, 59635592, 60981625)

require('colorspace')
seqcol <- heat_hcl(16, c=c(80,30), l = c(30,90), power = c(1/5, 1.5))
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

cond[(cond[,3]-cond[,2])>10000000,]
