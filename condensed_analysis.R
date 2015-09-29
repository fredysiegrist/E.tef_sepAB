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
cond <- read.delim(file="stretches_condensed.csv", header=FALSE, skip=0, comment.char = "")
colnames(cond) <- c("chr", "start", "end", "ori", "scaffold", "sstart", "send", "sori", "stretch", "number", "score", "reversion", "gir", "block", "position")
chrlen <- c(73840612, 77932577, 74440842, 68034345, 62352331, 62208772, 64342021, 55460251, 59635592, 60981625)


pdf(file=paste(getwd(),"/output/density_map_chromosomes_cond2.pdf", sep=""), paper="a4r", width = (2967/150)/2.54, height = (2099/150)/2.54)
for (chrno in 1:10) {
    n<-1
    plot(1:chrlen[chrno], (n:(chrlen[chrno])), type='n', sub=paste("Sorghum chromosome #",chrno), xlab="nt position on ", ylab="runif distributed E. tef scaffolds", main="Mapping of E. tef on Sorghum")
    apply(cond[cond$chr==chrno, c(2:3,15)], 1, function(z) {n<-sample(1:chrlen[chrno],1); x<- z[1:2]; y<-c(n,n); colr <- z[3]+1; lines(x, y, col=z, lwd=3)})
}

par(mfrow=c(2,5))
for (chrno in 1:10) {
    n<-1
    plot(1:chrlen[chrno], (n:(chrlen[chrno])), type='n', sub=paste("chr #",chrno), xlab="nt", ylab="scaffolds", main="E. tef on Sorghum")
    apply(cond[cond$chr==chrno, c(2:3,15)], 1, function(z) {n<-sample(1:chrlen[chrno],1); x<- z[1:2]; y<-c(n,n); colr <- z[3]+1; lines(x, y, col=z, lwd=2)})
}
dev.off()

# Find the stretches that cover more than 1 billion nucleotides:

cond[(cond[,3]-cond[,2])>10000000,]

"""
     chr    start      end ori     scaffold sstart   send sori stretch number
1732   2  7423978 71237703   1 scaffold6518  26669 132614   -1       3      1
4045   1 27532747 68908460   1 scaffold4243   9731  29813   -1       3      1
     score reversion gir block position
1732   144         r   3  3775        2
4045   150         r   3  2637        1
"""