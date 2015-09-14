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
if (!file.exists(statdir)){
    dir.create(statdir)
}


# Adressing the files
DAG <- read.delim(file="22790_24796.CDS-CDS.last.tdd10.cs0.filtered.dag.all.go", header=FALSE, skip=1)
head(DAG, 25)

# Number of mapps
n<-dim(DAG)[1]

# Reorder the factors of chromosome number (10 ordered between 1 and 2)
DAG$V5 = factor(DAG$V5,levels(DAG$V5)[c(1,3:10,2)])


#########################
#Graphical representation
#########################

# Leftover code: m<-sum(DAG$V5=="b24796_1")
DAG6 <- sapply(DAG$V6, function(n) unlist(strsplit(as.character(n), "||", fixed=TRUE)))
dim(DAG6) <- c(9, n)
DAG6 <- t(DAG6)

# Trunk to print out one chromosome profile in one pdf file.
#densChr1 <- apply(DAG6[DAG6[,1]=='1',2:3], 1, function(cord) seq(as.numeric(cord[1]), as.numeric(cord[2])) )
#densChr1.10 <- apply(DAG6[DAG6[,1]=='1',2:3], 1, function(cord) seq(as.numeric(cord[1])/10, as.numeric(cord[2])/10) )
#pdf(file="density_map_chromosom1.pdf", paper="a4r", width = (2967/150)/2.54, height = (2099/150)/2.54)
#hist(floor(unlist(densChr1.10)), freq=FALSE, breaks=diff(range(densChr1.10))/10)
#dev.off()
# 1. Error: unexpected numeric constant in "ulimit -t 600"
# 2. Killed

# Trunk to print out all chromosome profiles in one pdf file.
pdf(file=paste(getwd(),"/output/density_map_chromosomes.pdf", sep=""), paper="a4r", width = (2967/150)/2.54, height = (2099/150)/2.54)
par(mfrow=c(2,5))
for (n in as.character(1:10)) {
    chrDensPlot(n)
}
dev.off()


# do number

has for gene names / how many times gene name is covered /
build a test case to check code
read and summarize papaya/peach chapter in "Genetics and Genomics of Papaya"