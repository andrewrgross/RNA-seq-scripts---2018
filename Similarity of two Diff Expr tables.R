### Similarity of Two Diff. Expr. Tables - Andrew R Gross - 2018-09-12
### A standard exectution of the DESeq2 pipeline.  
### Requires the output data of "DESeq2 Pipeline"
### INPUT: This script requires differential expression tables generated using the DESeq2 pipeline.
### OUTPUT: This script generates list of shared genes and a Venn Diagram displaying the number of genes shared.

####################################################################################################################################################
### Header
library("pasilla")
library("DESeq2")
library("biomaRt")
library("VennDiagram")

##########################################################################
### Input - Desired input is the output of DESeq2 pipeline

setwd("C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E099 - RNAseq analysis of CHCHD10/Input files/")

results.als <- read.csv('DEGs - ALS 0.05 pv.csv', row.names= 1)
results.ko <- read.csv('DEGs - KO 0.05 pv.csv', row.names = 1)

### Find the overlap

shared.genes <- intersect(results.als$Gene, results.ko$Gene)

length(shared.genes)

### Draw Venn Diagram
draw.pairwise.venn(area1 = nrow(results.als), area2 = nrow(results.ko), cross.area = length(shared.genes), 
                   category = c("Mutant vs. Wildtype", "Knockout vs. Wildtype"),
                   lty = rep("blank", 2), fill = c("light blue", "pink"), 
                   alpha = rep(0.5, 2), cat.pos = c(0, 0), 
                   cat.dist = rep(0.025, 2))
                                                                                                                                                        0), cat.dist = rep(0.025, 2))

### Find expression values
shared.pos <- match(shared.genes, results.als$Gene)
results.als.shrd <- results.als[shared.pos,]

shared.pos <- match(shared.genes, results.ko$Gene)
results.ko.shrd <- results.ko[shared.pos,]

### Create data frame of shared genes
names(results.als.shrd)[1:3] <- c('Mut: Mean', 'Mut: lFC', 'Mut: pv')
names(results.ko.shrd)[1:3] <- c('KO: Mean', 'KO: lFC', 'KO: pv')

results.als.shrd <- results.als.shrd[order(row.names(results.als.shrd)),]
results.ko.shrd <- results.ko.shrd[order(row.names(results.ko.shrd)),]

results.shared.genes <- cbind(results.als.shrd, results.ko.shrd)
results.shared.genes <- results.shared.genes[c(12,13,1,2,3,14,15,16)]
shared.genes.order <- order(results.shared.genes$`Mut: pv`,decreasing = TRUE)
results.shared.genes <- results.shared.genes[shared.genes.order,]
dir.check <- sign(results.shared.genes$`Mut: lFC`*results.shared.genes$`KO: lFC`)==1
results.shared.genes$same.dir <- dir.check

setwd("C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E099 - RNAseq analysis of CHCHD10/DEG analyses/")

write.csv(results.shared.genes, 'Shared genes - 0.05 padj.csv')

### Sort by direction

up.reg <- results.ALS.shared$log2FoldChange > 0
als.up <- results.ALS.shared[up.reg,]
als.down <- results.ALS.shared[!up.reg,]

up.reg <- results.ko.shared$log2FoldChange > 0
ko.up <- results.ko.shared[up.reg,]
ko.down <- results.ko.shared[!up.reg,]

### Find Shared by direction

shared.up <- intersect(als.up$Gene, ko.up$Gene)
shared.down <- intersect(als.down$Gene, ko.down$Gene)

shared.up.pos <- match(shared.up, als.up$Gene)
shared.up.df <- als.up[shared.up.pos,]

shared.down.pos <- match(shared.down, als.down$Gene)
shared.down.df <- als.down[shared.down.pos,]

### Save results

setwd("C:/Users/grossar/Box/Sareen Lab Shared/Data/Prasanthi/")

write.csv(shared.up.df, "Shared UP reg - ALS-KO.csv")
write.csv(shared.down.df, "Shared DOWN reg - ALS-KO.csv")
write.csv(results.ALS.shared, 'Shared DEGs - ALS-KO.csv')

