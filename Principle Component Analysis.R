### Plotting DE in boxplots -- Andrew R Gross -- 2018-08-14
### This script will generate boxplots from differential expression data. 
### This script operates downstream of the standard DESeq2 Pipeline
### INPUT: Expression data; differential expression data is typical
### OUTPUT: Boxplot figures

############################################################################################
### Header
library(ggplot2)
library(ggbiplot)
library(biomaRt)
ensembl = useMart(host='www.ensembl.org',biomart='ENSEMBL_MART_ENSEMBL',dataset="hsapiens_gene_ensembl")
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)

############################################################################################
### Functions

add.description <- function(dataframe, identifier = (c('ensembl_gene_id', 'external_gene_name' ))) {
  descr <- getBM(attributes=c(identifier,'description'), filters= identifier, values=row.names(dataframe), mart=ensembl)
  descr <- descr[match(row.names(dataframe),descr[,1]),]
  descriptions <- c()
  for (rowNumber in 1:length(descr[,1])) {
    newDescr <- descr[rowNumber,][,2]
    newDescr <- strsplit(newDescr, " \\[")[[1]][1]
    descriptions <- c(descriptions, newDescr)
  }
  dataframe[length(dataframe)+1] <- descriptions
  names(dataframe)[ncol(dataframe)] <- "Description"
  return(dataframe)
}

############################################################################################
### Functions
#convert.ids <- function(dataframe) {
  ### This version is obsolte
#  ensemblIDs <- c()
 # gene.names <- c()
#  for (rowName in row.names(dataframe)) {
#    ensemblID <- strsplit(rowName,"\\.")[[1]][1]
 #   gene.name <- strsplit(rowName,"\\_")[[1]][2]
  #  ensemblIDs <- c(ensemblIDs, ensemblID)
   # gene.names <- c(gene.names, gene.name)
#  }
 # #row.names(dataframe) <- ensemblIDs
  #row.names(dataframe) <- make.unique(gene.names)
  #return(dataframe)
#}

convert.ids <- function(dataframe, add.gene.name.column = TRUE) {
  ### This function will convert a row name consisting of a contactenated ensembl ID and gene to one or the other,
  ### based on the users instruction (2018-10-04)
  ensemblIDs <- c()                                           # Empty lists are initialized to receive IDs as they're created
  gene.names <- c()
  for (rowName in row.names(dataframe)) {                     # Loops through all rows in the data frame
    ensemblID <- strsplit(rowName,"\\.")[[1]][1]                 # Splits the row name and declares the ensembl ID
    gene.name <- strsplit(rowName,"\\_")[[1]][2]                 # Splits the row name, declares the gene name
    ensemblIDs <- c(ensemblIDs, ensemblID)                       # Adds ensembl ID and gene name to appropriate lists
    gene.names <- c(gene.names, gene.name)
  }
  row.names(dataframe) <- ensemblIDs                          # assigns the new row names
  if(add.gene.name.column == TRUE) {
    dataframe$Gene <- gene.names
  }
  return(dataframe)                                           # Returns the data frame with new rows
}

new.names <- c()
for (names in test) {
  new.name <- strsplit(names, '\\_')[[1]][2]
  new.names <- c(new.names,new.name)
}


### Input

setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/RNAseq Data/Motor Neurons/E099 - PM - ALS v CTRL/')
tpm.als <- read.csv('PM-5119--07--02--2018_TPM.csv', row.names = 1)
sample.names.als <- c('02iCTR','03iCTR','159iALS','172iCTR','372iALS_n1','372iALS_n2','372iALS_n3','395iCTR')


setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/RNAseq Data/Motor Neurons/E0x - PP - D10 KO v CTRL/')
tpm.ko <- read.csv('TPM.csv', row.names = 1)
sample.names.ko <- read.csv('Sample names.txt', header = FALSE)[,1]







############################################################################################
### Format

### Rename columns
names(tpm.als) <- sample.names.als

names(tpm.ko) <- sample.names.ko

### Filter low expression genes
summary(tpm.als)
tpm.als.max <- apply(tpm.als, 1, max)
rows.to.keep <- which(tpm.als.max >10)
length(rows.to.keep)
### Reorder columns
results.als <- tpm.als[rows.to.keep,][c(1,2,4,8,3,5,6,7)]

tpm.ko.max <- apply(tpm.ko, 1, max)
rows.to.keep <- which(tpm.ko.max >5)
length(rows.to.keep)
results.ko <- tpm.ko[rows.to.keep,]

### Replace row names
#results.als <- convert.ids(results.als)

### Convert to matrix
results.als <- as.matrix(results.als)
results.als <- as.matrix(results.ko)

############################################################################################
### Calculate Principle Components

### Calculate the actual components
pca.als <- prcomp(t(results.als), scale = TRUE)

### Calculate the percent variation accounted for by each component
pca.als.var <- pca.als$sdev^2
pca.als.var.per <- round(pca.als.var/sum(pca.als.var)*100, 1)

### Identify the genes with the largest influence
# PC1
l.score.pc1 <- pca.als$rotation[,1]
l.score.pc1.ranked <- sort(abs(l.score.pc1), decreasing = TRUE)
l.score.pc1[names(l.score.pc1.ranked)][1:10]

# PC2
l.score.pc2 <- pca.als$rotation[,2]
l.score.pc2.ranked <- sort(abs(l.score.pc2), decreasing = TRUE)
l.score.pc2[names(l.score.pc2.ranked)][1:10]

############################################################################################
### Plot
plot(pca.als$x[,1], pca.als$x[,2])
plot(pca.als$x[,2], pca.als$x[,3])

barplot(pca.als.var.per, main = 'Scree Plot', xlab = 'Principle Component', ylab = 'Percent Variation')

### GGPlot

pca.als.to.plot <- data.frame(Sample = rownames(pca.als$x), 
                              X = pca.als$x[,1],
                              Y = pca.als$x[,2])

pca.als.to.plot
pca.als.to.plot

pca.plot <- ggplot(data = pca.als.to.plot, aes(x = X, y = Y, label = Sample)) +
  geom_text() + 
  xlab(paste('PC1 - ', pca.als.var.per[1], '%', sep = '')) +
  ylab(paste('PC2 - ', pca.als.var.per[2], '%', sep = '')) +
  theme_bw() +
  ggtitle('PCA of E099')

pca.plot

#ggbiplot(pca.als)

ggbiplot(pca.als, var.axes = FALSE, labels = rownames(pca.als$x))
#ggbiplot(mtcars.pca,ellipse=TRUE,obs.scale = 1, var.scale = 1,var.axes=FALSE,   labels=rownames(mtcars), groups=mtcars.country)

###############################################################################################
### Find the loadings that MATTER
### Convert all loadings into vectors of the first two components.  Find the axis of interest.
### Then, identify the loadings with the greatest magnitudes along this axis.

##############################################
### Find the slope and magnitude of all loadings
loadings <- data.frame(l.score.pc1, l.score.pc2)
loadings <- loadings*100
loadings$slope = loadings[,2]/loadings[,1]
loadings$degree = atan(loadings$slope)*180/pi
loadings$magnitude = abs((loadings[,1]^2+loadings[,2]^2)^0.5)
loadings$mag2 = abs(loadings[,1])+abs(loadings[,2])

### Check the distribution of loadings
plot(dist(loadings$magnitude[sample(1:nrow(loadings),100,replace = FALSE)]))
plot(dist(loadings$mag2[sample(1:nrow(loadings),100,replace = FALSE)]))

subsample <- loadings[sample(1:nrow(loadings),1000,replace = FALSE),]
plot(subsample$l.score.pc1, subsample$magnitude)
plot(subsample$l.score.pc2, subsample$magnitude)

plot(subsample$l.score.pc1, subsample$mag2)
plot(subsample$l.score.pc2, subsample$mag2)
hist(loadings$magnitude)
hist(loadings$mag2)

plot(subsample$degree, subsample$magnitude)
plot(subsample$degree, subsample$mag2)

subsample$mag3 <- (subsample$magnitude + subsample$mag2*0.5)
plot(subsample$degree, subsample$mag3)

##############################################
### Find the axis of interest
### Find center point of control samples:
mean.x.of.ctr <- mean(pca.als.to.plot$X[1:4])
mean.y.of.ctr <- mean(pca.als.to.plot$Y[1:4])
mean.ctr <- c('CTR mean', mean.x.of.ctr, mean.y.of.ctr)

mean.x.of.als <- mean(pca.als.to.plot$X[5:8])
mean.y.of.als <- mean(pca.als.to.plot$Y[5:8])
mean.als <- c('ALS mean', mean.x.of.als, mean.y.of.als)

#means.x <- c(mean.x.of.ctr, mean.x.of.als)
#means.y <- c(mean.y.of.ctr, mean.y.of.als)
#centers.of.sample.groups <- data.frame(means.x, means.y)

pca.als.to.plot <- data.frame(Sample = c(rownames(pca.als$x),'CTR mean','ALS mean'), 
                              X = c(pca.als$x[,1],mean.x.of.ctr, mean.x.of.als),
                              Y = c(pca.als$x[,2],mean.y.of.ctr, mean.y.of.als),
                              color = c('Pink','Pink','Pink','Pink','Blue','Blue','Blue','Blue','Red','Dark_blue'))

### Plot with means of each group
pca.plot <- ggplot(data = pca.als.to.plot, aes(x = X, y = Y, label = Sample)) +
    geom_text(aes(color = color)) + 
    #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9","Red")) +
    scale_color_manual(values = c('Blue','Dark Green','Red','Orange')) +
    xlab(paste('PC1 - ', pca.als.var.per[1], '%', sep = '')) +
    ylab(paste('PC2 - ', pca.als.var.per[2], '%', sep = '')) +
    theme_bw() + 
    ggtitle('PCA of E099') 

pca.plot


### Identify slope of interest
(slope.of.interest <- mean.y.of.als/mean.x.of.als)

### Visualize range of interest
range.1 = 1.2
range.2 = 2
range.3 = 10
pca.plot.w.slopes <- pca.plot + 
  geom_abline(slope=slope.of.interest*range.1, color = 'Gray20')+ geom_abline(slope=slope.of.interest/range.1, color = 'Gray20')+
  geom_abline(slope=slope.of.interest*range.2, color = 'Gray50')+ geom_abline(slope=slope.of.interest/range.2, color = 'Gray50')+
  geom_abline(slope=slope.of.interest*range.3, color = 'Gray70')+ geom_abline(slope=slope.of.interest/range.3, color = 'Gray70')
# geom_abline(slope = slope.of.interest)
pca.plot.w.slopes

##############################################
### Downselect loadings based on angle
rows.selected <- c()
### SET 1: Range 3
rows.of.interest <- intersect(which(loadings$slope < slope.of.interest/range.3),
                              which(loadings$slope > slope.of.interest*range.3))
#length(rows.of.interest)/nrow(loadings)
length(rows.of.interest)
loadings.of.interest <- loadings[rows.of.interest,]

subsample <- loadings.of.interest[sample(1:nrow(loadings),1000,replace = FALSE),]
plot(subsample$degree, subsample$magnitude)

### Find the nth percentile of magnitude
#(percentile.98 <- quantile(loadings.of.interest$magnitude, 0.9))
#rows.to.keep <- row.names(loadings.of.interest[which(loadings.of.interest$magnitude>percentile.98),])
#length(rows.to.keep)
#rows.selected <- c(rows.selected, rows.to.keep)
#length(rows.selected)
#loadings.of.interest <- loadings.of.interest[rows.to.keep,]

### SET 2: Range 2
rows.of.interest <- intersect(which(loadings$slope < slope.of.interest/range.2),
                              which(loadings$slope > slope.of.interest*range.2))
length(rows.of.interest)
loadings.of.interest <- loadings[rows.of.interest,]

### Find the nth percentile of magnitude
#(percentile.98 <- quantile(loadings.of.interest$magnitude, 0.99))
#rows.to.keep <- which(loadings.of.interest$magnitude>percentile.98)
#length(rows.to.keep)
#rows.selected <- c(rows.selected, rows.to.keep)
#length(rows.selected)

### SET 3: Range 1
rows.of.interest <- intersect(which(loadings$slope < slope.of.interest/range.1),
                              which(loadings$slope > slope.of.interest*range.1))
length(rows.of.interest)
loadings.of.interest <- loadings[rows.of.interest,]
### Find the nth percentile of magnitude
#(percentile.98 <- quantile(loadings.of.interest$magnitude, 0.98))
#rows.to.keep <- which(loadings.of.interest$magnitude>percentile.98)
#length(rows.to.keep)
#rows.selected <- c(rows.selected, rows.to.keep)
#length(rows.selected)

##############################################################################
#genes.to.plot <- loadings[1:2][unique(rows.selected),]*30
genes.to.plot <- loadings.of.interest
genes.to.plot[1:2] <- genes.to.plot[1:2]*30
genes.to.plot$x0 = 0
genes.to.plot$y0 = 0
genes.to.plot$Sample = 1
### Filter by level
genes.to.plot <- genes.to.plot[which(genes.to.plot$magnitude>1.4),]

#ggplot(data = genes.to.plot, aes(x = x0, y = y0, xend = l.score.pc1, yend = l.score.pc2, color = mag2)) +
#  geom_segment() +
#  scale_color_gradient(low = "white", high = "black")

pca.plot.w.slopes + geom_segment(data = genes.to.plot, aes(x = x0, y = y0, xend = l.score.pc1, yend = l.score.pc2))


##############################################################################
### Format row names
### Convert row names into gene names for easy reading or Entrez ids for GAGE analysis

loadings.null <- loadings[sample(1:nrow(loadings.of.interest),replace = FALSE),]

loadings.of.interest <- convert.ids(loadings.of.interest)
tpm.of.interest <- tpm.als[row.names(loadings.of.interest),]
tpm.of.interest <- convert.ids(tpm.of.interest)

tpm.o.i.ez <- convert.to.entrez(tpm.of.interest)
str(tpm.o.i.ez)
tpm.o.i.ez[4]
tpm.o.i.ez2 <- tpm.o.i.ez[[1]]
grep(8623, tpm.o.i.ez2$join)
tpm.o.i.ez2 <- tpm.o.i.ez2[-c(2202,2203),]
tpm.o.i.ez2 <- tpm.o.i.ez2[-581,]

### Assign new IDs to row names
row.names(tpm.o.i.ez2) <- tpm.o.i.ez2$join



loadings.null <- convert.ids(loadings.null)
tpm.null <- tpm.als[row.names(loadings.null),]
tpm.null <-convert.ids(tpm.null)
tpm.null <- convert.to.entrez(tpm.null)
tpm.null2 <- tpm.null[[1]]
tpm.null2 <- tpm.null2[-c(),]











### Filter TPM list based on loadings of interest

for.gage.loadings <- tpm.als[]

pway.loadings <- gage(for.gage.loadings, gsets = kegg.gs, ref = ctrl.index, samp = ko.index, compare = "unpaired")





### Annotate!

add.description(loadings.of.interest, 'external_gene_name')



load.scores.als <- pca.als$rotation[,1]
gene.scores <- abs(load.scores.als)
gene.scores.ranked <- sort(gene.scores, decreasing = TRUE)
gene.scores.ranked <- names(gene.scores.ranked[1:10])
pca.als$rotation[gene.scores.ranked,1]


setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E099 - RNAseq analysis of CHCHD10/DEG analyses/')
row.num.pos = 10
gene.name <- results.als$Gene[row.num.pos]
data.als <- results.als[row.num.pos,][-c(1,2,3,12,13)]
data.ko <- results.ko[row.num.pos,][-c(1,2,3,10,11)]
expression <- t(data.als)
expression.ko <- t(data.ko)

expression <- rbind(expression,expression.ko)
Disease <- c("CTR", "CTR", "CTR", "CTR", "ALS", "ALS", "ALS", "ALS", 'WT', 'WT', 'WT', 'KO', 'KO', 'KO')

expression <- data.frame(expression, Disease)
names(expression) <- c('tpm', 'dis')

expression$dis <- factor(expression$dis, c('CTR','ALS','WT','KO'))

boxplot(tpm~dis, data=expression, main=gene.name, xlab="Genome type", ylab="Expression [TPM]")

png(paste0(gene.name,'.png'))
boxplot(tpm~dis, data=expression, main=gene.name, xlab="Genome type", ylab="Expression [TPM]")
dev.off()


p <- ggplot(test, aes(Dis, TPM))
p + geom_boxplot()

c('Normalized Expression (TPM)', 'Genenome type')
test <- t(row.of.interest[-c(1,2,3,12,13)])
Dis  <- c("CTR", "CTR", "CTR", "CTR", "ALS", "ALS", "ALS", "ALS")
test <- cbind(test,Dis)
test <- data.frame(test)
names(test) <- c('TPM', 'Dis')

test[1] <- as.numeric(levels(test[,1]))[test[,1]]
boxplot(TPM~Dis, data=test, main=row.of.interest$gene.names, xlab="Expression", ylab="Diseases state")
