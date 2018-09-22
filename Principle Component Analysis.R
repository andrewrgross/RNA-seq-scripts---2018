### Plotting DE in boxplots -- Andrew R Gross -- 2018-08-14
### This script will generate boxplots from differential expression data. 
### This script operates downstream of the standard DESeq2 Pipeline
### INPUT: Expression data; differential expression data is typical
### OUTPUT: Boxplot figures

############################################################################################
### Header
library(ggplot2)

############################################################################################
### Functions
convert.ids <- function(dataframe) {
  ensemblIDs <- c()
  gene.names <- c()
  for (rowName in row.names(dataframe)) {
    ensemblID <- strsplit(rowName,"\\.")[[1]][1]
    gene.name <- strsplit(rowName,"\\_")[[1]][2]
    ensemblIDs <- c(ensemblIDs, ensemblID)
    gene.names <- c(gene.names, gene.name)
  }
  #row.names(dataframe) <- ensemblIDs
  row.names(dataframe) <- make.unique(gene.names)
  return(dataframe)
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
sample.names.ko <- read.csv('Sample names.txt', header = FALSE)







############################################################################################
### Format

### Rename columns
names(tpm.als) <- sample.names.als

### Filter low expression genes
summary(tpm.als)
tpm.als.max <- apply(tpm.als, 1, max)
rows.to.keep <- which(tpm.als.max >10)
length(rows.to.keep)
### Reorder columns
results.als <- tpm.als[rows.to.keep,][c(1,2,4,8,3,5,6,7)]

### Replace row names
results.als <- convert.ids(results.als)

### Convert to matrix
results.als <- as.matrix(results.als)

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

### Find the loadings that MATTER
### Convert all loadings into vectors of the first two components.  Find the axis of interest.
### Then, identify the loadings with the greatest magnitudes along this axis.

### Finding my axis of interest
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
                              Y = c(pca.als$x[,2],mean.y.of.ctr, mean.y.of.als))

### Plot with means of each group
ggplot(data = pca.als.to.plot, aes(x = X, y = Y, label = Sample)) +
    geom_text() + 
    xlab(paste('PC1 - ', pca.als.var.per[1], '%', sep = '')) +
    ylab(paste('PC2 - ', pca.als.var.per[2], '%', sep = '')) +
    theme_bw() +
    ggtitle('PCA of E099') 

### Identify slope of interest
slope.of.interest <- mean.y.of.als/mean.x.of.als


### Finding the slope of all loadings
loadings <- data.frame(l.score.pc1, l.score.pc2)
loadings$slope = loadings$l.score.pc2/loadings$l.score.pc1
loadings <- loadings[order(loadings$slope),]

### Create a range of interesting slopes
rows.of.interest <- intersect(which(loadings$slope > slope.of.interest/0.9),
                              which(loadings$slope < slope.of.interest*1.1))
loadings.of.interest <- loadings[rows.of.interest,]

### Create a range of interesting slopes
rows.of.interest <- intersect(which(loadings$slope < -slope.of.interest/0.9),
                              which(loadings$slope > -slope.of.interest*1.1))
loadings.of.interest <- loadings[rows.of.interest,]




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
