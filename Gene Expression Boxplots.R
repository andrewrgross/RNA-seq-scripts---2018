### Plotting DE in boxplots -- Andrew R Gross -- 2018-08-14
### This script will generate boxplots from differential expression data. 
### This script operates downstream of the standard DESeq2 Pipeline
### INPUT: Expression data; differential expression data is typical
### OUTPUT: Boxplot figures

### Header
library(ggplot2)

### Input

setwd("C:/Users/grossar/Box/Sareen Lab Shared/Data/Prasanthi/")

de.data <- read.csv("ALS DiffEx prelim data - Aug 14.csv", row.names = 1)

de.data <- results.als
### Format

### Break gene name into column
test <- row.names(head(de.data))

row.names.split <- strsplit(row.names(de.data), "_")

test.split <- strsplit(test, "_")

#for (item in test.split[1]) {
#  print(item[2])
#}

### Run through the rows and split

ensembl.id <- c()
gene.names <- c()
  
for (item.num in 1:length(row.names.split)) {
  #print(test.split[item.num][[1]][2])
  current.id <-   row.names.split[item.num][[1]][1]     # Assigns the ensembl ID of hte current row
  current.gene <- row.names.split[item.num][[1]][2]   # Assigns the gene name of the current row
  ensembl.id <- c(ensembl.id, current.id)             # Adds the current ID to the list
  gene.names <- c(gene.names, current.gene)           # Adds the current gene to the list
}

de.data <- cbind(de.data,gene.names)
row.names(de.data) <- ensembl.id

### Plot
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
