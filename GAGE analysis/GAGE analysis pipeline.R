### GAGE Pathway Analysis -- Andrew R Gross -- 2018-08-17
### Application of the GAGE pathway analysis pipeline
### This analysis is downstream of the DESeq2 Pipeline
### INPUT: Differential Expression table with ENTREZ IDs
### OUTPUT: Predicted pathways tables, predicted pathways bar graphs

####################################################################################################################################################
### Header
library(gage)
library(ggplot2)

####################################################################################################################################################
### Input
data(gse16873)
data(kegg.gs)
data(go.gs)

setwd("C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E099 - RNAseq analysis of CHCHD10/Input files/")

results.als.ez <- read.csv('DEGs - ALS 0.05 padj ENTREZ.csv', row.names= 1)
results.ko.ez <- read.csv('DEGs - KO 0.05 padj ENTREZ.csv', row.names = 1)

results.als.ez <- read.csv('DEGs - ALS 0.05 pv.csv', row.names= 1)
results.ko.ez <- read.csv('DEGs - KO 0.05 pv ENTREZ.csv', row.names = 1)

####################################################################################################################################################
### Formatting
#results.als.ez <- results.als.ez[results.als.ez$padj<0.05,]
for.gage.als <- as.matrix(results.als.ez[c(4,5,6,7,8,9,10,11)])

#results.ko.ez <- results.ko.ez[results.ko.ez$padj<0.05,]
for.gage.ko <- as.matrix(results.ko.ez[c(7,8,9,4,5,6)])

##########################################################################
### Confirm ID type for gene set
lapply(kegg.gs[1:3],head)
head(kegg.gs[1])

head(go.gs[1])    # They all appear to be Entrez from googling

##########################################################################
### Run analysis

### Example analysis
#cnts.kegg.p <- gage(cnts.norm, gsets = kegg.gs, ref = ref.idx, samp = samp.idx, compare ="unpaired")
#gse16873.kegg.p <- gage(gse16873, gsets = kegg.gs, ref = hn, samp = dcis)
##go.gs here only the first 1000 entries as a fast example.
#gse16873.go.p <- gage(gse16873, gsets = go.gs, ref = hn, samp = dcis)
##or we can do analysis on BP, MF, or CC subcategories if we've
##generated the data as above.
##gse16873.bp.p <- gage(gse16873, gsets = go.bp,
## ref = hn, samp = dcis)
#str(gse16873.kegg.p, strict.width='wrap')
#gse16873.kegg.2d.p <- gage(gse16873, gsets = kegg.gs, + ref = hn, samp = dcis, same.dir = F)

#results.kegg.2d.p <- gage(results.als.ez, gsets = kegg.gs, + ref = ctrl, samp = als, same.dir = F)

### Run real data
ctrl.index = 1:4
als.index = 5:8
pway.als.bd <- gage(for.gage.als, gsets = kegg.gs, ref = ctrl.index, samp = als.index, compare = "unpaired", same.dir=F)
pway.als.bd <- pway.als.bd$greater
pway.als <- gage(for.gage.als, gsets = kegg.gs, ref = ctrl.index, samp = als.index, compare = "unpaired")
pway.als.up <- pway.als$greater
pway.als.down <- pway.als$less

#head(pway.ALS$greater,20)
#head(pway.ALS$less)

ctrl.index = 1:3
ko.index = 4:6
pway.ko.bd <- gage(for.gage.ko, gsets = kegg.gs, ref = ctrl.index, samp = ko.index, compare = "unpaired", same.dir=F)
pway.ko.bd <- pway.ko.bd$greater
pway.ko <- gage(for.gage.ko, gsets = kegg.gs, ref = ctrl.index, samp = ko.index, compare = "unpaired")
pway.ko.up <- pway.ko$greater
pway.ko.down <- pway.ko$less

head(pway.ko.bd)
head(pway.ko.up)
head(pway.ko.down)
######################################################################################################################################
### Output data

############################################################
### Plot function
gen.pway.bar.chart <- function(data.frame, experiment, direction) {
  if(direction == 'Up') {color = 'Red'}
  if(direction == 'Down') {color = 'Blue'}
  if(direction == 'Bidir') {color = 'Purple'; direction = 'Bidirectional'}
  plot.df <- as.data.frame(data.frame)[1:20,]
  title = paste('Disreg. Pathways in', experiment, ' - ', direction)
  
  plot.df <- data.frame(Pathways = row.names(plot.df), log.p.val = -log(plot.df$p.val))
  plot.df$Pathways <- factor(plot.df$Pathways, levels = rev(plot.df$Pathways))
  
  plot <- ggplot(plot.df, aes(Pathways,log.p.val)) +
    geom_col(fill = color) + coord_flip() + theme_bw() + labs(title = title)
  return(plot)
}

experiment = 'Mutant'
direction = 'Up'
plot <- gen.pway.bar.chart(pway.als.up, paste(experiment, 'vs. WT'), direction)
plot

direction = 'Down'
plot <- gen.pway.bar.chart(pway.als.down, paste(experiment, 'vs. WT'), direction)
plot

direction = 'Bidir'
plot <- gen.pway.bar.chart(pway.als.bd, paste(experiment, 'vs. WT'), direction)
plot

experiment = 'KO'
direction = 'Up'
plot <- gen.pway.bar.chart(pway.ko.up, paste(experiment, 'vs. WT'), direction)
plot

direction = 'Down'
plot <- gen.pway.bar.chart(pway.ko.down, paste(experiment, 'vs. WT'), direction)
plot

direction = 'Bidir'
plot <- gen.pway.bar.chart(pway.ko.bd, paste(experiment, 'vs. WT'), direction)
plot

setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E099 - RNAseq analysis of CHCHD10/Pathway analyses/Pathway bar charts - 0.05 pv/')

output.name = paste0('Pway bar chrt - ', experiment, ' - ', direction, '.png')
png(filename = output.name,
    width = 1000, height = 800, units = "px", pointsize = 20,
    bg = "white",  res = 120)
plot
dev.off()


############################################################
### Manual plotting
plot.df <- as.data.frame(pway.als.up)[1:10,]
title = 'Disreg. Pathways in ALS lines - Up'; color = 'Red'
plot.df <- as.data.frame(pway.als.down)[1:10,]
title = 'Disreg. Pathways in ALS lines - Down'; color = 'Blue'

plot.df <- as.data.frame(pway.ko.up)[1:10,]
title = 'Disreg. Pathways in KO lines - Up'; color = 'Red'
plot.df <- as.data.frame(pway.ko.down)[1:10,]
title = 'Disreg. Pathways in KO lines - Down'; color = 'Blue'


plot.df <- data.frame(Pathways = row.names(plot.df), log.p.val = -log(plot.df$p.val))
### Set order
plot.df$Pathways <- factor(plot.df$Pathways, levels = rev(plot.df$Pathways))

ggplot(plot.df, aes(Pathways,log.p.val)) +
  geom_col(fill = color) + coord_flip() + theme_bw() + labs(title = title)

ggplot(plot.df, aes(Pathways,log.p.val)) +
  geom_col(fill = 'Purple') + coord_flip() + theme_bw() + labs(title = 'Disreg. Pathways in ALS lines - Bidirectional')


####################################
### Output
setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E099 - RNAseq analysis of CHCHD10/Pathway analyses/')

write.csv(pway.als.bd, 'pways - ALS 0.05pv - Bidir.csv')
write.csv(pway.als.up, 'pways - ALS 0.05pv - Up.csv')
write.csv(pway.als.down, 'pways - ALS 0.05pv - Down.csv')
write.csv(pway.ko.bd, 'pways - ko 0.05pv - Bidir.csv')
write.csv(pway.ko.up, 'pways - ko 0.05pv - Up.csv')
write.csv(pway.ko.down, 'pways - ko 0.05pv - Down.csv')





##########################################################################
### Pathview Analysis
setwd('C:/Users/grossar/Box/Sareen Lab Shared/Data/Andrew/E099 - RNAseq analysis of CHCHD10/Pathway analyses/Mapped pathways/')

fc.data <- as.matrix(results.als.ez)[, 2] ; suffix = 'als - 0.05pv' 
fc.data <- as.matrix(results.ko.ez)[, 2] ; suffix = 'ko - 0.05pv' 

pway.id = '04145'
pway.id = '00190'
pway.id = '04512'
pway.id = '04722'
pway.id = '04540'

pv.out <- pathview(gene.data = fc.data, pathway.id = pway.id, species = "hsa", out.suffix = suffix)

pv.out <- pathview(gene.data = fc.data, pathway.id = pway.id, species = "hsa", out.suffix = suffix)
