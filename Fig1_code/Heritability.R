#factor:  
options(stringsAsFactors = FALSE)
rm(list=ls())
gc()
library('lme4')
library('parallel')

################### Heritability for Tropical lines ############################
dfTrop <- read.csv('e340_trop_repUnmerged_Untransformed_fpkm.csv', header= TRUE, 
                   check.names = FALSE, row.names = 1)


###### only use duplicated genotypes #####
cn = colnames(dfTrop)
idx = colnames(dfTrop) %in% colnames(dfTrop)[duplicated(colnames(dfTrop))]
dfTrop <- dfTrop[, idx]
colnames(dfTrop) = cn[idx]
length(unique(colnames(dfTrop)))
group <- factor(colnames(dfTrop))
geneList = rownames(dfTrop)

############### Run mixed linear model ############################
HbTrop <- pbapply::pblapply(seq_along(geneList),  function(i){
  trait = geneList[i]
  gexp = as.numeric(dfTrop[trait, ])
  model = lmer(gexp ~ (1|group))
  varComp <- as.data.frame(VarCorr(model,comp='vcov'))
  sigmas <- varComp$vcov
  H2 <- sigmas[1]/sum(sigmas)
  H2
}, cl=32)

################### Heritability for Temperate lines ###########################
dfTemp <- read.csv('e340_temp_repUnmerged_Untransformed_fpkm.csv', header= TRUE, check.names = FALSE, row.names = 1)
#################### only use duplicated genotypes #############################
cn = colnames(dfTemp)
idx = colnames(dfTemp) %in% colnames(dfTemp)[duplicated(colnames(dfTemp))]
dfTemp <- dfTemp[, idx]
colnames(dfTemp) = cn[idx]
length(unique(colnames(dfTemp)))
group <- factor(colnames(dfTemp))
geneList = rownames(dfTemp)
################## Run the mixed linear model ##################################
HbTemp <- pbapply::pblapply(seq_along(geneList),  function(i){
  trait = geneList[i]
  gexp = as.numeric(dfTemp[trait, ])
  model = lmer(gexp ~ (1|group))
  varComp <- as.data.frame(VarCorr(model,comp='vcov'))
  sigmas <- varComp$vcov
  H2 <- sigmas[1]/sum(sigmas)
  H2
}, cl=32)

########################################### Combine ##################################
library(ggplot2)
library(grid)
library(MASS) 
library(reshape2) 
library(reshape)
HbTrop <- as.numeric(HbTrop)
HbTemp <- as.numeric(HbTemp)
dfHb_Comp = data.frame(HbTrop = HbTrop, HbTemp = HbTemp)
rownames(dfHb_Comp) <- geneList
write.csv(dfHb_Comp, file = 'Hb_Trop_VS_Temperate.csv')

dfHb_Comp <- read.csv('Hb_Trop_VS_Temperate.csv', row.names = 1)
dfHbTempLow <- subset(dfHb_Comp, dfHb_Comp$HbTrop*0.8 >= dfHb_Comp$HbTemp)
write.csv(dfHbTempLow, file = 'TempHbLow.csv')
write.csv(rownames(dfHbTempLow), file = 'HbTempLow.list', row.names = F,  quote = F)


#### plot the 2d density plot #######
dfHb_Comp <- read.csv('Hb_Trop_VS_Temperate.csv', row.names = 1)
dfHb_Comp <- subset(dfHb_Comp, HbTemp >0 & HbTrop > 0)

p1<- ggplot(data = dfHb_Comp, aes(HbTrop, HbTemp))+
  stat_density_2d(aes(fill = ..level..), geom='polygon')+
  geom_smooth(position = 'identity', method = lm)+
  scale_fill_distiller(palette=12, direction=1) +
  theme(legend.position='right')+
  xlim(0,1) + ylim(0,1) +
  xlab("H2 Tropical") +
  ylab("H2 Temperate") +
  theme_classic()

coef <- round(coef(lm(HbTemp~HbTrop, dfHb_Comp)), 2)

p1 = p1 + annotate("text", x=0.05,  y= 1, hjust=0, col="red", size=4, fontface="italic",
                   label=substitute(y==a*x+b, list(a=coef[1], b=coef[2])))
p1

########## distribution of the tropical H2 
dfHb_Comp <- read.csv('Hb_Trop_VS_Temperate.csv')
dfHb_trop_high <- subset(dfHb_Comp, HbTrop > HbTemp)
dfHb_trop_low <- subset(dfHb_Comp, HbTrop < HbTemp) 
dfHb_trop_high$comparison = 'TropicalHigh'
dfHb_trop_low$comparison = 'TropicalLow'
dfHb_Comp_merge = rbind(dfHb_trop_high, dfHb_trop_low)

ggplot(data = dfHb_Comp_merge, aes(HbTrop, fill = comparison)) +
  geom_density(alpha=0.4) +
  theme_classic()

######################################################
################### Heritability for all 340 lines ############################
dfexp <- read.csv('eGWAS340_repUnmerged_Untransformed_fpkm.csv', header= TRUE, check.names = FALSE, row.names = 1)
cn = colnames(dfexp)
idx = colnames(dfexp) %in% colnames(dfexp)[duplicated(colnames(dfexp))]
dfexp <- dfexp[, idx]
colnames(dfexp) = cn[idx]
length(unique(colnames(dfexp)))

group <- factor(colnames(dfexp))
gtList = unique(colnames(dfexp))
geneList = rownames(dfexp)
HbOutput <- pbapply::pblapply(seq_along(geneList),  function(i){
  trait = geneList[i]
  gexp = as.numeric(dfexp[trait, ])
  model = lmer(gexp ~ (1|group))
  varComp <- as.data.frame(VarCorr(model,comp='vcov'))
  sigmas <- varComp$vcov
  H2 <- sigmas[1]/sum(sigmas)
  H2
},cl=32)

Hb = unlist(HbOutput)
gexpMeanlog <- log2(rowMeans(dfexp)) #[sapply(dfexp,is.numeric)]))
dfHb = data.frame(geneID = geneList, Hb = Hb, log2GEMean = gexpMeanlog)
write.csv(dfHb, file = 'e340_Heritability.csv')

#########################################################################################
library(ggplot2)
library("ggpubr")
dfHb <- read.csv('e340_Heritability.csv', row.names = 1, header = T)
p1 <- ggplot(dfHb, aes(x=Hb)) + 
  geom_histogram(bins = 20, color='black', fill='grey') +
  xlab(bquote(~H^2))+ 
  ylab('Number of genes') +
  theme_bw()+ 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
p1
#p1 <- p1 + theme(axis.text.y = element_text(face = 'bold', size = 10, colour = 'black'), axis.text.x = element_text(face='bold', size = 10, colour = 'black'))
#p1 <- p1 + theme(axis.title.y = element_text(face = 'bold', size=10,  colour = 'black'), title = element_text(face='bold', size=10, colour = 'black'))
#p1 <- p1 + theme(axis.ticks.length = unit(.25,"cm"), axis.ticks = element_line(size=1,colour = 'black'))


p2 <- ggplot(dfHb, aes(x=log2GEMean, y=Hb) ) +
  stat_density_2d(aes(fill = ..level..), geom='polygon')+
  scale_fill_distiller(palette=12, direction=1) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(legend.position='right')+
  xlab(bquote(~log[2](FPKM))) +
  ylab(bquote(~H^2))
p2
#p2 <- p2 + theme(axis.text.y = element_text(face = 'bold', size = 10, colour = 'black'), axis.text.x = element_text(face='bold', size = 10, colour = 'black'))
#p2 <- p2 + theme(axis.title.y = element_text(face = 'bold', size=10,  colour = 'black'), title = element_text(face='bold', size=10, colour = 'black'))
#p2 <- p2 + theme(axis.ticks.length = unit(.15,"cm"), axis.ticks = element_line(size=0.5,colour = 'black'))
library(ggrepel)
library(gridExtra)
grid.arrange(p1,p2, nrow=1)


#################### Narrow sense heritability (use sommer) #################
rm(list=ls())
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
    library(data.table)
    library(pbapply)
    library(snpStats)
    library(sommer)
})

# get duplicated individuals
dfexp <- fread("eGWAS340_repUnmerged_Untransformed_fpkm.csv",nrow=1)
#group <- unique(colnames(dfexp))[-1]
group <- unique(colnames(dfexp)[duplicated(colnames(dfexp))])
length(group)

# calculate additive relationship matrix
plink <- read.plink("../../geno340/snps_prune.bed")
snp <- as(plink$genotypes[plink$fam$pedigree %in% group, ], 'numeric')
A <- A.mat(snp - 1)
dim(A)

# get gene expression matrix
f <- fread("eGWAS340_RepMerged_Untransformed_fpkm.csv", data.table = F)
f[1:5,1:5]
phe <- as.matrix(f[, match(group, colnames(f))])
rownames(phe) <- f$V1
traits <- rownames(phe)

# calculate h2
dfh2 <- pblapply(seq_along(traits), function(i){
    geneID <- traits[i]
    df2 <- data.frame(lineName=group, gexp=phe[geneID, ]) 
    q <-mmer(gexp ~ 1, random = ~vsr(lineName, Gu=A), 
             rcov = ~units, data = df2, tolParInv = 1e-3,
             verbose = FALSE)
    a <- vpredict(q, h2 ~ V1/(V1 + V2) )
    a$Estimate
}, cl = 32)

dfh2 <- data.frame(geneID = traits, h2 = unlist(dfh2))
write.csv(dfh2, file='e340_narrow_hb.csv')

######
dfh2 <- read.csv(file = "e340_narrow_hb.csv", row.names = 1)
dfHb <- read.csv(file = 'e340_Heritability.csv', row.names = 1)
dfHb <- subset(dfHb, Hb < 1 & Hb > 0)
dfh2 <- subset(dfh2, h2 < 1 & h2 > 0)

dfh = merge(dfh2, dfHb, by = 'geneID')
head(dfh)


library(ggplot2)

p1 <- ggplot(dfh, aes(x=Hb, y = h2))+
  stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white") +
  geom_smooth(method = lm, color = "red", fill = "#69b3a2", se=TRUE) +
  xlab('Broad sense heritability') +
  ylab('Narrow sense heritability') +
  theme_classic()
p1

p2 <- ggplot(dfh, aes(x=Hb, y = h2)) + 
  geom_point() + 
  geom_smooth(method = lm, color = "red", fill = "#69b3a2", se=TRUE) +
  xlab('Broad sense heritability') +
  ylab('Narrow sense heritability') +
  theme_classic()
p2






