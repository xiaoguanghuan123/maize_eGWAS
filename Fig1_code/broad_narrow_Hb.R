
################### Hb or all lines ############################
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
    library(data.table)
    library(pbapply)
    library('lme4')
})

dfexp <- fread('eGWAS340_repUnmerged_Untransformed_fpkm.csv')
dfexp[1:5,1:5]
group <- unique(colnames(dfexp))[-1]
#group <- unique(colnames(dfexp)[duplicated(colnames(dfexp))])
length(group)

geneList <- dfexp$V1
gt <- factor(colnames(dfexp)[colnames(dfexp) %in% group])
dfexp <- as.matrix(subset(dfexp, select = colnames(dfexp) %in% group))
rownames(dfexp) <- geneList
dim(dfexp)

HbOutput <- pblapply(seq_along(geneList),  function(i){
    df1 <- data.frame(gexp = dfexp[i, ], gt = gt)
    model = lmer(gexp ~ (1|gt), data = df1)
    sigmas <- as.data.frame(VarCorr(model))[,"vcov"]
    H2 <- sigmas[1]/sum(sigmas)
    H2
},cl=32)

Hb <- unlist(HbOutput)
gexpMeanlog <- log2(rowMeans(dfexp))
dfHb <- data.frame(geneID = geneList, Hb = Hb, log2GEMean = gexpMeanlog)
write.csv(dfHb, file = 'e219_Hb.csv')

#### tropical lines
trop <- scan("../../sweep/tropical.txt", what='')
trop_gt <- droplevels(gt[gt%in%trop])
trop_exp <- dfexp[, gt%in%trop]
dim(trop_exp)
tropHb <- pblapply(seq_along(geneList),  function(i){
    df1 <- data.frame(gexp = trop_exp[i, ], gt = trop_gt)
    model = lmer(gexp ~ (1|gt), data = df1)
    sigmas <- as.data.frame(VarCorr(model))[,"vcov"]
    H2 <- sigmas[1]/sum(sigmas)
    H2
},cl=32)

tropHb <- data.frame(geneID=geneList, Hb=unlist(tropHb))
write.csv(tropHb, file = 'trop22_Hb.csv')

#### tropical lines
temp <- scan("../../sweep/temperate.txt", what='')
temp_gt <- droplevels(gt[gt%in%temp])
temp_exp <- dfexp[, gt%in%temp]
tempHb <- pblapply(seq_along(geneList),  function(i){
    df1 <- data.frame(gexp = temp_exp[i, ], gt = temp_gt)
    model = lmer(gexp ~ (1|gt), data = df1)
    sigmas <- as.data.frame(VarCorr(model))[,"vcov"]
    H2 <- sigmas[1]/sum(sigmas)
    H2
},cl=32)
tempHb <- data.frame(geneID=geneList, Hb=unlist(tempHb))
write.csv(tempHb, file = 'temp39_Hb.csv')

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
write.csv(dfh2, file='e219_h2.csv')

###################### plot ############################
library(ggplot2)
x <- read.csv("e219_Hb.csv", row.names = 1)
y <- read.csv('e219_h2.csv', row.names = 1)

x <- subset(x, Hb > 0 & Hb < 1)
y <- subset(y, h2 > 0 & h2 < 1)

dfh = merge(x, y, by = 'geneID')
head(dfh)

p1 <- ggplot(dfh, aes(x= h2, y = Hb)) + ylim(0,1) +
    stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white") +
    geom_smooth(method = lm, color = "red", fill = "#69b3a2", se=TRUE) +
    ylab('Broad sense heritability') +
    xlab('Narrow sense heritability') +
    theme_classic()
p1

p2 <- ggplot(dfh, aes(x= h2, y = Hb)) + 
  geom_point() + ylim(0,1) +
  geom_smooth(method = lm, color = "red", fill = "#69b3a2", se=TRUE) +
  ylab('Broad sense heritability') +
  xlab('Narrow sense heritability') +
  theme_classic()
p2

pdf(width = 12, height = 6, file='../Fig_Hb_vs_h2.pdf')
cowplot::plot_grid(p1 + geom_abline(slope = 1,intercept = 0, lty=2, lwd=1.2, color=2), 
                   p2 + geom_abline(slope = 1,intercept = 0, lty=2, lwd=1.2, color=2), 
                   labels = 'auto')
dev.off()

