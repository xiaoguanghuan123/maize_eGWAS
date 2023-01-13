
setwd("/work3/maize_eGWAS/figures/Fig3_code/")
suppressPackageStartupMessages({
    library(fastICA)
    library(moments)
    library(matrixStats)
})

rawData<-read.csv("../../eGWAS340/eGWAS340_repUnmerged_Untransformed_fpkm.csv", 
                  row.names = 1, check.names = F)
dfHb <- read.csv('../Fig1_code/e340_Heritability.csv', row.names = 1)

Genekept <- dfHb$geneID[dfHb$Hb>=0.05]
geHb = rawData[Genekept, ] #subset genes with heritability higher than 0.05
ge1 <- subset(geHb, subset=matrixStats::rowSds(as.matrix(geHb))>1) #subset genes which expression sd higher than 1
dim(ge1)

# 16972   572

X <- t(scale(t(ge1)))

####### find optimal number of components #############
ev <- svd(X)$d^2  #perform singular value Decomposition with the raw dataset
Var <- ev/sum(ev) * 100  #calculate the variance explained by each of the components
head(Var)
#10.775900  6.789168  6.200510  3.487104  2.955675  2.463114

#### check how many ICs to calcualte that explaines 80 percent of the expression variation
n.comp1<-which(cumsum(ev/sum(ev))>0.8)[1]
n.comp1
#166

dfsvd <- data.frame('svd' = ev, 'varExplained' = Var, 'sumVarExplained' =cumsum(ev/sum(ev)) )
write.csv(dfsvd, file = 'svd.csv')

################ make screeplot ##################
pdf(file='screeplot.pdf',10,6)
#x11('',10,6)
par(mar=c(3,5,1,5))
plot(NA,xlim=c(1,200),ylim=log10(c(0.01,100)),xaxs='i',yaxs='i',xlab=NA,ylab=NA,yaxt='n')
axis(1,1)
axis(2,seq(-2,2,1),sprintf("%.2f%%",10^seq(-2,2,1)),las=2)
axis(4,seq(-2,2,0.4),paste0(seq(0,100,10),'%'),las=2)
text(-19,0,labels = 'Variance explained by each component',xpd=T,cex=1.2,font = 2,srt=90)
text(218,0,labels = 'Cumulative explained variance',xpd=T,cex=1.2,font = 2,srt=-90)
polygon(x=c(1:200,200,1),y=c(cumsum(Var)[1:200]/100*4-2,-2,-2),col=rgb(0.95,1,0.95),border = NA)
points(cumsum(Var)/100*4-2,col='green',lwd=2,type='l')
points(log10(Var),col='blue',lwd=2,type='l')
abline(v = 166, col="red", lwd=3, lty=2)
text(166,-1.9,labels = '166',cex=1.2,font = 2)
grid(0,5)
dev.off()

##################################
####### Perform ICA #############
set.seed(2021)

ica1 <- fastICA(X, n.comp1)

S = ica1$S
A = ica1$A

colnames(S) <- rownames(A) <- seq_len(n.comp1)
colnames(A) <- colnames(rawData)

write.csv(A, file = 'ICA_matrixA.csv', quote = F)

### Calculate kurtosis for each IC ####################
kurt_S <- apply(S, 2, kurtosis) # 1 means row wise, 2 means column wise
sum(kurt_S >= 6) #97 (only the ICs with kurtosis higher than 6 are kept for downstream)

################## Calculate fdr for each gene of each IC to identify the genes in the same cluster ######
fdr1<-apply(S, 2, function(x) {                                ## column of S is cluster, assign fdr column wise
  fdrtool::fdrtool(x, statistic = 'normal', plot=F, verbose = F)$qval
}) 

mod1<-t(apply(fdr1, 1, function(x) x < 0.01)) ## Get a binary matrix indicate whether a gene should be assigned to a cluster (column)
# colnames(mod1)<-1:n.comp1 ### each column is a IC ID (n.comp1 = 166)
colSums(mod1)[kurt_S>=6]  ### only keep the IC clusters with kurtosis higher than 6 and with more than 10 genes 

sum(apply(mod1, 1, any))

mod1_genes<-apply(mod1, 2, function(x) names(x[x]))
# names(mod1_genes)<-1:n.comp1

cluster1<-t(apply(A,1,function(x) table(kmeans(x,range(x))$cluster))) # using kmean to get two clusters and count number of individuals within each cluster
#cluster1<-t(apply(dfA,1,function(x) table(kmeans(x,range(x))$cluster))) # using kmean to get two clusters and count number of individuals within each cluster

sel <- kurt_S>=6 & rowAlls(cluster1>=10)
sum(sel)
A_sel <- subset(A, sel)
heatmap(A_sel,scale ='none')

write.csv(A_sel, 'ICA_matrixA_sel.csv')

srg<-mod1_genes[sel]  # take clusters that have kurtosis higher than 6 and kmean clusters with both having more than 10 genotypes
sgp<-ica1$A[sel, ] 
rownames(sgp)<- seq(n.comp1)[sel]
colnames(sgp)<-colnames(X)
save(ica1,srg,sgp,file='eGWAS340_ICA.rda')

write.csv(S[, sel], file='ICA_KeptMod.csv', quote=F)

########## plot kurtosis examples
library(ggplot2)
library(ggbreak)
library(patchwork)

dfexample = data.frame(geneID = row.names(S), IC1 = S[,1], IC38 = S[,38])
head(dfexample)
g1 <- ggplot(dfexample, aes(x= IC1)) + geom_histogram(color='black', fill = 'white',bins=50) # hist(S[,36], ylim = c(0,200))
g1 <- g1 + scale_y_break(c(50,200), scales = 0.5)
g1 <- g1 + labs(title = 'Kurtosis = 4.36')

g2 <- ggplot(dfexample, aes(x= IC38)) + geom_histogram(color='black', fill = 'white',bins=50) # hist(S[,36], ylim = c(0,200))
g2 <- g2 + scale_y_break(c(50,200), scales = 0.5)
g2 <- g2 + labs(title = 'Kurtosis = 11.8')

g1 + g2
