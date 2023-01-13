
####
setwd('/work3/maize_eGWAS/figures/Fig5_code/')
suppressPackageStartupMessages({
    library(data.table)
    library(GenomicRanges)
})

npg_pal<-ggsci::pal_npg()(10)

#### load data
res<-read.csv('GWAS_RES_SUMMARY.csv')
res<-res[res$eTrait=='Y' & res$SelectionSignal!='N',]
pheno<-read.csv('CandidateGene_exp_phe_subpop.csv')
ecis<-read.csv('top_CisEQTL_AF_Subpop.csv')
hb<-read.csv('../Fig1_code/Hb_Trop_VS_Temperate.csv')
colnames(hb)[1]<-'gene_id'

load('../../sweep/sel_genes_042121.rda')
load('../../sweep/snps_info.rda')
load('../../maize_genes.rda')
load('../../cis_eqtl_sqtl.rda')
source('../../manh.R')
df<-fread('Endosperm_Color.assoc.txt')
(cut<-0.05/nrow(df))
gr<-genes[genes$gene_id=='Zm00001d036345']
df<-df[df$chr==6 & df$ps>start(gr)-1e6 & df$ps<start(gr)+1e6,]
df<-data.frame(chr=df$chr,pos=df$ps,p=df$p_lrt)

#### gwas plot
pdf(width = 12,height = 3,file='gwas.pdf')
layout(matrix(rep(1:4,c(2,1,1,1)),nrow=1))
par(mar=c(2.5,4,1.5,4),yaxs='i',cex=0.95,cex.axis=0.95,cex.main=0.95,
    bty='l')
manh(df,col=0,ylim = c(0,72),chr=6,xlim=c(84.8,85.2),cut=cut,
     main='GWAS for endosperm color',font.main=1,ylab=NA)
mtext(expression(-log10(italic(P))),2,line = 2.2)
rect(start(gr)/1e6,0,end(gr)/1e6,71,col=5,border=5)
manh(df,col=1,ylim = c(0,72),chr=6,xlim=c(84.8,85.2),cut=cut,
     main='GWAS for endosperm color',font.main=1,ylab=NA, add=T)
mtext('chr6',1,at=84.75,line=0.5)
mtext('Mb',1,at=85.24,line=0.5)
x<- df$pos[which.min(df$p)]
y<- -log10(min(df$p))
points(x/1e6,y,pch=23, bg=2, cex=1.2)
#text(x/1e6,y,expression(QTL~peak%->%""),adj = c(1.1,0.5),xpd=T,col=2)

par(mar=c(2,3,1,0.5)+0.5)
p<-t.test(pheno$Endosperm_Color[pheno$subpop=='tropical'],
           pheno$Endosperm_Color[pheno$subpop=='temperate'])$p.value
p <- signif(p,3)
boxplot(pheno$Endosperm_Color[pheno$subpop=='tropical'],
        pheno$Endosperm_Color[pheno$subpop=='temperate'],outline=T,ylim=c(0.5,5.5),
        col=npg_pal[3],ylab='',las=1, names=NULL,yaxs='i',xaxt='n',boxwex=0.6,medlwd=1,
        xlab=NULL,main=bquote(italic(P)==.(p)))
axis(1,at=1:2,labels = c('Tropical','Temperate'),cex.axis=0.9)
mtext('Endosperm color',2,line=2.5)
p<-t.test(pheno$Peiffer_2014_GDD_DTA[pheno$subpop=='tropical'],
          pheno$Peiffer_2014_GDD_DTA[pheno$subpop=='temperate'])$p.value
p<-signif(p, 3)
boxplot(pheno$Peiffer_2014_GDD_DTA[pheno$subpop=='tropical']/24,
        pheno$Peiffer_2014_GDD_DTA[pheno$subpop=='temperate']/24,outline=T,
        ylim=c(50,110),medlwd=1,
        col=npg_pal[3],ylab=NULL,las=1, names=NULL,yaxs='i',xaxt='n',boxwex=0.6,
        xlab=NULL,main=bquote(italic(P)==.(p)) )
axis(1,at=1:2,labels = c('Tropical','Temperate'),cex.axis=0.9)
mtext('Flowering time (d)',2,line=2.5)
p<-t.test(pheno$Upper_Kernel_Shape[pheno$subpop=='tropical'],
          pheno$Upper_Kernel_Shape[pheno$subpop=='temperate'])$p.value
p <- signif(p,3)
boxplot(pheno$Upper_Kernel_Shape[pheno$subpop=='tropical'],
        pheno$Upper_Kernel_Shape[pheno$subpop=='temperate'],outline=T,ylim=c(0.5,6.5),
        col=npg_pal[3],ylab='',las=1, names=NULL,yaxs='i',xaxt='n',boxwex=0.6,medlwd=1,
        xlab=NULL,main=bquote(italic(P)==.(p)) )
axis(1,at=1:2,labels = c('Tropical','Temperate'),cex.axis=0.9)
mtext('Upper kernel shape',2,line=2.5)
dev.off()


#### y1
pdf(width = 12,height = 3,file = 'y1_sel.pdf')
gr<-genes[genes$gene_id=='Zm00001d036345']
af<-snps[snps$snp.nam==cis1$snps[cis1$gene=='Zm00001d036345'],]
layout(matrix(rep(1:4,c(2,1,1,1)),nrow=1))
par(mar=c(2.5,4,1.5,4),yaxs='i',cex=0.95,cex.axis=0.95,cex.main=0.95,bty='l')
manh(df1,log=F,cut=NULL,bty='u',ylab=NA,col='gray',chr=6,xlim=c(84.8,85.2),
     type='n',ylim=c(0,40))
rect(start(gr)/1e6,0,end(gr)/1e6,40,col=5,border=5)
manh(df1,log=F,add=T,type='l',chr=6,col='gray')
manh(df1,log=F,add=T,type='p',chr=6,col='gray')
points(df2$pos[df2$chr==6]/1e6,df2$pval[df2$chr==6]*100,col='gray',type='l')
points(df2$pos[df2$chr==6]/1e6,df2$pval[df2$chr==6]*100,col='gray',type='p',cex=0.7,pch=4)
axis(4,at=seq(0,120,10),seq(0,120,10)/100,las=1)
mtext('XP-CLR',2,line=2.5)
mtext('Fst',4,line=2.5)
legend(85,40,legend=c('XP-CLR   ','Fst'),pch=c(20,4),lty=1,
       col='gray',bty='n',ncol=2,xpd=T,xjust = 0.5, yjust = 0)
mtext('chr6',1,at=84.75,line=0.5)
mtext('Mb',1,at=85.24,line=0.5)
pe<-signif(cis1$pvalue[cis1$gene=='Zm00001d036345'],3)
points(af$pos/1e6,1,pch=25,bg=2,col=2,xpd=T)
legend(84.8,40,legend=substitute(eQTL~(italic(P)==x),list(x=pe)),pch=25,
       col=2,bty='n',pt.bg=2, text.col=2)

par(mar=c(2,3,1,0.5)+0.5)
bp<-barplot(as.numeric(hb[hb$gene_id=='Zm00001d036345',c('HbTrop','HbTemp')]),
            col=npg_pal[3],las=1,ylim=c(0,1),width = 0.6, space = 4/6, 
            xlim = c(0.2,2.2))
axis(1,at=bp,labels = c('Tropical','Temperate'),cex.axis=0.9)
mtext(expression(H^2~of~italic('y1')~expression),2,line=2.5)
axis(1,pos=0,c(0,3))
p<-t.test(pheno$Zm00001d036345[pheno$subpop=='tropical'],
           pheno$Zm00001d036345[pheno$subpop=='temperate'])$p.value
p<-signif(p,3)
boxplot(pheno$Zm00001d036345[pheno$subpop=='tropical'],
        pheno$Zm00001d036345[pheno$subpop=='temperate'],outline=T,ylim=c(0,8),
        col=npg_pal[3],las=1, names=NULL,yaxs='i',xaxt='n',boxwex=0.6,medlwd=1,
        xlab=NULL,ylab=NULL,main=bquote(italic(P)==.(p)) )
axis(1,at=1:2,labels = c('Tropical','Temperate'),cex.axis=0.9)
mtext(expression(italic('y1')~root~expression),2,line=2.2)
x<-rbind(1-af[,c('RAF.trop', 'RAF.temp')],
         af[,c('RAF.trop', 'RAF.temp')])
rownames(x)<-c(af[,'a2'],af[,'a1'])
bp<-barplot( as.matrix(x),col=npg_pal[2:1],xaxt='s',legend.text =T,
             args.legend = list(x=1.3,y=1,xpd=T,bty='n',ncol=2,xjust=0.5,yjust=0),
             xaxt='n',las=1,width = 0.6, space = 4/6, xlim = c(0.2,2.2))
axis(1,at=bp,labels = c('Tropical','Temperate'),cex.axis=0.9,pos=0)
mtext(substitute(eQTL~(chr*x:y),list(x=af$chr,y=format(af$pos,big.mark=","))),2,line=2.5)
axis(1,pos=0,c(0,3))
dev.off()


#### sweet2
pdf(width = 12,height = 3,file = 'sweet2_sel.pdf')
gr<-genes[genes$gene_id=='Zm00001d009365']
af<-snps[snps$snp.nam==cis1$snps[cis1$gene=='Zm00001d009365'],]
layout(matrix(rep(1:4,c(2,1,1,1)),nrow=1))
par(mar=c(2.5,4,1.5,4),yaxs='i',cex=0.95,cex.axis=0.95,cex.main=0.95,bty='l')
manh(df1,log=F,cut=NULL,bty='u',ylab=NA,col="gray",
     chr=8,xlim=c(59.2,60.2),type='n',ylim=c(0,15))
rect(start(gr)/1e6,0,end(gr)/1e6,15,col=5,border = 5)
manh(df1,log=F,add=T,type='l',chr=8,col="gray")
manh(df1,log=F,add=T,type='p',chr=8,col="gray")
points(df2$pos[df2$chr==8]/1e6,df2$pval[df2$chr==8]*30,col="gray",type='l')
points(df2$pos[df2$chr==8]/1e6,df2$pval[df2$chr==8]*30,col="gray",type='p',cex=0.7, pch=4)
axis(4,at=seq(0,15,3),seq(0,15,3)/30,las=1)
mtext('XP-CLR',2,line=2.5)
mtext('Fst',4,line=2.5)
legend(59.7,15,legend=c('XP-CLR   ','Fst'),pch=c(20,4),lty=1,
       col='gray',bty='n',ncol=2,xpd=T,xjust = 0.5, yjust = 0)
mtext('chr8',1,at=59.05,line=0.5)
mtext('Mb',1,at=60.3,line=0.5)
pe<-signif(cis1$pvalue[cis1$gene=='Zm00001d009365'],3)
ps<-signif(cis$pvalue.sqtl[cis$gene=='Zm00001d009365'],3)
points(af$pos/1e6,0.5,pch=25,bg=2,col=2,xpd=T)
points(cis$pos.sqtl[cis$gene=='Zm00001d009365']/1e6,0.5,pch=25,bg=1,col=1)
legend(59.2,15,legend = c(as.expression(substitute(eQTL~(italic(P)==x),list(x=pe))),
                        as.expression(substitute(sQTL~(italic(P)==x),list(x=ps)))),
       pch=25,col=c(2,1),bty='n',pt.bg=c(2,1), text.col=c(2,1),y.intersp=1.2)

par(mar=c(2,3,1,0.5)+0.5)
bp<-barplot(as.numeric(hb[hb$gene_id=='Zm00001d009365',c('HbTrop','HbTemp')]),
            col=npg_pal[3],las=1,ylim=c(0,1),width = 0.6, space = 4/6, xlim = c(0.2,2.2))
axis(1,at=bp,labels = c('Tropical','Temperate'),cex.axis=0.9)
mtext(expression(H^2~of~italic('sweet2')~expression),2,line=2.5)
axis(1,pos=0,c(0,3))
p<-t.test(pheno$Zm00001d009365[pheno$subpop=='tropical'],
           pheno$Zm00001d009365[pheno$subpop=='temperate'])$p.value
p<-signif(p, 3)
boxplot(pheno$Zm00001d009365[pheno$subpop=='tropical'],
        pheno$Zm00001d009365[pheno$subpop=='temperate'],outline=T,ylim=c(0,30),
        col=npg_pal[3],las=1, names=NULL,yaxs='i',xaxt='n', boxwex=0.6,medlwd=1,
        xlab=NULL,ylab=NULL,main=bquote(italic(P)==.(p)) )
axis(1,at=1:2,labels = c('Tropical','Temperate'),cex.axis=0.9)
mtext(expression(italic('sweet2')~root~expression),2,line=2.2)
x<-rbind(1-af[,c('RAF.trop', 'RAF.temp')],
         af[,c('RAF.trop', 'RAF.temp')])
rownames(x)<-c(af[,'a2'],af[,'a1'])
bp<-barplot( as.matrix(x),col=npg_pal[2:1],xaxt='s',legend.text =T,
             args.legend = list(x=1.3,y=1,xpd=T,bty='n',ncol=2,xjust=0.5,yjust=0),
             xaxt='n',las=1,width = 0.6, space = 4/6, xlim = c(0.2,2.2))
axis(1,at=bp,labels = c('Tropical','Temperate'),cex.axis=0.9,pos=0)
mtext(substitute(eQTL~(chr*x:y),list(x=af$chr,y=format(af$pos,big.mark=","))),2,line=2.5)
axis(1,pos=0,c(0,3))
dev.off()


#### mcm4
pdf(width = 12,height = 3,file = 'mcm4_sel.pdf')
gr<-genes[genes$gene_id=='Zm00001d009374']
af<-snps[snps$snp.nam==cis1$snps[cis1$gene=='Zm00001d009374'],]
layout(matrix(rep(1:4,c(2,1,1,1)),nrow=1))
par(mar=c(2.5,4,1.5,4),yaxs='i',cex=0.95,cex.axis=0.95,cex.main=0.95,bty='l')
manh(df1,log=F,cut=NULL,bty='u',ylab=NA,col="gray",
     chr=8,xlim=c(60.6,61),type='n',ylim=c(0,20))
rect(start(gr)/1e6,0,end(gr)/1e6,40,col=5,border =5)
manh(df1,log=F,add=T,type='l',chr=8,col="gray")
manh(df1,log=F,add=T,type='p',chr=8,col="gray")
points(df2$pos[df2$chr==8]/1e6,df2$pval[df2$chr==8]*40,col="gray",type='l')
points(df2$pos[df2$chr==8]/1e6,df2$pval[df2$chr==8]*40,col="gray",type='p',cex=0.7, pch=4)
axis(4,at=seq(0,20,4),seq(0,20,4)/40,las=1)
mtext('XP-CLR',2,line=2.5)
mtext('Fst',4,line=2.5)
legend(60.8,20,legend=c('XP-CLR   ','Fst'),pch=c(20,4),lty=1,
       col='gray',bty='n',ncol=2,xpd=T,xjust = 0.5, yjust = 0)
mtext('chr8',1,at=60.55,line=0.5)
mtext('Mb',1,at=61.04,line=0.5)
pe<-signif(cis1$pvalue[cis1$gene=='Zm00001d009374'],3)
points(af$pos/1e6,0.5,pch=25,bg=2,col=2)
legend(60.6,20,legend=substitute(eQTL~(italic(P)==x),list(x=pe)),pch=25,
       col=2,bty='n',pt.bg=2, text.col=2)

par(mar=c(2,3,1,0.5)+0.5)
bp<-barplot(as.numeric(hb[hb$gene_id=='Zm00001d009374',c('HbTrop','HbTemp')]),
            col=npg_pal[3],las=1,ylim=c(0,1),width = 0.6, space = 4/6, xlim = c(0.2,2.2))
axis(1,at=bp,labels = c('Tropical','Temperate'),cex.axis=0.9)
mtext(expression(H^2~of~italic('mcm4')~expression),2,line=2.5)
axis(1,pos=0,c(0,3))
p<-t.test(pheno$Zm00001d009374[pheno$subpop=='tropical'],
           pheno$Zm00001d009374[pheno$subpop=='temperate'])$p.value
p<-signif(p,3)
boxplot(pheno$Zm00001d009374[pheno$subpop=='tropical'],
        pheno$Zm00001d009374[pheno$subpop=='temperate'],outline=T,ylim=c(0,10),
        col=npg_pal[3],las=1, names=NULL,yaxs='i',xaxt='n',boxwex=0.6,medlwd=1,
        xlab=NULL,ylab=NULL,main=bquote(italic(P)==.(p)) )
axis(1,at=1:2,labels = c('Tropical','Temperate'),cex.axis=0.9)
mtext(expression(italic('mcm4')~root~expression),2,line=2.2)
x<-rbind(1-af[,c('RAF.trop', 'RAF.temp')],
         af[,c('RAF.trop', 'RAF.temp')])
rownames(x)<-c(af[,'a2'],af[,'a1'])
bp<-barplot( as.matrix(x),col=npg_pal[2:1],xaxt='s',legend.text =T,las=1,
             args.legend = list(x=1.3,y=1,xpd=T,bty='n',ncol=2,xjust=0.5,yjust=0),
             xaxt='n',width = 0.6, space = 4/6, xlim = c(0.2,2.2))
axis(1,at=bp,labels = c('Tropical','Temperate'),cex.axis=0.9,pos=0)
mtext(substitute(eQTL~(chr*x:y),list(x=af$chr,y=format(af$pos,big.mark=","))),2,line=2.5)
axis(1,pos=0,c(0,3))
dev.off()


#### Aldolase
pdf(width = 12,height = 3,file = 'aldolase_sel.pdf')
gr<-genes[genes$gene_id=='Zm00001d049559']
af<-snps[snps$snp.nam==cis1$snps[cis1$gene=='Zm00001d049559'],]
layout(matrix(rep(1:4,c(2,1,1,1)),nrow=1))
par(mar=c(2.5,4,1.5,4),yaxs='i',cex=0.95,cex.axis=0.95,cex.main=0.95,bty='l')
manh(df1,log=F,cut=NULL,bty='u',ylab=NA,col="gray",
     chr=4,xlim=c(34.8,35.4),type='l',ylim=c(0,15))
rect(start(gr)/1e6,0,end(gr)/1e6,40,col=5,border = 5)
manh(df1,log=F,add=T,type='l',chr=4,col="gray")
manh(df1,log=F,add=T,type='p',chr=4,col="gray")
points(df2$pos[df2$chr==4]/1e6,df2$pval[df2$chr==4]*60,col="gray",type='l')
points(df2$pos[df2$chr==4]/1e6,df2$pval[df2$chr==4]*60,col="gray",type='p',cex=0.7, pch=4)
axis(4,at=seq(0,15,3),seq(0,15,3)/60,las=1)
mtext('XP-CLR',2,line=2.5)
mtext('Fst',4,line=2.5)
legend(35.1,15,legend=c('XP-CLR   ','Fst'),pch=c(20,4),lty=1,
       col='gray',bty='n',ncol=2,xpd=T,xjust = 0.5, yjust = 0)
mtext('chr4',1,at=34.73,line=0.5)
mtext('Mb',1,at=35.47,line=0.5)
pe<-signif(cis1$pvalue[cis1$gene=='Zm00001d049559'],3)
points(af$pos/1e6,0.5,pch=25,bg=2,col=2)
legend(34.8,15,legend=substitute(eQTL~(italic(P)==x),list(x=pe)),pch=25,
       col=2,bty='n',pt.bg=2, text.col=2)

par(mar=c(2,3,1,0.5)+0.5)
bp<-barplot(as.numeric(hb[hb$gene_id=='Zm00001d049559',c('HbTrop','HbTemp')]),
            col=npg_pal[3],las=1,ylim=c(0,1),width = 0.6, space = 4/6, 
            las=1,xlim = c(0.2,2.2))
axis(1,at=bp,labels = c('Tropical','Temperate'),cex.axis=0.9)
mtext(expression(H^2~of~italic('aldolase')~expression),2,line=2.5)
axis(1,pos=0,c(0,3))
p<-t.test(pheno$Zm00001d049559[pheno$subpop=='tropical'],
           pheno$Zm00001d049559[pheno$subpop=='temperate'])$p.value
p<-signif(p, 3)
boxplot(pheno$Zm00001d049559[pheno$subpop=='tropical'],
        pheno$Zm00001d049559[pheno$subpop=='temperate'],outline=T,ylim=c(0,5),
        col=npg_pal[3],las=1, names=NULL,yaxs='i',xaxt='n',boxwex=0.6,medlwd=1,
        xlab=NULL,ylab=NULL,main=bquote(italic(P)==.(p)) )
axis(1,at=1:2,labels = c('Tropical','Temperate'),cex.axis=0.9)
mtext(expression(italic('aldolase')~root~expression),2,line=2.2)
x<-rbind(1-af[,c('RAF.trop', 'RAF.temp')],
         af[,c('RAF.trop', 'RAF.temp')])
rownames(x)<-c(af[,'a2'],af[,'a1'])
bp<-barplot( as.matrix(x),col=npg_pal[2:1],xaxt='s',legend.text =T,
             args.legend = list(x=1.3,y=1,xpd=T,bty='n',ncol=2,xjust=0.5,yjust=0),
             xaxt='n',width = 0.6, space = 4/6, xlim = c(0.2,2.2),las=1)
axis(1,at=bp,labels = c('Tropical','Temperate'),cex.axis=0.9,pos=0)
mtext(substitute(eQTL~(chr*x:y),list(x=af$chr,y=format(af$pos,big.mark=","))),2,line=2.5)
axis(1,pos=0,c(0,3))
dev.off()


#### dsc1
pdf(width = 12,height = 3,file = 'dsc1_sel.pdf')
gr<-genes[genes$gene_id=='Zm00001d049872']
af<-snps[snps$snp.nam==cis1$snps[cis1$gene=='Zm00001d049872'],]
layout(matrix(rep(1:4,c(2,1,1,1)),nrow=1))
par(mar=c(2.5,4,1.5,4),yaxs='i',cex=0.95,cex.axis=0.95,cex.main=0.95,bty='l')
manh(df1,log=F,cut=NULL,bty='u',ylab=NA,col="gray",
     chr=4,xlim=c(49.4,50.0),type='l',ylim=c(0,15))
rect(start(gr)/1e6,0,end(gr)/1e6,40,col=5,border = 5)
manh(df1,log=F,add=T,type='l',chr=4,col="gray")
manh(df1,log=F,add=T,type='p',chr=4,col="gray")
points(df2$pos[df2$chr==4]/1e6,df2$pval[df2$chr==4]*50,col="gray",type='l')
points(df2$pos[df2$chr==4]/1e6,df2$pval[df2$chr==4]*50,col="gray",type='p',cex=0.7, pch=4)
axis(4,at=seq(0,15,3),seq(0,15,3)/50,las=1)
mtext('XP-CLR',2,line=2.5)
mtext('Fst',4,line=2.5)
legend(49.7,15,legend=c('XP-CLR   ','Fst'),pch=c(20,4),lty=1,
       col='gray',bty='n',ncol=2,xpd=T,xjust = 0.5, yjust = 0)
mtext('chr4',1,at=49.33,line=0.5)
mtext('Mb',1,at=50.07,line=0.5)
pe<-signif(cis1$pvalue[cis1$gene=='Zm00001d049872'],3)
points(af$pos/1e6,0.5,pch=25,bg=2,col=2)
legend(49.4,13,legend=substitute(eQTL~(italic(P)==x),list(x=pe)),pch=25,
       col=2,bty='n',pt.bg=2, text.col=2)

par(mar=c(2,3,1,0.5)+0.5)
bp<-barplot(as.numeric(hb[hb$gene_id=='Zm00001d049872',c('HbTrop','HbTemp')]),
            col=npg_pal[3],las=1,ylim=c(0,1),width = 0.6, space = 4/6, 
            las=1,xlim = c(0.2,2.2))
axis(1,at=bp,labels = c('Tropical','Temperate'),cex.axis=0.9)
mtext(expression(H^2~of~italic('dsc1')~expression),2,line=2.5)
axis(1,pos=0,c(0,3))
p<-t.test(pheno$Zm00001d049872[pheno$subpop=='tropical'],
          pheno$Zm00001d049872[pheno$subpop=='temperate'])$p.value
p <-signif(p, 3)
boxplot(pheno$Zm00001d049872[pheno$subpop=='tropical'],
        pheno$Zm00001d049872[pheno$subpop=='temperate'],outline=T,ylim=c(0,5),
        col=npg_pal[3],las=1, names=NULL,yaxs='i',xaxt='n',boxwex=0.6,medlwd=1,
        xlab=NULL,ylab=NULL,main=bquote(italic(P)==.(p)) )
axis(1,at=1:2,labels = c('Tropical','Temperate'),cex.axis=0.9)
mtext(expression(italic('dsc1')~root~expression),2,line=2.2)
x<-rbind(1-af[,c('RAF.trop', 'RAF.temp')],
         af[,c('RAF.trop', 'RAF.temp')])
rownames(x)<-c(af[,'a2'],af[,'a1'])
bp<-barplot( as.matrix(x),col=npg_pal[2:1],xaxt='s',legend.text =T,
             args.legend = list(x=1.3,y=1,xpd=T,bty='n',ncol=2,xjust=0.5,yjust=0),
             xaxt='n',width = 0.6, space = 4/6, xlim = c(0.2,2.2),las=1)
axis(1,at=bp,labels = c('Tropical','Temperate'),cex.axis=0.9,pos=0)
mtext(substitute(eQTL~(chr*x:y),list(x=af$chr,y=format(af$pos,big.mark=","))),2,line=2.5)
axis(1,pos=0,c(0,3))
dev.off()

