
#######
setwd('/work3/maize_eGWAS/figures/Fig4_code/')
suppressPackageStartupMessages({
    library(data.table)
    library(GenomicRanges)
    library(ggplot2)
    library(ggplotify)
    library(cowplot)
    library(gaston)
})

source('manh.R')
load('../../sweep/sel_genes_042121.rda')
npg_pal<-ggsci::pal_npg()(10)

#### manhattan for region selection
hb <- read.csv("../Fig1_code/Hb_Trop_VS_Temperate.csv")
colnames(hb)[1] <- "gene_id"
ann <- sel_genes
ann$HbTrop <- hb$HbTrop[match(ann$gene_id, hb$gene_id)]
ann$HbTemp <- hb$HbTemp[match(ann$gene_id, hb$gene_id)]

idx <-  ann$HbTemp < ann$HbTrop/2 & !is.na(ann$HbTemp) & ann$symbol!='' & ann$characterized
ann<-ann[idx | ann$gene_id=='Zm00001d041715'| ann$gene_id=='Zm00001d007490']
ann$symbol[ann$gene_id=='Zm00001d041715']<-'vil1'
ann$symbol[ann$gene_id=='Zm00001d007490']<-'pfd2'
ann<-data.frame(ann)[,c('seqnames','start','xpclr','symbol')]

pdf(width = 12,height = 5,file='region_sel_manh.pdf')
par(mfrow=c(2,1))
par(mar=c(0,4,0.5,0.5)+0.5)
manh(df1,log=F,cut=cut_xpclr,bty='l',ylab='XP-CLR',cex=0.2,x.axis = F,yaxs='i',
     ylim=c(0,120),ann=ann,ann.col=1,col=ggsci::pal_npg()(10),
     xaxs='i',xlim=c(-0.03e9,2.24e9),ann.cex = 0.8,cut.col=1)
par(mar=c(2,4,0,0.5)+0.5)
manh(df2,log=F,cut=cut_fst,bty='l',ylab = 'Fst',cex=0.2,yaxs='i',ylim=c(0,0.6),
     col=ggsci::pal_npg()(10),xaxs='i',xlim=c(-0.03e9,2.24e9),cut.col=1)
mtext('Chromosome',1,at=0,line=0.5)
dev.off()


#### pfd2 selection
load('../../cis_eqtl_sqtl.rda')
load('../../sweep/snps_info.rda')

pfd2_eqtl<-cis1[gene=='Zm00001d007490']
pfd2_af<-snps[pfd2_eqtl$snps,]
pfd2_af<-rbind(1-pfd2_af[,c('RAF.trop', 'RAF.temp')],
               pfd2_af[,c('RAF.trop', 'RAF.temp')])
rownames(pfd2_af)<-c(pfd2_eqtl$major,pfd2_eqtl$minor)

vil1_sqtl<-cis2[gene=='Zm00001d041715'][order(pvalue)][1]
vil1_af<-snps[vil1_sqtl$snps,]
vil1_af<-rbind(1-vil1_af[,c('RAF.trop', 'RAF.temp')],
               vil1_af[,c('RAF.trop', 'RAF.temp')])
rownames(vil1_af)<-c(vil1_sqtl$major,vil1_sqtl$minor)

pdf(width = 4.5,height = 3,file = 'pfd2_sel.pdf')
par(mfrow=c(1,1),mar=c(2.5,4,0.5,4),cex.axis=0.9,yaxs='i')
manh(df1,log=F,cut=NULL,bty='u',ylab=NA,col=ggsci::pal_npg()(10)[2],
     chr=2,xlim=c(232.2,233.2),type='l',ylim=c(0,120))
manh(df1,log=F,add=T,type='p',chr=2,col=ggsci::pal_npg()(10)[2])
points(df2$pos[df2$chr==2]/1e6,df2$pval[df2$chr==2]*200,col=ggsci::pal_npg()(10)[2],type='l')
points(df2$pos[df2$chr==2]/1e6,df2$pval[df2$chr==2]*200,col=ggsci::pal_npg()(10)[2],type='p',cex=0.6,pch=4)
abline(v=pfd2_eqtl$pos/1e6,col=2,lwd=1.5,lty=2)
axis(4,at=seq(0,120,20),seq(0,120,20)/200,las=1)
mtext('XP-CLR',2,line=2.5)
mtext('Fst',4,line=2.5)
legend(232.2,120,legend=c('XP-CLR','Fst'),pch=c(20,4),lty=1,
       col=ggsci::pal_npg()(10)[2],bty='n')
dev.off()

#### pfd2 ld map
pdf(width = 3,height = 1.8,file = 'pfd2_ld.pdf')
pfd2<-read.bed.matrix('../../sweep/pfd2_snps')
set.seed(2021)
pfd2<-pfd2[,sort(c(sample(ncol(pfd2),199),which(pfd2@snps$id==pfd2_eqtl$snps)))]
pfd2_ld<-LD(pfd2, c(1,ncol(pfd2)))
pfd2_pos<-pfd2@snps$pos
LD.plot(pfd2_ld,pfd2_pos,write.snp.id = F,write.ld = NULL, draw.chr = F,
        polygon.par = list(border = NA))
segments(1:200,0,(pfd2_pos-232.8e6)/0.2e6*200,10,col='gray',lwd=0.5)
idx<-which(pfd2@snps$id==pfd2_eqtl$snps)
segments(c(1:200)[idx],0,(pfd2@snps$pos[idx]-232.8e6)/0.2e6*200,10,col='red')
axis(3,at=seq(0,200,20),pos=10,labels = FALSE,xpd=T,tck = -0.02)
dev.off()

#### vil1 selection
pdf(width = 4.5,height = 3,file = 'vil1_sel.pdf')
par(mfrow=c(1,1),mar=c(2.5,4,0.5,4),cex.axis=0.9,yaxs='i')
manh(df1,log=F,cut=NULL,bty='u',ylab=NA,col=ggsci::pal_npg()(10)[3],
     chr=3,xlim=c(134,135),type='l',ylim=c(0,70))
manh(df1,log=F,add=T,type='p',chr=3,col=ggsci::pal_npg()(10)[3])
points(df2$pos[df2$chr==3]/1e6,df2$pval[df2$chr==3]*100,col=ggsci::pal_npg()(10)[3],type='l')
points(df2$pos[df2$chr==3]/1e6,df2$pval[df2$chr==3]*100,col=ggsci::pal_npg()(10)[3],type='p',cex=0.6,pch=4)
abline(v=vil1_sqtl$pos/1e6,col=2,lwd=1.5,lty=2)
axis(4,at=seq(0,70,10),seq(0,70,10)/100,las=1)
mtext('XP-CLR',2,line=2.5)
mtext('Fst',4,line=2.5)
legend(134,70,legend=c('XP-CLR','Fst'),pch=c(20,4),lty=1,
       col=ggsci::pal_npg()(10)[3],bty='n')
dev.off()

#### vil1 ld map
pdf(width = 3,height = 1.8,file = 'vil1_ld.pdf')
vil1<-read.bed.matrix('../../sweep/vil1_snps')
set.seed(2021)
vil1<-vil1[,sort(c(sample(ncol(vil1),199),which(vil1@snps$id==vil1_sqtl$snps)))]
vil1_ld<-LD(vil1, c(1,ncol(vil1)))
vil1_pos<-vil1@snps$pos
LD.plot(vil1_ld,vil1_pos,write.snp.id = F,write.ld = NULL, draw.chr = F,
        polygon.par = list(border = NA))
segments(1:200,0,(vil1_pos-134.6e6)/0.2e6*200,10,col='gray',lwd=0.5)
idx<-which(vil1@snps$id==vil1_sqtl$snps)
segments(c(1:200)[idx],0,(vil1_pos[idx]-134.6e6)/0.2e6*200,10,col='red')
axis(3,at=seq(0,200,20),pos=10,labels = FALSE,xpd=T,tck = -0.02)
dev.off()

#### Allele frequency
hb<-read.csv('../Fig1_code/Hb_Trop_VS_Temperate.csv')
colnames(hb)[1] <- "gene_id"
pfd2_hb<-round(as.numeric(hb[hb$gene_id=='Zm00001d007490',c('HbTrop','HbTemp')]),2)
vil1_hb<-round(as.numeric(hb[hb$gene_id=='Zm00001d041715',c('HbTrop','HbTemp')]),2)
pdf(width = 1.6,height = 5,file='allele_frequency.pdf')
par(mfrow=c(2,1),mar=c(3,3.5,1,0.5)+0.5,cex=0.9,cex.axis=0.9)
bp<-barplot(pfd2_hb,col=npg_pal[3],ylim=c(0,1),las=1)
text(bp,-0.05,labels = c('Tropical','Temperate'),srt=45,xpd=T,adj = c(1,0.5),cex=0.9)
abline(h=0)
mtext(bquote(H^2),2,line=2.5,cex=0.9)
bp<-barplot( as.matrix(pfd2_af),col=npg_pal[2:1],xaxt='n',legend.text =T,las=1,
             args.legend = list(x=1.3,y=1,xpd=T,bty='n',ncol=2,xjust=0.5,yjust=0))
text(bp,-0.05,labels = c('Tropical','Temperate'),srt=45,xpd=T,adj = c(1,0.5),cex=0.9)
abline(h=0)
mtext('Allele frequency',2.5,line=2.5,cex=0.9)
bp<-barplot(vil1_hb,col=npg_pal[3],ylim=c(0,1),las=1)
text(bp,-0.05,labels = c('Tropical','Temperate'),srt=45,xpd=T,adj = c(1,0.5),cex=0.9)
abline(h=0)
mtext(bquote(H^2),2,line=2.5,cex=0.9)
bp<-barplot( as.matrix(vil1_af),col=npg_pal[2:1],xaxt='n',legend.text =T,las=1,
             args.legend = list(x=1.3,y=1,xpd=T,bty='n',ncol=2,xjust=0.5,yjust=0))
text(bp,-0.05,labels = c('Tropical','Temperate'),srt=45,xpd=T,adj = c(1,0.5),cex=0.9)
abline(h=0)
mtext('Allele frequency',2,line=2.5,cex=0.9)
dev.off()






