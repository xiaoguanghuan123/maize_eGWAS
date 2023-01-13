
####
suppressPackageStartupMessages({
    library(data.table)
    library(GenomicRanges)
    library(parallel)
})
source('../manh.R')
load('../maize_genes.rda')
load('../cis_eqtl_sqtl.rda')
load('../psi340_in_expr_genes.rda')
pp<-c('All','Tropical','Temperate','Before1980','After1990')
col<-c(2:6,RColorBrewer::brewer.pal(5,'Dark2'))

#### pi
pi1<-fread('tropical.windowed.pi')
colnames(pi1)<-c('chr','start','end','snps.trop','pi.trop')
pi1<-GRanges(pi1)
pi2<-fread('temperate.windowed.pi')
colnames(pi2)<-c('chr','start','end','snps.temp','pi.temp')
pi2<-GRanges(pi2)
pi<-pi1[pi1%in%pi2]
pi$snps.temp<-pi2$snps.temp[match(pi,pi2)]
pi$pi.temp<-pi2$pi.temp[match(pi,pi2)]
pi$ratio<-pi$pi.trop/pi$pi.temp
cut_pi<-median(pi$ratio)  # 0.9464664
pi<-pi[pi$ratio>=cut_pi]
rm(pi1,pi2)

#### xpclr
load('sel_res.rda')
xpclr_sel<-res[!is.na(res$xpclr)]
cut_xpclr<-quantile(xpclr_sel$xpclr,0.9) # 6.150673
xpclr_sel$xpclr_norm<- (xpclr_sel$xpclr - mean(xpclr_sel$xpclr))/sd(xpclr_sel$xpclr)
xpclr_sel<-xpclr_sel[xpclr_sel%in%pi & xpclr_sel$xpclr>cut_xpclr]

df1<-data.frame(chr=seqnames(res),pos=start(res),pval=res$xpclr)
manh(df1, cut = cut_xpclr,log=F,ylab='XP-CLR')

xpclr_genes<-subsetByOverlaps(genes,xpclr_sel,type='within')
xpclr_genes$xpclr<-sapply(seq_along(xpclr_genes),function(i) {
    mean(xpclr_sel$xpclr[xpclr_sel%over%xpclr_genes[i]])
})
    
#### fst
fst_sel<-fread('region_fst.windowed.weir.fst')
colnames(fst_sel)<-c('chr','start','end','snps','w_fst','m_fst')
fst_sel<-GRanges(fst_sel)
cut_fst<-quantile(fst_sel$w_fst,0.9)     # 0.1774624  

df2<-data.frame(chr=seqnames(fst_sel),pos=start(fst_sel),pval=fst_sel$w_fst)
manh(df2, cut = cut_fst,log=F,ylab='Fst')

fst_sel<-fst_sel[fst_sel%in%pi & fst_sel$w_fst>=cut_fst]
fst_genes<-subsetByOverlaps(genes,fst_sel,type='within')
fst_genes$fst<-sapply(seq_along(fst_genes),function(i) {
    mean(fst_sel$w_fst[fst_sel%over%fst_genes[i]])
})

#### selection genes
sel_genes<-genes[genes$gene_id%in%xpclr_genes$gene_id & 
                 genes$gene_id%in%fst_genes$gene_id]  
sel_genes$cis_eqtl<-cis1$snps[match(sel_genes$gene_id,cis1$gene)]
sel_genes$pval_eqtl<-cis1$pvalue[match(sel_genes$gene_id,cis1$gene)]
sel_genes$cis_sqtl<-cis2$snps[match(sel_genes$gene_id,cis2$gene)]
sel_genes$pval_sqtl<-cis2$pvalue[match(sel_genes$gene_id,cis2$gene)]
sel_genes$ld<-cis$ld[match(sel_genes$gene_id,cis$gene)]
sel_genes$xpclr<-xpclr_genes$xpclr[match(sel_genes$gene_id,xpclr_genes$gene_id)]
sel_genes$fst<-fst_genes$fst[match(sel_genes$gene_id,fst_genes$gene_id)]

write.csv(sel_genes,'sel_genes_042121.csv')
save(pi,cut_pi,xpclr_sel,xpclr_genes, cut_xpclr, fst_sel,fst_genes,cut_fst,
     sel_genes,df1, df2, file='sel_genes_042121.rda')
