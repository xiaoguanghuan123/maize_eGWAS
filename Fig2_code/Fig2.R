#### 
setwd("/work3/maize_eGWAS/figures/Fig2_code/")
suppressPackageStartupMessages({
    library(ggplotify)
    library(ggplot2)
    library(cowplot)
    library(ggpubr)
    library(data.table)
    library(GenomicRanges)
})
source('locuscompare.R')
npg<-ggsci::pal_npg('nrc')(10) 
theme_set(theme_classic(12) + 
          theme( axis.ticks = element_line(colour = 1),
                 axis.text = element_text(colour = 1,size=rel(0.9)),
                 legend.text = element_text(size=rel(0.9)),
                 line=element_line(size = 0.5),
                 legend.background = element_blank()))

#### Hb vs PVE
hb <- read.csv("../Fig1_code/e340_Heritability.csv", row.names = 1)
load("../../eGWAS340/eqtl_peak.rda")
freq <- table(eqtls$gene)
eqtls <- eqtls[gene%in%names(freq[freq<=10])]
eqtls$Hb <- hb$Hb[match(eqtls$gene, hb$geneID)]

p1<-ggplot(eqtls,aes(x=Hb,y=PVE)) + geom_point(col='gray',cex=0.2) +
    xlab(bquote(H^2)) + ylab("eQTL PVE") +
    geom_abline(slope = 1,intercept = 0,col=npg[8],lty=2) +
    scale_x_continuous(expand = expansion(c(0,0.05))) +
    scale_y_continuous(expand = expansion(c(0,0.05)))

#### PVE/Hb
hb <- hb[hb$geneID %in% eqtls$gene, ]
hb$PVE <- tapply(eqtls$PVE, 
                 factor(eqtls$gene, levels = hb$geneID),
                 sum)
hb$prop <- hb$PVE/hb$Hb

p2<-ggplot(hb,aes(prop)) + 
    geom_density(col=1,fill=1,lwd=1,alpha=0.1) + 
    xlab(bquote("eQTL_PVE/"*H^2)) + ylab("Density") +
    scale_x_continuous(expand = expansion(c(0,0.05)), limits = c(0,1)) +
    scale_y_continuous(expand = expansion(c(0,0.05)))

#### sqtl_loc
load('../../psi340_sd0.01_miss0.05.rda')
cis2_int<-int[match(cis2$intron,int$id)]
d3<-ifelse(as.character(strand(cis2_int))=='+',cis2$pos-start(cis2_int),
           end(cis2_int)-cis2$pos)
names(d3)<-cis2_int$id
d3[d3>=width(cis2_int)]<- (d3-width(cis2_int))[d3>=width(cis2_int)]+2000
d3[d3>0 & d3<width(cis2_int)]<- (d3*2000/width(cis2_int))[d3>0 & d3<width(cis2_int)]
d3<-d3[d3>=-2e3 & d3<=4e3]
y3<-table(findInterval(d3,seq(-2e3,4e3,100)))
y3<-y3[match(seq(60),names(y3))]
y3[is.na(y3)]<-0

p3<-ggplot(data.frame(x=seq(60),y=as.numeric(y3))) +
    geom_line(aes(x=x,y=y),color=npg[2],lwd=0.8) +
    geom_vline(xintercept = c(20,41),col=npg[4],lty=2) + 
    xlab('Distance to affected introns (kb)') +
    ylab(bquote(Number~of~italic("cis")*"-sQTLs")) +
    scale_x_continuous(breaks=c(0,10,20,41,51,61),
                       labels=c(-2,-1,"5'SS","3'SS",1,2)) +
    annotate('text',x=30.5,y=0, label='<-intron->',col=npg[4],vjust=-0.5) +
    theme(axis.text.x = element_text(color=c(1,1,npg[4],npg[4],1,1))) +
    annotate('rect',xmin=20,xmax = 41,ymin=0,ymax = 350,fill=npg[4],alpha=0.1) +
    scale_y_continuous(expand = expansion(c(0,0.05)))

#### tss_tes
load('../../cis_eqtl_sqtl.rda')
y1<-table(findInterval(d1,seq(-2000,4000,100)))
y1<-y1[match(seq(60),names(y1))]
y1[is.na(y1)]<-0; names(y1)<-seq(60)
y2<-table(findInterval(d2,seq(-2000,4000,100)))
y2<-y2[match(seq(60),names(y2))]
y2[is.na(y2)]<-0; names(y2)<-seq(60)

p4<-ggplot(data.frame(x=seq(60),eQTLs=as.numeric(y1),sQTLs=as.numeric(y2))) +
    geom_line(aes(x=x,y=eQTLs,color='eQTLs'),lwd=0.8) + 
    geom_line(aes(x=x,y=sQTLs,color='sQTLs'),lwd=0.8) +
    theme(legend.position = c(1/6,0.85)) +
    labs(x="Distance to affected genes (kb)",
         y=bquote(Number~of~italic("cis")*"-QTLs"),color='') +
    scale_color_manual(values = c('eQTLs'=npg[1],'sQTLs'=npg[2])) + 
    geom_vline(xintercept =c(21,40),color=npg[4],lty=2) +
    scale_x_continuous(breaks=c(1,11,21,40,50,60),
                       labels=c(-2,-1, "TSS","TES",1,2)) +
    annotate('text',x=30.5,y=0,label='<-gene->',color=npg[4],vjust=-0.5) +
    theme(axis.text.x = element_text(color=c(1,1,npg[4],npg[4],1,1))) +
    annotate('rect',xmin=21,xmax = 40,ymin=0,ymax = 350,fill=npg[4],alpha=0.1) +
    scale_y_continuous(expand = expansion(c(0,0.05)))


#### ld_hist
load('../../cis_eqtl_sqtl.rda')
lab<-paste0('italic(r^2)',c('<','>='),0.6,': ',
            c(paste0('"',format(sum(cis$ld<0.6),big.mark = ','),'"'),
              sum(cis$ld>=0.6)))
p5<-ggplot(cis,aes(ld)) + 
    labs(x=bquote(italic(r^2)~(ld)),y='Number of genes') +
    stat_bin(breaks=seq(0,1,0.01),lwd=0.3,
             fill=c(rep(NA,60),rep(npg[5],40)),color=1) +
    geom_vline(xintercept = 0.6, col='red',lty=2) + 
    annotate('text',x=c(0.35,0.82),y=200,label = parse(text=lab),size=4)  +
    scale_x_continuous(expand = expansion(c(0,0.05)),breaks = seq(0,1,0.2)) +
    scale_y_continuous(expand = expansion(c(0,0.05)))

## co-localization
com<-cis[ld>=0.6]
load('../../maize_genes.rda')
genes<-genes[com$gene]
d1<-d1[names(d1)%in%com$gene]
d2<-d2[names(d2)%in%com$intron]

y1<-table(findInterval(d1,seq(-2000,4000,100)))
y1<-y1[match(seq(60),names(y1))]
y1[is.na(y1)]<-0; names(y1)<-seq(60)
y2<-table(findInterval(d2,seq(-2000,4000,100)))
y2<-y2[match(seq(60),names(y2))]
y2[is.na(y2)]<-0; names(y2)<-seq(60)

p6<-ggplot(data.frame(x=seq(60),eQTLs=as.numeric(y1),sQTLs=as.numeric(y2))) +
    geom_line(aes(x=x,y=eQTLs,color='eQTLs'),lwd=0.8) + 
    geom_line(aes(x=x,y=sQTLs,color='sQTLs'),lwd=0.8) +
    theme(legend.position = c(1/6,0.85)) +
    labs(x="Distance to affected genes (kb)",
         y=bquote(Number~of~italic("cis")*"-QTLs"),color='') +
    scale_color_manual(values = c('eQTLs'=npg[1],'sQTLs'=npg[2])) + 
    geom_vline(xintercept =c(21,40),color=npg[4],lty=2) +
    scale_x_continuous(breaks=c(1,11,21,40,50,60),
                       labels=c(-2,-1, "TSS","TES",1,2)) +
    annotate('text',x=30.5,y=0,label='<-gene->',color=npg[4],vjust=-0.5) +
    theme(axis.text.x = element_text(color=c(1,1,npg[4],npg[4],1,1))) +
    annotate('rect',xmin=21,xmax = 40,ymin=0,ymax = 30,fill=npg[4],alpha=0.1) +
    scale_y_continuous(expand = expansion(c(0,0.05)))


#### CT2
load('CT2_eqtl_sqtl.rda')
eqtl<-eqtl[eqtl$pvalue<0.01,]
sqtl<-sqtl[sqtl$pvalue<0.01,]
eqtl$pos<-bed$map$position[match(eqtl$snps,bed$map$snp.name)]/1e6
sqtl$pos<-bed$map$position[match(sqtl$snps,bed$map$snp.name)]/1e6
snp<-eqtl$snps[which.min(eqtl$pvalue)]

ct2_zoom1<-ggplot(eqtl,aes(x=pos,y=-log10(pvalue))) +
    geom_point(fill='gray',color='gray',pch=21,size=1) + 
    ylab(bquote(-log[10](italic(P)))) +
    annotate('text',16,20,label='eQTL') +
    theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
    scale_y_continuous(expand = expansion(c(0,0.05)), limits = c(0,25)) + 
    geom_point(inherit.aes = FALSE,x=eqtl$pos[which.min(eqtl$pvalue)],
               y=max(-log10(eqtl$pvalue)),pch=23,fill='purple',size=2)  + 
    annotate('text',x=eqtl$pos[which.min(eqtl$pvalue)],y=max(-log10(eqtl$pvalue)),
             label=snp,size=3.5,hjust=-0.1) 

ct2_zoom2<-ggplot(sqtl,aes(x=pos,y=-log10(pvalue)))  +
    geom_point(fill='gray',color='gray',pch=21,size=1) + 
    xlab('chr1 (Mb)') + ylab(bquote(-log[10](italic(P)))) +
    annotate('text',16,35,label='sQTL') +
    scale_y_continuous(expand = expansion(c(0,0.05)), limits = c(0,40)) +
    geom_point(inherit.aes = FALSE,x=sqtl$pos[sqtl$snps==snp],
               y=-log10(sqtl$pvalue[sqtl$snps==snp]),pch=23,fill='purple',size=2)  + 
    annotate('text',x=sqtl$pos[sqtl$snps==snp],y=-log10(sqtl$pvalue)[sqtl$snps==snp],
             label=snp,size=3.5,hjust=-0.1)

ct2_zoom<-plot_grid(ct2_zoom1,ct2_zoom2,align = 'v',ncol=1,rel_heights = c(0.8,1)) 

ct2<-locuscompare(eqtl,sqtl,bed=bed, combine = F,region = 1e6,legend_position = "topleft")
ct2_comp<-ct2$locuscompare +
    annotate('text',x=-log10(eqtl$pvalue)[eqtl$snps==snp],
             y=-log10(sqtl$pvalue)[sqtl$snps==snp],
             label=snp,size=3.5,vjust=3,hjust=1)
p7<-plot_grid(ct2_zoom,ct2_comp)


#### Cys2
load('Cys2_eqtl_sqtl.rda')
eqtl<-eqtl[eqtl$pvalue<0.01,]
sqtl<-sqtl[sqtl$pvalue<0.01,]
eqtl$pos<-bed$map$position[match(eqtl$snps,bed$map$snp.name)]/1e6
sqtl$pos<-bed$map$position[match(sqtl$snps,bed$map$snp.name)]/1e6
snp<-eqtl$snps[which.min(eqtl$pvalue)]

cys2_zoom1<-ggplot(eqtl,aes(x=pos,y=-log10(pvalue))) +
    geom_point(fill='gray',color='gray',pch=21,size=1) + 
    ylab(bquote(-log[10](italic(P)))) +
    annotate('text',178.25,20,label='eQTL') +
    theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
    scale_y_continuous(expand = expansion(c(0,0.05)), limits = c(0,25)) + 
    scale_x_continuous(limits = c(178,180.02)) +
    geom_point(inherit.aes = FALSE,x=eqtl$pos[eqtl$snps=='1-177050943'],
               y=-log10(eqtl$pvalue)[eqtl$snps=='1-177050943'],pch=23,fill='purple',size=2)  + 
    annotate('text',x=eqtl$pos[which.max(-log10(eqtl$pvalue))],size=3.5,
             y=max(-log10(eqtl$pvalue)),label='1-177050943',hjust=-0.1)

cys2_zoom2<-ggplot(sqtl,aes(x=pos,y=-log10(pvalue)))  +
    geom_point(fill='gray',color='gray',pch=21,size=1) + 
    xlab('chr1 (Mb)') +
    ylab(bquote(-log[10](italic(P)))) +
    annotate('text',178.25,100,label='sQTL') +
    scale_y_continuous(expand = expansion(c(0,0.05)), limits = c(0,120)) +
    scale_x_continuous(limits = c(178,180)) +
    geom_point(inherit.aes = FALSE,x=sqtl$pos[sqtl$snps=='1-177050943'],
               y=-log10(sqtl$pvalue)[sqtl$snps=='1-177050943'],pch=23,fill='purple',size=2)  + 
    annotate('text',x=sqtl$pos[sqtl$snps=='1-177050943'],size=3.5,
             y=-log10(sqtl$pvalue)[sqtl$snps=='1-177050943'],label='1-177050943',hjust=-0.1)

cys2_zoom<-plot_grid(cys2_zoom1,cys2_zoom2,align = 'v',ncol=1,rel_heights = c(0.8,1)) 
cys2<-locuscompare(eqtl,sqtl,bed=bed, combine = F,region = 1e6,legend_position = "topleft")
cys2_comp<-cys2$locuscompare +
    annotate('text',x=-log10(eqtl$pvalue)[eqtl$snps==snp],
             y=-log10(sqtl$pvalue)[sqtl$snps==snp],
             label=snp,size=3.5,vjust=0.5,hjust=1.1)
p8<-plot_grid(cys2_zoom,cys2_comp)

#### arrange plots by row
load("ct2_cys2_plist.rda")
p5 <- p5+theme(axis.title.x = element_text(margin = margin(0,0,-2.5,0)))
pr1 <- plot_grid(p1,p2,p3,labels = 'auto',nrow=1)
pr2 <- plot_grid(p4, p5, p6,labels = c('d','e','f'),nrow=1)

p7 <- annotate_figure(p7, 
                      top = text_grob("Compact plant2 (ct2)", 
                                      face="bold.italic", color='darkblue'))
pr3 <- plot_grid(p7,
                 plot_grid(plotlist = plist[c(1,3)], nrow=2),
                 labels = c("g", "h"), nrow=1, 
                 rel_widths = c(2,1))

p8 <- annotate_figure(p8, 
                      top = text_grob("Cysteine synthase2 (cys2)", 
                                      face="bold.italic", color='darkblue'))
pr4 <- plot_grid(p8,
                 plot_grid(plotlist = plist[c(2,4)], nrow=2),
                 labels = c("i", "j"), nrow=1, 
                 rel_widths = c(2,1))

#### output
pdf(width = 10,height = 12,file='../Fig2.pdf')
plot_grid(pr1,pr2,pr3,pr4,ncol=1,rel_heights = c(1,1,1.5,1.5))
dev.off()


