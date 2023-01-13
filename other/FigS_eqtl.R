#### function ggpie
ggpie<-function(value,start=pi/2,direction=1.1,rawNumber=TRUE,label.pos=1,
                show.legend=T,legend.position='right',legend.title=NA,
                legend.text=NA,label.size=4) {
    require(ggplot2)
    if (is.null(names(value))) {
        group<-paste0('group',seq_along(value))
    }else {
        group<-names(value)
    }
    value<-as.numeric(value)
    data<-data.frame(value=value,group=group)
    labels<- paste0(round(data$value/sum(data$value) * 100, 2), "%")
    if (rawNumber) {
        labels<- paste0(format(data$value,big.mark = ','),'\n(',labels,')')
    }
    p<- ggplot(data,aes(x=factor(1),y=value,fill=group),size=label.size*11/4) + 
        geom_bar(stat="identity",show.legend = show.legend,color='white') +
        geom_text(label=labels,aes(x=label.pos),size=label.size,
                  position = position_stack(vjust = 0.5)) +
        theme_void() + 
        coord_polar("y", start=start,direction=direction) +
        theme(legend.position=legend.position,
              legend.text = element_text(size=11*label.size/4),
              legend.title = element_text(size=11*label.size/4),
              legend.key.size = unit(label.size*1.2/4,'lines'),
              panel.border = element_rect(fill=NA),
              plot.margin = unit(rep(0.2,4), "lines")) 
    if (!missing(legend.text)) {
        p<- p + scale_fill_discrete(name='group',labels=legend.text)
    }
    if (!missing(legend.title)) {
        p<- p + guides(fill=guide_legend(title=legend.title)) 
    }
    p
}

####
suppressPackageStartupMessages({
    library(ggplot2)
    library(ggpubr)
    library(cowplot)
    library(data.table)
    library(GenomicRanges)
})
load('../../eGWAS340/eqtl_peak.rda')

#### frequency
freq<-table(eqtls$gene)
eqtls<-eqtls[eqtls$gene%in%names(freq[freq<=10])]
eqtl1<-eqtls[eqtls$gene%in%names(freq[freq==1])]

freq<-table(freq)
freq[11]<-sum(freq[-seq(10)])
freq<-freq[seq(11)]
names(freq)[11]<-'>10'

p1<-ggplot(data.frame(freq),aes(x=Var1,y=Freq))+
    geom_bar(stat="identity",width = 0.6,fill=c(rep("steelblue",10),'gray')) + 
    xlab('Number of eQTLs') +
    ylab(expression("Number of genes")) +
    theme_minimal(12)  + 
    theme(axis.line = element_line(),panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank(), 
          axis.text = element_text(colour = 'black',size=12),
          axis.ticks = element_line(colour = 'black')) +
    geom_text(aes(label=format(Freq,big.mark = ',')), vjust=-0.5,size=3.5) +
    scale_y_continuous(n.breaks =8,limits = c(NA,max(freq)*1.05)) 
    

p1.1<-ggpie(table(eqtl1$type),rawNumber = T,
            legend.position = 'top',
            legend.title = 'Single eQTL gene:',
            label.pos = 1.1,label.size = 3,show.legend = T) +
    theme(legend.text = element_text(face = 3))

#### scatter plot
chr<-read.table('chrLength.txt',header = T)
chr$cumsum<-cumsum(c(0,chr$length[-10]))
load('../../maize_genes.rda')
gene_pos<-genes[eqtls$gene]
eqtl_pos<-GRanges(eqtls$chr,eqtls$pos)
p<- -log10(eqtls$pvalue)
p[p>101]<-101

gene_pos<-start(gene_pos)+chr$cumsum[match(as.character(seqnames(gene_pos)),chr$chr)]
eqtl_pos<-start(eqtl_pos)+chr$cumsum[match(as.character(seqnames(eqtl_pos)),chr$chr)]
df<-data.frame(x=gene_pos,y=eqtl_pos,p=p)

p2<-ggplot(df)+geom_point(aes(x=x,y=y),size=0.2,color='darkblue') + theme_classic() +
    labs(x='Gene position (chromosome)',y='eQTL position (chromosome)',
         color=bquote("-"*log[10]*(italic(P)))) +
    scale_x_continuous(breaks = chr$length/2+chr$cumsum,label=1:10) +
    scale_y_continuous(breaks = chr$length/2+chr$cumsum,label=1:10)

#### pie plot
p3<-ggpie(table(eqtls$type),rawNumber = T,legend.position = c(0.5,0),
          legend.title = NULL,label.pos = 1) + 
    ggtitle("eQTL types") +
    theme(panel.border = element_blank(),
          legend.text = element_text(face = 3),
          legend.direction = 'horizontal', 
          plot.title = element_text(hjust = 0.5))

qtl_type<-table(eqtls$gene,eqtls$type)
qtl1<- sum(qtl_type[,1]>0 & qtl_type[,2]==0)
qtl2<- sum(qtl_type[,1]==0 & qtl_type[,2]>0)
qtl3<- sum(qtl_type[,1]>0 & qtl_type[,2]>0)
p4<-ggpie(c(qtl1,qtl2,qtl3),legend.position = c(0.5,0),
          legend.title = NULL, label.pos = 1.1,
          legend.text = c('cis','trans','cis+trans')) +
    ggtitle("Gene (e-trait) types") +
    theme(panel.border = element_blank(),
          legend.text = element_text(face = 3),
          legend.direction = 'horizontal',
          plot.title = element_text(hjust = 0.5))

#### gene and eQTL distribution
p5<-paste0('chr',1:10)
for(i in 1:10) {
    x<-c(start(genes[seqnames(genes)==i]),eqtls$pos[eqtls$chr==i])
    type<-rep(c('Genes','eQTLs'),c(sum(seqnames(genes)==i),sum(eqtls$chr==i)))
    assign(p5[i],ggplot(data.frame(x=x,type=type),aes(x=x,y=after_stat(density))) +
        geom_density(aes(fill=type,color=type),size=1,alpha=0.3) +
        theme_cowplot() + 
        xlab(paste0('chr',i)) +
        theme(axis.title.y = element_blank(),
              axis.text = element_blank(),
              legend.position = 'none')
    )
}


p6<-plot_grid(get(p5[1]),
          get(p5[2]),get(p5[3]),get(p5[4]),get(p5[5]),
          get(p5[6]),get(p5[7]),get(p5[8]),get(p5[9]),get(p5[10]),
          nrow=2,ncol=5)
####
p7 <- ggdensity(eqtls, x='MAF', color = "type", fill='type',alpha=0.5,
          xlab = "Minor allele frequency", ylab='Density') +
    theme(legend.text = element_text(face =3), legend.title = element_blank(),
          legend.position = c(0.5, 0.9), legend.direction = 'horizontal')
p8 <- ggdensity(eqtls, x='PVE', color = "type", fill='type',alpha=0.5,
                xlab = "eQTL PVE", ylab='Density') +
    theme(legend.text = element_text(face =3), legend.title = element_blank(),
          legend.position = c(0.5, 0.9), legend.direction = 'horizontal')

### permutation
load('../../eGWAS340/SliceData_fpkm.rda')
load('../../eGWAS340/perm_res.rda')

f <- factor(res$gene, levels = unlist(fpkm$rowNameSlices))
pval <- tapply(res$pvalue, f, min) 
pval[is.na(pval)] <- 1e-5

quantile(pval ,c(0.01,0.02,0.05))

snps <- data.table::fread("../../geno340/snps.bim")
perm0.05 <- -log10(quantile(pval, 0.05))
cutoff <- -log10(0.05/nrow(snps))

p9 <- as_grob(function () {
    par(mar=c(4, 4, 0.5, 0.5))
    plot(density(-log10(pval),from=5), xlim=c(5, 15), col=2, lwd=2,
         xlab=bquote(-log[10](italic(P))), main=NA, las=1)
    abline(v=perm0.05, col=2, lty=2, lwd=2)
    lines(density(-log10(eqtls$pvalue), from=cutoff), col=4, lwd=2)
    abline(v=cutoff, col=4, lty=2, lwd=2)
    legend(10, 0.6, legend = c("Permutation", "eQTLs"), col=c(2,4), lwd=2,bty='n')
    mtext(round(perm0.05,2), 1, col=2, at=perm0.05, line = 0.5, adj = 1, font=2)
    mtext(round(cutoff,2), 1, col=4, at=cutoff, line = 0.5, adj = 0, font=2)
})

####
pdf(width = 12,height = 8,file='../FigS_eqtl.pdf')
pr1 <- plot_grid(ggdraw(p1)+draw_grob(as_grob(p1.1),scale = 0.7,x=0.2,y=0.15),
                 p7, p8, rel_widths = c(1, 0.5, 0.5), labels = 'auto', 
                 scale = 0.98, nrow=1)
pr2 <- plot_grid(p3, p4, p9, rel_widths = c(0.5, 0.5, 1), nrow=1,
                 labels = c('d','e','f'), scale = 0.98)
plot_grid(pr1, pr2, ncol=1)
dev.off()

