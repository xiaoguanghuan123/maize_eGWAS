#######
setwd('/work3/maize_eGWAS/figures/Fig3_code/')
suppressPackageStartupMessages({
    library(ggplotify)
    library(ggplot2)
    library(cowplot)
    library(GenomicRanges)
    library(data.table)
})
source('manh.R')

theme_set(theme_bw(12) + 
              theme( axis.ticks = element_line(colour = 1),
                     axis.text = element_text(colour = 1,size=rel(0.9)),
                     legend.title = element_text(size=rel(0.9)),
                     legend.text = element_text(size=rel(0.9)),
                     plot.title = element_text(hjust = 0.5,size = 12)))

#### IC_scatter
mod<-read.csv('moduleHb_Size.csv')
test<-cor.test(log2(mod$ModelSize),mod$Heritability)
lab <-  mod[mod$ModelName %in% c(39, 79, 101), ]
lab$Func <- c("(Oxidation reduction;\nResponse to water deprivation)",
              "(Translation)",
              "(Protein biogenesis/protein folding;\n Plant-pathogen interaction)")
p1<-ggplot(mod,aes(x=log2(ModelSize),y=Heritability)) + geom_point() +
    geom_smooth(method = 'lm') +
    labs(x=bquote(log[2]*'(Model size)'),
         y=bquote(Model~mean~H^2), 
         title = substitute(italic(R)==x*','~italic(P)==y,
                            list(x=round(test$estimate,2),y=signif(test$p.value,3)))) +
    annotate('text',x=log2(lab$ModelSize),y=lab$Heritability,
             label=paste0('IC',lab$ModelName, " ", lab$Func),
             vjust=c(0.2,-1, 0.2), hjust=c(-0.05,0.2,-0.05), color=2, size=3) +
    annotate('point',x=log2(lab$ModelSize),y=lab$Heritability,
             shape=23, fill="blue", color="darkred", size=3)
pdf(width = 6, height = 5, file='IC_scatter.pdf')
print(p1)
dev.off()

### IC39
load('../../sweep/sel_genes_042121.rda')
load('../../maize_genes.rda')
load('../../eGWAS340/eqtl_peak.rda')

f1 <- fread("../../gemma/model_res/mod39.assoc.txt")
cut <- 0.05/nrow(f1)
df1$pval <- df1$pval * -log10(cut)/cut_xpclr
df2$pval <- df2$pval * -log10(cut)/cut_fst
x1 <- sort(genes[c('Zm00001d047757','Zm00001d047755','Zm00001d047758')],
           ignore.strand=T)

colnames(f1)[colnames(f1) == 'p_lrt'] <- 'pval'
f1 <- subset(f1, pval < 0.001 )

p2 <- as_grob(function() {
    par(layout(matrix(rep(1:2, c(3,4)))), mar=c(2, 4.5, 0.5, 0.5), cex=1)
    manh(f1, col=ggsci::pal_igv()(10), cut=cut,cut.col = 1, cex=0.5)
    par(mar=c(4.5, 4.5, 2.5, 0.5))
    manh(f1, col=ggsci::pal_igv()(10), cut=cut,cut.col = 1, chr=9, 
         xlim=c(140.8,141.2), cex=0.5)
    points(x=141083854/1e6, y=-log10(f1$pval[f1$ps==141083854]),pch=23, 
           col=1, bg=2, cex=1)
    #text(x=141083854/1e6, y=-log10(f1$pval[f1$ps==141083854]), label=141083854,
    #     adj = c(-0.1, 0.5), col=2)
    rect(start(x1)/1e6, 43, end(x1)/1e6, 44, col=2, xpd=T)
    text(x=start(x1[1])/1e6, y=44, xpd=T, labels=x1$symbol[1], font=3,
         adj=c(0.5, -0.5))
    text(x=start(x1[2])/1e6, y=44, xpd=T, labels=x1$symbol[2], font=3,
         adj=c(1, -0.5))
    text(x=end(x1[3])/1e6, y=44, xpd=T, labels=x1$symbol[3], font=3,
         adj=c(0, -0.5))
    #manh(df1, col=2, chr=9, type='l', log=F, lwd=2, xlim=c(140.6,141.4), add=T)
    #abline(v=141083854/1e6)
    #manh(df2, col=4, chr=9, type='l', log=F, lwd=2, xlim=c(140.6,141.4), add=T)
})
pdf (width = 6, height = 5, file='IC39.pdf')
plot_grid(p2)
dev.off()

#### IC79
f2 <- fread("../../gemma/model_res/mod79.assoc.txt")
colnames(f2)[colnames(f2) == 'p_lrt'] <- 'pval'
f2 <- subset(f2, pval < 0.001 )
x2 <- genes['Zm00001d014664']

p3 <- as_grob(function() {
    par(layout(matrix(rep(1:2, c(3,4)))), mar=c(2, 4.5, 0.5, 0.5), cex=1)
    manh(f2, col=ggsci::pal_igv()(10), cut=cut,cut.col = 1, cex=0.5)
    par(mar=c(4.5, 4.5, 2.5, 0.5))
    manh(f2, col=ggsci::pal_igv()(10), cut=cut,cut.col = 1, chr=5, xlim=c(57.5,58.5), cex=0.5)
    lead <- f2[which.min(f2$pval)]
    points(x=lead$ps/1e6, y=-log10(lead$pval), pch=23, cex=1, bg=2)
    rect(start(x2)/1e6, 13.5, end(x2)/1e6, 14, col=2, xpd=T)
    text(x=start(x2)/1e6, y=14, labels = 'cia2', font=3, adj=c(0.5, -0.5),xpd=T)
    #manh(df1, col=2, chr=1, type='l', log=F, lwd=2, xlim=c(57,59), add=T)
    #manh(df2, col=4, chr=1, type='l', log=F, lwd=2, xlim=c(57,59), add=T)
})
pdf (width = 6, height = 5, file='IC79.pdf')
plot_grid(p3)
dev.off()

#### IC101
f3 <- fread("../../gemma/model_res/mod101.assoc.txt")
colnames(f3)[colnames(f3) == 'p_lrt'] <- 'pval'
f3 <- subset(f3, pval < 0.001 )
x3 <- genes['Zm00001d042533']

p4 <- as_grob(function() {
    par(layout(matrix(rep(1:2, c(3,4)))), mar=c(2, 4.5, 0.5, 0.5), cex=1)
    manh(f3, col=ggsci::pal_igv()(10), cut=cut,cut.col = 1)
    par(mar=c(4.5, 4.5, 2.5, 0.5))
    manh(f3, col=ggsci::pal_igv()(10), cut=cut,cut.col = 1, chr=3, xlim=c(170.6,171.3))
    lead <- f3[which.min(f3$pval)]
    points(x=lead$ps/1e6, y=-log10(lead$pval), pch=23, cex=1, bg=2)
    rect(start(x3)/1e6, 16, end(x3)/1e6, 16.5, col=2, xpd=T)
    text(x=start(x3)/1e6, y=16.5, labels = x3$symbol, font=3, adj=c(0.5, -0.5),xpd=T)
    #manh(df1, col=2, log=F, type='l', cut=cut,cut.col = 1, chr=3, xlim=c(170.5,171.5),add=T, lwd=2)
    #manh(df2, col=4, log=F, type='l', cut=cut,cut.col = 1, chr=3, xlim=c(170.5,171.5),add=T, lwd=2)
    yaxis2 <- c(0, 15) * cut_xpclr/(-log10(cut))
    #axis(4, at=pretty(yaxis2) * -log10(cut)/cut_xpclr, labels = pretty(yaxis2), las=1)
    #mtext('XPCLR', 4, line=2, cex=1)
})

pdf (width = 6, height = 5, file='IC101.pdf')
plot_grid(p4)
dev.off()


#### combine
pdf(width = 12,height = 10,file='../Fig3.pdf')
plot_grid(p1,p2,p3,p4, labels = 'auto')
dev.off()
