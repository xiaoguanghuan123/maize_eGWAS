
#######
setwd('/work3/maize_eGWAS/figures/Fig1_code/')
suppressPackageStartupMessages({
    library(ggpubr)
    library(cowplot)
})

npg<-ggsci::pal_npg('nrc')(10)
hb<-read.csv('e340_Heritability.csv',row.names = 1)

#### Hb histgram
p1<-ggplot(hb,aes(Hb)) + 
    stat_bin(breaks=seq(0,1,0.05),fill='gray',color=1) + 
    theme_classic(12) + 
    theme(axis.text = element_text(colour = 'black'),
          axis.ticks = element_line(colour = 'black')) +
    xlab(bquote(H^2)) + ylab('Number of genes') +
    scale_x_continuous(expand = expansion(c(0,0.05)),limits = c(0,1)) +
    scale_y_continuous(expand = expansion(c(0,0.05))) 

#### Hb 2D density
p2<-ggplot(hb,aes(log2GEMean,Hb)) + 
    stat_density_2d(geom='polygon',aes(fill = after_stat(level)),
                    breaks=seq(0,1,0.05)) +
    scale_fill_distiller(palette=12, direction=1) +
    xlab(bquote(log[2](FPKM))) + ylab(bquote(H^2)) +
    theme_classic(12) +
    theme(legend.position = c(0.9,0.8),legend.key.size = unit(12,'point'),
          axis.text = element_text(colour = 'black'),
          axis.ticks = element_line(colour = 'black')) +
    scale_x_continuous(expand = expansion(c(0,0.05)),limits = c(0,NA)) +
    scale_y_continuous(expand = expansion(c(0,0.05)),limits = c(0,1)) 


#### Volcano plot
library(EnhancedVolcano)
vol <- read.csv('Hb_GO_redundancyReduced.csv',row.names=1)
head(vol)
max(vol[which(vol$FDR<=0.05),"pValue"])

BP = vol$desc[vol$Category =='Biological process']
CC = vol$desc[vol$Category =='Cellular component']
MF = vol$desc[vol$Category =='Molecular function']
celltype1 <- BP
celltype2 <- MF
keyvals.shape <- ifelse(vol$desc %in% celltype1, 17,
                        ifelse(vol$desc %in% celltype2, 19, 18))

keyvals.shape[is.na(keyvals.shape)] <-18
names(keyvals.shape)[keyvals.shape == 18] <- 'Cellular component'
names(keyvals.shape)[keyvals.shape == 17] <- 'Biological process'
names(keyvals.shape)[keyvals.shape == 19] <- 'Molecular function'
volSig = subset(vol, abs(log2FoldChange)>log2(1.2) & FDR <=0.05) # & Category != 'Cellular component')
rownames(volSig)

p3 <- EnhancedVolcano(vol,lab = vol$desc,x = 'log2FoldChange',y = 'pValue',
                     shapeCustom = keyvals.shape,
                     pCutoff = 4.94e-05, FCcutoff = log2(1.1), colAlpha = .7,
                     ylim = c(0, 12), xlim = c(-0.8,0.8),
                     pointSize = 2,
                     labSize = 2,
                     caption = 'total = 981 GO terms',
                     gridlines.major = T,gridlines.minor = F,
                     title = NULL, subtitle = NULL,
                     xlab = bquote(log[2]("fold change")),
                     ylab = bquote(-log[10](italic(P))),
                     borderWidth = 1.5,
                     legendLabSize=9,
                     legendIconSize = 2.5,
                     drawConnectors = F) +
    theme_classic(12)   +
    scale_y_continuous(expand = expansion(c(0,0.05)),limits = c(0,12)) +
    guides(color = FALSE) + xlim(-1,0.6) +
    theme(legend.position = 'top', legend.title = element_blank(),
          legend.text = element_text(size=9),
          axis.text = element_text(colour = 'black'),
          axis.ticks = element_line(colour = 'black'))

#### Hb_Trop_VS_Temperate
df <- read.csv('Hb_Trop_VS_Temperate.csv',row.names = 1)
df <- subset(df, HbTemp >0 & HbTrop > 0)

p4<-ggplot(data = df, aes(HbTemp, HbTrop)) +
    stat_density_2d(aes(fill = ..level..), geom='polygon') +
    # geom_smooth(position = 'identity', method = lm) +
    scale_fill_distiller(palette=12, direction=1) +
    theme(legend.position='right')+
    xlim(0,1) + ylim(0,1) +
    xlab(bquote(H^2~"Tropical")) +
    ylab(bquote(H^2~"Temperate")) +
    theme_classic(12)
# coef <- round(coef(lm(HbTemp~HbTrop, df)), 2)
# p4 <- p4 + annotate("text", x=0.05,  y= 0.95, hjust=0, col=2, size=4, fontface="italic",
#                    label=substitute(y==a*x+b, list(a=coef[2], b=coef[1]))) +
#     theme(legend.position = c(0.9,0.3),legend.key.size = unit(12,'point'),
#           axis.text = element_text(colour = 'black'),
#           axis.ticks = element_line(colour = 'black')) +
#     scale_x_continuous(expand = expansion(c(0,0.05)),limits = c(0,1)) +
#     scale_y_continuous(expand = expansion(c(0,0.05)),limits = c(0,1)) 
p4

p5<-ggplot(data.frame(x=c(df$HbTemp, df$HbTrop), 
                      y=rep(c("Temperate","Tropical"),each=nrow(df))),
           aes(x, fill=y, color=y)) +
    geom_density(alpha=0.3,lwd=1,show.legend = T) +
    theme_classic(12)  +
    labs(x=bquote(H^2), y='Density') +
    scale_x_continuous(expand = expansion(c(0,0.05)),limits = c(0,1)) +
    scale_y_continuous(expand = expansion(c(0,0.05)),limits = c(0,NA)) +
    theme(legend.position = c(0.8,0.9),legend.title = element_blank(),
          legend.key.size = unit(12,'point'), legend.text = element_text(size=8),
          axis.text = element_text(colour = 'black'),
          axis.ticks = element_line(colour = 'black')) 

####
pc1 <- plot_grid(p1,p2,align = 'v',ncol=1,labels = letters[1:2])
pc2 <- plot_grid(p3,labels=letters[3])
pc3 <- plot_grid(p4,p5,align = 'v',ncol=1,labels = letters[4:5])

pdf(width = 10,height = 6, file='../Fig1.pdf')
plot_grid(pc1,pc2,pc3,nrow=1,rel_widths = c(1,1.5))
dev.off()

    
    