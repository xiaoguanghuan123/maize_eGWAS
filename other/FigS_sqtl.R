####
ggpie<-function(value,start=pi/2,direction=1.1,rawNumber=TRUE,label.pos=1,
                show.legend=T,legend.position='right',legend.title=NA,
                legend.text=NA) {
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
    p<- ggplot(data,aes(x=factor(1),y=value,fill=group)) + 
        geom_bar(stat="identity",show.legend = show.legend,color='white') +
        geom_text(label=labels,aes(x=label.pos),
                  position = position_stack(vjust = 0.5)) +
        theme_void() + 
        coord_polar("y", start=start,direction=direction) +
        theme(legend.position=legend.position,
              legend.text = element_text(size=12),
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
load('../../sqtl_peak.rda')
suppressPackageStartupMessages({
    library(ggplot2)
    library(ggpubr)
    library(data.table)  
})


freq<-table(sqtls$intron)
sqtls<-sqtls[sqtls$intron%in%names(freq[freq<=10])]

freq<-table(freq)
freq[11]<-sum(freq[-seq(10)])
freq<-freq[seq(11)]
names(freq)[11]<-'>10'

p1<-ggplot(data.frame(freq),aes(x=Var1,y=Freq))+
    geom_bar(stat="identity",width = 0.6,fill=c(rep("steelblue",10),'gray')) + 
    xlab('Number of sQTLs') +
    ylab(expression("Number of introns")) +
    theme_minimal(12)  + 
    theme(axis.line = element_line(),panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank(), 
          axis.text = element_text(color=1, size=11),
          axis.ticks = element_line(color=1)) +
    geom_text(aes(label=format(Freq,big.mark = ',')), vjust=-0.5,size=2.8) +
    scale_y_continuous(n.breaks =8,limits = c(NA,max(freq)*1.05))

#### 
p2<-ggpie(table(sqtls$type),rawNumber = T,legend.position = c(0.5,0.05),
          legend.title = NULL) + 
    ggtitle("sQTL types") +
    theme(panel.border = element_blank(),
          legend.text = element_text(face = 3),
          legend.direction = "horizontal",
          plot.title = element_text(hjust=0.5, vjust=0))

####
qtl_type<-table(sqtls$intron,sqtls$type)
qtl1<- sum(qtl_type[,1]>0 & qtl_type[,2]==0)
qtl2<- sum(qtl_type[,1]==0 & qtl_type[,2]>0)
qtl3<- sum(qtl_type[,1]>0 & qtl_type[,2]>0)
p3<-ggpie(c(qtl1,qtl2,qtl3),legend.position = c(0.5,.05),
          legend.title = NULL, legend.text = c('cis','trans','cis+trans')) +
    ggtitle("Intron (splicing event) types") +
    theme(panel.border = element_blank(),
          legend.direction = "horizontal",
          plot.title = element_text(hjust=0.5, vjust = 0),
          legend.text = element_text(face = 3))

####
p4 <- ggdensity(sqtls, x='PVE', color='type', fill='type',
                xlab = 'sQTL PVE', ylab='Density') +
    theme(legend.text = element_text(face =3), legend.title = element_blank(),
          legend.position = c(0.5, 0.9), legend.direction = 'horizontal')

df<-data.frame(sqtls[,c('type','PVE','R2')])
lab<-paste0('italic(',c('cis','trans'),')*','"-sQTLs"')
p4.1<-ggplot(df,aes(x=type,y=PVE))+
    geom_violin(aes(fill=type),show.legend = FALSE,draw_quantiles=0.5) +
    theme_light() + xlab(NULL) + ylab('PVE') +
    scale_x_discrete(labels = parse(text = lab)) +
    theme(axis.text.x = element_text(angle = 30,hjust = 1))

p4.2<- ggplot(df,aes(x=PVE*100,y=after_stat(density)))+
    geom_density(aes(fill=type,color=type),size=1,alpha=0.3) +
    theme_classic() + xlab('PVE (%)') + labs(fill='sQTL type:',col='sQTL type:') +
    scale_x_continuous(expand = expansion(0,0.05),limits = c(0,105)) +
    scale_y_continuous(name='Density',expand = expansion(0,0)) +
    theme(legend.position = c(0.7,0.75),
          legend.text = element_text(face=3))

p5<- cowplot::plot_grid(p1, p4.2, align = "h", nrow = 1,labels = c(NA,'b'),
                        rel_widths = c(0.7,0.3), axis = 'bt')
p6<- cowplot::plot_grid(p2, p3, align = "h", nrow = 1,labels = c(NA,'d'))

p7<- cowplot::plot_grid(p5, p6, nrow=2,labels =c('a','c'),scale=0.98)

####
pdf(width = 10,height = 7,file='../FigS_sqtl.pdf')
plot_grid(p1, p2, p4, p3, labels = 'auto')
dev.off()
