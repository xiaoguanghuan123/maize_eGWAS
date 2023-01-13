chr_sizes.maize<-c(307041717,
                   244442276,
                   235667834,
                   246994605,
                   223902240,
                   174033170,
                   182381542,
                   181122637,
                   159769782,
                   150982314)
names(chr_sizes.maize)<-seq(10)
chr_sizes.rice<-c(43270923,
                  35937250,
                  36413819,
                  35502694,
                  29958434,
                  31248787,
                  29697621,
                  28443022,
                  23012720,
                  23207287,
                  29021106,
                  27531856)
names(chr_sizes.rice)<-seq(12)
# col<-RColorBrewer::brewer.pal(5,'Dark2') 
# col<-RColorBrewer::brewer.pal(5,'Paired') 
manh<-function(x,chr=NULL,xlim=NULL,ylim=NULL,cut=1e-5,chr_sizes=chr_sizes.maize,
               col=c('darkblue','yellowgreen'),pch=20,cex=0.8,chr.split=TRUE,
               x.axis=T,ann=NULL,ann.col=1,ann.cex=0.8,ann.font=3,cut.col=2,
               log=TRUE,add=FALSE,ylab=expression(-log[10](italic(P))),...) {
    x<-as.data.frame(x)
    chr_sizes<-chr_sizes+chr_sizes/20
    chrom<-grep('^chr',colnames(x),ignore.case =T)
    pos<-grep('^bp|^pos|^ps',colnames(x),ignore.case =T)
    pval<-grep('^p$|^pval',colnames(x),ignore.case =T)
    if (length(chrom)!=1) {
        stop('there is no unique chromosome column in "x"')
    }
    if (length(pos)!=1) {
        stop('there is no unique position column in "x"')
    }
    if (length(pval)!=1) {
        stop('there is no unique pvalue column in "x"')
    }
    x<-x[order(x[,chrom],x[,pos]),]
    chrom<-x[,chrom]
    pos<-x[,pos]
    pval<-x[,pval]
    if (log) {
        pval<- -log10(pval)
        cut<- -log10(cut)
    }
    sizes<-c(0,cumsum(chr_sizes[-length(chr_sizes)]))
    names(sizes)<-names(chr_sizes)
    pos<-pos+sizes[match(chrom,names(sizes))]
    color<-rep_len(col,length(chr_sizes))
    names(color)<-names(chr_sizes)
    color<-color[match(chrom,names(sizes))]
    if (is.null(chr)) {
        if (add) {
            points(x=pos,y=pval,pch=pch,cex=cex,col=color,...)
        }else{
            plot(x=pos,y=pval,pch=pch,cex=cex,col=color,xlab='Chromosome',
                 xlim=xlim,ylim=ylim,xaxt='n',las=1,ylab=ylab, ...)
            if (x.axis) axis(1,cumsum(chr_sizes)-chr_sizes/2,names(chr_sizes))
            abline(h=cut,col=cut.col,lwd=1,lty=2)
            if (chr.split) {
                abline(v=(cumsum(chr_sizes)-chr_sizes/40)[-length(chr_sizes)],
                       col='gray') 
            }
        }
        if (!missing(ann)) {
            h<-max(ylim)
            x0<-ann[,2]+sizes[match(ann[,1],names(sizes))]
            y0<-ann[,3]
            y1<-rep(seq(h,h*0.4,-h*0.1),length=length(y0))
            y1[y1<y0]<-y1+h/5
            segments(x0,y0+h/30,x0,y1-h/30)
            points(x0,y0,pch=20,col=ann.col)
            text(x0,y1,labels= ann[,4],cex=ann.cex,col=ann.col,font=ann.font,xpd=T)
        }
    }else{
        if (add) {
            points(x=x[chrom==chr,grep('^bp|^pos|^ps',colnames(x),ignore.case =T)]/1e6,
                   y=pval[chrom==chr],pch=pch,cex=cex,col=color[chrom==chr],...)
        }else{
            plot(x=x[chrom==chr,grep('^bp|^pos|^ps',colnames(x),ignore.case =T)]/1e6,
                 y=pval[chrom==chr],pch=pch,cex=cex,col=color[chrom==chr],xlim=xlim,
                 xlab='Chromosome position (Mb)',ylim=ylim,ylab=ylab,las=1,...)
            abline(h=cut,col=cut.col,lwd=1,lty=2)
        }
    }
}