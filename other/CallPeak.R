#### call peak from output of MatrixEQTL
CallPeak<-function(qtl,snps, wind=1e6, cutoff=0.05/nrow(snps),snpN=1,cores=1){
    require(parallel)
    require(data.table)
    qtl<-subset(qtl,pvalue<=cutoff)
    if (nrow(qtl)==0) return(NULL)
    res<-merge(data.table(qtl),data.table(snps)[,c('snps','chr','pos')],by='snps')
    res<-split(res,res$gene)
    res<-mclapply(res,function(x) {
        chr<-split(x,x$chr)
        chr<-lapply(chr, function(y) {
            y<-y[order(pos)]
            idx<-which(diff(y$pos)>wind)
            idx.start<-c(1,idx+1)
            idx.end<-c(idx,nrow(y))
            n<-idx.end-idx.start+1
            if (all(n<snpN)) return(NULL)
            rbindlist(lapply(seq_along(idx.start)[n>=snpN], function(i) {
                z<-y[seq(idx.start[i],idx.end[i])]
                start<-min(z$pos); end<-max(z$pos)
                distance<-end-start
                snp_number<-nrow(z)
                z<-z[which.min(pvalue)]
                data.table(z,snp_number,start,end,distance)
            }))
        })
        chr<-rbindlist(chr)
        if (nrow(chr)==0) return(NULL)
        chr<-chr[order(pvalue)]
        chr
    },mc.cores = cores)
    rbindlist(res)
}