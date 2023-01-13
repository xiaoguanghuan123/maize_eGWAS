
####
suppressPackageStartupMessages({
    library(data.table)
    library(GenomicRanges)
    library(parallel)
    load('gmap.rda')
})

snp<-fread('region_snp.bim')
colnames(snp)<-c('chr','snp','gpos','ppos','a1','a2')
snp<-split(snp,snp$chr)
res<-mclapply(seq(10),function(i) {
    x<-snp[[i]]
    x$gpos<-x$ppos*gmap$rate[i]/100
    ir<-IRanges(seq(1,max(x$ppos),0.9e6),width = 1e6)
    end(ir)[length(ir)]<-max(x$ppos)
    if (width(ir)[length(ir)]<0.1e6) {
        ir<-ir[-length(ir)]
        end(ir)[length(ir)]<-max(x$ppos)
    }
    for (j in seq_along(ir)) {
        slice<-sprintf('chr%d_%s_%s',i,start(ir[j]),end(ir[j]))
        cmd<-paste('plink --bfile region_snp --keep-allele-order --recode 01 transpose --output-missing-genotype 9',
                   sprintf("--keep-fam temperate.txt --chr %d --from-bp %s --to-bp %s --out geno_slices/%s",
                           i,start(ir[j]),end(ir[j]),slice))
        system(cmd)
        system(sprintf("cut -d' ' -f 5- geno_slices/%s.tped > geno_slices/%s.geno1",slice,slice))
        
        cmd<-paste('plink --bfile region_snp --keep-allele-order --recode 01 transpose --output-missing-genotype 9',
                   sprintf("--keep-fam tropical.txt --chr %d --from-bp %s --to-bp %s --out geno_slices/%s",
                           i,start(ir[j]),end(ir[j]),slice))
        system(cmd)
        system(sprintf("cut -d' ' -f 5- geno_slices/%s.tped > geno_slices/%s.geno2",slice,slice))
        
        write.table(x[x$ppos %within% ir[j],c(2:1,3:6)],sprintf('geno_slices/%s.snp',slice),
                    row.names = F,col.names = F,quote = F,sep='\t')
        
        cmd<-paste(sprintf("XPCLR -xpclr geno_slices/%s.geno1 geno_slices/%s.geno2 geno_slices/%s.snp out/%s",
                           slice, slice, slice, slice),
                   sprintf("-w1 0.0005 100 100 %d -p1 0.95",i))
        system(cmd)
    }
},mc.cores=10)




