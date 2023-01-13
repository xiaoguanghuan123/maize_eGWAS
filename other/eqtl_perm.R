
#######
suppressPackageStartupMessages({
    library(parallel)
    library(MatrixEQTL)
    library(data.table)
    library(GenomicRanges)  
})

load('SliceData_fpkm.rda')
load('../SlicedData_snps.rda')

ls()
snps
m <- nrow(fpkm)
n <- ncol(fpkm)

set.seed(2023)
snps$ColumnSubsample(sample(n))
snps

idx <- split(seq_len(m), cut(seq_len(m), 36, labels=F))

res <- pbapply::pblapply(idx, function(i) {
    fpkm$RowReorder(i)
    suppressMessages(
        me <- Matrix_eQTL_engine(
        snps = snps, 
        gene = fpkm, 
        cvrt = cvrt, 
        output_file_name = NULL, 
        pvOutputThreshold = 1e-5, 
        useModel = modelLINEAR, 
        errorCovariance = numeric(), 
        verbose = FALSE,
        pvalue.hist = FALSE,
        min.pv.by.genesnp = FALSE,
        noFDRsaveMemory = FALSE
    ))
    me$all$eqtls
}, cl = 36)

res <- do.call(rbind, res)

save(res,file='perm_res.rda')

####
load('SliceData_fpkm.rda')
load('perm_res.rda')
load('eqtl_peak.rda')

f <- factor(res$gene, levels = unlist(fpkm$rowNameSlices))
pval <- tapply(res$pvalue, f, min) 
pval[is.na(pval)] <- 1e-5

quantile(pval ,c(0.01,0.02,0.05))

snps <- data.table::fread("../geno340/snps.bim")
perm0.05 <- -log10(quantile(pval, 0.05))
cutoff <- -log10(0.05/nrow(snps))

pdf(width = 7, height = 4, file='permutation.pdf')
par(mar=c(4, 4, 0.5, 0.5))
plot(density(-log10(pval),from=5), xlim=c(5, 15), col=2, lwd=2,
     xlab=bquote(-log[10](italic(P))), main=NA, las=1)
abline(v=perm0.05, col=2, lty=2, lwd=2)
lines(density(-log10(eqtls$pvalue), from=cutoff), col=4, lwd=2)
abline(v=cutoff, col=4, lty=2, lwd=2)
legend(12, 0.6, legend = c("Permutation", "eQTLs"), col=c(2,4), lwd=2,bty='n')
mtext(round(perm0.05,2), 1, col=2, at=perm0.05, line = 0.5, adj = 1, font=2)
mtext(round(cutoff,2), 1, col=4, at=cutoff, line = 0.5, adj = 0, font=2)
dev.off()

