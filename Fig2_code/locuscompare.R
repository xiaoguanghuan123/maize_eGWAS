#### rewrite locuscompare
locuscompare<-function(in_fn1, in_fn2, bed, snp=NULL, region=500e3,  
                       marker_col1 = "snps", pval_col1 = "pvalue", title1 = "eQTL: ", 
                       marker_col2 = "snps", pval_col2 = "pvalue", title2 = "sQTL: ",  
                       combine = TRUE, legend = TRUE, 
                       legend_position = c("bottomright", "topright", "topleft"), 
                       lz_ylab_linebreak = FALSE ){
    
    #### loading libraries 
    library(snpStats)
    library(ggplot2)
    library(locuscomparer)
    options('ggrepel.max.overlaps'=Inf)
    theme_set(theme_classic(12) + 
                  theme( axis.ticks = element_line(colour = 1),
                         axis.text = element_text(colour = 1),
                         legend.text = element_text(size=10),
                         legend.background = element_blank()))
    
    #### add_label
    add_label = function(merged, snp){
        merged$label = ifelse(merged$rsid %in% snp, merged$rsid, '')
        return(merged)
    }
    
    #### read_metal
    read_metal<-function (in_fn, marker_col = "rsid", pval_col = "pval") {
        if (is.character(in_fn)) {
            d = read.table(in_fn, header = TRUE)
        }
        else if (is.data.frame(in_fn)) {
            d = in_fn
        }
        else {
            stop("The argument \"in_fn\" must be a string or a data.frame")
        }
        colnames(d)[which(colnames(d) == marker_col)] = "rsid"
        colnames(d)[which(colnames(d) == pval_col)] = "pval"
        d$logp = -log10(d$pval)
        return(d[, c("rsid", "pval", "logp")])
    }
    
    #### rewrite make_scatterplot 
    make_scatterplot<-function (merged, title1, title2, color, shape, size, legend = TRUE, 
                                legend_position = c("bottomright", "topright", "topleft") ) {
        require(grid)
        p = ggplot(merged, aes(logp1, logp2)) + 
            geom_point(aes(fill = rsid, size = rsid, shape = rsid), alpha = 0.8) + 
            geom_point(data = merged[merged$label !=  "", ], 
                       aes(logp1, logp2, fill = rsid, size = rsid, shape = rsid)) + 
            xlab(bquote(.(title1)*-log[10](italic(P)))) + 
            ylab(bquote(.(title2)*-log[10](italic(P)))) + 
            scale_fill_manual(values = color, guide = "none") + 
            scale_shape_manual(values = shape, guide = "none") + 
            scale_size_manual(values = size, guide = "none") 
        #ggrepel::geom_text_repel(aes(label = label),size=3.5)
        if (legend == TRUE) {
            legend_position = match.arg(legend_position)
            if (legend_position == "bottomright") {
                legend_box = data.frame(x = 0.8, y = seq(0.4, 0.2, -0.05))
            }else if (legend_position == "topright") {
                legend_box = data.frame(x = 0.8, y = seq(0.8, 0.6, -0.05))
            }else {
                legend_box = data.frame(x = 0.1, y = seq(0.8, 0.6, -0.05))
            }
            p = p + 
                annotation_custom(rectGrob(legend_box$x,legend_box$y,0.05,0.05,
                                           just = c(0,0),
                                           gp=gpar(fill = rev(c("blue4", "skyblue", "darkgreen", "orange", "red"))))) +
                annotation_custom(textGrob(seq(0.8,0.2,-0.2), x = legend_box$x[-1] + 0.05, 
                                           y = legend_box$y[-1]+0.05,just = c(-0.3,0.5),gp=gpar(fontsize=10,lineheight=0.9))) + 
                annotation_custom(textGrob(parse(text = "italic(r^2)"), x = legend_box$x[1] + 0.05, 
                                           y = legend_box$y[1]+0.05,just = c(0.5,-1),gp=gpar(fontsize=10,lineheight=0.9)))
        }
        return(p + theme(plot.margin = theme_get()$plot.margin-margin(0,0,4.5,0)))
    }
    
    #### make_locuszoom
    make_locuszoom<-function (metal, title, chr, color, shape, size, ylab_linebreak = FALSE) {
        p = ggplot(metal, aes(x = pos, logp)) + 
            geom_point(aes(fill = rsid, size = rsid, shape = rsid), alpha = 0.8) + 
            geom_point(data = metal[metal$label != "", ], aes(x = pos, logp, fill = rsid, size = rsid, shape = rsid)) + 
            scale_fill_manual(values = color, guide = "none") + 
            scale_shape_manual(values = shape, guide = "none") + 
            scale_size_manual(values = size, guide = "none") + 
            ggrepel::geom_text_repel(aes(label = label),size=3.5) + 
            xlab(paste0("chr", chr, " (Mb)")) + 
            ylab(bquote(.(title)*-log[10](italic(P))))
        if (ylab_linebreak == TRUE) {
            p = p + ylab(bquote(atop(.(title), -log[10](italic(P)))))
        }
        return(p)
    }
    
    #### make_combined_plot
    make_combined_plot <-function (merged, title1, title2, ld, chr, snp = NULL, 
                                   combine = TRUE, legend = TRUE, 
                                   legend_position = c("bottomright", "topright",  "topleft"), 
                                   lz_ylab_linebreak = FALSE) {
        snp = get_lead_snp(merged, snp)
        color = assign_color(merged$rsid, snp, ld)
        shape = ifelse(merged$rsid == snp, 23, 21)
        names(shape) = merged$rsid
        size = ifelse(merged$rsid == snp, 3, 2)
        names(size) = merged$rsid
        merged = add_label(merged, snp)
        p1 = make_scatterplot(merged, title1, title2, color, shape, 
                              size, legend, legend_position)
        metal1 = merged[, c("rsid", "logp1", "chr", "pos", "label")]
        colnames(metal1)[which(colnames(metal1) == "logp1")] = "logp"
        p2 = make_locuszoom(metal1, title1, chr, color, shape, size, 
                            lz_ylab_linebreak)
        metal2 = merged[, c("rsid", "logp2", "chr", "pos", "label")]
        colnames(metal2)[which(colnames(metal2) == "logp2")] = "logp"
        p3 = make_locuszoom(metal2, title2, chr, color, shape, size, 
                            lz_ylab_linebreak)
        if (combine) {
            p2 = p2 + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
            p4 = cowplot::plot_grid(p2, p3, align = "v", nrow = 2, 
                                    rel_heights = c(0.8, 1))
            p5 = cowplot::plot_grid(p1, p4)
            return(p5)
        }
        else {
            return(list(locuscompare = p1, locuszoom1 = p2, locuszoom2 = p3))
        }
    }
    
    #### run
    cat('---read input files:\t')
    d1 = read_metal(in_fn1, marker_col1, pval_col1)
    d2 = read_metal(in_fn2, marker_col2, pval_col2)
    merged = merge(d1, d2, by = "rsid", suffixes = c("1", "2"), all = FALSE)
    cat('OK!\n---read plink files:\t')
    if (is.character(bed)){
        plink<-read.plink(bed=bed, select.snps=merged$rsid)
    }else if (is.list(bed) & all(names(bed) == c('genotypes','fam','map'))){
        plink<-bed
    }else{
        stop("The argument \"bed\" must be a string or a list of result from read.plink")
    }
    cat('OK!\n---make locuscompare:\t')
    merged = data.frame(merged,
                        plink$map[match(merged$rsid,plink$map$snp.name),c("chromosome","position")])
    colnames(merged)[colnames(merged) == "chromosome"] = "chr"
    colnames(merged)[colnames(merged) == "position"] = "pos"
    merged$pos<-merged$pos/1e6
    snp<-get_lead_snp(merged,snp=snp)
    snp.index<-which(merged$rsid==snp)
    merged<-subset(merged,chr==merged$chr[snp.index] & 
                       pos>=merged$pos[snp.index]-region/2 &
                       pos<=merged$pos[snp.index]+region/2)
    chr<-unique(merged$chr)
    R2 <- ld(plink$genotypes[,merged$rsid],plink$genotypes[,snp],stats="R.squared")
    ld<-data.frame(SNP_A=colnames(R2),SNP_B=rownames(R2),R2=R2[,1])
    cat('OK!\n---make combined plot:\t')
    p = make_combined_plot(merged, title1, title2, ld, chr, snp, combine=FALSE, 
                           legend, legend_position, lz_ylab_linebreak)
    p1<-p$locuscompare
    p2<-p$locuszoom1 
    p3<-p$locuszoom2
    if (combine) {
        p2 = p2 + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
        p4 = cowplot::plot_grid(p2, p3, align = "v", nrow = 2, 
                                rel_heights = c(0.8, 1))
        p5 = cowplot::plot_grid(p1, p4)
        cat('done!\n')
        return(p5)
    }
    else {
        cat('done!\n')
        return(list(locuscompare = p1, locuszoom1 = p2, locuszoom2 = p3))
    }
}