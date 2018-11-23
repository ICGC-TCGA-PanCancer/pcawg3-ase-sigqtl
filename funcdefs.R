library(ggplot2)
library(colorspace)

type2text = function(x) {
	txt = gsub("log\\((.*)\\)", "\\1 [log]", x)
	txt = gsub("\\(", "", txt)
	txt = gsub("\\)", "", txt)
	txt = gsub("isase", "observed AEI", txt)
	txt = gsub("isimprinted", "imprinted", txt)
	txt = gsub("asefrac", "observed AEI recurrence", txt)
	txt = gsub("miss_stop", "LoF / GoF mutations", txt)
	txt = gsub("total_counts", "read depth (gene)", txt)
	txt = gsub('avg_counts', 'read depth (sample)', txt)
	txt = gsub("sites_per_gene", "no. ASE sites", txt)
	txt = gsub("pred_ase_germline", "predicted germline AEI", txt)
	txt = gsub("pred_ase_selected", "predicted AEI (SCNA, 5'UTR, 3'UTR, missense, PTC)", txt)
	txt = gsub("pred_ase_somatic", "predicted somatic AEI", txt)
	txt = gsub("pred_ase_global", "predicted AEI (germline + somatic)", txt)
	txt = gsub("pred_ase_scnas", "predicted somatic AEI (SCNAs only)", txt)
	txt = gsub("pred_ase_snvs", "predicted somatic AEI (SNVs only)", txt)
	txt = gsub("pred_cancer", "pred. cancer gene", txt)
	txt = gsub("type", "", txt)
	txt = gsub("mark", "", txt)
	txt = gsub("_prime", "'", txt)
	txt = gsub("is_loh", 'LOH', txt)
	txt = gsub("frac_normal", "fraction AEI in normal tissue", txt)
	txt = gsub("_plus", "+", txt)
	txt = gsub("_", " ", txt)
	txt = gsub("cnratio", "copy-number ratio", txt)
	txt = gsub("cnlogR", "copy-number ratio", txt)
	txt = gsub("cn", "copy-number", txt)
	txt = gsub("totalcn", "total copy-number", txt)
	txt = gsub("abs\\(", "", txt)
	txt = gsub("ishet", "lead eQTL het.", txt)
	txt = gsub("ishetTRUE", "lead eQTL het.", txt)
	txt = gsub("isconserved", "site conserved", txt)
	txt = gsub("timingcode", "mutation is early", txt)
	txt = gsub("iscn", "gene has SCNA", txt)
	txt = gsub("ismutated", "gene has SNV", txt)
	txt = gsub("distcatdist", "distance to TSS <", txt)
	txt = gsub("distcat", "", txt)
	txt = gsub("rel", "\\(relative\\)", txt)
	txt = gsub("diff", "\\(difference\\)", txt)
	txt = gsub("ref", "", txt)
	txt = gsub("gene variant", "", txt)
	txt = gsub("variant", "", txt)
	txt = gsub("[[:space:]]+", " ", txt)
	txt = gsub("upstream", "-", txt)
	txt = gsub("downstream", "\\+", txt)
	return(txt)
}

plot_linear_model = function(model, colors = NULL, fixAspect = FALSE, font.size = 12, coef.to.plot = NULL, sort.coef = FALSE, fdr = 0.05, coef.thresh=0.1) {
  
	# Create result dataframe
	#result = as.data.frame(summary(model)$coefficients[-1,, drop = FALSE])
	result = as.data.frame(summary(model)$coefficients)
	colnames(result) = c("value", "ci", "zval", "pval")
	result$padj = p.adjust(result$pval, method = 'BH')
	tmp = sapply(rownames(result), function(x) { gsub("type_", "", x) })
	result$type = factor(tmp, levels = tmp)
  
	# Limit coefficients to plot
	if (!is.null(coef.to.plot)) {
		if (length(coef.to.plot) > nrow(result)) {
			coef.to.plot = coef.to.plot[1:nrow(result)]
		} else if (length(coef.to.plot) < nrow(result)) {
			coef.to.plot = c(coef.to.plot, rep(coef.to.plot[length(coef.to.plot)], nrow(result)-length(coef.to.plot)))
		}
		result  = result[coef.to.plot,, drop = FALSE]
	}
  
	# Order results
	result = result[nrow(result):1,]
	tmp = rownames(result)
	if (sort.coef) {
		result = result[order(result$type, decreasing = TRUE),]
	}
  
	# Plot results
	return(plot_linear_model_result(result, colors, fixAspect, font.size, fdr, coef.thresh))
}

plot_linear_model_result = function(result, colors = NULL, fixAspect = FALSE, font.size = 12, fdr = 0.05, coef.thresh=0.1) {

	thelabels = unique(type2text(result$type))
	style = rep('plain', length(thelabels))
	style[result$padj < fdr & abs(result$value) >= coef.thresh] = 'bold'

	p = ggplot(result, aes(x = value, y = type, color = type)) + geom_point() +
	geom_errorbarh(aes(xmin = value - ci, xmax = value + ci, height = 0.2)) +
	xlab("") + ylab("") + theme_bw() + theme(legend.position = "none", axis.text = element_text(size = font.size, face = style)) +
	scale_y_discrete(limits = result$type, labels = thelabels) +
	geom_vline(aes(xintercept = 0), linetype = 'dotted', colour = 'grey50')
	if (fixAspect) {
		p = p + theme(aspect.ratio = 1)
	}
	if (!is.null(colors)) {
		p = p + scale_colour_manual(name = "type", values = colors)
	}
	return(p)
}

plot_mixed_model = function(model, colors = NULL, fixAspect = TRUE, font.size = 12, coef.to.plot = NULL, sort.coef = FALSE) {
	if (is.null(coef.to.plot)) {
		result = as.data.frame(summary(model)$coefficients[-1,, drop = FALSE])
	} else {
		result = as.data.frame(summary(model)$coefficients[coef.to.plot,, drop = FALSE])
	}
	result = result[nrow(result):1,]
	colnames(result) = c("value", "ci", "zval")
	result$type = sapply(rownames(result), function(x) { gsub("type_", "", x) })
	if (sort.coef) {
		result = result[order(result$type, decreasing = TRUE),]
	}
	p = ggplot(result, aes(x = value, y = type, color = type)) + geom_point() +
	geom_errorbarh(aes(xmin = value - ci, xmax = value + ci, height = 0.2)) +
	xlab("size of effect on ASE ratio") + ylab("") + theme_bw() + theme(legend.position = "none", axis.text = element_text(size = font.size)) +
	scale_y_discrete(limits = result$type, labels = unique(type2text(result$type)))
	if (fixAspect) {
		p = p + theme(aspect.ratio = 1)
	}
	if (!is.null(colors)) {
		p = p + scale_colour_manual(name = "type", values = colors)
	}
	return(p)
}

## boxplots for germline/somatic effects
plotgene = function(gene, data, residual = FALSE) {

	tmp = data[which(data$geneid == gene),]
	if (nrow(tmp) == 0) {
		return(NULL)
	}
	genename = getBM(attributes = c('external_gene_name', 'description'), filters = 'ensembl_gene_id', values = gene, mart = ensembl)
	title = sprintf('%s / %s / %s', gene, genename$external_gene_name, gsub(" ?\\[.*\\]", '', genename$description))

	grp1 = tmp$nmutselect == 0
	grp2 = !grp1
	plotdat = data.frame(ase = c(tmp$aseflipped[grp1], tmp$aseflipped[grp2]))
	plotdat$type = c(rep('wild type', sum(grp1)), rep('mutated', sum(grp2)))
	p = ggplot(plotdat, aes(x = type, y = ase, fill = type)) + ylab('absolute allelic imbalance') +
	geom_boxplot(notch = TRUE) + theme_bw() + scale_fill_brewer(palette = "Set1") + xlab(NULL) + theme(legend.position = 'none') +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) +
	#geom_jitter(height=0, width=0.2) +
	ggtitle(title)
	return(p)
}

standardize = function(x, center = TRUE) {
	if (center) {
		y = apply(x, 2, function(x) {(x - mean(x)) / sd(x) })
	} else {
		y = apply(x, 2, function(x) { x / sd(x) })
	}
	if (is.data.frame(x)) {
		y = as.data.frame(y)
	}
	return(y)
}

center = function(x, ...) {
	y = apply(x, 2, function(x) { x - mean(x, ...) })
	if (is.data.frame(x)) {
		y = as.data.frame(y)
	}
	return(y)
}

randomize = function(x) {
	return(x[sample(1:nrow(x)),])
}

print.geneids = function(x) {
	cat(paste(x, collapse = '\n'))
}

write.geneids = function(x, filename) {
	writeLines(paste(unique(x), collapse = '\n'), filename)
}

gidc = function(x) {
	writeClipboard(paste(x, collapse = '\n'))
}

plot.enrichment = function(x, geneset, pc = 0) {
	xord = as.data.frame(apply(x, 2, function(z) {(cumsum(as.numeric(geneset[order(z, decreasing = TRUE)])) + pc) / ((1:length(z)) + pc) }))
	xord$x = 1:nrow(xord)
	plotdat = melt(xord, id.vars = 'x')
	levels(plotdat$variable) = type2text(levels(plotdat$variable))
	p = ggplot(plotdat, aes(col = variable, y = value, x = x)) + theme_bw() +
	geom_point() +
	ylab('fraction cancer genes') + xlab('rank') +
	theme(legend.title = element_blank(), legend.box.margin = margin(c(5, 5, 5, 5)), legend.justification = c(1, 1), legend.position = c(1, 1), legend.key = element_blank()) +
	scale_colour_brewer(palette = "Set1")
	return(p)

}


plot.enrichment.bar = function(x, geneset, n = 500) {
	xord = as.data.frame(apply(x, 2, function(z) { as.numeric(geneset[order(z, decreasing = TRUE)]) }))
	plotdat = melt(apply(xord[1:n,], 2, sum))
	plotdat$variable = factor(type2text(rownames(plotdat)))
	plotdat$variable = factor(plotdat$variable, levels = levels(plotdat$variable)[c(1, 4, 2, 3)])
	p = ggplot(plotdat, aes(fill = variable, y = value, x = variable)) + theme_bw() +
	geom_bar(stat = 'identity') +
	ylab('number of cancer genes') + xlab('') +
	theme(legend.position = 'none') +
	scale_fill_brewer(palette = "Accent") +
	geom_hline(yintercept = sum(geneset) / nrow(x) * n, colour = 'black', linetype = 2)
	return(p)
}

reactomePathwaysENSEMBL = function(genes) {
	stopifnot(requireNamespace("reactome.db"))
	stopifnot(requireNamespace("AnnotationDbi"))
	pathways <- na.omit(AnnotationDbi::select(reactome.db::reactome.db,
		keys = genes, c("PATHID"), keytype = "ENSEMBL"))
	pathways <- split(pathways$ENTREZID, pathways$PATHID)
	pathway2name <- as.data.table(na.omit(AnnotationDbi::select(reactome.db::reactome.db,
		names(pathways), c("PATHNAME"), "PATHID")))
	PATHNAME = NULL
	pathway2name[, `:=`(PATHNAME, sub("^[^:]*: ", "", PATHNAME))]
	names(pathways) <- pathway2name$PATHNAME
	pathways
}


annotate.go = function(terms) {
	terms$nHits = as.numeric(sapply(strsplit(terms$GeneRatio, '/'), function(x) { x[1] }))
	terms$nGenes = as.numeric(sapply(strsplit(terms$GeneRatio, '/'), function(x) { x[2] }))
	terms$nSet = as.numeric(sapply(strsplit(terms$BgRatio, '/'), function(x) { x[1] }))
	terms$nUniverse = as.numeric(sapply(strsplit(terms$BgRatio, '/'), function(x) { x[2] }))
	terms$Expect = terms$nGenes * (terms$nSet / terms$nUniverse)
	terms$logFC = log2(terms$nHits / terms$Expect)
	return(terms)
}
