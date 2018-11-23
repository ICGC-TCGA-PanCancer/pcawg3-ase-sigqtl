library(lme4)
library(ape)
library(ggplot2)
library(reshape2)
library(RevoScaleR)
rxOptions(numCoresToUse = -1)
rxSetComputeContext('localpar')
library(colorspace)
library(biomaRt)
ensembl = useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", GRCh = 37)
library(viridis)
library(hexbin)
library(car)
library(gplots)
library(clusterProfiler)
library(vegan)
library(fgsea)
library(stringr)
library(boot)
library(xlsx)


options(max.print = 2000)

figfolder = './'
datafolder = './'
metadatafolder = './'
srcfolder = './'
suppfolder = './'

source(paste(srcfolder, 'funcdefs.R', sep=''))
source(paste(srcfolder, 'import_data.R', sep=''))

type_cols = readLines(paste(suppfolder, 'type_cols.txt', sep = ''))
type_cols_ht1 = paste(type_cols, 'ht1', sep = '_')
type_cols_ht2 = paste(type_cols, 'ht2', sep = '_')

## annotations
sample.annotation = read.table(paste(metadatafolder, '1188_wl_donors_annotations.tsv', sep='/'), sep='\t', quote=NULL, header=TRUE)
rownames(sample.annotation) = sample.annotation$wgs_aliquot_id

## read table of cancer genes
cancer.genes = read.table(paste(suppfolder, 'cosmic.tsv', sep = '/'), sep='\t', header = TRUE)
cancer.genes = str_extract(cancer.genes$Synonyms, 'ENSG[0-9]*')
cancer.genes = cancer.genes[!is.na(cancer.genes)]

imprinting.table = read.table(paste(suppfolder, 'imprinted_genes.txt', sep = '/'), sep = '\t', header = TRUE)
imprinted.genes = imprinting.table$Ensembl.ID[imprinting.table$Imprinting.status=='I']
imprinted.genes = as.character(unique(imprinted.genes[imprinted.genes!='-']))

#seqtls = read.xlsx2(paste(datafolder, 'SF_cis_somatic_eGenes_FDR005.xlsx', sep='/'), 1, colClasses = c(rep('character',3), 'integer', 'character', rep('numeric',4),'character','character'))
seqtls = read.table(paste(datafolder, 'allsamples.sum.vaf.cis.FDR005.Gencode.annotations.tsv', sep='/'), sep='\t', quote='', header=TRUE, comment.char='')
seqtls$ensembl = gsub('\\.[0-9]*', '', seqtls$eGeneID)
seqtls$Ensemble_Reg_build_annotation = as.character(seqtls$Ensemble_Reg_build_annotation)
seqtls$Ensemble_Reg_build_annotation[is.na(seqtls$Ensemble_Reg_build_annotation)] = 'not annotated'

seqtl.genes.tested = gsub('\\.[0-9]*', '', readLines(paste(datafolder, 'all.genes.tested.txt', sep = '/')))

sigqtl = read.xlsx2(paste(datafolder, 'ST2.xls', sep = '/'), 1)
sigqtl$geneid = gsub('\\.[0-9]*', '', sigqtl$Ensembl.ID)
sig.genes = split(sigqtl$geneid, list(sigqtl$Signature))

tsgs = read.table(paste(suppfolder, 'tsgs.txt', sep = '/'), header = TRUE, row.names = NULL, quote = '', sep = '\t')
ts.genes = as.character(tsgs$geneid[(tsgs$Classification_mRNA == 'C1' | tsgs$Classification_mRNA == 'C2')])

mark_cols_ref = readLines(paste(suppfolder, 'mark_cols_ref.txt', sep = ''))

## BEGIN READ
ase = rxImport(paste(datafolder, "asesomatic_with_normal_dlvalues.tsv", sep = ''), delimiter = "\t", firstRowIsColNames = TRUE, stringsAsFactors = FALSE)
ase = post_process_data(ase, type_cols)

mark_cols = colnames(ase)[grepl('mark_', colnames(ase))]

ase = ase[, !(colnames(ase) %in% setdiff(mark_cols, mark_cols_ref))]
ase[,type_cols] = ase[,type_cols_ht1] + ase[,type_cols_ht2]

#ase = ase[, !(colnames(ase) %in% union(type_cols_ht1, type_cols_ht2))]
ase = ase[, !(colnames(ase) %in% c('type_start_lost_ht1', 'type_start_lost_ht2', 'type_stop_lost_ht1', 'type_stop_lost_ht2'))]
#ase = ase[, !(colnames(ase) %in% c('avg_counts','avg_counts_normal', 'AFphased', 'cnlogR', 'isconserved', 'ismutated', 'iscn', 'burden_ht1', 'burden_ht2'))]

## filter out HLA genes
ann = getBM(attributes = c('ensembl_gene_id', 'external_gene_name','chromosome_name', 'start_position', 'end_position', 'transcript_length'), filters = 'ensembl_gene_id', values = unique(ase$geneid), mart = ensembl)
hla_genes = ann[grepl('^HLA', ann$external_gene_name),]$ensembl_gene_id
ase = ase[!(ase$geneid %in% hla_genes),]
ase$gene_length = (ann$end_position - ann$start_position)[match(ase$geneid, ann$ensembl_gene_id)]
ase$transcript_length = ann$transcript_length[match(ase$geneid, ann$ensembl_gene_id)]
#rm(ann)

ase$iscancergene = ase$geneid %in% cancer.genes
ase$isimprinted = ase$geneid %in% imprinted.genes
ase$has_seqtl = ase$geneid %in% seqtls$ensembl
ase$has_sigqtl = ase$geneid %in% sigqtl$geneid
ase$istsg = ase$geneid %in% ts.genes

ase$seqtl_reg = seqtls$Ensemble_Reg_build_annotation[match(ase$geneid, seqtls$ensembl)]
ase$seqtl_reg[is.na(ase$seqtl_reg)] = 'none'
ase$seqtl_reg = factor(ase$seqtl_reg, levels = c('none', 'not annotated', "CTCF Binding Site", "Enhancer", "Open chromatin", "Promoter", "Promoter Flanking Region", "TF binding site"))

ase$seqtl_cat = as.character(seqtls$type[match(ase$geneid, seqtls$ensembl)])
ase$seqtl_cat[is.na(ase$seqtl_cat)] = 'none'
ase$seqtl_cat = factor(ase$seqtl_cat, levels = c('none', 'flanking', 'introns', 'exons'))

p = ase$pval_normal
p = ase$pval_normal + min(ase$pval_normal[ase$pval_normal!=0], na.rm=TRUE)
normal = aggregate(data.frame(logp = -log10(p[!is.na(p)]), isase = as.numeric(ase$padj_normal[!is.na(ase$padj_normal)] < 0.05), n = rep(1, sum(!is.na(ase$padj_normal)))), by = list(geneid = ase$geneid[!is.na(ase$padj_normal)]), sum)
normal$fracreg = (normal$isase + 1) / (normal$n + 2)
normal$logp = normal$logp / normal$n
rownames(normal) = normal$geneid
ase$frac_normal = normal[ase$geneid, 'fracreg']
ase$frac_normal[is.na(ase$frac_normal)] = mean(ase$frac_normal[!is.na(ase$frac_normal)])
ase$logp_normal = normal[ase$geneid, 'logp']
ase$logp_normal[is.na(ase$logp_normal)] = -log10(0.5)
rm(p)

ase$logp = - log10(ase$pval + min(ase$pval[ase$pval != 0]))
ase$logp_std = qqnorm(ase$logp, plot.it = FALSE)$x
ase$logp_normal_std = qqnorm(ase$logp_normal, plot.it = FALSE)$x
ase$logp_normal = -log10(ase$pval_normal + min(ase$pval_normal[ase$pval_normal != 0], na.rm=TRUE))

ase = ase[(ase$ht1_count + ase$ht2_count) >= 15,]
ase[, mark_cols_ref] = abs(ase[, mark_cols_ref])
#ase = ase[complete.cases(ase[,mark_cols_ref]),]

outliers = detect_outlier_samples(ase)

ase = ase[!(ase$rna_seq_aliquot_id %in% outliers),]

ase$histology_abbreviation = sample.annotation[ase$wgs_aliquot_id, 'histology_abbreviation']
gc()

saveRDS(ase, paste(datafolder, 'ase.revision.Rsave', sep='/'))
##ase = readRDS(paste(datafolder, 'ase.revision.Rsave', sep='/'))
# sort(sapply(ls(), function(x) { object.size(get(x)) }))

## END READ

## check ASE of APOBEC genes
#apobec.genes = c('ENSG00000128383','ENSG00000179750','ENSG00000111732','ENSG00000239713', 'ENSG00000243811', 'ENSG00000100298','ENSG00000178896')
#apobec = aggregate(ase[ase$geneid %in% apobec.genes,c('isase','ase_flipped','ht1_cn','ht2_cn')], by=list(geneid=ase$geneid[ase$geneid %in% apobec.genes]), mean)

## plot AEI per cancer type
#plotdat = aggregate(data.frame(isase = ase$isase, iscn = as.numeric((ase$ht1_cn != 1) | (ase$ht2_cn != 1))), by = list(rna_seq_aliquot_id = ase$rna_seq_aliquot_id, hist = ase$histology_tier2), mean)
plotdat = aggregate(data.frame(isase = ase$isase, iscn = as.numeric(ase$cnratio != 0.5)), by = list(rna_seq_aliquot_id = ase$rna_seq_aliquot_id, hist = ase$histology_abbreviation), mean)
plotdat$hist = as.character(plotdat$hist)
plotdat$hist = gsub("/", "\n", plotdat$hist)
plotdat.agg = aggregate(plotdat$isase, by = list(hist = plotdat$hist), median)
plotdat$hist = factor(plotdat$hist, levels = plotdat.agg$hist[order(plotdat.agg$x)])
plotdat = melt(plotdat, id.vars = c('rna_seq_aliquot_id', 'hist'))
plotdat$variable = as.character(plotdat$variable)
plotdat$variable = gsub('iscn', 'SCNA', plotdat$variable)
plotdat$variable = gsub('isase', 'AEI', plotdat$variable)
p = ggplot(plotdat, aes(x = hist, y = value, fill = variable)) + ylab('rel. frequency') +
geom_boxplot(notch = TRUE) + theme_bw() + xlab(NULL) + theme(legend.position = c(0.05, 0.9), legend.title = element_blank(), legend.background = element_blank()) +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
	panel.background = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)) +
scale_fill_brewer(palette = 'Set1')
print(p)
ggsave(paste(figfolder, "ase_dist_per_histology.pdf", sep = ''), plot = p, width = 15, height = 5, scale = 1, useDingbats = FALSE, bg = 'transparent')

## plot tumour normal comparisons
plotdat = aggregate(ase[!is.na(ase$isase_normal), c('isase', 'isase_normal')], by = list(rna_seq_aliquot_id = ase$rna_seq_aliquot_id[!is.na(ase$isase_normal)]), mean)
plotdat = melt(plotdat, id.vars='rna_seq_aliquot_id', value.name='asefreq', variable.name='cat')
levels(plotdat$cat)[levels(plotdat$cat) == "isase"] = "tumour"
levels(plotdat$cat)[levels(plotdat$cat) == "isase_normal"] = "normal"
plotdat$cat = factor(plotdat$cat, levels(plotdat$cat)[c(2, 1)])
p = ggplot(plotdat, aes(x = cat, y = asefreq, fill = cat)) + ylab('AEI frequency per patient (5% FDR)') +
geom_boxplot(notch = TRUE) + theme_bw() + xlab(NULL) + theme(legend.position = 'none') +
scale_fill_brewer(palette = "Set1", direction = -1) +
geom_point(alpha=0.3) +
geom_line(aes(group = rna_seq_aliquot_id), size = 0.05, alpha = 0.15) +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
	panel.background = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0))
print(p)
ggsave(paste(figfolder, "ase_dist_normal_vs_tumour.pdf", sep = ''), plot = p, width = 5, height = 5, scale = 1, useDingbats = FALSE, bg = 'transparent')


## show correlation CN with ASE
plotdat = aggregate(ase[, c('iscn', 'isase')], by = list(ase$wgs_aliquot_id), sum)
p = ggplot(plotdat, aes(x = iscn, y = isase)) +
theme_bw() +
geom_point() +
#stat_smooth(method='lm', formula='y~x + I(x^2) + I(x^3)') +
stat_smooth(method = 'lm') +
theme(legend.position = c(0.118, 1), legend.justification = c(1, 1)) +
labs(x = "number of genes with SCNA", y = "number of genes with ASE")
print(p)
ggsave(paste(figfolder, "regression_ase_vs_scna.pdf", sep = ''), plot = p, width = 6, height = 6, scale = 1, bg = 'transparent', useDingbats = FALSE)

## quantify relative effect sizes of SCNA, SNVs (timed) and eQTL
GLMmut.df = data.frame(
		wgs_aliquot_id = ase$wgs_aliquot_id,
		baf_phased = ase$baf_phased,
		burden_total = ase$burden_total,
		iscn = ase$iscn,
		cnratio = ase$cnratio,
		ishet = ase$ishet,
		aseratio = ase$ase_phased,
		frac_normal = ase$frac_normal,
		timingearly = as.numeric(ase$timing == 'early'),
		timinglate = as.numeric(ase$timing == 'late'),
		isimprinted = as.integer(ase$isimprinted),
		total_counts = ase$ht1_count + ase$ht2_count,
		sites_per_gene = ase$sites_per_gene
)
GLMmut.df$timingearly[is.na(GLMmut.df$timingearly)] = 0
GLMmut.df$timinglate[is.na(GLMmut.df$timinglate)] = 0
GLMmut.df$timingearly = GLMmut.df$timingearly * GLMmut.df$burden_total
GLMmut.df$timinglate = GLMmut.df$timinglate * GLMmut.df$burden_total
GLMmut.df$burden_coding = sign(ase$type_missense_variant + ase$type_stop_gained + ase$type_synonymous_variant)
GLMmut.df$burden_noncoding = sign(GLMmut.df$burden_total - GLMmut.df$burden_coding)

GLMmut.df = cbind(isase = ase$isase, standardize(GLMmut.df[,-1]))
GLMmut = glm(isase ~ cnratio + ishet + burden_coding + burden_noncoding + isimprinted, family = 'binomial', data = GLMmut.df)
summary(GLMmut)


plotdat = as.data.frame(summary(GLMmut)$coefficients)[2:6,]
colnames(plotdat) = c('est', 'se', 'z', 'p')
plotdat$pos = plotdat$est + (0.02 * max(plotdat$est))
#labels = c('SCNA', 'eQTL het.', 'fraction normal', 'SNV\n(late)', 'SNV\n(early)')
labels = c('SCNA', 'eQTL het.', 'SNV\n(coding)', 'SNV\n(non-coding)', 'imprinting')
#plotdat$type = factor(labels, levels = labels[c(1, 3, 2, 5, 4)])
plotdat$type = factor(labels, levels = labels[c(1, 2, 4, 3, 5)])
relvarexplained = setNames(plotdat$est / sum(plotdat$est), rownames(plotdat))
limits = aes(ymax = est + se, ymin = est - se)
write.table(plotdat, file = paste(figfolder, 'global_regression.tsv', sep = ''), sep = '\t')
p = ggplot(plotdat, aes(fill = type, y = est, x = type, ymax = est + se, ymin = est - se)) + theme_bw() +
geom_bar(stat = 'identity') +
geom_text(data = plotdat, aes(x = type, y = pos, label = sprintf('%.2f', est)), colour = "black") +
ylab('effect size') +
xlab('') +
theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.title = element_blank(), legend.position = 'bottom', legend.key = element_blank()) +
scale_fill_brewer(palette = "Set1")
print(p)
ggsave(paste(figfolder, "global_regression.pdf", sep = ''), plot = p, width = 5, height = 5, scale = 1)
rm(GLMmut.df, GLMmut)
gc()
## MUTATION CATEGORIES

## model of cancergene status
#GLMbinom.form = sprintf('iscancergene ~ cnratio + purity + ishet + seqtl_cat + %s', paste(type_cols, collapse = ' + '))
#fixed = GLMbinom.form
#GLMbinom = glm(as.formula(fixed), data = ase , family = binomial, y = FALSE, x = FALSE, model = FALSE)
#pbinom = plot_linear_model(GLMbinom, fdr = 0.05, colors = NULL, fixAspect = FALSE, coef.to.plot = c(FALSE, TRUE, rep(TRUE, sum(!is.na(GLMbinom$coefficients)) - 2)))
#print(pbinom)
#summary(GLMbinom)
#ggsave(paste(figfolder, "logistic_regression_cosmic.pdf", sep = ''), plot = pbinom + xlab('change in log-odds'), width = 10, height = 10, scale = 1, useDingbats = FALSE, bg = 'transparent')
#ase$pred_cancer = predict(GLMbinom, type='response')

## global model of ASE including germline
GLMglobal.df = ase[, c('ht1_cn', 'ht2_cn', 'avg_counts', 'isimprinted', 'wgs_aliquot_id', 'geneid', 'isase', 'burden_total', 'burden_diff', 'sites_per_gene', 'transcript_length', 'gene_length', 'cnratio', 'baf_phased', 'purity', 'ishet', 'has_sigqtl', 'has_seqtl', 'ht1_count', 'ht2_count', 'histology_tier2', type_cols)]
GLMglobal.df$type_downstream_gene_variant_50kb_plus = GLMglobal.df$type_downstream_gene_variant_50kb + GLMglobal.df$type_downstream_gene_variant_60kb + GLMglobal.df$type_downstream_gene_variant_70kb + GLMglobal.df$type_downstream_gene_variant_80kb + GLMglobal.df$type_downstream_gene_variant_90kb + GLMglobal.df$type_downstream_gene_variant_100kb
GLMglobal.df$type_upstream_gene_variant_50kb_plus = GLMglobal.df$type_upstream_gene_variant_50kb + GLMglobal.df$type_upstream_gene_variant_60kb + GLMglobal.df$type_upstream_gene_variant_70kb + GLMglobal.df$type_upstream_gene_variant_80kb + GLMglobal.df$type_upstream_gene_variant_90kb + GLMglobal.df$type_upstream_gene_variant_100kb
type_cols2 = c('type_upstream_gene_variant_50kb_plus', 'type_upstream_gene_variant_40kb', 'type_upstream_gene_variant_30kb', 'type_upstream_gene_variant_20kb', 'type_upstream_gene_variant_10kb', 'type_promoter_variant', 'type_5_prime_UTR_variant', 'type_intron_variant', 'type_splice_region_variant', 'type_missense_variant', 'type_synonymous_variant', 'type_stop_gained', 'type_3_prime_UTR_variant', 'type_downstream_gene_variant_10kb', 'type_downstream_gene_variant_20kb', 'type_downstream_gene_variant_30kb', 'type_downstream_gene_variant_40kb', 'type_downstream_gene_variant_50kb_plus')
GLMglobal.df = GLMglobal.df[, ! (colnames(GLMglobal.df) %in% setdiff(type_cols, type_cols2))]
GLMglobal.df$purity = GLMglobal.df$purity
GLMglobal.df$cnratio = GLMglobal.df$cnratio - mean(GLMglobal.df$cnratio)
GLMglobal.df$baf_phased = GLMglobal.df$baf_phased - mean(GLMglobal.df$baf_phased)
GLMglobal.df$isimprinted = as.integer(GLMglobal.df$isimprinted)
GLMglobal.df$total_counts = (GLMglobal.df$ht1_count + GLMglobal.df$ht2_count)
#GLMglobal.df$total_count = pmin(GLMglobal.df$ht1_count, GLMglobal.df$ht2_count)
#GLMglobal.df[, type_cols2] = (GLMglobal.df[, type_cols2] * 1000 )/ (GLMglobal.df$gene_length)
type_cols_cds = c('type_synonymous_variant', 'type_missense_variant', 'type_stop_gained')
#GLMglobal.df[, type_cols_cds] = (GLMglobal.df[, type_cols_cds] * 1000) / (GLMglobal.df$transcript_length)
#GLMglobal.df$type_intron_variant = (GLMglobal.df$type_intron_variant * 1000) / pmax(1,(GLMglobal.df$gene_length-GLMglobal.df$transcript_length))

GLMglobal.form = sprintf('isase ~ cnratio + isimprinted + log(sites_per_gene) + log(avg_counts) + log(total_counts) + purity + ishet + %s', paste(type_cols2, collapse = ' + '))
GLMglobal = glm(as.formula(GLMglobal.form), data = GLMglobal.df, family = binomial)
summary(GLMglobal)
write.table(as.data.frame(summary(GLMglobal)$coefficients), file=paste(figfolder, 'logistic_regression_ase_types.tsv', sep=''), sep='\t')

#ase$pred_ase_global = predict(GLMglobal, type = 'response')
pglobal = plot_linear_model(GLMglobal, fdr = 0.15, coef.thresh=0.1, colors = NULL, fixAspect = FALSE, coef.to.plot = c(FALSE, FALSE, TRUE))
print(pglobal)
summary(GLMglobal)
ggsave(paste(figfolder, "logistic_regression_ase_types.pdf", sep = ''), plot = pglobal + xlab('change in log-odds'), width = 7.5, height = 10, scale = 1, useDingbats = FALSE, bg = 'transparent')

## predict globally
terms = predict(GLMglobal, type='terms')
terms = as.data.frame(terms)
colnames(terms) = paste('term_', colnames(terms), sep = '')
pred_ase_snvs = apply(terms[, grepl('term_type', colnames(terms))], 1, sum)
#pred_ase_somatic = pred_ase_snvs + apply(terms[, c('term_cnratio','term_baf_phased')], 1, sum)
pred_ase_somatic = pred_ase_snvs + apply(terms[, 'term_cnratio',drop=FALSE], 1, sum)
predictions = data.frame(wgs_aliquot_id=ase$wgs_aliquot_id, geneid=ase$geneid, pred_ase_snvs, pred_ase_somatic)
write.table(predictions, paste(datafolder,'ase_predictions.revision.tsv', sep='/'), sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)

## cross validation for model
#cost = function(r, pi = 0) mean(abs(r - pi) > 0.5)
#CV2 = cv.glm(ase, GLMglobal, K=10, cost=cost)

## test if our prediction enriches for cancer genes

noncoding = seqtls$ensembl[((seqtls$type == 'flanking') & ((seqtls$Gencode_annotation == 'intergenic') | (seqtls$Gencode_annotation == 'intronic')))]

aggregate.ase = function(xterms, x) {
	agg = aggregate(cbind(xterms, x[, c('ht1_count', 'ht2_count', 'gene_length', 'sites_per_gene', 'iscn', 'ase_phased', 'ase_flipped', 'cnratio', 'burden_total', 'burden_diff', type_cols)]), by = list(geneid = x$geneid), mean)
	rownames(agg) = agg$geneid
	agg$asecount = aggregate(x$isase, by = list(geneid = x$geneid), sum)$x
	agg$n = aggregate(rep(1, nrow(x)), by = list(geneid = x$geneid), sum)$x
	agg$asefrac = agg$asecount / agg$n
	agg$iscancergene = agg$geneid %in% cancer.genes
	agg$miss_stop = agg$type_missense_variant + agg$type_stop_gained
	agg$isimprinted = agg$geneid %in% imprinted.genes
	agg$has_seqtl = agg$geneid %in% seqtls$ensembl
	agg$has_seqtl_noncoding = agg$geneid %in% noncoding
	agg$has_sigqtl = agg$geneid %in% sigqtl$geneid
	agg$sigqtl = setNames(sigqtl$Signature, sigqtl$geneid)[agg$geneid]
	agg$total_count = agg$ht1_count + agg$ht2_count
	agg$istsg = agg$geneid %in% ts.genes

	agg$pred_ase_germline = agg$term_ishet
	agg$pred_ase_global = apply(agg[, colnames(agg)[grepl('term_', colnames(agg))]], 1, sum)
	agg$pred_ase_snvs = apply(agg[, colnames(agg)[grepl('term_type', colnames(agg))]], 1, sum)
	agg$pred_ase_scnas = agg[, 'term_cnratio']
	agg$pred_ase_somatic = agg$pred_ase_snvs + agg$pred_ase_scnas
	agg$pred_ase_test = agg$pred_ase_germline + agg$pred_ase_somatic


	agg = agg[agg$n>=5, ]
	agg = agg[order(agg$pred_ase_snvs, decreasing = TRUE),]

	agg.ann = getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'chromosome_name', 'start_position', 'end_position'), filters = 'ensembl_gene_id', values = agg$geneid, mart = ensembl)
	agg$contig = agg.ann$chromosome_name[match(agg$geneid, agg.ann$ensembl_gene_id)]

	return(agg)
}

agg = aggregate.ase(terms, ase)
write.table(agg, paste(datafolder, 'ase_aggregate.revision.tsv', sep = '/'), sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

mf=5
p = plot.enrichment(agg[agg$n >= mf, c('miss_stop', 'pred_ase_somatic', 'asefrac', 'pred_ase_germline')], agg$iscancergene[agg$n >= mf], 0)
p = p + ylim(c(0, 0.08)) + xlim(c(500, 5000))
print(p)
ggsave(paste(figfolder, "cancer_gene_enrichment.pdf", sep = ''), plot = p, width = 5, height = 5, scale = 1, useDingbats = FALSE, bg = 'transparent')

p = plot.enrichment(agg[, c('miss_stop', 'pred_ase_snvs', 'pred_ase_scnas', 'asefrac', 'pred_ase_germline', 'pred_ase_global')], agg$iscancergene, 0)
p = p + ylim(c(0, 0.08)) + xlim(c(500, 5000)) + theme(legend.background=element_rect(fill='transparent'))
print(p)
ggsave(paste(figfolder, "cancer_gene_enrichment_detailed.pdf", sep = ''), plot = p, width = 5, height = 5, scale = 1, useDingbats = FALSE, bg = 'transparent')

n=round(nrow(agg)*0.1)
genes.somatic = agg$geneid[order(agg$pred_ase_somatic, decreasing = TRUE)][1:n]
genes.snvs = agg$geneid[order(agg$pred_ase_snvs, decreasing = TRUE)][1:n]
genes.scnas = agg$geneid[order(agg$pred_ase_scnas, decreasing = TRUE)][1:n]
genes.mutated = agg$geneid[order(agg$miss_stop, decreasing = TRUE)][1:n]
genes.burden = agg$geneid[order(agg$burden_total, decreasing = TRUE)][1:n]
genes.lofgof = agg$geneid[order(agg$miss_stop, decreasing = TRUE)][1:n]
genes.observed = agg$geneid[order(agg$asefrac, decreasing = TRUE)][1:n]
genes.germline = agg$geneid[order(agg$pred_ase_germline, decreasing = TRUE)][1:n]

venn(list(lofgof = genes.lofgof, burden = genes.burden, snvs = genes.snvs), universe = unique(agg$geneid))
venn(list(somatic=genes.somatic, snvs=genes.snvs, scnas=genes.scnas), universe = unique(agg$geneid))
boxplot(log(agg$gene_length) ~ I(agg$geneid %in% genes.snvs), notch=TRUE, xlab='is gene in top genes?', ylab='log(gene length)')
boxplot(log(agg$gene_length) ~ I(agg$geneid %in% genes.scnas), notch = TRUE, xlab = 'is gene in top genes?', ylab = 'log(gene length)')
barplot(table(agg$contig[agg$geneid %in% genes.snvs]))
barplot(table(agg$contig[agg$geneid %in% genes.scnas]))

gidc(genes.somatic)

test = chisq.test(agg$geneid %in% genes.somatic, agg$geneid %in% cancer.genes)
test$observed
test$expected
test

test = chisq.test(agg$geneid %in% genes.snvs, agg$geneid %in% cancer.genes)
test$observed
test$expected
test

test = chisq.test(agg$geneid %in% genes.scnas, agg$geneid %in% cancer.genes)
test$observed
test$expected
test

observed.table = table(gene_has_AEI=agg$geneid %in% genes.observed, gene_in_COSMIC=agg$geneid %in% cancer.genes)
test = chisq.test(observed.table)
test$observed
test$expected
test

test = chisq.test(agg$geneid %in% genes.germline, agg$geneid %in% cancer.genes)
test$observed
test$expected
test

test = chisq.test(agg$geneid %in% genes.mutated, agg$geneid %in% cancer.genes)
test$observed
test$expected
test


## somatic ASE
GO = enrichGO(genes.somatic, 'org.Hs.eg.db', universe = agg$geneid, keyType = 'ENSEMBL', ont = "BP", pvalueCutoff = 1, pAdjustMethod = "fdr", qvalueCutoff = 1)
goterms = as.data.frame(GO)
goterms = annotate.go(goterms)
goterms = goterms[order(goterms$pvalue, decreasing=FALSE),]

GO.burden = enrichGO(genes.lofgof, 'org.Hs.eg.db', universe = agg$geneid, keyType = 'ENSEMBL', ont = "BP", pvalueCutoff = 1, pAdjustMethod = "fdr", qvalueCutoff = 1)
goterms.burden = as.data.frame(GO.burden)
goterms.burden = annotate.go(goterms.burden)
goterms.burden = goterms.burden[order(goterms.burden$pvalue, decreasing = FALSE),]

#goterms.chosen = rownames(goterms)[goterms$p.adjust < 0.05]
#goterms.chosen = union(rownames(goterms)[goterms$p.adjust < 0.05], rownames(goterms.burden)[goterms.burden$p.adjust < 0.05])
#goterms.chosen = intersect(goterms.chosen, rownames(goterms.burden))
goterms.chosen = union(rownames(goterms)[1:20], rownames(goterms.burden)[1:20])


plotdat = data.frame(
	ID = goterms.chosen,
	Description = goterms[goterms.chosen, 'Description'],
	logFC = goterms[goterms.chosen, 'logFC'],
	padj = -log10(goterms[goterms.chosen, 'pvalue'])
)
plotdat$type = 'somatic AEI'
plotdat2 = data.frame(
	ID = goterms.chosen,
	Description = goterms[goterms.chosen, 'Description'],
	logFC = goterms.burden[goterms.chosen, 'logFC'],
	padj = -log10(goterms.burden[goterms.chosen, 'pvalue'])
)
plotdat2$type = 'LoF/GoF mutations'
plotdat$logFCdiff = plotdat$logFC - plotdat2$logFC
plotdat2$logFCdiff = plotdat$logFC - plotdat2$logFC
plotdat = rbind(plotdat, plotdat2)
#plotdat$GOterm = paste(plotdat$ID, plotdat$Description, sep=': ')

plotdat$GOterm = plotdat$Description
plotdat$GOterm = reorder(plotdat$GOterm, plotdat$logFCdiff, function(x)x[1])

p = ggplot(plotdat, aes(fill = type, y = logFC, x = GOterm)) + theme_bw() +
geom_bar(stat = 'identity', position = 'dodge') +
ylab('log fold change over expectation') +
xlab('') +
theme(axis.text.x = element_text(angle = -90, hjust=0), legend.title = element_blank(), legend.position = 'bottom', legend.key = element_blank()) +
scale_fill_brewer(palette='Set1') + coord_flip()
print(p)
ggsave(paste(figfolder, "ASE_somatic_GO_enrichment.pdf", sep = ''), plot = p, width = 7.5, height = 10, scale = 1, useDingbats = FALSE, bg = 'transparent')

## cumulative plot of how many genes are measured at least x times in the cohort
cum = aggregate(cbind(n = rep(1, nrow(ase)), isase = ase$isase), by = list(geneid = ase$geneid), sum)
cum = cum[order(cum$n, decreasing = TRUE),]
cum$index = 1:nrow(cum)
p = ggplot(cum, aes(y = n, x = index)) + theme_bw() + geom_step() + xlab('number of genes') + ylab('number of samples')
print(p)

## TSGs
length(intersect(agg$geneid, ts.genes))

test = chisq.test(seqtl.genes.tested %in% ts.genes, seqtl.genes.tested %in% seqtls$ensembl)
test$observed
test$expected
test

test = chisq.test(agg$geneid %in% genes.somatic, agg$geneid %in% ts.genes)
test$observed
test$expected
test

snvs2tsgs = table(snv=agg$geneid %in% genes.snvs, tsg=agg$geneid %in% ts.genes)
test = chisq.test(snvs2tsgs)
test$observed
test$expected
test

plotdat = as.data.frame(prop.table(snvs2tsgs, 2))
p = ggplot(plotdat, aes(fill = snv, x = tsg, y = Freq)) + theme_bw() +
	geom_bar(stat = 'identity', position = 'stack')
print(p)

ts.snv = intersect(genes.snvs, ts.genes)

ts2timing = table(timing = ase$timing, snv = (ase$geneid %in% ts.snv))
test = chisq.test(ts2timing)
test$observed
test$expected
test

plotdat = as.data.frame(prop.table(ts2timing, 2))
p = ggplot(plotdat, aes(fill = timing, x = snv, y=Freq)) + theme_bw() +
	geom_bar(stat = 'identity', position = 'stack')
print(p)

plotdat = ase[, c('timing', 'geneid')]
plotdat$tsg = plotdat$geneid %in% ts.genes
p = ggplot(plotdat, aes(fill = timing, x = tsg)) + theme_bw() +
	geom_bar(position = 'stack')
print(p)

## GSEA
gseaParam = 1
values = setNames(qqnorm(agg$pred_ase_snvs[agg$geneid %in% genes.snvs], plot.it = FALSE)$x, agg$geneid[agg$geneid %in% genes.snvs])
GSEA = fgsea(list(cancertestisgenes = ts.genes), values, nperm = 2000, gseaParam = gseaParam)
GSEA
p = plotEnrichment(ts.genes, values, gseaParam = gseaParam) + labs(title = 'Cancer testis gene enrichment in somatic ASE (SNVs)', subtitle = sprintf('NES=%.2f, p=%.4f', GSEA$NES, GSEA$padj))
plot(p)
ggsave(paste(figfolder, "GSEA_ct_genes_in_snvs.pdf", sep = ''), plot = p, width = 8, height = 4, scale = 1, useDingbats = FALSE, bg = 'transparent')

gseaParam = 1
values = setNames(qqnorm(agg$pred_ase_scnas[agg$geneid %in% genes.scnas], plot.it = FALSE)$x, agg$geneid[agg$geneid %in% genes.scnas])
GSEA = fgsea(list(cancertestisgenes = ts.genes), values, nperm = 2000, gseaParam = gseaParam)
GSEA
p = plotEnrichment(ts.genes, values, gseaParam = gseaParam) + labs(title = 'Cancer testis gene enrichment in somatic ASE (SCNAs)', subtitle = sprintf('NES=%.2f, p=%.4f', GSEA$NES, GSEA$padj))
plot(p)
ggsave(paste(figfolder, "GSEA_ct_genes_in_scnas.pdf", sep = ''), plot = p, width = 8, height = 4, scale = 1, useDingbats = FALSE, bg = 'transparent')

gseaParam = 1
values = setNames(qqnorm(agg$asefrac, plot.it = FALSE)$x, agg$geneid)
GSEA = fgsea(list(cancergenes = cancer.genes), values, nperm = 2000, gseaParam = gseaParam)
GSEA
p = plotEnrichment(cancer.genes, values, gseaParam = gseaParam) + labs(title = 'Cancer gene enrichment in observed ASE', subtitle = sprintf('NES=%.2f, p=%.4f', GSEA$NES, GSEA$padj))
plot(p)
ggsave(paste(figfolder, "GSEA_is_ase.pdf", sep = ''), plot = p, width = 8, height = 4, scale = 1, useDingbats = FALSE, bg = 'transparent')

gseaParam = 1
values = setNames(qqnorm(agg$pred_ase_global, plot.it = FALSE)$x, agg$geneid)
GSEA = fgsea(list(cancergenes = cancer.genes), values, nperm = 2000, gseaParam = gseaParam)
GSEA
p = plotEnrichment(cancer.genes, values, gseaParam = gseaParam) + labs(title = 'Cancer gene enrichment in predicted total ASE', subtitle = sprintf('NES=%.2f, p=%.4f', GSEA$NES, GSEA$padj))
plot(p)
ggsave(paste(figfolder, "GSEA_predicted_ase_global.pdf", sep = ''), plot = p, width = 8, height = 4, scale = 1, useDingbats = FALSE, bg = 'transparent')

gseaParam = 1
values = setNames(qqnorm(agg$pred_ase_germline, plot.it = FALSE)$x, agg$geneid)
GSEA = fgsea(list(cancergenes = cancer.genes), values, nperm = 2000, gseaParam = gseaParam)
GSEA
p = plotEnrichment(cancer.genes, values, gseaParam = gseaParam) + labs(title = 'Cancer gene enrichment in predicted germline ASE', subtitle = sprintf('NES=%.2f, p=%.4f', GSEA$NES, GSEA$padj))
plot(p)
ggsave(paste(figfolder, "GSEA_predicted_ase_germline.pdf", sep = ''), plot = p, width = 8, height = 4, scale = 1, useDingbats = FALSE, bg = 'transparent')

gseaParam = 1
values = setNames(qqnorm(agg$pred_ase_somatic, plot.it = FALSE)$x, agg$geneid)
GSEA = fgsea(list(cancergenes = cancer.genes), values, nperm = 2000, gseaParam = gseaParam)
GSEA
p = plotEnrichment(cancer.genes, values, gseaParam = gseaParam) + labs(title = 'Cancer gene enrichment in predicted somatic ASE (SNVs + SCNAs)', subtitle = sprintf('NES=%.2f, p=%.4f', GSEA$NES, GSEA$padj))
plot(p)
ggsave(paste(figfolder, "GSEA_predicted_ase_somatic.pdf", sep = ''), plot = p, width = 8, height = 4, scale = 1, useDingbats = FALSE, bg = 'transparent')
leadingEdge_pred_ase_somatic = GSEA$leadingEdge[[1]]

gseaParam = 1
values = setNames(qqnorm(agg$pred_ase_snvs, plot.it=FALSE)$x, agg$geneid)
GSEA = fgsea(list(cancergenes = cancer.genes), values, nperm = 2000, gseaParam = gseaParam)
GSEA
p = plotEnrichment(cancer.genes, values, gseaParam = gseaParam) + labs(title = 'Cancer gene enrichment in predicted somatic ASE (SNVs)', subtitle = sprintf('NES=%.2f, p=%.4f', GSEA$NES, GSEA$padj))
plot(p)
ggsave(paste(figfolder, "GSEA_predicted_ase_snvs.pdf", sep = ''), plot = p, width = 8, height = 4, scale = 1, useDingbats = FALSE, bg = 'transparent')
leadingEdge_pred_ase_selected = GSEA$leadingEdge[[1]]

gseaParam = 1
values = setNames(qqnorm(agg$pred_ase_scnas, plot.it = FALSE)$x, agg$geneid)
GSEA = fgsea(list(cancergenes = cancer.genes), values, nperm = 2000, gseaParam = gseaParam)
GSEA
p = plotEnrichment(cancer.genes, values, gseaParam = gseaParam) + labs(title = 'Cancer gene enrichment in predicted somatic ASE (SCNAs)', subtitle = sprintf('NES=%.2f, p=%.4f', GSEA$NES, GSEA$padj))
plot(p)
ggsave(paste(figfolder, "GSEA_predicted_ase_scnas.pdf", sep = ''), plot = p, width = 8, height = 4, scale = 1, useDingbats = FALSE, bg = 'transparent')

values = setNames(qqnorm(agg$miss_stop, plot.it=FALSE)$x, agg$geneid)
p = plotEnrichment(cancer.genes, values, gseaParam = gseaParam) + labs(title = 'Cancer gene enrichment in LoF/GoF mutations', subtitle = sprintf('NES=%.2f, p=%.4f', GSEA$NES, GSEA$padj))
plot(p)
ggsave(paste(figfolder, "GSEA_lof_gof_mutations.pdf", sep = ''), plot = p, width = 8, height = 4, scale = 1, useDingbats = FALSE, bg = 'transparent')
GSEA = fgsea(list(cancergenes = cancer.genes), values, nperm = 2000, gseaParam = gseaParam)
GSEA
leadingEdge_totalBurden = GSEA$leadingEdge[[1]]


## directed model
GLMlin.df = data.frame(
	geneid = ase$geneid,
	wgs_aliquot_id = ase$wgs_aliquot_id,
	ishet = ase$ishet,
	purity = ase$purity,
	nmut = ase$nmut,
	ht1_cn = ase$ht1_cn,
	ht2_cn = ase$ht2_cn,
	baf_phased = ase$baf_phased,
	burden_diff = ase$burden_diff,
	aselogR = ase$aselogR,
	ase_phased = ase$ase_phased,
	histology_tier1 = ase$histology_tier1,
	histology_tier2 = ase$histology_tier2,
	histology_tier4 = ase$histology_tier4,
	transcript_length = ase$transcript_length,
	gene_length = ase$gene_length
)

GLMlin.df[,type_cols] = ase[, type_cols_ht1] - ase[,type_cols_ht2]
GLMlin.df = GLMlin.df[(GLMlin.df$nmut > 0),]

form.types = sprintf('aselogR ~ ht1_cn + ht2_cn + purity + log(gene_length) + log(transcript_length) + ishet + burden_diff + %s', paste(type_cols, collapse = ' + '))
GLMlin = lm(as.formula(form.types), data = GLMlin.df)
summary(GLMlin)

plin = plot_linear_model(GLMlin, colors = NULL, fixAspect = FALSE, coef.to.plot = c(FALSE, FALSE, FALSE, TRUE))
print(plin)
ggsave(paste(figfolder, "linear_regression_ase_types.pdf", sep = ''), plot = plin + xlab('directed effect size'), width = 7.5, height = 10, scale = 1, useDingbats = FALSE, bg = 'transparent')

## normal validation
filter = (!is.na(ase$ase_phased_normal)) & ((ase$ht1_count_normal + ase$ht2_count_normal) >= 5) & (ase$burden_diff!=0) #&
#(ase$type_3_prime_UTR_variant != 0 | ase$type_5_prime_UTR_variant != 0 | ase$type_stop_gained != 0 | ase$type_splice_region_variant!=0)

asenorm = ase[filter,]
dim(asenorm)

asenorm$ishet = 0
asenorm$gene_length = 1
asenorm$transcript_length = 1

pred = predict(GLMlin, asenorm)
r = cor(asenorm$ase_phased_normal + pred, asenorm$ase_phased)

asenorm.nomut = asenorm
asenorm.nomut[, type_cols] = 0
pred.nomut = predict(GLMlin, asenorm.nomut)
r.nomut = cor(asenorm.nomut$ase_phased_normal + pred.nomut, asenorm.nomut$ase_phased)

asenorm.nocnv = asenorm
asenorm.nocnv[, c("ht1_cn","ht2_cn")] = 0
pred.nocnv = predict(GLMlin, asenorm.nocnv)
r.nocnv = cor(asenorm.nocnv$ase_phased_normal + pred.nocnv, asenorm.nocnv$ase_phased)

bardata = c(r, r.nomut, r.nocnv)
bardata

plotdat = data.frame(type = c('full model', 'w/o SNV', 'w/o SCNA'), r = bardata, pos = bardata + 0.01)
plotdat$type = factor(plotdat$type, level = plotdat$type[c(1, 2, 3)])
p = ggplot(plotdat, aes(fill = type, y = r, x = type)) + theme_bw() +
	geom_bar(stat = 'identity') +
	geom_text(data = plotdat, aes(x = type, y = pos, label = sprintf('%.4f', r)), colour = "black") +
	ylab('r') +
	xlab('') +
	theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.title = element_blank(), legend.position = 'bottom', legend.key = element_blank()) +
	scale_fill_brewer(palette = "Set1")
print(p)
ggsave(paste(figfolder, "normal_validation.pdf", sep = ''), plot = p, width = 4, height = 6, scale = 1)

## NMD plots
grp1 = (ase$cnratio == 0.5) & (ase$type_stop_gained != 0) & (ase$type_synonymous_variant == 0)
grp2 = (ase$cnratio == 0.5) & (ase$type_stop_gained == 0) & (ase$type_synonymous_variant != 0)

plotdat = data.frame(ase = c(ase$ase_flipped[grp1], ase$ase_flipped[grp2]))
plotdat$type = c(rep('stop gained', sum(grp1)), rep('synonymous', sum(grp2)))
p = ggplot(plotdat, aes(x = type, y = ase, fill = type)) + ylab('absolute allelic imbalance') +
geom_boxplot(notch = TRUE) + theme_bw() + scale_fill_brewer(palette = "Set1") + xlab(NULL) + theme(legend.position = 'none') +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
print(p)
ggsave(paste(figfolder, "stop_gain_vs_syn.pdf", sep = ''), plot = p, width = 4, height = 4, scale = 0.7, useDingbats = FALSE, bg = 'transparent')
