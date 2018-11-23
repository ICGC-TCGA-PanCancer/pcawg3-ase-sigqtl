
post_process_data = function(ase, type_cols) {
	ase$isase = as.numeric(ase$padj < 0.05)
	ase$isase_normal = as.numeric(ase$padj_normal < 0.05)
	ase$aselogR = log2(ase$ht1_count / ase$ht2_count)
	ase$gt_lead_eqtl[is.na(ase$gt_lead_eqtl)] = 0
	ase$isconserved = as.numeric(ase$phylop > 0)
	ase$ishet = as.numeric(ase$gt_lead_eqtl == 1)
	timing = as.character(ase$timing)
	timing[is.na(timing)] = 'none'
	ase$timing = factor(timing, levels = c('none', 'late', 'early'))

	type_cols_ht1 = paste(type_cols, 'ht1', sep = '_')
	type_cols_ht2 = paste(type_cols, 'ht2', sep = '_')

	ase$burden_ht1 = apply(ase[, type_cols_ht1], 1, sum)
	ase$burden_ht2 = apply(ase[, type_cols_ht2], 1, sum)
	ase$burden_diff = ase$burden_ht1 - ase$burden_ht2
	ase$burden_total = ase$burden_ht1 + ase$burden_ht2
	ase$cnratio = ase$ht1_cn / (ase$ht2_cn + ase$ht1_cn)
	ase$iscn = as.numeric((ase$ht1_cn!=1) | (ase$ht2_cn!=1))
	ase$ase_flipped = pmax(ase$ht1_count, ase$ht2_count) / (ase$ht1_count + ase$ht2_count)
	return(ase)
}

detect_outlier_samples = function(ase) {
	agg = aggregate(rep(1, nrow(ase)), by = list(rna_seq_aliquot_id = ase$rna_seq_aliquot_id, hist = ase$histology_tier2), sum)

	no_histology = agg$rna_seq_aliquot_id[agg$hist == '<unknown>']

	tmp = (unique(ase[, c('rna_seq_aliquot_id', 'purity')]))
	tmp = setNames(tmp$purity, tmp$rna_seq_aliquot_id)
	agg$purity = tmp[agg$rna_seq_aliquot_id]

	less_than_200_genes = agg$rna_seq_aliquot_id[agg$x < 200]

	return(union(no_histology, less_than_200_genes))
}
