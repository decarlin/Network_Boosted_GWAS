import pandas as pd

# Load SNP p-value association table from file
def load_SNP_pvals(snp_pval_file, delimiter='\t', header=False, cols='0,1,2,3'):
	# Check for valid 'cols' parameter
	try:
		cols_idx = [int(c) for c in cols.split(',')]
	except:
		raise ValueError('Invalid column index string')
	# Load gene_pos_file
	if header:
		SNP_summary = pd.read_csv(snp_pval_file, delimiter=delimiter)
	else:
		SNP_summary = pd.read_csv(snp_pval_file, delimiter=delimiter, header=-1)
	# Check gene positions table format
	if (SNP_summary.shape[1] < 4) | (max(cols_idx) >  SNP_summary.shape[1]-1):
		raise ValueError('Not enough columns in SNP Summary File')
	# Construct gene position table
	SNP_summary = SNP_summary[cols_idx]
	SNP_summary.columns = ['Marker', 'Chr', 'Pos', 'P-Value']
	return SNP_summary


# Load gene positions from file
def load_gene_pos(gene_pos_file, delimiter='\t', header=False, cols='0,1,2,3'):
	# Check for valid 'cols' parameter
	try:
		cols_idx = [int(c) for c in cols.split(',')]
	except:
		raise ValueError('Invalid column index string')
	# Load gene_pos_file
	if header:
		gene_positions = pd.read_csv(gene_pos_file, delimiter=delimiter)
	else:
		gene_positions = pd.read_csv(gene_pos_file, delimiter=delimiter, header=-1)
	# Check gene positions table format
	if (gene_positions.shape[1] < 4) | (max(cols_idx) >  gene_positions.shape[1]-1):
		raise ValueError('Not enough columns in Gene Positions File')
	# Construct gene position table
	gene_positions = gene_positions[cols_idx]
	gene_positions.columns = ['Gene', 'Chr', 'Start', 'End']
	return gene_positions.set_index('Gene')

# Returns the closest gene to each SNP on each chromosome
def closest_gene(SNP_summary, gene_positions, window):
	snplist = list(SNP_summary.index)
	snp_closest_genes = []
	for SNP in snplist:
		SNP_info = SNP_summary.ix[SNP]
		marker = SNP_info['Marker']
		chrom = SNP_info['Chr']
		pos = SNP_info['Pos']
		# Get all genes on same chromosome
		gene_positions_filt1 = gene_positions[gene_positions['Chr']==chrom]
		# Get distance of all genes on same chromosome to SNP
		gene_distances = abs(gene_positions_filt1['Start']-pos)
		# Get closest gene by absolute genomic distance to SNP
		closest_gene = gene_distances.argmin()
		closest_gene_dist = gene_distances.ix[closest_gene]
		snp_closest_genes.append([SNP, chrom, pos, SNP_info['P-Value'], closest_gene, closest_gene_dist])
	snp_closest_gene_table = pd.DataFrame(snp_closest_genes, columns = ['rsID', 'Chr', 'SNP Pos', 'SNP P-Value', 'Closest Gene', 'Closest Gene Distance'])
	# Filter all SNPs where the closest gene is not within the specified window size
	dist = window*1000
	snp_closest_gene_table_filt = snp_closest_gene_table[snp_closest_gene_table['Closest Gene Distance'] < dist]
	# Assign the best SNP p-value for each gene
	genelist = list(set(snp_closest_gene_table_filt['Closest Gene']))
	closest_gene_min_p = []
	for gene in genelist:
		gene_SNPs = snp_closest_gene_table_filt[snp_closest_gene_table_filt['Closest Gene']==gene]
		min_p_data = gene_SNPs.ix[gene_SNPs['SNP P-Value'].argmin()]
		gene_info = gene_positions.ix[gene]
		chrom = gene_info['Chr']
		start = gene_info['Start']
		stop = gene_info['End']
		closest_gene_min_p.append([gene, chrom, start, stop, gene_SNPs.shape[0], min_p_data['rsID'], int(min_p_data['SNP Pos']), min_p_data['SNP P-Value']])
	min_p_table = pd.DataFrame(closest_gene_min_p, columns = ['Gene', 'Chr', 'Gene Start', 'Gene End', 'nSNPs', 'TopSNP', 'TopSNP Pos', 'TopSNP P-Value'])
	min_p_table = min_p_table.set_index('Gene').dropna().sort_values(by=['TopSNP P-Value', 'Chr', 'Gene Start'])
	return min_p_table		


# Filters list of gene distances to SNPs to only SNPs where the transcription start site
# is within a specified window of the SNP
def closest_gene_filter(gene_distance_table, window):
	dist = window*1000

# Assign p-values to genes by the minimum p-value of all SNPs near a single gene
def min_p(SNP_summary, gene_positions, window):
	dist = window*1000
	genelist = list(gene_positions.index)
	min_p_list = []
	for gene in genelist:
		gene_info = gene_positions.ix[gene]
		chrom = gene_info['Chr']
		start = gene_info['Start']
		stop = gene_info['End']
		# Get all SNPs on same chromosome
		SNP_summary_filt1 = SNP_summary[SNP_summary['Chr']==chrom]
		# Get all SNPs after window start position
		SNP_summary_filt2 = SNP_summary_filt1[SNP_summary_filt1['Pos'] >= start-dist]
		# Get all SNPs before window end position
		SNP_summary_filt3 = SNP_summary_filt2[SNP_summary_filt2['Pos'] <= stop+dist]
		# Get min_p statistics for this gene
		if len(SNP_summary_filt3) >= 1:
			min_p_data = SNP_summary_filt3.ix[SNP_summary_filt3['P-Value'].argmin()]
			min_p_list.append([gene, chrom, start, stop, SNP_summary_filt3.shape[0], min_p_data['Marker'], int(min_p_data['Pos']), min_p_data['P-Value']])
		else:
			min_p_list.append([gene, chrom, start, stop, 0, None, None, None])
	min_p_table = pd.DataFrame(min_p_list, columns = ['Gene', 'Chr', 'Gene Start', 'Gene End', 'nSNPs', 'TopSNP', 'TopSNP Pos', 'TopSNP P-Value'])
	min_p_table = min_p_table.set_index('Gene').dropna().sort_values(by=['TopSNP P-Value', 'Chr', 'Gene Start'])
	return min_p_table