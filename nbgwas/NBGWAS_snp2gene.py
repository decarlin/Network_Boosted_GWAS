import pandas as pandas

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
	return SNP_summary.set_index('Marker')


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

