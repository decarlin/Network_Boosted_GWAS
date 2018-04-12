import pandas as pd

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



##### Old plotting functions that can be repurposed later ############


# results_table (output of random_walk function from network_GWAS_functions module) loaded as a pandas DataFrame
def PRC_plots(results_table, gold_standard_genes, p_thresh='5e-7', outdir=None, plot_title=None, file_prefix=None):
    # Identify propagated seeds from p_thresh (set as the same value as in random_walk function)
    prop_seeds = results_table[results_table['GWAS P'] <= p_thresh].index
    # Identify propagated seeds that are also in the gold standard
    prop_seeds_gs = list(set(prop_seeds).intersection(set(gold_standard_genes)))
    # Filter results table to remove all network seeds
    network_nodes = list(results_table.index)
    network_non_seed = list(set(network_nodes)-set(prop_seeds))
    network_non_seed_gs = list(set(network_non_seed).intersection(set(gold_standard_genes)))
    results_table_filt = results_table.ix[network_non_seed]
    # Get precision, recall, and AUPRC from GWAS, propagation, and combined GWAS x propagation results
    # Note: THE RESULTS TABLE WILL HAVE THE CORRECT COLUMN NAMES, DO NOT CHANGE THEM AFTER LOADING
    GWAS_pre, GWAS_rec, GWAS_AUPRC = PRC(results_table_filt['GWAS P'], network_non_seed_gs, True)
    prop_pre, prop_rec, prop_AUPRC = PRC(results_table_filt['Prop-Value'], network_non_seed_gs, False)
    adjp_pre, adjp_rec, adjp_AUPRC = PRC(results_table_filt['Adjusted Association'], network_non_seed_gs, True)    

    # Plot recovery results
    plt.subplots(2,2,figsize=(12,8))

    ax1 = plt.subplot(2,2,1)
    ax1.plot(GWAS_rec, GWAS_pre)
    ax1.set_yscale('log')
    ax1.set_title('GWAS Precision Recall Curve')
    ax1.set_xlabel('Recall')
    ax1.set_ylabel('Precision')

    ax2 = plt.subplot(2,2,2)
    ax2.plot(prop_rec, prop_pre)
    ax2.set_yscale('log')
    ax2.set_title('Propagation Precision Recall Curve')
    ax2.set_xlabel('Recall')
    ax2.set_ylabel('Precision')

    ax3 = plt.subplot(2,2,3)
    ax3.plot(adjp_rec, adjp_pre)
    ax3.set_yscale('log')
    ax3.set_title('GWAS x Propagation Precision Recall Curve')
    ax3.set_xlabel('Recall')
    ax3.set_ylabel('Precision')

    ax4 = plt.subplot(2,2,4)
    index, bar_width = np.arange(3), 0.7
    ax4.bar(index, [GWAS_AUPRC, prop_AUPRC, adjp_AUPRC], bar_width)
    ax4.set_xlim((bar_width-1,3))
    ax4.set_xticks(index + bar_width / 2)
    ax4.set_xticklabels(('GWAS', 'Prop', 'GWAS x Prop'))
    ax4.set_xlabel('Values')
    ax4.set_ylabel('AUPRC')

    plt.tight_layout()
    plt.subplots_adjust(top=0.9)
    plt.suptitle(plot_title, fontsize=16)
    plt.savefig(outdir+file_prefix+'_Gold_Standard_Recovery_Performance_Plots.pdf')
    plt.close()    

    return GWAS_AUPRC, prop_AUPRC, adjp_AUPRC

# Measure the distribution of the propagation scores from random_walk()
def prop_score_distributions(results_table, gs_sets, p_thresh=5e-7, weak_p_thresh=1e-4, outdir=None, plot_title=None):
    # Get union of all gold standard genes
    all_gs_genes = []
    for gs_set in gs_sets:
        all_gs_genes = all_gs_genes + list(gs_sets[gs_set])
    all_gs_genes = list(set(all_gs_genes))

    # Network gene results
    network_nodes = list(results_table.index)
    # Identify propagated seeds from p_thresh (set as the same value as in random_walk function)
    prop_seeds = results_table[results_table['GWAS P'] <= p_thresh].index
    # Identify all genes that were not network
    network_non_seed = list(set(network_nodes)-set(prop_seeds))
    results_table_filt = results_table.ix[network_non_seed]

    f, (ax1, ax2) = plt.subplots(2,1,figsize=(10,8))
    # Plot propagation score/rank distributions for all gold standard gene sets
    for gs_set in gs_sets:
        # Filter results table to get only results of gold standard genes that were not seed genes
        network_non_seed_gs = list(set(network_non_seed).intersection(set(gs_sets[gs_set])))
        # Plot propagation score/rank distributions for all gold standard gene sets
        network_non_seed_gs_values = results_table_filt.ix[network_non_seed_gs]['log10(Prop-Value)']
        network_non_seed_gs_ranks = results_table_filt.ix[network_non_seed_gs]['Prop Rank']
        sns.kdeplot(network_non_seed_gs_values, ax=ax1, label=gs_set+' Genes n=('+str(len(network_non_seed_gs))+')', shade=True)
        sns.kdeplot(network_non_seed_gs_ranks, ax=ax2, label=gs_set+' Genes n=('+str(len(network_non_seed_gs_ranks))+')', shade=True)

    # Plot propagation score/rank distribution for non-gold standard, non-seed, "weakly associated" genes
    network_weak_assoc = results_table[results_table['GWAS P'] <= weak_p_thresh].index
    network_weak_assoc_non_seed = list(set(network_non_seed).intersection(set(network_weak_assoc)))
    network_weak_assoc_non_seed_non_gs = list(set(network_weak_assoc_non_seed)-set(all_gs_genes))
    network_non_seed_assoc_values = results_table_filt.ix[network_weak_assoc_non_seed_non_gs]['log10(Prop-Value)']
    network_non_seed_assoc_ranks = results_table_filt.ix[network_weak_assoc_non_seed_non_gs]['Prop Rank']
    sns.kdeplot(network_non_seed_assoc_values, ax=ax1, label='Weakly Associated Genes n=('+str(len(network_weak_assoc_non_seed_non_gs))+')', shade=True)
    sns.kdeplot(network_non_seed_assoc_ranks, ax=ax2, label='Weakly Associated Genes n=('+str(len(network_weak_assoc_non_seed_non_gs))+')', shade=True)        
    
    # Plot propagation score/rank for all non gold standard, non-associated genes, non-propagated genes
    network_bg = list(set(network_non_seed) - set(network_weak_assoc).union(set(all_gs_genes)))
    network_bg_values = results_table_filt.ix[network_bg]['log10(Prop-Value)']
    network_bg_ranks = results_table_filt.ix[network_bg]['Prop Rank']
    sns.kdeplot(network_bg_values, ax=ax1, label='Background Genes n=('+str(len(network_bg))+')', shade=True)
    sns.kdeplot(network_bg_ranks, ax=ax2, label='Background Genes n=('+str(len(network_bg))+')', shade=True)        

    ax1.set_xlabel(r'$\mathregular{log_{10}(Propagation-Value)}$')
    ax1.set_ylabel('Density')

    ax2.set_xlim((0, len(results_table_filt)))
    ax2.set_xlabel('Propagation Value Rank')
    ax2.set_ylabel('Density')

    plt.tight_layout()
    plt.subplots_adjust(top=0.92)
    plt.suptitle(plot_title, fontsize=16)
    plt.savefig(outdir+'Propagation_Distributions.pdf')
    plt.close()

    return



#### Iterative random walk propagation that needs to be debugged ######

# Random walk function (assign propagated values to nodes)
# network is NetworkX Graph object to perform network propagation over
# GWAS_results is pandas Series (NOT DataFrame) to assign "significance threshold" for determining which genes to propagate
def random_walk(network, GWAS_results, p_thresh=5e-7, alpha=None, binary_prop=True, 
	max_iter=250, prop_min_rtol=1e-8, prop_min_atol=1e-12, save_path=None, verbose=True):

	network_nodes = network.nodes()
	# Determine list of genes to propagate over network
	# If a network node does not have a significance value given, assign that value as p=1
	network_GWAS_results = GWAS_results.ix[network_nodes].fillna(1)
	network_GWAS_results_rank = network_GWAS_results.rank()
	prop_nodes = network_GWAS_results[network_GWAS_results <= p_thresh].index

	# If there are no nodes to propagate, give empty binary result
	if len(prop_nodes)==0:
		if verbose:
			print 'No propagation seed nodes found in network'
		propagated_values = np.zeros(len(network.nodes()))
	# If there are nodes to propagate, determine propagation coefficient and perform propagation
	else:
		if verbose:
			print 'Significant genes to be propagated:', prop_nodes
		if alpha is None:
			a = calculate_alpha(network, verbose=verbose)
		else:
			a = alpha
		norm_adj_mat = create_norm_adj_mat(network, verbose=verbose)
		# Currently set as binary 1/0 for propagation level, but will be edited to include "effect size" propagation
		init_array = prop_seed(network, prop_nodes)
		propagated_values = propagate(norm_adj_mat, init_array, a, max_iter=max_iter, prop_min_rtol=prop_min_rtol, prop_min_atol=prop_min_atol, verbose=verbose)
	log_prop_values = np.log10(propagated_values)
	prop_value_table = pd.DataFrame([network.nodes(), list(propagated_values), list(log_prop_values)], 
									index=['Gene', 'Prop-Value', 'log10(Prop-Value)']).T.sort_values(by='log10(Prop-Value)', ascending=False).set_index('Gene')
	#prop_value_table = prop_value_table.reset_index(drop=True)
	prop_value_table['Prop Rank'] = prop_value_table['Prop-Value'].rank(ascending=False)
	prop_value_table['Empirical P'] = prop_value_table['Prop Rank'] / (prop_value_table.shape[0]+1)

	# Combine GWAS results and propagation results together
	prop_results_table = pd.concat([network_GWAS_results, network_GWAS_results_rank, prop_value_table], axis=1)
	prop_results_table.columns = ['GWAS P', 'GWAS Rank']+list(prop_value_table.columns)
	prop_results_table['Adjusted Association'] = prop_results_table['GWAS P'].multiply(prop_results_table['Empirical P'])
	if verbose:
		print 'Combined network propagation results constructed'

	if save_path is not None:
		prop_results_table.to_csv(save_path, sep='\t')

	return prop_results_table

# Calculate propagated values of network nodes from seed nodes
# matrix is the degree-normalized adjacency matrix
# score currently blows up for large networks like PCNet...Need to investigate further
# ********* BUG ALERT *************
# Using network propagation methods from Network_Evaluation_Tools for now
# def propagate(matrix, init_array, alpha, max_iter=250, prop_min_rtol=1e-6, prop_min_atol=1e-08, verbose=False):
# 	F = init_array.copy()
# 	F_0 = (1-alpha)*init_array
# 	F_new = alpha*np.dot(F, matrix) + F_0
# 	for i in range(max_iter):
# 		if not(np.allclose(F, F_new, rtol=prop_min_rtol, atol=prop_min_atol)):
# 			F = F_new
# 			F_new = alpha*np.dot(F, matrix) + F_0
# 		else:
# 			break
# 	if verbose:
# 		print 'Network propagation complete:', i+1, 'iterations'
# 	return F_new		
	    