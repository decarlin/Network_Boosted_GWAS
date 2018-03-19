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