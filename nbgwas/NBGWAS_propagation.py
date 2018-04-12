import networkx as nx
import numpy as np
import os
import pandas as pd
import time

# Load 2-column network from edge list file 
def load_network(network_file, delimiter='\t'):
	load_start = time.time()
	network = nx.read_edgelist(network_file, delimiter='\t', data=False)
	network_nodes = list(network.nodes())
	network_node_degree = pd.Series(dict(network.degree()), name='Degree')
	print 'Network File Loaded:', time.time()-load_start, 'seconds'
	print 'Number of network nodes:', len(network_nodes)
	print 'Number of network edges:', len(network.edges())
	return network, network_nodes, network_node_degree

# Shuffle network by preserving node-degree
def degree_shuffNet(network):
	shuff_time = time.time()
	edge_len=len(network.edges())
	shuff_net=network.copy()
	try:
		nx.double_edge_swap(shuff_net, nswap=edge_len, max_tries=edge_len*10)
	except:
		print 'Note: Maximum number of swap attempts ('+repr(edge_len*10)+') exceeded before desired swaps achieved ('+repr(edge_len)+').'
	# Evaluate Network Similarity
	shared_edges = len(set(network.edges()).intersection(set(shuff_net.edges())))
	print 'Network shuffled:', time.time()-shuff_time, 'seconds. Edge similarity:', shared_edges/float(edge_len)
	return shuff_net	

# Calculate optimal propagation coefficient
# Model from Huang and Carlin et al 2018
def calculate_alpha(network, m=-0.02935302, b=0.74842057):
    log_edge_count = np.log10(len(network.edges()))
    alpha_val = round(m*log_edge_count+b,3)
    if alpha_val <=0:
        raise ValueError('Alpha <= 0 - Network Edge Count is too high')
        # There should never be a case where Alpha >= 1, as avg node degree will never be negative
    else:
        return alpha_val

# Create degree-normalized adjacency matrix
def create_norm_adj_mat(network, symmetric_norm=False):
	adj_mat = nx.adjacency_matrix(network)
	adj_array = np.array(adj_mat.todense())
	if symmetric_norm:
		D = np.diag(1/np.sqrt(sum(adj_array)))
		adj_array_norm = np.dot(np.dot(D, adj_array), D)
		print 'Adjacency matrix normalized (symmetric)'
	else:
		degree = sum(adj_array)
		adj_array_norm = (adj_array*1.0/degree).T
		print 'Adjacency matrix normalized (asymmetric - default)'		
	return adj_array_norm		

# Closed form random-walk propagation (as seen in HotNet2) for each subgraph: Ft = (1-alpha)*Fo * (I-alpha*norm_adj_mat)^-1
# Concatenate to previous set of subgraphs
def fast_random_walk(alpha, binary_mat, subgraph_norm, prop_data):
    term1=(1-alpha)*binary_mat
    term2=np.identity(binary_mat.shape[1])-alpha*subgraph_norm
    term2_inv = np.linalg.inv(term2)
    subgraph_prop = np.dot(term1, term2_inv)
    return np.concatenate((prop_data, subgraph_prop), axis=1)

# Construct influence matrix of each network node propagated across network to use as kernel in AUPRC analysis
# Input: NetowkrkX object. No propagation constant or alpha model required, can be calculated
def construct_network_kernel(network, alpha=None, outdir=None):
    # Calculate propagation coefficient alpha
    if alpha is None:
        alpha_val = calculate_alpha(network)
    else:
        alpha_val = alpha
    print "Alpha = ", alpha_val

    # Calculate the network kernel by performing closed form random-walk network propagation using an identity matrix
    # Construct identity matrix
    starttime = time.time()
    network_nodes = network.nodes()
    I = pd.DataFrame(data=np.identity(len(network_nodes)), index=network_nodes, columns=network_nodes)
    # Separate network into connected components and calculate propagation values of each sub-sample on each connected component
    subgraphs = list(nx.connected_component_subgraphs(network))
    # Initialize propagation results by propagating first subgraph
    subgraph = subgraphs[0]
    subgraph_nodes = list(subgraph.nodes)
    prop_data_node_order = list(subgraph_nodes)
    I_filt = np.array(I.T.ix[subgraph_nodes].fillna(0).T)
    subgraph_norm = normalize_network(subgraph, symmetric_norm=symmetric_norm)
    prop_data_empty = np.zeros((I_filt.shape[0], 1))
    prop_data = fast_random_walk(network_alpha, I_filt, subgraph_norm, prop_data_empty)
    # Get propagated results for remaining subgraphs
    for subgraph in subgraphs[1:]:
        subgraph_nodes = list(subgraph.nodes)
        prop_data_node_order = prop_data_node_order + subgraph_nodes
        I_filt = np.array(I.T.ix[subgraph_nodes].fillna(0).T)
        subgraph_norm = normalize_network(subgraph, symmetric_norm=symmetric_norm)
        prop_data = fast_random_walk(network_alpha, I_filt, subgraph_norm, prop_data)
    # Return propagated result as dataframe
    prop_data_df = pd.DataFrame(data=prop_data[:,1:], index = binary_matrix.index, columns=prop_data_node_order)
    # Save network kernel to file
    prop_data_df = prop_data_df.ix[prop_data_df.columns]
    print 'Propagated network kernel constructed', time.time()-starttime, 'seconds'
    if outdir is not None:
        prop_data_df.to_hdf(outdir+'network_kernel.hdf', key='Kernel', mode='w')
        print 'Propagated network kernel saved', outdir+'network_kernel.hdf'
    return prop_data_df    

# Construct binary propagation vector with top % threshold or p-value.
# Note: If using a p_thresh, pct_thresh must be set to None, otherwise having a p_thresh and pct_thresh both defined will throw error
def set_network_seeds(gene_pval_table, network_genes, pct_thresh = 0.01, p_thresh = None):
    if (pct_thresh is not None) & (p_thresh is None):
        threshold_genes = gene_pval_table.ix[:int(np.ceil(pct_thresh*gene_pval_table.shape[0]))]
        prop_vector = (gene_pval_table.set_index('Gene', drop=False).ix[network_genes]['Gene'].isin(threshold_genes['Gene'])).astype(float)
        threshold_name = ' '.join(['Top', str(pct_thresh*100)+'%'])
    elif (pct_thresh is None) & (p_thresh is not None):
        threshold_genes = gene_pval_table[gene_pval_table['TopSNP P-Value'] < p_thresh]
        prop_vector = (gene_pval_table.set_index('Gene').ix[network_genes]['TopSNP P-Value'] < p_thresh).astype(float)
        threshold_name = ' '.join(['p', '<', str(p_thresh)])
    else:
        raise ValueError("Exactly one of 'pct_thresh' or 'p_thresh' must be defined.")
    
    prop_vector.name = threshold_name
    print threshold_genes.shape[0], 'genes above threshold ('+threshold_name+')'
    print int(prop_vector.sum()), 'threshold genes in network as seeds'
        
    return threshold_genes, pd.DataFrame(prop_vector).ix[network_genes].T

# Construct a vector of "null seeds" that is degree matched from the actual network seeds
# Uses output of set_network_seeds()
def set_null_network_seeds(network_seeds_vector, network_node_degree):
    network_seeds = list(network_seeds_vector.ix[0][network_seeds_vector.ix[0]==1].index)
    null_seeds = []
    fixed_seeds = []
    for seed_gene in network_seeds:
        seed_gene_degree = network_node_degree.ix[seed_gene]
        degree_matched_genes = list(set(network_node_degree[network_node_degree==seed_gene_degree].index)-set(null_seeds))
        if len(degree_matched_genes) == 1:
            fixed_seeds.append(seed_gene)
            null_seeds.append(seed_gene)
        else:
            degree_matched_null_gene = np.random.choice(degree_matched_genes, 1)[0]
            null_seeds.append(degree_matched_null_gene)
        null_seed_vector = pd.Series(0.0, index=network_seeds_vector.columns, name='Null Seeds')
        null_seed_vector.ix[null_seeds]=1.0
    print len(fixed_seeds), 'seed genes could not be swapped:', fixed_seeds
    return null_seeds, pd.DataFrame(null_seed_vector).ix[network_seeds_vector.columns].T

# Perform network propagation using the network kernel on seed genes with the given threshold
# Also perform network propagation on a generated degree-matched null seed gene set
# Combine propagation results with restuls of other network-based gene re-orderings
def kernel_propagation(kernel, gene_pval_table, network_node_degree, pct_thresh=0.01, p_thresh=None, outdir=None):
    # Get network genes from kernel data frame
    network_genes = list(kernel.index)

    # Determine seed genes
    threshold_genes, network_seeds_vector = set_network_seeds(gene_pval_table, network_genes, pct_thresh = pct_thresh, p_thresh = p_thresh)

    # Create vector of degree-matched random seeds
    null_seeds, null_seeds_vector = set_null_network_seeds(network_seeds_vector, network_node_degree)
    null_genes = gene_pval_table.ix[null_seeds]
    prop_seeds = pd.concat([network_seeds_vector, null_seeds_vector])[network_genes]

    # Propagate seed genes
    starttime = time.time()
    prop_val_matrix = np.dot(prop_seeds, kernel)
    prop_val_table = pd.DataFrame(prop_val_matrix, index = prop_seeds.index, columns = prop_seeds.columns)
    print 'Seed genes propagated:', time.time()-starttime, 'seconds'

    # Construct propagation result summary tables
    thresh_prop_results_table = kernel_prop_results_table(gene_pval_table, prop_val_table.ix[0], threshold_genes, network_node_degree, outdir=outdir)
    null_prop_results_table = kernel_prop_results_table(gene_pval_table, prop_val_table.ix[1], null_genes, network_node_degree, outdir=outdir)

    return thresh_prop_results_table, null_prop_results_table
    
# Construct a table summarizing the propagation results of a network propagation from seed genes
# Also re-order genes based on results from other network-based methods
def kernel_prop_results_table(gene_pval_table, prop_results, seed_genes, network_node_degree, outdir=None):
    # Get all non-seed genes in network
    network_non_seed = list(set(prop_results.index) - set(seed_genes['Gene']))

    # Combine P-values with propagation values for non-seed genes
    non_seed_gwas_p = gene_pval_table.set_index('Gene').ix[network_non_seed]['TopSNP P-Value'].fillna(1.0)
    non_seed_snp_dist = gene_pval_table.set_index('Gene').ix[network_non_seed]['SNP Distance']
    non_seed_prop_results_table = pd.concat([non_seed_gwas_p, non_seed_snp_dist, network_node_degree.ix[network_non_seed], 
                                             prop_results.T.ix[network_non_seed]], axis=1)
    non_seed_prop_results_table.columns = ['GWAS P-Value', 'GWAS SNP Dist', 'Degree', 'Prop Value']

    # Update non-seed gene propagation table with non-seed gene ranks/empirical p-values, and adjusted ranks
    non_seed_prop_results_table['Prop Rank'] = non_seed_prop_results_table['Prop Value'].rank(ascending=False)
    non_seed_prop_results_table['Prop Empirical P-Value'] = non_seed_prop_results_table['Prop Rank'] / float(non_seed_prop_results_table.shape[0])
    non_seed_prop_results_table['Adjusted P-Value'] = non_seed_prop_results_table['GWAS P-Value'].multiply(non_seed_prop_results_table['Prop Empirical P-Value'])
    non_seed_prop_results_table['Adjusted Rank'] = non_seed_prop_results_table['Adjusted P-Value'].rank()
    non_seed_prop_results_table['Seed Gene'] = False
    print 'Non-seed gene propagation results compiled'

    # Construct table for threshold genes (seed genes are a subset of these genes)
    top_genes = seed_genes[['Gene','TopSNP P-Value','SNP Distance']].set_index('Gene')
    seed_prop_results_table = pd.concat([top_genes, network_node_degree.ix[top_genes.index], 
                                         prop_results.T.ix[top_genes.index]], axis=1)
    seed_prop_results_table.columns = ['GWAS P-Value', 'GWAS SNP Dist', 'Degree', 'Prop Value']
    seed_prop_results_table['Prop Rank'] = 0
    seed_prop_results_table['Prop Empirical P-Value'] = None
    seed_prop_results_table['Adjusted P-Value'] = None
    seed_prop_results_table['Adjusted Rank'] = 0  
    seed_prop_results_table['Seed Gene'] = True
    print 'Seed gene propagation results compiled'

    # Combine seed genes table and new propagation ranking table
    combined_prop_gwas_table = pd.concat([seed_prop_results_table, non_seed_prop_results_table])
    combined_prop_gwas_table_sorted = combined_prop_gwas_table[seed_prop_results_table.columns].sort_values(by=['Seed Gene', 'Prop Rank', 'Adjusted Rank', 
                                                                                                                'GWAS P-Value','GWAS SNP Dist','Degree'],
                                                                                                         ascending=[False,True,True,True,True,False])
    print 'Propagation results constructed and sorted'
    if outdir is not None:
        fn = '_'.join(prop_results.name.split(' '))+'_prop_results.csv'
        combined_prop_gwas_table_sorted.to_csv(outdir+fn)
        print 'Propagation results saved:', outdir+fn    
    return combined_prop_gwas_table_sorted




