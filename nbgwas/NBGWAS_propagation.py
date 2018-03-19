import networkx as nx
import numpy as np
import os
import pandas as pd
import time

# Load 2-column network from edge list file 
def load_network(network_file, delimiter='\t', verbose=False):
	load_start = time.time()
	network = nx.read_edgelist(network_file, delimiter='\t', data=False)
	if verbose:
		print 'Network File Loaded:', time.time()-load_start, 'seconds'
		print 'Number of network nodes:', len(network.nodes())
		print 'Number of network edges:', len(network.edges())
	return network

# Shuffle network by preserving node-degree
def degree_shuffNet(network, verbose=False):
	shuff_time = time.time()
	edge_len=len(network.edges())
	shuff_net=network.copy()
	try:
		nx.double_edge_swap(shuff_net, nswap=edge_len, max_tries=edge_len*10)
	except:
		if verbose:
			print 'Note: Maximum number of swap attempts ('+repr(edge_len*10)+') exceeded before desired swaps achieved ('+repr(edge_len)+').'
	if verbose:
		# Evaluate Network Similarity
		shared_edges = len(set(network.edges()).intersection(set(shuff_net.edges())))
		print 'Network shuffled:', time.time()-shuff_time, 'seconds. Edge similarity:', shared_edges/float(edge_len)
	return shuff_net	

# Calculate optimal propagation coefficient
# Model from Huang and Carlin et al 2018
def calculate_alpha(network, m=-0.02935302, b=0.74842057, verbose=False):
	log_edge_count = np.log10(len(network.edges()))
	alpha_val = round(m*log_edge_count+b,3)
	if verbose:
		print 'Calculated Alpha:', alpha_val
	if alpha_val <=0:
		raise ValueError('Alpha <= 0 - Network Edge Count is too high')
		# There should never be a case where Alpha >= 1, as avg node count will never be negative
	else:
		return alpha_val

# Create degree-normalized adjacency matrix
def create_norm_adj_mat(network, symmetric_norm=False, verbose=False):
    adj_mat = nx.adjacency_matrix(network)
    adj_array = np.array(adj_mat.todense())
    if symmetric_norm:
        D = np.diag(1/np.sqrt(sum(adj_array)))
        adj_array_norm = np.dot(np.dot(D, adj_array), D)
    else:
        degree = sum(adj_array)
        adj_array_norm = (adj_array*1.0/degree).T
    if verbose:
    	if symmetric_norm:
    		print 'Adjacency matrix normalized (symmetric)'
    	else:
    		print 'Adjacency matrix normalized (asymmetric - default)'
    return adj_array_norm		

# Find initial seed nodes for random walk propagation
def prop_seed(network, prop_nodes):
    network_nodelist = network.nodes()
    seed_list, network_seeds = [], []
    for node in network_nodelist:
        if node not in prop_nodes:
            seed_list.append(0)
        else:
            seed_list.append(1)
            network_seeds.append(node)
    return np.array(seed_list)    

# Calculate propagated values of network nodes from seed nodes
def propagate(matrix, init_array, alpha, max_iter=250, prop_min_rtol=1e-6, prop_min_atol=1e-08, verbose=False):
    F = init_array.copy()
    F_0 = (1-alpha)*init_array
    F_new = alpha*np.dot(F, matrix) + F_0
    for i in range(max_iter):
        if not(np.allclose(F, F_new, rtol=prop_min_rtol, atol=prop_min_atol)):
            F = F_new
            F_new = alpha*np.dot(F, matrix) + F_0
        else:
            break
    if verbose:
    	print 'Network propagation complete:', i+1, 'iterations'
    return F_new    	



