import matplotlib
matplotlib.use('agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('white')

import numpy as np
import pandas as pd
import sklearn.metrics as metrics

# Construct dictionary of gold standard gene sets from input text file to perform NBGWAS results analysis
# File format: Each line is a delimited list with the first item in the list is the name of the gene set
# All other genes in the list follow the gene set name
def load_GS_sets(gs_sets_file, delimiter='\t', verbose=False):
    f = open(gs_sets_file)
    gs_set_lines = f.read().splitlines()
    gs_set_lines_split = [line.split(delimiter) for line in gs_set_lines]
    f.close()
    gs_sets = {gs_set[0]:set(gs_set[1:]) for gs_set in gs_set_lines_split}
    if verbose:
        print len(gs_sets), 'gold standard gene sets loaded from:', gs_sets_file
    return gs_sets

# Construct precision recall curve
def PRC(values, gold_standard_set, presorted=True, sort_value_ascending=True):
    # Set up sorted values table
    y_actual = pd.Series(0, index=values.index, dtype=int)
    y_actual.ix[gold_standard_set]+=1
    y_actual.name = 'Gold Standard Gene'
    if not presorted:
        if sort_value_ascending:
            sorted_table = pd.concat([values, y_actual], axis=1).sort_values(by=[values.name, y_actual.name], ascending=[True, True])
        else:
            sorted_table = pd.concat([values, y_actual], axis=1).sort_values(by=[values.name, y_actual.name], ascending=[False, True])
        sorted_gs_genes = sorted_table[sorted_table[y_actual.name]==1].index

    # Calculate precision recall curve
    TP, FN, count = 0, len(gold_standard_set), 1
    precision, recall = [1.0], [0.0]
    for node in sorted_gs_genes:
        TP += 1.0
        FN -= 1.0
        precision.append(TP/float(sorted_table.ix[:node].shape[0]))
        recall.append(TP/float(TP+FN))
        count += 1

    AUPRC = metrics.auc(recall, precision)
    
    return precision, recall, AUPRC

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