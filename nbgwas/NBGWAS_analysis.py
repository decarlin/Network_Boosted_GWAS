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
    else:
        sorted_table = pd.concat([values, y_actual], axis=1)
        sorted_gs_genes = sorted_table[sorted_table[y_actual.name]==1].index
        
    # Calculate precision recall curve
    TP, FN, count = 0, len(gold_standard_set), 1
    precision, recall = [1.0], [0.0]
    for gene in sorted_gs_genes:
        TP += 1.0
        FN -= 1.0
        precision.append(TP/float(sorted_table.ix[:gene].shape[0]))
        recall.append(TP/float(TP+FN))
        count += 1

    AUPRC = metrics.auc(recall, precision)
    
    return precision, recall, AUPRC


# Can be updated to compare against degree-matched null gene sets
def ranksum(ranked_gene_list, gene_list):
    # Calculate ranksum of actual top genes 
    GS_prop_results = ranked_gene_list[ranked_gene_list['Gene'].isin(gene_list)]
    GS_ranksum = GS_prop_results['Rank'].sum()

    # Generate null ranksums (creates normal null distribution)
    n = GS_prop_results.shape[0]
    genes = list(ranked_gene_list['Gene'])
    null_ranksums = []
    for i in range(1000):
        null_ranksums.append(ranked_gene_list[ranked_gene_list['Gene'].isin(np.random.choice(genes, n, replace=False))]['Rank'].sum())

    # Calculate Z-score and p-value of GS_ranksum
    Z = (GS_ranksum - np.mean(null_ranksums)) / np.std(null_ranksums)
    p = stats.norm.sf(abs(Z))
    print 'Z-score and p-value of TopProp Only Rank Sum Test:'
    print 'Z =', Z
    print 'p =', p

    return Z, p