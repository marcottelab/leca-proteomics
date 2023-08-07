#!/usr/bin/env python3
#Date: Jan 19, 2023 #Edit July 2, 2023
#Authors: Muyoung Lee & Rachael Cox
#Description:
# 1) Propagate scores based on the input network and nodes.
# 2) Rank nodes based on the new scores.
# 3) Draw ROC curves and calculate AUROC.
# 4) Compare to random.
# 5) Provide list of input IDs and top hits.
# Usage: python3 [THIS SCRIPT] --ppi_network [NETWORK FILE] --annotations [ANNOTATION FILE] --outfile_name [OUTFILE]

import random
import numpy as np
import pandas as pd
import os
import re
import argparse
import matplotlib.pyplot as plt

''' Functions '''
def get_ntwk_scores(network_file):
    network = {}
    node_scores = {}
    with open(network_file) as NETWORK:
        NETWORK.readline()
        for line in NETWORK:
            edge, edge_score = line.strip().split(',')
            node1, node2 = edge.split(' ')	
            node1, node2 = sorted([node1, node2])
            network[(node1, node2)] = float(edge_score)
            if node1 not in node_scores:
                node_scores[node1] = 0.0
            if node2 not in node_scores:
                node_scores[node2] = 0.0
    #print("How many nodes?", len(node_scores), sep="\t")
    #print("How many edges?", len(network), sep="\t")
    return network, node_scores

def propagation(network, node_scores, node_set):
	new_scores = node_scores.copy()
	for edge in network:
		node1, node2 = edge
		weight = network[edge]
		if node1 in node_set:
			new_scores[node2] += weight
		if node2 in node_set:
			new_scores[node1] += weight
	return new_scores

def get_ranks(node_scores, ascend):
    series = pd.Series(node_scores.values()) 
    series.index = node_scores.keys()
    # default .rank() params = highest scores will have the highest rank
    # ties = average rank across group (e.g., a tie between rank 2 and rank 3 = both are assigned 2.5)
    ranks = series.rank(ascending=ascend)
    return ranks

def ranks_from_leave_one_out(network, node_scores, node_set):
	new_node_scores = propagation(network, node_scores, node_set)
	ranks = get_ranks(new_node_scores, False)
	return ranks.sort_values(ascending=True)

def cal_TPR_FPR_AUC(nodes_of_interest, network, ranks_of_nodes):
	values = []
	num_total_nodes = len(ranks_of_nodes)
	dTPRy = 1 / len(nodes_of_interest)
	dFPRx = 1 / (num_total_nodes - len(nodes_of_interest))
	x, y, AUC = 0, 0, 0
	for node in ranks_of_nodes.index:
		dAUC = 0
		if node in nodes_of_interest:
			y += 1
		else:
			x += 1
			dAUC = dFPRx * y * dTPRy
		AUC += dAUC
		values.append((x * dFPRx, y * dTPRy, AUC))
	remainder_AUC = (1 - x * dFPRx) * (y * dTPRy + 1) * 0.5
	total_AUC = AUC + remainder_AUC
	values.append((1.0, 1.0, total_AUC))
	return values

def cal_TPR_FPR_AUC_for_RANDOM(nodes_of_interest, network, ranks_of_nodes):
	num_total_nodes = len(ranks_of_nodes)
	dTPRy = 1 / len(nodes_of_interest)
	dFPRx = 1 / (num_total_nodes - len(nodes_of_interest))
	x, y, AUC = 0, 0, 0
	for node in ranks_of_nodes.index:
		dAUC = 0
		if node in nodes_of_interest:
			y += 1
		else:
			x += 1
			dAUC = dFPRx * y * dTPRy
		AUC += dAUC
	remainder_AUC = (1 - x * dFPRx) * (y * dTPRy + 1) * 0.5
	total_AUC = AUC + remainder_AUC
	return total_AUC

def draw_curve_without_types(coordinates_for_graph):
	x = []
	y = []
	for value in coordinates_for_graph:
		x.append(value[0])
		y.append(value[1])
	total_AUC = coordinates_for_graph[-1][-1]
	if total_AUC > 0.5:
		color = "#072A6C"
	else:
		color = "#CC5500"
	plt.plot(x, y, color=color, linewidth=1, alpha=0.6)
	plt.xlim([0, 1])
	plt.ylim([0, 1])


def get_disease_scores(disease_file, node_scores):
    node_dict = {}
    disease_dict = {}
    with open(disease_file) as NODE_SETS:
        for line in NODE_SETS:
            words = line.strip("\n").split("\t")
            diseaseID = words[0]
            node_set = set(words[2].split(",")) & set(node_scores)
            if len(node_set) >= 5:
                node_dict[diseaseID] = node_set
                disease_dict[diseaseID] = words[1]
    return node_dict, disease_dict

def convert_ids(mapping_file, disease_node_dict, disease_auc_dict):
    egg_uniacc = {}
    egg_symbol = {}
    with open(mapping_file) as CONVERT_TABLE:
        for line in CONVERT_TABLE:
            eggnog, uniacc, symbol = line.strip().split("\t")
            egg_uniacc[eggnog] = uniacc
            egg_symbol[eggnog] = symbol
    egg_col = []
    gene_col = []
    up_col = []
    for diseaseID in disease_auc_dict:
        human_uniacc_list = []
        human_symbol_list = []
        egg_list = sorted(list(disease_node_dict[diseaseID]))
        for eggnog in egg_list:
            if eggnog in egg_uniacc:
                human_uniacc_list.append(egg_uniacc[eggnog])
                human_symbol_list.append(egg_symbol[eggnog])
            else:
                human_uniacc_list.append("")
                human_symbol_list.append("")
        egg_col.append(";".join(egg_list))
        gene_col.append(";".join(human_symbol_list))
        up_col.append(";".join(human_uniacc_list))
    return egg_col, gene_col, up_col

def format_out(outfile_name=None):
    # format outfile paths/names
    if outfile_name:
        path, filename = os.path.split(os.path.realpath(outfile_name))
        outpath = path+'/'
        if '.' in filename:
            filename = re.search('.*(?=\.)', filename)[0]
        file_out = outpath+filename
    else:
        path = os.getcwd()
        outpath = path+'/'
        file_out = outpath+'disease_propagation'
    return(file_out)
    
''' Main '''
def main():
    
    # get file name for results
    out_prefix = format_out(args.outfile_name)
    print(f'Results will be saved to: {out_prefix}')
    
    # network and node_scores
    network, node_scores = get_ntwk_scores(args.ppi_network)

    # disease node sets
    disease_nodes, disease_ids = get_disease_scores(args.disease_network, node_scores)
    
    # calc auc, roc plot
    diseaseID_AUC = {}
    for diseaseID in disease_nodes:
        ranks = ranks_from_leave_one_out(network, node_scores, disease_nodes[diseaseID])
        # TODO: write this to file to plot nicely in R later
        values = cal_TPR_FPR_AUC(disease_nodes[diseaseID], network, ranks)
        draw_curve_without_types(values)
        diseaseID_AUC[diseaseID] = values[-1][-1]

    plt.plot([0,1], [0,1], color="gray", linewidth=1, linestyle="dashed")
    plt.savefig(f"{out_prefix}_ROC.png", dpi=300)
    plt.clf()

    sorted_diseaseID_AUC = dict(sorted(diseaseID_AUC.items(), key=lambda item: item[1], reverse=True))
    
    egg_ids, gene_ids, up_ids = convert_ids(args.annotations, disease_nodes, sorted_diseaseID_AUC)
    
    cols = ['disease_id', 'disease', 'disease_size', 'leave1out_auc', 'random_auc', 'top_hits']
    rows = []
    
    # format df output
    for diseaseID in sorted_diseaseID_AUC:
        
        # for random AUC calc
        nodes = random.sample(list(node_scores.keys()), len(disease_nodes[diseaseID]))
        random_ranks = ranks_from_leave_one_out(network, node_scores, nodes)
        random_auc = cal_TPR_FPR_AUC_for_RANDOM(nodes, network, random_ranks)
        
        # actual ranks
        ranks = ranks_from_leave_one_out(network, node_scores, disease_nodes[diseaseID])
        
        egg_list = set(disease_nodes[diseaseID])
        hits = []
        rank_list = []
        
        for n, i in enumerate(ranks.index):
            rank_of_i = ranks[n]
            if set(i) not in egg_list:
                print(i, '\t', rank_of_i, '\t', disease_ids[diseaseID])
                hits.append(i)
                rank_list.append(rank_of_i)
            if args.num_hits:
                if len(hits) == args.num_hits:
                    break
                    
        if args.one_file_per_disease:
            df = ranks.to_frame().reset_index()
            df.columns = ['ID', 'rank']
            df['disease'] = disease_ids[diseaseID]
            df['status'] = ['known' if i in egg_list else 'candidate' for i in df['ID']]
            dis_fmt = disease_ids[diseaseID].replace(' ', '_')
            df.to_csv(f'{out_prefix}_{dis_fmt}_ranked_ids.csv', index=False)
            print(df)
                
        rows.append([diseaseID, disease_ids[diseaseID], len(disease_nodes[diseaseID]), sorted_diseaseID_AUC[diseaseID], random_auc, ";".join(hits)])
        
    df = pd.DataFrame(rows, columns=cols)
    df['input_nogs'] = egg_ids
    df['gene_names'] = gene_ids
    df['uniprot_ids'] = up_ids
    print(df)
    df.to_csv(f'{out_prefix}.csv', index=False)
        
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()

    # specify scored PPI file
    parser.add_argument("-i", "--ppi_network", help="(Required) Path to scored interactions file with a column called 'ID' containing space-separated PPIs and a column called 'ppi_score' containing PPI scores in descending order (i.e., as output by ppi_predict.py).")

    # specify disease file
    parser.add_argument("-d", "--disease_network", action="store", help="(Required) Path to tab-separated file with columns for disease ID, disease name, and a comma-separated list of associated eggNOG groups.")
    
    # specify number of candidates to return (optional)
    parser.add_argument("--num_hits", action="store", type=int, default=None, help="(Optional) Specify number of top candidates to return for each disease (default = all, ranked descending).")
    
    # specify ID annotation file (optional)
    parser.add_argument("-a", "--annotations", action="store", default=None, help="(Optional) Specify a tab-separated file containing a table that maps eggNOG groups to human gene names.")
    
    # specify outfile name (optional)
    parser.add_argument("-o", "--outfile_name", action="store", default=None, help="(Optional) Specify the outfile path/name. Default is 'disease_propagation.csv' in the same directory as the input PPI network file.")
    
    # specify if you want a file with a ranked list for each disease (optional)
    parser.add_argument("--one_file_per_disease", action="store_true", default=False, help="(Optional) Specify if you want a file with a ranked list for each disease.")
    
    args = parser.parse_args()
    main()


    
