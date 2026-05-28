import os
import click
import logging
from rich.logging import RichHandler
import string
import itertools
import pandas as pd
import numpy as np
import networkx as nx
import scikit_posthocs as spo
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from collections import defaultdict
from Bio import AlignIO
from Bio import Phylo
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from scipy.stats import kruskal, mannwhitneyu
from statsmodels.stats.multitest import multipletests

def calculadjustp(data):
    outgroup_data = data[0]
    stats_ps = [mannwhitneyu(outgroup_data, ingroup_data, alternative='two-sided') for ingroup_data in data[1:]]
    raw_p_values = [i[1] for i in stats_ps]
    reject, pvals_corrected, _, _ = multipletests(raw_p_values, alpha=0.05, method='bonferroni')
    return pvals_corrected

def assign_letters(p_matrix, group_names, df_work, col, alpha=0.05):
    G = nx.Graph()
    G.add_nodes_from(group_names)
    for i in range(len(group_names)):
        for j in range(i + 1, len(group_names)):
            g1 = group_names[i]
            g2 = group_names[j]
            if p_matrix.loc[g1, g2] >= alpha: G.add_edge(g1, g2)
    cliques = list(nx.find_cliques(G))
    alphabet = string.ascii_lowercase
    letters_dict = {name: [] for name in group_names}
    for i, clique in enumerate(cliques):
        letter = alphabet[i]
        for node in clique:
            letters_dict[node].append(letter)
    for node in letters_dict:
        letters_dict[node] = "".join(sorted(letters_dict[node]))
    medians = df_work.groupby('Root')[col].median().sort_values(ascending=False)
    return letters_dict

def addstats(ax,groups_ordered,data_group_ordered,col,df_work):
    outgroup_data = data_group_ordered[0]
    stats_ps = [mannwhitneyu(outgroup_data, ingroup_data, alternative='two-sided') for ingroup_data in data_group_ordered[1:]]
    raw_p_values = [i[1] for i in stats_ps]
    reject, pvals_corrected, _, _ = multipletests(raw_p_values, alpha=0.05, method='bonferroni')
    group_names = df_work["Group"].unique()
    p_matrix = spo.posthoc_dunn(df_work, val_col=col, group_col='Group', p_adjust='bonferroni')
    letters_dict = assign_letters(p_matrix, group_names, df_work, col)
    return letters_dict

def addletters(ax,data_group_ordered,letters_dict,groups_ordered,whiskers):
    ordered_letters = [letters_dict[name] for name in groups_ordered]
    for i in range(len(data_group_ordered)):
        top_whisker_y = whiskers[2*i + 1].get_ydata()[1]
        ax.text(i+1, top_whisker_y, ordered_letters[i], ha='center', va='bottom', fontsize=10, fontweight='bold')

def cldanalyze(dic,dic2,dic3,dic4,dic5,dic6):
    ordered_root = ['Bryophyte','Setaphyte','Hornwort','Moss','Liverwort']
    data1 = [dic[root] for root in ordered_root]
    data2 = [dic2[root] for root in ordered_root]
    data3 = [dic3[root] for root in ordered_root]
    data4 = [dic4[root] for root in ordered_root]
    data5 = [dic5[root] for root in ordered_root]
    data6 = [dic6[root] for root in ordered_root]
    cols = ["Alignment Length","Parsimony-informative sites","GRAVY (mean)","GRAVY (cv)","Tree clock-likeness","Substitution rate"]
    fig, axes = plt.subplots(2,3,figsize=(18, 12))
    colors = cm.viridis(np.linspace(0, 1, len(ordered_root)))
    for col,data,ax in zip(cols,(data1,data2,data3,data4,data5,data6),axes.flatten()):
        col1 = [root for root,da in zip(ordered_root,data) for _ in range(len(da))]
        col2 = list(itertools.chain.from_iterable(data))
        df = pd.DataFrame.from_dict({col:col2,"Root":col1},orient="columns")
        bplot = ax.boxplot(data,patch_artist=True,sym='')
        for patch, cr in zip(bplot['boxes'], colors):
            patch.set_facecolor(cr)
        ax.set_xticklabels(ordered_root)
        if col == "Tree clock-likeness": ax.set_ylabel("Coefficient of variance of root-to-tip distance",fontsize=10)
        elif col == "GRAVY (mean)": ax.set_ylabel("Mean",fontsize=10)
        elif col == "GRAVY (cv)": ax.set_ylabel("Coefficient of variance",fontsize=10)
        elif col == "Alignment Length": ax.set_ylabel("Length",fontsize=10)
        elif col == "Parsimony-informative sites": ax.set_ylabel("Number of sites",fontsize=10)
        elif col == "Substitution rate": ax.set_ylabel("Number of substitutions per site per branch",fontsize=10)
        else: ax.set_ylabel(col,fontsize=10)
        ax.set_title(col,fontsize=10,fontweight='bold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        stat, p_value = kruskal(*data)
        if p_value < 0.05:
            adjust_ps = calculadjustp(data)
            letters_dict = addstats(ax,ordered_root,data,col,df.copy().sort_values(by=["Root"])[[col,"Root"]])
            addletters(ax,data,letters_dict,ordered_root,bplot['whiskers'])
        ax.text(0.01,0.97,"$P$ = {:.3f}; Kruskal-Wallis test".format(p_value),transform=ax.transAxes)
    return fig,ax

def getparsimony(fn):
    with open(fn,"r") as f:
        for line in f:
            if " parsimony-informative" in line:
                num = int(line.split(" parsimony-informative")[0])
                return num

def calculaa(seq):
    seq_clean = "".join([s for s in seq if s != "-" and s!= "?" and s!= "*" and s!= "X"])
    analysis = ProteinAnalysis(seq_clean)
    return analysis

def calculaacv(msa):
    seqs = [str(sq.seq) for sq in msa]
    gravy = []
    for seq in seqs:
        seq_clean = "".join([s for s in seq if s != "-" and s!= "?" and s!= "*" and s!= "X"])
        analysis = ProteinAnalysis(seq_clean)
        gravy.append(analysis.gravy())
    u = np.mean(gravy)
    sigma = np.std(gravy)
    cv = sigma/abs(u) if u != 0 else float('inf')
    return cv

def getcvroot(fn):
    tree = Phylo.read(fn, "newick")
    Depths = tree.depths()
    distances = [Depths[tip] for tip in tree.get_terminals()]
    mean_dist = np.mean(distances)
    std_dist = np.std(distances)
    cv = (std_dist / mean_dist)
    return cv

def getrateroot(fn):
    tree = Phylo.read(fn, "newick")
    Depths = tree.depths()
    distances = [Depths[tip] for tip in tree.get_terminals()]
    mean_dist = np.mean(distances)
    return mean_dist

@click.group(context_settings={'help_option_names': ['-h', '--help']})
@click.option('--verbosity', '-v', type=click.Choice(['info', 'debug']), default='info', help="Verbosity level, default = info.")
def cli(verbosity):
    """
    :This script is to conduct CLD analysis
    :Usage python cld.py cldanalysis gfpathlist -p msa_path
    """
    logging.basicConfig(
        format='%(message)s',
        handlers=[RichHandler()],
        datefmt='%H:%M:%S',
        level=verbosity.upper())
    logging.info("Proper Initiation")
    pass

@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.argument('gfpathlist', type=click.Path(exists=True))
@click.option('--msa_path', '-p', default=None, show_default=True, help="path to msa")
def cldanalysis(gfpathlist,msa_path):
    """
    :Post-hoc Multiple Comparison Test with Compact Letter Display (CLD)
    :gfpathlist is a file indicating the path to the list of gene families supporting different scenarios
    :msa_path is the path containing msa, iqtree log, and tree files for different gene families
    """
    dic,dic2,dic3,dic4,dic5,dic6 = [defaultdict(list) for _ in range(6)]
    with open(gfpathlist,"r") as f:
        for line in f:
            root = os.path.basename(line.strip())
            with open(line.strip(),"r") as ff:
                for gf in ff:
                    msafn = os.path.join(msa_path,gf.strip())
                    logfn = os.path.join(msa_path,gf.strip())
                    treefn = os.path.join(msa_path,gf.strip())
                    msa = AlignIO.read(msafn,"fasta")
                    dic[root] += [msa.get_alignment_length()]
                    dic2[root] += [getparsimony(logfn)]
                    dic3[root] += [calculaa(msa)]
                    dic4[root] += [calculaacv(msa)]
                    dic5[root] += [getcvroot(treefn)]
                    dic6[root] += [getrateroot(treefn)]
    fig,_ = cldanalyze(dic,dic2,dic3,dic4,dic5,dic6)
    fig.tight_layout()
    fig.savefig("cld.pdf")

if __name__ == '__main__':
	cli()
