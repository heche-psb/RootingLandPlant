import pandas as pd
import click
import logging
import os
from rich.logging import RichHandler

@click.group(context_settings={'help_option_names': ['-h', '--help']})
@click.option('--verbosity', '-v', type=click.Choice(['info', 'debug']), default='info', help="Verbosity level, default = info.")
def cli(verbosity):
    """
    :This script is to get rSOGs from multi-copy orthogroups
    :Usage python rSOG.py multi2singlecopy orthotsv dmdresultpath
    """
    logging.basicConfig(
        format='%(message)s',
        handlers=[RichHandler()],
        datefmt='%H:%M:%S',
        level=verbosity.upper())
    logging.info("Proper Initiation")
    pass

@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.argument('orthotsv', type=click.Path(exists=True))
@click.argument('dmdresultpath', type=click.Path(exists=True))
def multi2singlecopy(orthotsv,dmdresultpath):
    """
    :Assuming your have already finished diamond analysis and all the diamond results are in the directory "dmdresultpath"
    :Make sure that your diamond result files are named such as OG0004110.fa.dmd.tsv, OG0001064.fa.dmd.tsv
    :"orthotsv" is the Orthogroups.tsv file from OrthoFinder
    """
    df = pd.read_csv(orthotsv,header=0,index_col=0,sep='\t')
    df_new = df.copy()
    for gf in df.index:
        dmdfilepath = os.path.join(dmdresultpath,gf+".fa.dmd.tsv")
        df_dmd = pd.read_csv(dmdfilepath,header=None,index_col=None,sep='\t')
        bitscore_name = df_dmd.columns[-1]
        df_dmd = df_dmd.rename(columns={bitscore_name:"Bit-Score"})
        df_dmd["Pair"] = ["__".join(sorted([i,j])) for i,j in zip(df_dmd[0],df_dmd[1])]
        pair_scores = {i:j for i,j in zip(df_dmd["Pair"],df_dmd["Bit-Score"])}
        for sp in df.columns:
            multi_genes = df.loc[gf,sp].split(", ")
            multi_genes_scores = {gene:0 for gene in multi_genes}
            other_sps = [left_sp for left_sp in df.columns if left_sp!=sp]
            other_gene_pools = []
            for left_sp in other_sps:
                other_gene_pools += df.loc[gf,left_sp].split(", ")
            for gene in multi_genes:
                gene = gene
                gene_pairs = ["__".join(sorted([gene,other_gene])) for other_gene in other_gene_pools]
                for gene_pair in gene_pairs:
                    multi_genes_scores[gene]+= pair_scores.get(gene_pair,0)
            best_gene,best_score = sorted(multi_genes_scores.items(),key=lambda x:x[1],reverse=True)[0]
            df_new.loc[gf,sp] = best_gene
    df_new.to_csv("rSOGs.tsv",header=True,index=True,sep='\t')

if __name__ == '__main__':
	cli()
