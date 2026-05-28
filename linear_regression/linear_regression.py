import os
import click
import logging
from rich.logging import RichHandler
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
from Bio import Phylo
from scipy import stats

def univariate_linear_regression(x,y,titles=None,vname=None,fig_=None,axs_=None,responsvari=None,panellabel=None):
    x_y = [(i,j) for i,j in zip(x,y)]
    x_y_sorted = sorted(x_y, key=lambda x:x[0])
    x,y = np.array([i[0] for i in x_y_sorted]),np.array([i[1] for i in x_y_sorted])
    slope, intercept, r, p, se = stats.linregress(x, y)
    y_fit = intercept + slope * x
    norder = titles.index(vname)
    ax = axs_.flatten()[norder]
    ax.plot(x, y, 'ok', alpha = 0.3,markersize=3)
    ax.plot(x, intercept + slope*x, 'b', label='Linear fit')
    se = np.sqrt(np.sum((y - y_fit)**2) / (len(y) - 2))
    t_value = stats.t.ppf(1 - 0.025, df=len(x) - 2)
    ci = t_value * se * np.sqrt(1/len(x) + (x - np.mean(x))**2 / np.sum((x - np.mean(x))**2))
    ax.fill_between(x, y_fit-ci, y_fit+ci, color='gray', alpha=0.5, label='95% confidence band',zorder=3)
    ax.set_ylabel(responsvari,fontsize = 20)
    ax.set_xlabel(vname,fontsize = 20)
    ax.plot([], [], ' ', label="R-squared: {:4.4f}".format(r**2))
    ax.plot([], [], ' ', label="$P$: {:4.4f}".format(p))
    ax.plot([], [], ' ', label="Slope: {:4.4f}".format(slope))
    ax.legend(fontsize=15,frameon=False)
    ax.set_title("{} against {}".format(responsvari,vname),fontdict={'fontsize':20})
    if panellabel is not None: ax.text(-0.115, 1.05, "{})".format(panellabel), transform=ax.transAxes, fontsize=20, weight='bold')
    if norder == len(titles) - 1:
        fig_.tight_layout()
        fig_.savefig("linear_regression.pdf")
        plt.close()

def covariate_linear_regression(x,y,z):
    data = pd.DataFrame({"y":y,"x":x,"z":z})
    model = smf.ols(formula='y ~ x + z', data=data).fit()
    model.summary()

def getr2tdis(Tree):
    Root = Tree.root
    r2t_dis = []
    for tip in Tree.get_terminals():
        r2t_dis.append(Root.distance(tip))
    r2t_dis_narray = np.array(r2t_dis)
    cv = np.std(r2t_dis_narray)/np.mean(r2t_dis_narray)
    return cv

def getroot2tipdistancefrompath(treepath):
    cvs = []
    fns = []
    with open(treepath,"r") as f:
        for li in f:
            Tree = Phylo.read(li.strip(),'newick')
            cv = getr2tdis(Tree)
            cvs += [cv]
            fns += [os.path.basename(li.strip())]
    df = pd.DataFrame.from_dict({"MSA":fns,"CV":cvs},orient="columns")
    return df

@click.group(context_settings={'help_option_names': ['-h', '--help']})
@click.option('--verbosity', '-v', type=click.Choice(['info', 'debug']), default='info', help="Verbosity level, default = info.")
def cli(verbosity):
    """
    :This script is to conduct linear regression analysis
    :Usage python linear_regression.py linearregression ritsv entropytsv divscoretsv treepath sntsv quartettsv
    """
    logging.basicConfig(
        format='%(message)s',
        handlers=[RichHandler()],
        datefmt='%H:%M:%S',
        level=verbosity.upper())
    logging.info("Proper Initiation")
    pass

@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.argument('ritsv', type=click.Path(exists=True))
@click.argument('entropytsv', type=click.Path(exists=True))
@click.argument('divscoretsv', type=click.Path(exists=True))
@click.argument('treepath', type=click.Path(exists=True))
@click.argument('sntsv', type=click.Path(exists=True))
@click.argument('quartettsv', type=click.Path(exists=True))
def linearregression(ritsv,entropytsv,divscoretsv,treepath,sntsv,quartettsv):
    """
    :Linear regression analysis
    :ritsv is the file with Retention Index (RI) information
    :entropytsv is the file with Shannon Entropy information
    :divscoretsv is the file with Div score information
    :treepath is the file indicating the path to tree files
    :sntsv is the file with sCF, sDF, and sN information
    :quartettsv is the file with q1, q2, q3, and LPP information
    """
    df_ri = pd.read_csv(ritsv,header=0,index_col=None,sep="\t")
    df_entropy = pd.read_csv(entropytsv,header=None,index_col=None,sep="\t")
    df_div = pd.read_csv(divscoretsv,header=None,index_col=None,sep="\t")
    df_ccr2tdis = getroot2tipdistancefrompath(treepath)
    df_sn = pd.read_csv(sntsv,header=0,index_col=None,sep=r"\s+")
    df_q = pd.read_csv(quartettsv,header=0,index_col=None,sep=r"\s+")
    ri_array = df_ri["RI"].to_numpy()
    entropy_array = df_entropy.iloc[:,1].to_numpy()
    div_array = df_div.iloc[:,1].to_numpy()
    cv_array = df_ccr2tdis["CV"].to_numpy()
    sN_array = df_sn["sN"].to_numpy()
    sCF_sDF_ratio_array = df_sn["sCF"].to_numpy() / (100 - df_sn["sCF"].to_numpy())
    lpp = df_q["lpp"].to_numpy()
    q1_q23_ratio_array = df_q["q1"].to_numpy() / (1 - df_q["q1"].to_numpy())
    y_array = sCF_sDF_ratio_array
    predictive_variable_names = ["RI","Shannon Entropy","Div score","CV","q1/(q2+q3)","LPP"]
    predictive_variable_arrays = [ri_array,entropy_array,div_array,cv_array,q1_q23_ratio_array,lpp]
    ncol = int(np.ceil(len(predictive_variable_names)/2))
    fig, axs = plt.subplots(2,ncol,figsize=(6*ncol, 12))
    for x_array,vname in zip(predictive_variable_arrays,predictive_variable_names):
        univariate_linear_regression(x_array,y_array,titles=predictive_variable_names,vname=vname,fig_=fig,axs_=axs,responsvari="sCF/sDF")
        covariate_linear_regression(x_array,y_array,div_array)

if __name__ == '__main__':
    cli()
