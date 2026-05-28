import os
import click
import logging
from rich.logging import RichHandler
import statsmodels.api as sm
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

@click.group(context_settings={'help_option_names': ['-h', '--help']})
@click.option('--verbosity', '-v', type=click.Choice(['info', 'debug']), default='info', help="Verbosity level, default = info.")
def cli(verbosity):
    """
    :This script is to conduct logistic regression analysis
    :Usage python logistic_regression.py logisticregression
    """
    logging.basicConfig(
        format='%(message)s',
        handlers=[RichHandler()],
        datefmt='%H:%M:%S',
        level=verbosity.upper())
    logging.info("Proper Initiation")
    pass

def plot_logistic(df, model, odds_ratio, ci, p_value):
    score_range = np.linspace(df["score"].min() - 0.5, df["score"].max() + 0.5, 300)
    X_pred = sm.add_constant(pd.DataFrame({"score": score_range}))
    probs = model.predict(X_pred)
    pred = model.get_prediction(X_pred)
    ci_frame = pred.summary_frame(alpha=0.05)
    if "mean_ci_lower" in ci_frame.columns: ci_lo, ci_hi = ci_frame["mean_ci_lower"], ci_frame["mean_ci_upper"]
    else: ci_lo, ci_hi = ci_frame["ci_lower"], ci_frame["ci_upper"]
    obs = df.groupby("score")["group"].mean().reset_index()
    obs.columns = ["score", "obs_prob"]
    counts = df.groupby("score")["group"].count().reset_index()
    counts.columns = ["score", "n"]
    p_label = "$P$ < 0.001" if p_value < 0.001 else f"$P$ = {p_value:.3f}"
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(score_range, probs, color="blue", linewidth=2.5, label="Fitted logistic curve")
    ax.fill_between(score_range, ci_lo, ci_hi, alpha=0.5, color="gray", label="95% confidence band")
    ax.scatter(obs["score"], obs["obs_prob"], color="k", s=80, zorder=5, label="Observed proportion")
    for _, row in counts.iterrows(): ax.text(row["score"], -0.06, f'$N$ = {int(row["n"])}', ha="center", va="top", fontsize=11, color="gray")
    ax.axhline(0.5, color="gray", linestyle="--", linewidth=1)
    or_lo, or_hi = np.exp(ci[0]), np.exp(ci[1])
    llr_pvalue = model.llr_pvalue
    pseudo_r2  = model.prsquared
    ax.plot([],[],label=f"LLR $P$ = {llr_pvalue:.3e}, Pseudo $R²$ = {pseudo_r2:.3f}",alpha=0)
    ax.plot([],[],label=f"OR = {odds_ratio:.3f}, 95% CI [{or_lo:.3f}, {or_hi:.3f}], {p_label}",alpha=0)
    ax.set_xlabel("Number of retained Zygnematophyceae representative species", fontsize=11)
    ax.set_ylabel("Probability of recovering bryophyte monophyly", fontsize=11)
    ax.set_title("Logistic regression", fontsize=11)
    ax.set_xticks(sorted(df["score"].unique()))
    ax.set_ylim(-0.1, 1.1)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.legend(fontsize=10, loc=0,frameon=False)
    plt.tight_layout()
    plt.savefig("logistic_regression.pdf")

@cli.command(context_settings={'help_option_names': ['-h', '--help']})
def logisticregression():
    """
    :Logistic regression analysis
    :Data from the paper
    """
    N_out = [1]*7 + [2]*21 + [3]*35 + [4]*35 + [5]*21 + [6]*7
    group = [1]*2 + [0]*5 + [1]*7 + [0]*14 + [1]*17 + [0]*18 + [1]*22 + [0]*13 + [1]*18 + [0]*3 + [1]*7
    df = pd.DataFrame({"score": N_out,"group": group})
    X = df[["score"]]
    y = df["group"]
    X_const = sm.add_constant(X)
    model = sm.Logit(y, X_const).fit()
    coef = model.params["score"]
    p_value = model.pvalues["score"]
    ci = model.conf_int().loc["score"]
    odds_ratio = np.exp(coef)
    plot_logistic(df, model, odds_ratio, ci, p_value)

if __name__ == '__main__':
    cli()
