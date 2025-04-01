"""
File containing helper plotting functions for working with ATAC/RNA AnnData and PG-SCOTT scores
"""
import matplotlib.pyplot as plt
import scanpy as sc
import anndata as ad
import seaborn as sns
import pandas as pd

palette = {'Unstim': '#2ca02c', 'TGFb': '#d62728', 'TGFb+JQ1': '#1f77b4'}
markers = {'Lung': 'o', 'Heart': '^', 'Liver': 's'}

def plot_organ_pca(adata: ad.AnnData, **kwargs) -> plt.Figure:
  """
  PCA plot for all samples as well as grouped per-organ. Colored by condition.
  """
  fig, axs = plt.subplots(2, 2)
  for i, (organ, marker) in enumerate(markers.items()):
    sub = adata[adata.obs['organ'] == organ].copy()
    sc.pl.embedding(sub, 'pca', color='condition', palette=palette, s=50, title=organ, marker=marker, ax=axs[i % 2][int(i / 2)], show=False, use_raw=False, legend_loc='best' if i == 2 else None, **kwargs)
    sc.pl.embedding(sub, 'pca', color='condition', palette=palette, s=50, title='All', marker=marker, ax=axs[1][1], legend_loc=None, use_raw=False, show=False, **kwargs)

  plt.tight_layout()
  return fig


def plot_expression(adata: ad.AnnData, gene: str, layer: str = None, points = False, **kwargs) -> plt.Axes:
  """
  Plot expression of the given gene as a boxtplot across conditions and organs.
  """
  gene_rna_adata = adata[:, adata.var['gene'] == gene]
  ax = sns.boxplot(gene_rna_adata.obs, x = 'organ', y = gene_rna_adata.to_df(layer = layer).sum(axis = 1), order=markers.keys(), hue = 'condition', palette=palette, hue_order = palette.keys(), fliersize = 0 if points else None, **kwargs)
  if 'ax' not in kwargs:
    kwargs['ax'] = ax
  if points:
    sns.stripplot(gene_rna_adata.obs, x = 'organ', y = gene_rna_adata.to_df(layer = None).sum(axis = 1), order=markers.keys(), hue = 'condition', palette=palette, edgecolor='black', linewidth=1, hue_order = palette.keys(), dodge=True, **kwargs)
  ax.set_ylabel(f'{gene} expression {layer if layer is not None else ""}')
  return ax


def plot_peak(adata: ad.AnnData, peak: str, layer: str = None, points = False, **kwargs) -> plt.Axes:
  """
  Plot accessibility of the given peak as a boxtplot across conditions and organs.
  """
  ax = sns.boxplot(adata.obs, x = 'organ', y = adata[:, str(peak)].to_df(layer = layer).sum(axis = 1), order=markers.keys(), hue = 'condition', palette=palette, hue_order = palette.keys(), fliersize = 0 if points else None, **kwargs)
  if 'ax' not in kwargs:
    kwargs['ax'] = ax
  if points:
    sns.stripplot(adata.obs, x = 'organ', y = adata[:, str(peak)].to_df(layer = None).sum(axis = 1), order=markers.keys(), hue = 'condition', palette=palette, edgecolor='black', linewidth=1, hue_order = palette.keys(), dodge=True, **kwargs)
  ax.set_ylabel(f'{peak} peaks {layer if layer is not None else ""}')
  return ax

# TODO: figure out peak range plot which requires gene_to_peaks map

def plot_pg_pairs(gene_score: pd.Series, peak_score: pd.Series, best_peaks: list[tuple[str, str]] = [], **kwargs) -> plt.Axes: 
  """
  Plot peak-gene pairs based on the provided scores. Optionally hightlight specific pairs.
  """
  plot_args = {
    **kwargs
  }
  if len(best_peaks) > 0:
    plot_args['hue'] = gene_score.index.isin(best_peaks)
  ax = sns.scatterplot(x = gene_score, y = peak_score.loc[gene_score.index], legend = None, **plot_args)
  for (gene, peak) in best_peaks:
    point = gene_score.loc[(gene, peak)], peak_score.loc[(gene, peak)]
    ax.text(point[0], point[1], f'{gene}-{peak}', size='x-small', va='bottom')
  return ax

