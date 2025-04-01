"""
Script to generate PG-SCOTT scores and candidate lists for the given input data.
Individual functions can also be imported for use in other notebooks/scripts.
"""
from argparse import ArgumentParser
from pathlib import Path
import pandas as pd
import anndata as ad
import numpy as np


def gene_to_peaks_map(loc_file: Path) -> dict[str, list[str]]:
  """Map gene to peaks within 150kb"""
  gene_locations = pd.read_table(loc_file, index_col = 'Gene ID', usecols = ['Gene ID', 'Begin', 'End', 'Chromosome', 'Symbol', 'Gene Type', 'Accession'])
  gene_locations.drop(index = [idx for idx, accession in gene_locations['Accession'].items() if not accession.startswith('NC')], inplace = True)
  new_gene_locs = gene_locations.copy()
  new_gene_locs.index = new_gene_locs['Symbol']
  gene_to_idx = {
    gene: new_gene_locs.loc[gene] for gene in rna_adata.var['gene'].unique() if gene in new_gene_locs.index and len(new_gene_locs.loc[gene].shape) == 1
  }
  gene_to_peaks_map = {
    gene: atac_adata.var_names[(atac_adata.var['Chr'] == f'chr{row['Chromosome']}') & 
                              (atac_adata.var['End'] >= row['Begin'] - 1.5 * 10 ** 5) & 
                              (atac_adata.var['Start'] <= row['End'] + 1.5 * 10 ** 5)]
    for gene, row in gene_to_idx.items()
  }
  return gene_to_peaks_map


def calc_condition_gene_averages(adata: ad.AnnData, organ: str, layer: str = None) -> pd.DataFrame:
  """Average expression counts across replicates, summing across isoforms"""
  organ_adata = adata[adata.obs['organ'] == organ, :]
  organ_df = organ_adata.to_df(layer = layer)
  organ_df['condition'] = organ_adata.obs['condition']
  condition_avg_df = organ_df.groupby('condition', observed = False).aggregate('mean').T  # Average across replicates
  condition_avg_df['gene'] = organ_adata.var['gene']
  return condition_avg_df.groupby('gene', observed = False).aggregate('sum')  # Sum across gene variants


def calc_condition_peak_averages(adata: ad.AnnData, organ: str, layer: str = None) -> pd.DataFrame:
  """Average peak counts across replicates"""
  organ_adata = adata[adata.obs['organ'] == organ, :]
  organ_df = organ_adata.to_df(layer = layer)
  organ_df['condition'] = organ_adata.obs['condition']
  return organ_df.groupby('condition', observed = False).aggregate('mean').T  # Average across replicates


def calc_lfc_scoring(organ_df: pd.DataFrame) -> pd.DataFrame:
  """Uses log fold change to determine which genes are up- and down-regulated in each condition"""
  odf_copy = organ_df + 1
  lfc_df = pd.DataFrame({
    'up': np.sqrt(np.log2(odf_copy.mean(axis = 1))) * np.log2(odf_copy['TGFb'] / odf_copy['Unstim']),
    'down': np.sqrt(np.log2(odf_copy.mean(axis = 1))) * np.log2(odf_copy['TGFb+JQ1'] / odf_copy['TGFb']),
  })
  t = np.sign(np.sign(odf_copy['Unstim'] - odf_copy['TGFb']) - np.sign(odf_copy['TGFb'] - odf_copy['TGFb+JQ1']))
  lfc_df['score'] = lfc_df['up'] * lfc_df['down'] * t
  return lfc_df


def merge_gene_peak_lfc(gene_lfc_df: pd.DataFrame, peak_lfc_df: pd.DataFrame, gene_to_peaks: dict) -> pd.DataFrame:
  """Merge gene and peak scores based on 150kb range from gene coding region"""
  gene_peak_pair_df = pd.DataFrame({(gene, peak): {'gene_score': gene_lfc_df['score'][gene], 'peak_score':peak_lfc_df['score'][peak]} for gene, peaks in gene_to_peaks.items() for peak in peaks}).T
  gene_peak_pair_df.index.names = ['Gene', 'Peak']
  gene_peak_pair_df = gene_peak_pair_df.sort_values(['gene_score', 'peak_score'], ascending=False)
  return gene_peak_pair_df


def calc_condition_sensitivity(rna_adata: ad.AnnData, atac_adata: ad.AnnData, organ: str, layer: str = None) -> pd.DataFrame:
  """Calculate peak and gene condition sensitivity scores for a given organ"""
  organ_rna_df = calc_condition_gene_averages(rna_adata, organ, layer = layer)
  organ_gene_lfc_df = calc_lfc_scoring(organ_rna_df)
  organ_atac_df = calc_condition_peak_averages(atac_adata, organ, layer = layer)
  organ_peak_lfc_df = calc_lfc_scoring(organ_atac_df)
  return organ_gene_lfc_df, organ_peak_lfc_df


def get_best_pg_pairs(gene_score: pd.Series, peak_score: pd.Series, n = 10) -> list[tuple[str, str]]:
  """Select n top peak-gene combinations based on multiplicative scores"""
  reordered_peaks = peak_score.loc[gene_score.index]
  ordered_pairs = (gene_score * reordered_peaks * (np.sign(reordered_peaks) + np.sign(gene_score)) - 1).sort_values(ascending=False)
  top_genes = ordered_pairs.index.get_level_values(0).unique()[0:n]
  top_pairs = ordered_pairs.groupby(level=0).nth(0).loc[top_genes].index
  return top_pairs


def main(rna_adata: ad.AnnData, atac_adata: ad.AnnData, gene_to_peaks: dict[str, list[str]], n_candidates: int) -> list[tuple[str, pd.DataFrame]]:
  """
  Calculate PG-SCOTT scores for the given data. Returns a list of label-dataframe tuples.
  """
  outdata = []
  organs = rna_adata.obs['organ'].unique()

  # Calculate sensitivity scores per-organ
  organ_to_pg_lfc_scores = {
    organ: calc_condition_sensitivity(rna_adata, atac_adata, organ) for organ in organs
  }
  organ_to_merged_sensitivity_scores = {
    organ: merge_gene_peak_lfc(gene_lfc_df, peak_lfc_df, gene_to_peaks) for organ, (gene_lfc_df, peak_lfc_df) in organ_to_pg_lfc_scores.items()
  }

  # Merge sensitivity scores across organs, normalize, and save
  sensitivity_df = pd.concat(organ_to_merged_sensitivity_scores.values(), axis = 'columns', keys = organ_to_merged_sensitivity_scores.keys())
  sensitivity_df = sensitivity_df.fillna(0.)
  sensitivity_df = sensitivity_df / sensitivity_df.max()
  sensitivity_df = sensitivity_df.where(sensitivity_df > 0., 0.)
  outdata.append(('pg_pairs_sensitivity', sensitivity_df))

  # Calculate specificity scores per-organ
  sums = sensitivity_df.T.groupby(level = 1).agg('sum').T
  n_organs = organs.shape[0]
  specificity_df = sensitivity_df.apply(lambda col: col - (sums[col.name[1]] - col) / (n_organs - 1), axis = 0)
  specificity_df = specificity_df.sort_values(list(specificity_df.columns), ascending=False)
  outdata.append(('pg_pairs_specificity', specificity_df))
  
  # Calculate universality scores (not per-organ)
  universality_df = (1 - specificity_df.abs().T.groupby(level = 1).agg('mean').T) * (sensitivity_df.T.groupby(level = 1).agg('min').T)
  outdata.append(('pg_pairs_universality', universality_df))

  # Get top PG pairs where peak is organ-specific or gene is universal
  organ_to_best_specificity_pairs = {
    organ: get_best_pg_pairs(specificity_df[organ, 'gene_score'], specificity_df[organ, 'peak_score'], n = n_candidates) for organ in organs
  }
  organ_to_universal_gene_specific_peaks = {
    organ: get_best_pg_pairs(universality_df['gene_score'], specificity_df[organ, 'peak_score'], n = n_candidates) for organ in organs
  }
  universal_pg_pairs = get_best_pg_pairs(universality_df['gene_score'], universality_df['peak_score'], n_candidates)

  # Generate candidate list based on interesting conditions
  candidate_dict = {}
  candidate_dict.update({('organ-specific peak-gene', organ): pairs for organ, pairs in organ_to_best_specificity_pairs.items()})
  candidate_dict.update({('organ-specific peak', organ): pairs for organ, pairs in organ_to_universal_gene_specific_peaks.items()})
  candidate_dict.update({('universal peak-gene', ''): universal_pg_pairs})
  candidate_df = pd.DataFrame(candidate_dict)
  outdata.append(('candidate_pg_pairs', candidate_df))


if __name__ == '__main__':
  parser = ArgumentParser()
  parser.add_argument("--rna", help="Filepath to an AnnData object containing the bulk RNA data", type=Path, required=False, default='/projectnb/paxlab/dillon/OrganFibrosisConditionalGRNs/data/rna.h5ad')
  parser.add_argument("--atac", help="Filepath to an AnnData object containing the bulk ATAC data", type=Path, required=False, default='/projectnb/paxlab/dillon/OrganFibrosisConditionalGRNs/data/atac.h5ad')
  parser.add_argument("--gene-locs", help="Filepath to table mapping genes to chromosome locations", type=Path, required=False, default='/projectnb/paxlab/dillon/OrganFibrosisConditionalGRNs/data/gene_locations_hg38.tsv')
  parser.add_argument("--outfolder", help="Folderpath to output score and candidate csvs in", type=Path, required=False, default='/projectnb/paxlab/dillon/OrganFibrosisConditionalGRNs/data/')
  parser.add_argument("--n-candidates", help="Number of top peak-gene pairs to extract for candidate list", type=int, required=False, default=10)

  args = parser.parse_args()

  # Read in data
  rna_adata = ad.read_h5ad(args.rna)
  atac_adata = ad.read_h5ad(args.atac)
  gene_to_peaks = gene_to_peaks_map(args.gene_locs)

  outfiles = main(rna_adata, atac_adata, gene_to_peaks, args.n_candidates)

  for name, df in outfiles:
    df.to_csv(f'{args.outfolder}/{name}.csv')
