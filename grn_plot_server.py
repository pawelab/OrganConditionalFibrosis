"""
Sets up a server that takes in nodes from the grn in `visualize_grn.js`, plots their behavior, and places figures in a temp directory.
"""

from flask import Flask, request, send_from_directory
from tempfile import TemporaryDirectory
import os
from pathlib import Path
from plot_utils import plot_expression, plot_peak
import argparse
import anndata as ad

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument("--rna_adata", help="Path to the bulk RNA AnnData file", type=Path, default='./data/rna.h5ad')
  parser.add_argument("--atac_adata", help="Path to the bulk ATAC AnnData file", type=Path, default='./data/atac.h5ad')

  args = parser.parse_args()

  rna_ad = ad.read_h5ad(args.rna_adata)
  atac_ad = ad.read_h5ad(args.atac_adata)

  with TemporaryDirectory(dir='./') as tmpdir:
    app = Flask(__name__, static_url_path='', static_folder='./')
    app.config['tmpdir'] = tmpdir

    @app.route('/plot', methods = ['POST'])
    def plot():
      data = request.get_json()
      node_id = data['id']
      is_gene = data['node_type'] == 'gene'
      temp_file_path = os.path.join(os.path.basename(tmpdir), f"{node_id}.png")
      if not os.path.exists(temp_file_path):
        ax = plot_expression(rna_ad, node_id) if is_gene else plot_peak(atac_ad, f'Interval_{node_id}')
        ax.figure.savefig(temp_file_path)
        ax.figure.clear()
      return { 'path': temp_file_path }

    app.run()