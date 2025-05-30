"""
Sets up a server that takes in nodes from the grn in `visualize_grn.js`, plots their behavior, and places figures in a temp directory.
"""

from flask import Flask, request, render_template
from tempfile import TemporaryDirectory
import os
from dotenv import load_dotenv
from plot_utils import plot_expression, plot_peak
import anndata as ad
import requests
import shutil

app = Flask(__name__, template_folder='./', static_url_path='', static_folder='./')

load_dotenv()

if os.environ['ENV_TYPE'] == 'production':
  with requests.get(os.environ['RNA_ADATA'], stream=True) as response:
    with open('rna_adata.h5ad', 'wb') as out_file:
      shutil.copyfileobj(response.raw, out_file)

  rna_ad = ad.read_h5ad('rna_adata.h5ad')

  with requests.get(os.environ['ATAC_ADATA'], stream=True) as response:
    with open('atac_adata.h5ad', 'wb') as out_file:
      shutil.copyfileobj(response.raw, out_file)

  atac_ad = ad.read_h5ad('atac_adata.h5ad')
else:
  rna_ad = ad.read_h5ad(os.environ['RNA_ADATA'])
  atac_ad = ad.read_h5ad(os.environ['ATAC_ADATA'])

grn_path = os.environ['GRN_DATA']

@app.route('/')
def grn():
  return render_template('grn.html', grn_data = grn_path)

@app.route('/plot', methods = ['POST'])
def plot():
  data = request.get_json()
  node_id = data['id']
  is_gene = data['node_type'] == 'gene'
  temp_file_path = os.path.join(app.config['tmpdir'], f"{node_id}.png")
  if not os.path.exists(temp_file_path):
    ax = plot_expression(rna_ad, node_id) if is_gene else plot_peak(atac_ad, f'Interval_{node_id}')
    ax.figure.savefig(temp_file_path)
    ax.figure.clear()
  return { 'path': temp_file_path }

if __name__ == '__main__':  # run locally
  with TemporaryDirectory(dir='./') as tmpdir:
    app.config['tmpdir'] = os.path.basename(tmpdir)
    app.run()
else:  # run on server
  os.mkdir('tmp')
  app.config['tmpdir'] = 'tmp'