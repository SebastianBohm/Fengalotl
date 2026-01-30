from pathlib import Path
import numpy as np
import pandas as pd

BASE_DIR = Path(__file__).resolve().parent.parent.parent
DATA_DIR = BASE_DIR / 'data'

# Clustering options - leiden clustering and structure annotation
CLUSTERING_OPTIONS = {
    'spatial_leiden_e30_s8': 'Leiden clustering',
    'structure': 'Structure annotation'
}

GENES = np.load(DATA_DIR / 'genes.npy', allow_pickle=True).tolist()
GENES_LABEL = GENES  # No dash format in this dataset
DATA = np.load(DATA_DIR / 'samples.npy', allow_pickle=True).tolist()

# Load annotation metadata for differential gene expression markers
ANNOTATION_FILE = BASE_DIR.parent / 'Result' / 'Adult_meta_DGE_markers.csv'
try:
    ANNOTATION_DF = pd.read_csv(ANNOTATION_FILE, index_col=0)
    
    # Create gene annotation lookup (gene_id -> human gene symbol)
    GENE_ANNOTATION = {}
    for gene_id in ANNOTATION_DF.index.unique():
        row = ANNOTATION_DF.loc[gene_id]
        if isinstance(row, pd.DataFrame):
            row = row.iloc[0]
        # Try to get a meaningful gene name (prefer Axolotl annotation, then human)
        gene_name = row.get('Axolotl_tanaka_annotated_gene', '-')
        if pd.isna(gene_name) or gene_name == '-' or gene_name == 'N/A':
            gene_name = row.get('hs_gene', '-')
        if pd.isna(gene_name) or gene_name == '-':
            gene_name = gene_id
        GENE_ANNOTATION[gene_id] = gene_name

    # Get unique cell types from the DGE markers
    CELLTYPES = sorted(ANNOTATION_DF['Celltype'].unique().tolist())
except FileNotFoundError:
    ANNOTATION_DF = None
    GENE_ANNOTATION = {}
    CELLTYPES = []
