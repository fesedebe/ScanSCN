import scanpy as sc
import pandas as pd
from sklearn.model_selection import train_test_split

def load_adata_for_modeling(h5ad_path, gene_names_path=None):
    adata = sc.read_h5ad(h5ad_path)
    if gene_names_path:
        gene_names = pd.read_csv(gene_names_path, sep="\t")['gene_names'].values
        adata.var.index = gene_names
    return adata

def get_data_and_target(adata, target_col="SCNBalanis_AUC"):
    X = adata.X
    y = adata.obs[target_col].values
    return X, y

def split_train_test(X, y, test_size=0.2, random_state=824):
    return train_test_split(X, y, test_size=test_size, random_state=random_state)