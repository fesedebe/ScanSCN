
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from pathlib import Path

try:
    from factor_analyzer import Rotator  # Optional: for varimax
except ImportError:
    Rotator = None

DEFAULT_SCN_PATH = Path("data/internal/balanisSCN_log2UQ.parquet")
ALTERNATE_CSV_PATH = Path("data/internal/balanisSCN_log2UQ.csv.gz")

def load_default_scn(path: Path = DEFAULT_SCN_PATH) -> pd.DataFrame:
    if path.suffix == ".parquet":
        return pd.read_parquet(path)
    elif path.suffix in [".csv", ".gz"]:
        return pd.read_csv(path, index_col=0)
    else:
        raise ValueError(f"Unsupported file format for: {path}")

def match_genes(scn_df: pd.DataFrame, query_df: pd.DataFrame) -> tuple:
    common = scn_df.index.intersection(query_df.index)
    return scn_df.loc[common], query_df.loc[common]

def run_pca_projection(
    query_df: pd.DataFrame,
    scn_df: pd.DataFrame = None,
    n_components: int = 3,
    varimax: bool = True,
    rotate_pc1: bool = True,
) -> pd.DataFrame:
    if scn_df is None:
        scn_df = load_default_scn()

    # Match genes
    scn_df, query_df = match_genes(scn_df, query_df)

    # Transpose to shape: samples x genes
    scn_mat = scn_df.T
    query_mat = query_df.T

    # PCA
    pca = PCA(n_components=n_components)
    scn_scaled = StandardScaler().fit_transform(scn_mat)
    pca.fit(scn_scaled)

    # Project query into SCN PCA space
    query_scaled = StandardScaler().fit_transform(query_mat)
    scores = query_scaled @ pca.components_.T

    # Optional: varimax rotation
    if varimax:
        if Rotator is None:
            raise ImportError("Install `factor_analyzer` to use varimax.")
        rotator = Rotator()
        scores = rotator.fit_transform(scores)

    # Optional: flip PC1
    if rotate_pc1:
        scores[:, 0] *= -1

    projected_df = pd.DataFrame(
        scores, index=query_mat.index,
        columns=[f"PC{i+1}" for i in range(scores.shape[1])]
    )
    projected_df.insert(0, "Sample", projected_df.index)
    return projected_df
