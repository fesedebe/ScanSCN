
import pandas as pd
from scan_scn.bulk.scoring.scn_scoring import run_pca_projection

def run_scn_projection(query_path, output_path):
    query_df = pd.read_csv(query_path, index_col=0)
    scores = run_pca_projection(query_df)
    scores.to_csv(output_path, sep="\t", index=False)

if __name__ == "__main__":
    query_path = "data/input/ovarian_log2UQ.txt"
    output_path = "data/output/UCLA_projected_onto_SCN_PCA_predicted.scores.txt"
    run_scn_projection(query_path, output_path)
    print(f"Saved projected scores to {output_path}")
