#!/usr/bin/env python3

import argparse
from scan_scn.singlecell.scn_gene_importance.ml_utils import (
    load_adata_for_modeling, get_data_and_target, split_train_test
)
from scan_scn.singlecell.scn_gene_importance import train_random_forest, train_xgboost

def run_training(
    model_type, h5ad_path, gene_path, output_path, top_n=250
):
    adata = load_adata_for_modeling(h5ad_path, gene_path)
    X, y = get_data_and_target(adata)
    gene_names = adata.var.index

    X_train, X_test, y_train, y_test = split_train_test(X, y)

    if model_type == "random_forest":
        model, params = train_random_forest.train_random_forest(X_train, y_train, param_grid={
            'n_estimators': [200, 300],
            'max_depth': [20, 30, None],
            'min_samples_split': [2, 5],
            'min_samples_leaf': [1, 2],
            'max_features': ['sqrt', 'log2']
        })
        metrics = train_random_forest.evaluate_model(model, X_test, y_test)
        top_genes = train_random_forest.get_feature_importance(model, gene_names, top_n)

    elif model_type == "xgboost":
        model, params = train_xgboost.train_xgboost(X_train, y_train, param_grid={
            'n_estimators': [200, 300],
            'max_depth': [20, 30, None],
            'learning_rate': [0.01, 0.1],
            'reg_alpha': [0, 1],
            'reg_lambda': [1, 10],
            'colsample_bytree': [0.8, 1]
        })
        metrics = train_xgboost.evaluate_model(model, X_test, y_test)
        top_genes = train_xgboost.get_feature_importance(model, gene_names, top_n)

    else:
        raise ValueError("Model type must be 'random_forest' or 'xgboost'")

    print("Best parameters:", params)
    print("Test RMSE:", metrics["rmse"])
    print("Test R2:", metrics["r2"])
    top_genes.to_csv(output_path, sep="\t", index=False)
    print(f"Saved top {top_n} feature importance to {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", choices=["random_forest", "xgboost"], required=True)
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--gene_names", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--top_n", type=int, default=250)
    args = parser.parse_args()

    run_training(
        model_type=args.model,
        h5ad_path=args.h5ad,
        gene_path=args.gene_names,
        output_path=args.output,
        top_n=args.top_n
    )