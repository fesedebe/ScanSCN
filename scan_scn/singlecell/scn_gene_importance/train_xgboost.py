import xgboost as xgb
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import mean_squared_error, r2_score
import pandas as pd

def train_xgboost(X_train, y_train, param_grid=None, cv=5, n_jobs=4, verbose=1):
    model = xgb.XGBRegressor(random_state=824)
    if param_grid:
        search = GridSearchCV(model, param_grid, cv=cv, n_jobs=n_jobs, scoring="r2", verbose=verbose)
        search.fit(X_train, y_train)
        return search.best_estimator_, search.best_params_
    else:
        model.fit(X_train, y_train)
        return model, {}

def evaluate_model(model, X_test, y_test):
    y_pred = model.predict(X_test)
    return {
        "rmse": mean_squared_error(y_test, y_pred, squared=False),
        "r2": r2_score(y_test, y_pred)
    }

def get_feature_importance(model, gene_names, top_n=250):
    booster = model.get_booster()
    importances = booster.get_score(importance_type='gain')

    importance_df = pd.DataFrame(list(importances.items()), columns=["feature", "importance"])
    importance_df["feature_index"] = importance_df["feature"].str.extract(r"f(\d+)").astype(int)
    importance_df["gene"] = gene_names[importance_df["feature_index"]].values
    importance_df = importance_df.sort_values("importance", ascending=False)
    importance_df["weight"] = importance_df["importance"] / importance_df["importance"].sum()
    return importance_df[["gene", "importance", "weight"]].head(top_n)