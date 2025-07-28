import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import GridSearchCV

def train_random_forest(X_train, y_train, param_grid=None, cv=5, n_jobs=4, verbose=1):
    model = RandomForestRegressor(random_state=824)
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
    importances = model.feature_importances_
    df = pd.DataFrame({"gene": gene_names, "importance": importances})
    df = df.sort_values("importance", ascending=False).reset_index(drop=True)
    df["weight"] = df["importance"] / df["importance"].sum()
    return df.head(top_n)