import pandas as pd
import dask.dataframe as dd
import json
from tfsage import utils
from snakemake.script import snakemake


def compute_metrics(inputs, output_file, params, wildcards):
    predictions = load_predictions(inputs, params, wildcards)
    test_set = pd.read_parquet(inputs.test_set)

    n_samples = len(test_set)
    p_positive = test_set["positive"].mean()

    if predictions is not None:
        # Intersect predictions with test set
        df = utils.intersect_predictions_with_test_set(predictions, test_set)
        df["score"] = utils.sanitize_scores(df["score"])

        # Get true labels and scores
        y_true = df["positive"]
        y_score = df["score"]
        y_pred = y_score > 0

        # Compute metrics
        metrics = utils.compute_classification_metrics(y_true, y_score, y_pred)
    else:
        metrics = None

    # Combine everything into one dictionary
    results = {
        "method_class": params.method_class,
        "method_name": params.method_name,
        "query_id": wildcards.query_id,
        "factor": wildcards.factor,
        "target_id": wildcards.target_id,
        "n_samples": n_samples,
        "p_positive": p_positive,
        "metrics": metrics,
    }

    # Save results to file
    with open(output_file, "w") as f:
        json.dump(results, f, indent=4)


def load_predictions(inputs, params, wildcards) -> pd.DataFrame | None:
    if params.method_class == "sage":
        predictions = pd.read_parquet(inputs.predictions).assign(
            score=lambda x: x["sum"]
        )
    elif params.method_class == "scan":
        predictions = load_predictions_motif_scan(
            inputs.predictions, wildcards.factor, inputs.m2f
        )
    else:
        raise ValueError(f"Unknown method class: {params.method_class}")
    return predictions


def load_predictions_motif_scan(
    predictions_path, factor, m2f_path
) -> pd.DataFrame | None:
    factor = factor.replace("-human", "")
    motif_list = (
        pd.read_csv(m2f_path, sep="\t").query("Factor == @factor")["Motif"].tolist()
    )
    if len(motif_list) == 0:
        return None

    ddf = dd.read_csv(predictions_path, sep="\t", header=None, comment="#")
    ddf = ddf[ddf[3].isin(motif_list)]
    predictions = ddf.compute()
    predictions.columns = ["chrom", "start", "end", "motif", "score", "strand"]
    return predictions


compute_metrics(
    snakemake.input, snakemake.output[0], snakemake.params, snakemake.wildcards
)
