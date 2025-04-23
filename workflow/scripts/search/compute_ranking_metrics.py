import pandas as pd
from sklearn import metrics
from tqdm import tqdm
from snakemake.script import snakemake


def compute_ranking_metrics(inputs, output_file, wildcards):
    unperturbed_ids = load_unperturbed_ids(inputs.metadata)
    df = (
        pd.read_parquet(inputs.search)
        .query("id in @unperturbed_ids")
        .assign(hit=lambda x: x["cell_type"] == x["cell_type_query"])
    )
    benchmark = (
        pd.read_csv(inputs.benchmark)
        .filter(["id_query", "target_factor"], axis=1)
        .drop_duplicates()
    )

    records = []
    for id_query, target_factor in tqdm(
        benchmark.itertuples(index=False), total=len(benchmark)
    ):
        ranked_list = df.query("id_query == @id_query & target == @target_factor")
        record = ranked_list_metrics(ranked_list)
        record["id_query"] = id_query
        record["target_factor"] = target_factor
        records.append(record)

    results = (
        pd.DataFrame.from_records(records)
        .merge(
            df[["id_query", "assay_query", "cell_type_query"]].drop_duplicates(),
            on=["id_query"],
        )
        .assign(
            features=wildcards.features,
            embeddings=wildcards.embeddings,
            metric=wildcards.metric,
        )
    )

    results.to_parquet(output_file, index=False)


def ranked_list_metrics(ranked_list: pd.DataFrame) -> dict:
    ranked_list = ranked_list.assign(
        rank=lambda x: x["score"].rank(ascending=False),
        precision=lambda x: x["hit"].cumsum() / x["rank"],
    )
    reciprocal_rank = 1 / ranked_list.query("hit", engine="python")["rank"].min()
    average_precision = ranked_list.query("hit", engine="python")["precision"].mean()
    ndcg_score = metrics.ndcg_score([ranked_list["hit"]], [ranked_list["score"]])
    results = {
        "reciprocal_rank": reciprocal_rank,
        "average_precision": average_precision,
        "ndcg_score": ndcg_score,
    }
    return results


def load_unperturbed_ids(metadata_file):
    df = (
        pd.read_parquet(metadata_file)
        .query("`Biosample treatments`.isna()")
        .query("`Biosample treatments amount`.isna()")
        .query("`Biosample treatments duration`.isna()")
        .query("`Biosample genetic modifications methods`.isna()")
        .query("`Biosample genetic modifications categories`.isna()")
        .query("`Biosample genetic modifications targets`.isna()")
        .query("`Biosample genetic modifications site coordinates`.isna()")
        .filter(
            [
                "File accession",
                "Assay",
                "Biosample term name",
                "Experiment target",
            ],
            axis=1,
        )
        .rename(
            columns={
                "File accession": "id",
                "Assay": "assay",
                "Biosample term name": "cell_type",
                "Experiment target": "target",
            }
        )
    )

    unperturbed_ids = df["id"].unique().tolist()
    return unperturbed_ids


compute_ranking_metrics(snakemake.input, snakemake.output[0], snakemake.wildcards)
