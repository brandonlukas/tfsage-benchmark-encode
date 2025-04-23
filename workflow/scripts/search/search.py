import pandas as pd
import numpy as np
from scipy.spatial.distance import cdist, pdist
from snakemake.script import snakemake


def search(inputs, output_file, params, wildcards):
    embeddings = pd.read_parquet(inputs.embeddings).set_index("__index_level_0__")
    metadata = (
        pd.read_parquet(inputs.metadata)
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
    p = pdist(embeddings.values, metric=wildcards.metric).var()
    scoring_function = lambda x: np.exp(-(x**2) / p)

    factor_ids = metadata.query("target in @params.factor_list")["id"].values
    embeddings_factor = embeddings.query("index in @factor_ids")
    embeddings_query = embeddings.query("index in @params.query_ids")

    df = (
        pd.DataFrame(
            cdist(embeddings_factor, embeddings_query, wildcards.metric),
            index=embeddings_factor.index,
            columns=embeddings_query.index,
        )
        .merge(
            metadata,
            left_index=True,
            right_on="id",
            how="left",
        )
        .melt(
            id_vars=metadata.columns,
            var_name="query_id",
            value_name="distance",
        )
        .merge(
            metadata,
            left_on="query_id",
            right_on="id",
            how="left",
            suffixes=("", "_query"),
        )
        .drop(columns=["query_id"])
        .assign(score=lambda df: scoring_function(df.distance))
        .reset_index(drop=True)
    )

    df.to_parquet(output_file, index=False)


search(snakemake.input, snakemake.output[0], snakemake.params, snakemake.wildcards)
