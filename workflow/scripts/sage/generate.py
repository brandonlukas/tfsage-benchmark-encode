import os
import pandas as pd
import dask.dataframe as dd
from tfsage.generate import synthesize_experiments
from snakemake.script import snakemake


def generate(input_file, output_file, params, wildcards):
    # Fetch ranked list from search data
    ranked_list = (
        dd.read_parquet(input_file)
        .query("id_query == @wildcards.query_id", local_dict={"wildcards": wildcards})
        .query("target == @wildcards.factor", local_dict={"wildcards": wildcards})
        .query("cell_type != cell_type_query")
        .compute()
        .nlargest(int(wildcards.n), "score", "all")
    )

    bed_files = [
        os.path.join(params.downloads_dir, f"{x}.bed") for x in ranked_list["id"]
    ]
    weights = ranked_list["score"].tolist()
    result = generate_result(bed_files, weights)
    result.to_parquet(output_file)


def generate_result(bed_files, weights) -> pd.DataFrame:
    if len(bed_files) == 1:
        # "Hack" to make the function work with only one bed file
        bed_files = [bed_files[0], bed_files[0]]
        df = synthesize_experiments(bed_files, weights).drop("file_1", axis=1)
    else:
        df = synthesize_experiments(bed_files, weights)
    return df


generate(snakemake.input[0], snakemake.output[0], snakemake.params, snakemake.wildcards)
