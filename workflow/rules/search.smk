import os


rule search:
    input:
        embeddings=os.path.join(
            config["index"]["data_dir"],
            "embeddings/{features}/{embeddings}.parquet",
        ),
        metadata=config["index"]["metadata"],
    output:
        os.path.join(
            config["output_dir"],
            "search/search/{features}/{embeddings}/{metric}.parquet",
        ),
    params:
        query_ids=query_ids,
        factor_list=factor_list,
    script:
        "../scripts/search/search.py"


rule compute_ranking_metrics:
    input:
        search=rules.search.output[0],
        metadata=config["index"]["metadata"],
        benchmark=config["inputs"]["benchmark"],
    output:
        os.path.join(
            config["output_dir"],
            "search/ranking_metrics/{features}/{embeddings}/{metric}.parquet",
        ),
    script:
        "../scripts/search/compute_ranking_metrics.py"


rule collect_ranking_metrics:
    input:
        collect(
            rules.compute_ranking_metrics.output[0],
            features=config["search"]["features"],
            embeddings=config["search"]["embeddings"],
            metric=config["search"]["metrics"],
        ),
    output:
        os.path.join(config["output_dir"], "search/collect_results.csv"),
    run:
        import pandas as pd

        df_list = [pd.read_parquet(f) for f in input]
        df = pd.concat(df_list, ignore_index=True)
        df.to_csv(output[0], index=False)
