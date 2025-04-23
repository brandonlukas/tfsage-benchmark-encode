import os


rule generate:
    input:
        rules.search.output[0].format(
            features=config["sage"]["features"],
            embeddings=config["sage"]["embeddings"],
            metric=config["sage"]["metric"],
        ),
    output:
        os.path.join(
            config["output_dir"],
            "sage/{n}/{query_id}/{factor}.parquet",
        ),
    params:
        downloads_dir=config["index"]["downloads_dir"],
    script:
        "../scripts/generate/generate.py"
