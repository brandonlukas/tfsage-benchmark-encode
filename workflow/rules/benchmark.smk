import os


rule test_set:
    input:
        query_file=os.path.join(config["index"]["downloads_dir"], "{query_id}.bed"),
        target_file=os.path.join(config["index"]["downloads_dir"], "{target_id}.bed"),
    output:
        os.path.join(
            config["output_dir"], "benchmark/test_set/{query_id}/{target_id}.parquet"
        ),
    params:
        cfg=config["benchmark"]["test_set"],
    script:
        "../scripts/benchmark/test_set.py"


rule compute_metrics_sage:
    input:
        predictions=rules.generate.output[0],
        test_set=rules.test_set.output[0],
    output:
        os.path.join(
            config["output_dir"],
            "benchmark/metrics/sage/{n}/{query_id}/{factor}/{target_id}.json",
        ),
    params:
        method_class="sage",
        method_name="TFSage-{n}",
    script:
        "../scripts/benchmark/compute_metrics.py"


rule compute_metrics_scan:
    input:
        predictions=rules.motif_scan.output[0],
        test_set=rules.test_set.output[0],
        m2f=rules.filter_motifs.output.m2f,
    output:
        os.path.join(
            config["output_dir"],
            "benchmark/metrics/scan/{motif_db}/{query_id}/{factor}/{target_id}.json",
        ),
    params:
        method_class="scan",
        method_name="{motif_db}",
    script:
        "../scripts/benchmark/compute_metrics.py"


rule collect_metrics:
    input:
        expand(
            expand(
                rules.compute_metrics_sage.output[0],
                zip,
                query_id=benchmark["id_query"],
                factor=benchmark["target_factor"],
                target_id=benchmark["id_factor"],
                allow_missing=True,
            ),
            n=config["sage"]["n"],
        ),
        expand(
            expand(
                rules.compute_metrics_scan.output[0],
                zip,
                query_id=benchmark["id_query"],
                factor=benchmark["target_factor"],
                target_id=benchmark["id_factor"],
                allow_missing=True,
            ),
            motif_db=config["scan"]["motif_dbs"],
        ),
    output:
        os.path.join(config["output_dir"], "benchmark/collect_results.csv"),
    script:
        "../scripts/benchmark/collect_metrics.py"
