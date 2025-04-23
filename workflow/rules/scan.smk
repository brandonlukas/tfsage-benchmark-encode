import os


rule filter_motifs:
    input:
        m2f=os.path.join(
            config["scan"]["motif_dbs_dir"], "{motif_db}.motif2factors.txt"
        ),
        pfm=os.path.join(config["scan"]["motif_dbs_dir"], "{motif_db}.pfm"),
    output:
        m2f="resources/motif_databases_filtered/{motif_db}_filtered.motif2factors.txt",
        pfm="resources/motif_databases_filtered/{motif_db}_filtered.pfm",
    params:
        factors=factor_list,
    script:
        "../scripts/scan/filter_motifs.py"


rule motif_scan:
    input:
        os.path.join(config["index"]["downloads_dir"], "{query_id}.bed"),
        pfm=rules.filter_motifs.output.pfm,
    output:
        os.path.join(config["output_dir"], "scan/{motif_db}/{query_id}.bed"),
    params:
        opts=lambda w, input: f"-g GRCh38 -p {input.pfm}",
    threads: 24
    resources:
        mem_mb=24000,
    conda:
        "gimme_env"
    shell:
        """
        gimme scan {input[0]} {params.opts} -N {threads} -b > {output}
        """
