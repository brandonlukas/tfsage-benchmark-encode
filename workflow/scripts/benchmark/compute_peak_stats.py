import json
import pybedtools
from snakemake.script import snakemake


def compute_peak_stats(inputs, output_file, wildcards):
    bedtools = pybedtools.BedTool()
    query_bed = pybedtools.BedTool(inputs.query_file).sort()
    target_bed = pybedtools.BedTool(inputs.target_file).sort()

    query_in_target_bed = bedtools.intersect(
        a=query_bed,
        b=target_bed,
        u=True,
        sorted=True,
    )
    target_in_query_bed = bedtools.intersect(
        a=target_bed,
        b=query_bed,
        u=True,
        sorted=True,
    )

    results = {
        "query_id": wildcards.query_id,
        "target_id": wildcards.target_id,
        "n_peaks_query": query_bed.count(),
        "n_peaks_target": target_bed.count(),
        "n_peaks_query_in_target": query_in_target_bed.count(),
        "n_peaks_target_in_query": target_in_query_bed.count(),
    }

    # Save results to file
    with open(output_file, "w") as f:
        json.dump(results, f, indent=4)


compute_peak_stats(snakemake.input, snakemake.output[0], snakemake.wildcards)
