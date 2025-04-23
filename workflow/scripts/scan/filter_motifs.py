import pandas as pd
from tqdm import tqdm
from snakemake.script import snakemake


def filter_motifs(inputs, outputs, params):
    # Filter motif2factors
    factors = [f.replace("-human", "") for f in params.factors]
    m2f_filtered = pd.read_csv(inputs.m2f, sep="\t").query("Factor in @factors")
    motifs = m2f_filtered["Motif"].unique()

    # Filter PFM file
    pfm_filtered = filter_pfm(inputs.pfm, motifs)

    # Save filtered motif2factors and PFM
    m2f_filtered.to_csv(outputs.m2f, sep="\t", index=False)
    with open(outputs.pfm, "w") as f:
        f.writelines(pfm_filtered)


def filter_pfm(pfm_file, motifs):
    with open(pfm_file, "r") as f:
        pfm_lines = f.readlines()

    include = False
    lines = []
    for line in tqdm(pfm_lines):
        if line.startswith(">"):
            motif = line.replace(">", "").strip()
            include = motif in motifs
        if include or line.startswith("#"):
            lines.append(line)

    return lines


filter_motifs(snakemake.input, snakemake.output, snakemake.params)
