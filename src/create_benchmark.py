import os
import pandas as pd


def all_human_transcription_factors(filepath="resources/mmc2.xlsx"):
    # The Human Transcription Factors - Lambert et al. 2018
    # https://www.cell.com/cell/fulltext/S0092-8674(18)30106-5
    df = pd.read_excel(filepath, sheet_name="Table S1. Related to Figure 1B")
    human_tfs = df[df["Is TF?"].str.strip().str.lower() == "yes"]["Unnamed: 1"]
    return human_tfs


def get_filtered_transcription_factors(df, human_tfs, n_cell_types=5):
    # TFs present in at least n_cell_types
    # and are in the list of human transcription factors
    filtered_tfs = (
        df.query("Assay == 'TF ChIP-seq'")
        .groupby(["Experiment target"])["Biosample term name"]
        .nunique()
        .to_frame("n")
        .query("n >= @n_cell_types")
        .query("index.str.replace('-human', '') in @human_tfs")
        .index.tolist()
    )
    return filtered_tfs


df = pd.read_parquet("/mnt/wdssd/tfsage/encode/resources/index.parquet")
human_tfs = all_human_transcription_factors()
filtered_tfs = get_filtered_transcription_factors(df, human_tfs, n_cell_types=5)

# Remove all perturbed experiments
df_filtered = (
    df.query("`Biosample treatments`.isna()")
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

# Match TFs with ATAC or H3K27ac in same cell type
benchmark = pd.merge(
    df_filtered.query("target in @filtered_tfs"),
    df_filtered.query("assay != 'TF ChIP-seq'"),
    how="inner",
    on=["cell_type"],
    suffixes=("_factor", "_query"),
)

os.makedirs("data", exist_ok=True)
benchmark.to_csv("data/benchmark.csv", index=False)
