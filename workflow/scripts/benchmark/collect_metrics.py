import pandas as pd
import json
from tqdm import tqdm
from snakemake.script import snakemake


def collect_metrics(input_files, output_file):
    df = aggregate_json(input_files)
    df.to_csv(output_file, index=False)


def aggregate_json(json_files):
    data = [read_json(f) for f in tqdm(json_files)]
    df = pd.json_normalize(data)
    return df


def read_json(json_file):
    with open(json_file, "r") as f:
        data = json.load(f)
    return data


collect_metrics(snakemake.input, snakemake.output[0])
