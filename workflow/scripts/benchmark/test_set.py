from tfsage import utils
from types import SimpleNamespace
from snakemake.script import snakemake


def test_set(inputs, output_file, params):
    params = SimpleNamespace(**params.cfg)
    df = utils.generate_test_set(
        inputs.query_file,
        inputs.target_file,
        params.peak_width,
    )
    df = df.query("chrom not in @params.exclude_chroms")
    df = utils.stratified_sample(
        df,
        params.n_samples,
        params.p_positive,
        params.random_state,
    )
    df.to_parquet(output_file, index=False)


test_set(snakemake.input, snakemake.output[0], snakemake.params)
