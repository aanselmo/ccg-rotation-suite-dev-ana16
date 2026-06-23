#!/usr/local/miniforge3/bin/python

"""
Collect the data for the "FP_concordance_matrix" table in the aggregate report.
"""

__docformat__ = "reStructuredText"
__credits__ = ["Brandon Sie"]


from os.path import basename
import pandas as pd
import polars as pl
import narwhals.stable.v2 as nw


def debug_df(name, df, n=10):
    print(f"\n===== {name} (top {n}) =====")
    print(df.head(n))
    print("Shape:", df.shape)


def create_concordance_matrix(somalier_pairs: str):
    """
    Transform the somalier paired data output into a wide matrix.  
    """

    fingerprint_table = nw.from_native(pl.read_csv(
        somalier_pairs,
        separator="\t",
    ))

    # to make a complete matrix, combine (1) original paired data (2) inverted data (swap sample_a and sample_b
    # and (3) self-self pairs with relatedness value set to 1.0.

    pl_fingerprint_original = (fingerprint_table.to_native()
                               .select("#sample_a", "sample_b", "relatedness").rename({"#sample_a": "sample_a"}))
    pl_fingerprint_inverted = (pl_fingerprint_original.rename({"sample_a": "sample_b", "sample_b": "sample_a"})
        .select("sample_a", "sample_b", "relatedness"))
    fingerprint_samples_unique = (list(set(fingerprint_table["#sample_a"].to_list())
                                       .union(set(fingerprint_table["sample_b"].to_list()))))
    pl_fingerprint_self = pl.DataFrame({
        "sample_a": fingerprint_samples_unique,
        "sample_b": fingerprint_samples_unique,
        "relatedness": 1.0
    })

    pl_fingerprint_tall = pl.concat([pl_fingerprint_original, pl_fingerprint_inverted, pl_fingerprint_self], 
                                    how = "vertical").sort(by=["sample_a", "sample_b"])

    pl_fingerprint_matrix = (pl_fingerprint_tall.pivot(index = "sample_a", on = "sample_b", values = "relatedness")
                             .rename({"sample_a": "Sarek Relatedness"}))

    concordance_matrix = nw.from_native(pl_fingerprint_matrix)
    debug_df("concordance_matrix", concordance_matrix.to_native(), n=10)

    return concordance_matrix
