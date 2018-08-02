This repo will contain all the ccgd rotation scripts, rewriten to be more robust, reuable, and extendable.

Scripts were originally started 3 years ago (curr date 2017/11) but were written as single files containing all aspects, including aspects
that are reused. New versions will be rewritten in python instead of R, and will be written so they are more friendly to other analysts
to modify and extend. These core goals necessitate implementing the scripts using proper software engineering practices.

First scripts will be the non essential scripts in order to get started. This will include scripts for:
1) Running the ConPair pipeline for contamination testing of tumor-normal pairs
2) Comparing Picard hsInterval metrics to Agilent Design files to very the interval files were properly created

After these are done and implemented we will move onto moving the existing working scripts into python:
1) Lane Level report
2) Aggregation report
3) Variant report
4) Breakmer report
5) CNV report
6) Lego Plot
7) Fingerprint report

Notes: When Mahesh was here briefly he attempted a rewrite into python using snakemake. This approach is a little onerous in that snakemake
may be esoteric for some new analysts (implied requirement being that we want to keep everything in common libraries/languages) and makes some
added steps for running scripts that was actually simpler in R. Additionally some of his code mixed hardcoded aspects with abstracted functions,
we want to abstract and document as much as possible. That said, his project Metriq will be used as a starting point for the code/structure
we want to replicate.