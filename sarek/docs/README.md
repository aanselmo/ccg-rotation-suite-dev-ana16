# Sarek Directions
These are directions for running the Sarek pipeline and its reporting scripts.

## Basic Workflow
Steps of running sarek are:
1. [Setup for Sarek Run](./setup_for_sarek_run.md#rename-fastqs-and-create-samplesheet)
2. [Run the sarek pipeline](./run_sarek_pipeline.md#basic-kickoff)
3. [Run somalier/Fingerprinting](./run_sarek_pipeline.md#run-somalierfingerprinting)
4. [Run Picard hsmetrics](./run_sarek_pipeline.md#run-picard-hsmetrics)
5. [Run reporting scripts](./run_reporting.md)

## Edge Cases
### Demultiplex
Please see [Demultiplexing Docs](./demultiplex.md). MBCF usually does the demultiplexing for us so you don't have to.

### Migrate Old Bait and PON
Please see [Migrate Old Bait and PON Docs](./migrate_old_bait_and_pon.md)