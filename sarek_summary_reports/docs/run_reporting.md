# Sarek Reporting Scripts

For Sarek, there are various reporting scripts. I would like to incorporate these into the Sarek pipeline, someday.
However, that is a lot of work that remains to be done.
1. [Multi QC (FastQC) report](#1-multiqc-fastqc-report)
2. [Aggregation Metrics Report](#2-aggregation-metrics-report)
3. [Copy Number Variation (CNV) report](#3-copy-number-variation-cnv-report)
4. [Single nucleotide variant and Insertion-Deletion (SNV & Indel) Report](#4-single-nucleotide-variant-and-insertion-deletion-snv--indel-report)
5. 🚧 **Under contruction** Structure Variant (SV) Report 🚧
6. [Download and copy to the writable](#6-download-and-copy-to-the-writable)

## 1. MultiQC (FastQC) Report
1. This is run automatically by Sarek. You should see it in the multiqc folder that Sarek generates. No additional steps
are needed.

## 2. Aggregation Metrics report
Remember to run [somalier](./run_sarek_pipeline.md#run-somalierfingerprinting) before this.
1. Run the following commands.
```bash
gcloud compute ssh --zone "us-central1-a" "nextflow4" --project "ccgd-174313"
PROJECT_CODE=[YOUR PROJECT CODE]
cd /mnt/storage/centers/ccgd/Illumina/pipeline_sarek/${PROJECT_CODE}
/opt/ccgd-rotation-suite/sarek_summary_reports/aggregate_report.py -p ${PROJECT_CODE} --bucket .
```

## 3. Copy Number Variation (CNV) report
1. Run the following commands.
```bash
gcloud compute ssh --zone "us-central1-a" "nextflow4" --project "ccgd-174313"
PROJECT_CODE=[YOUR PROJECT CODE]
cd /mnt/storage/centers/ccgd/Illumina/pipeline_sarek/${PROJECT_CODE}
/opt/ccgd-rotation-suite/sarek_summary_reports/cnv_report.py -p ${PROJECT_CODE} --bucket .
```

## 4. Single nucleotide variant and Insertion-Deletion (SNV & Indel) Report
1. Run the following commands. Most of the reporting commands are pretty instant, but the SNV & Indel report takes
about 30 seconds - 1 minute to run.
```bash
gcloud compute ssh --zone "us-central1-a" "nextflow4" --project "ccgd-174313"
PROJECT_CODE=[YOUR PROJECT CODE]
cd /mnt/storage/centers/ccgd/Illumina/pipeline_sarek/${PROJECT_CODE}
/opt/ccgd-rotation-suite/sarek_summary_reports/snv_indel_report.py -p ${PROJECT_CODE} --bucket . \
    --bed-file /mnt/storage/centers/ccgd/reference_files/baits/GRCh37/[Your bed file]
```


## 6. Download and copy to the writable
2. Download the files to your desktop to view:
```bash
gcloud storage cp -r "gs://ccg-bioinformatics/mnt/storage/centers/ccgd/Illumina/pipeline_sarek/${PROJECT_CODE}/reports/ccg_*_reports" .
```
3. Copy the files to the writable:
```bash
cd /media/rcsm/rc-stor8/writable/Illumina/Sarek_Qcs_and_Reports/
mkdir -p [parent_project]/[sub_project]
cp -r [report_dir] [parent_project]/[sub_project]
```

## Edge Cases
### Regenerate MultiQC
Let's say you ran Sarek in pieces (e.g. `--tools manta,tiddit` first; then ran `--tools mutect2,cnvkit` after
you created the PONs). You can recreate the multiqc report like so to include info from all Sarek runs:
```bash
PROJECT_CODE=[your project code]
docker run -v /mnt/storage/centers/ccgd/Illumina/pipeline_sarek/${PROJECT_CODE}:/mnt/outdir \
    -it community.wave.seqera.io/library/multiqc:1.33--ee7739d47738383b bash
cd /mnt/outdir
ls  # make sure you have the right dir
multiqc .
```
NOTE: the docker image url can be found in `modules/nf-core/multiqc/main.nf` in the sarek repo.
Move the files to mirror the sarek structure.
```
gcloud storage mv gs://ccg-bioinformatics/mnt/storage/centers/ccgd/Illumina/pipeline_sarek/${PROJET_CODE}/multiqc_data gs://ccg-bioinformatics/mnt/storage/centers/ccgd/Illumina/pipeline_sarek/${PROJET_CODE}/multiqc/
```