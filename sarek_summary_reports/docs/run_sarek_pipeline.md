# Running Sarek
Directions on how to run Sarek

## Basic Kickoff
We use GATK.GRCh37 by default as the reference. You might have to change that in the future when CCG migrates
to GRCh38. I recommend using tmux instead of screen with nextflow. It looks much better in the
terminal (i.e. no screen tearing). For a basic tutorial of tmux: https://www.redhat.com/en/blog/introduction-tmux-linux
```bash
gcloud compute ssh --zone "us-central1-a" "nextflow4" --project "ccgd-174313"
tmux  # start new session
# Make sure to add the bait profile
/opt/nextflow/25.10.2/nextflow-25.10.2-dist run  /opt/sarek/main.nf \
    -profile docker,google,ccg,[bait_profile] \
    --input gs://ccg-bioinformatics/mnt/storage/centers/ccgd/nap24/pipeline_rnd/pipeline_sarek/SE81226a-ULPWGS/csv/SE81226a_sarek_samplesheet.csv \
    --outdir gs://ccg-bioinformatics/mnt/storage/centers/ccgd/nap24/pipeline_rnd/pipeline_sarek/SE81226a-ULPWGS \
    -w gs://ccg-bioinformatics/mnt/scratch/ccgd/global/nextflow_execution/SE81226a-ULPWGS
```

## Generate PON
Pretty much the same except use the /opt/sarek-create-pon directory, create_pon profile, and make sure no PON is set.
Use the create_pon flag as the last profile to make sure the PON is not set. I set Sarek to use the Broaod PON by
default so that's why you have to do this.
```bash
/opt/nextflow/25.10.2/nextflow-25.10.2-dist run  /opt/sarek-create-pon/main.nf \
    -profile docker,google,ccg,create_pon \
    --input gs://ccg-bioinformatics/mnt/storage/centers/ccgd/nap24/pipeline_rnd/pipeline_sarek/SCB1046d_pon_popv3/csv/recalibrated.csv \
    --outdir gs://ccg-bioinformatics/mnt/storage/centers/ccgd/nap24/pipeline_rnd/pipeline_sarek/SCB1046d_pon_popv3 \
    -w gs://ccg-bioinformatics/mnt/scratch/ccgd/global/nextflow_execution/SCB1046d_pon_popv3
```


## Troubleshooting
### Restart from middle of run
It's pretty common for nextflow to trip over itself, you still have to check on it every few hours while it's running.
It also has a e-mail feature, but I haven't figured out how to use it, yet.
An example error that just needs a restart:
```
ERROR ~ Error executing process > 'NFCORE_SAREK:SAREK:FASTQ_PREPROCESS_GATK:BAM_BASERECALIBRATOR:GATK4_BASERECALIBRATOR (G06-65026-A9-P2-B5_S31)'

Caused by:
  ALREADY_EXISTS: Resource "projects/738703343162/locations/us-central1/jobs/nf-90474ae9-1768526411909" already exists
```

Use the `-resume` flag to resume. Example Resume command:
```bash
/opt/nextflow/25.10.2/nextflow-25.10.2-dist run /opt/sarek/main.nf \
  -resume \
  -profile docker,google,ccg,POPv3_1_S3050142 \
  --input gs://ccg-bioinformatics/mnt/storage/centers/ccgd/nap24/pipeline_rnd/pipeline_sarek/20260114_CCG1200c_update_bed_hugo_and_filter_false/csv/CCG1200c_sarek_3.7.1_samplesheet.csv \
  --outdir gs://ccg-bioinformatics/mnt/storage/centers/ccgd/nap24/pipeline_rnd/pipeline_sarek/20260114_CCG1200c_update_bed_hugo_and_filter_false \
  -w gs://ccg-scratch/mnt/scratch/ccgd/global/nextflow_execution/20260114_CCG1200c_update_bed_hugo_and_filter_false
```

⚠️ **NOTE:** Nextflow does take about 1 hour to resume towards the end of the pipeline. Unlike cromwell's resume
feature, which just saves the last failure in the database, nextflow has to retrace it's steps through the run's
execution graph. On a positive note, Cromwell's resume only works on clusters (not the cloud) so nextflow has the
upper hand there.

### Rerun a tool
Let's say you when you first ran a pipeline run, you forgot to make a PON. You can just rerun Mutect2 or CNVKit. Instead
of the starting samplesheet, pass the `recalibrated.csv` samplesheet created automatically by Sarek. Add
`--step variant_calling` and `--tools cnvkit,mutect2` to the command.
Putting it all together, here is an example rerun command:
```bash
/opt/nextflow/25.10.2/nextflow-25.10.2-dist run /opt/sarek/main.nf \
  -profile docker,google,ccg,wgs    \
  --step variant_calling \
  --tools cnvkit,mutect2 \
  --input gs://ccg-bioinformatics/mnt/storage/centers/ccgd/Illumina/pipeline_sarek/SE81226a-ULPWGS/csv/recalibrated.csv \
  --outdir gs://ccg-bioinformatics/mnt/storage/centers/ccgd/Illumina/pipeline_sarek/SE81226a-ULPWGS \
  -w gs://ccg-scratch/mnt/scratch/ccgd/global/nextflow_execution/SE81226a-ULPWGS
```