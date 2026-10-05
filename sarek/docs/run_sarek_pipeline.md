# Running Sarek
Directions on how to run Sarek

## Basic Kickoff
We use GATK.GRCh37 by default as the reference. You might have to change that in the future when CCG migrates
to GRCh38. I recommend using tmux instead of screen with nextflow. It looks much better in the
terminal (i.e. no screen tearing). For a basic tutorial of tmux: https://www.redhat.com/en/blog/introduction-tmux-linux
This typically takes 5-8 hours to run.
Nextflow creates a .nextflow directory that allows for `-resume`
    * This .nextflow directory you run in must support [file system locks](https://github.com/ncbi/egapx/issues/150)
    * If you want a new simultaneous nextflow run, you must
      [create a new .nextflow directory](https://github.com/nextflow-io/nextflow/discussions/3883) otherwise the locks
      will interfere.
```bash
gcloud compute ssh --zone "us-central1-a" "nextflow4" --project "ccgd-174313"
tmux  # start new session
PROJECT_CODE=[your_project_code]
mkdir /data/${PROJECT_CODE}
cd /data/${PROJECT_CODE}  # When you run nextflow, it adds a hidden .nextflow folder.
# Make sure to add the bait profile
/opt/nextflow/25.10.2/nextflow-25.10.2-dist run  /opt/sarek/main.nf \
    -profile docker,google,ccg,[bait_profile] \
    --input gs://ccg-bioinformatics/mnt/storage/centers/ccgd/Illumina/pipeline_sarek/${PROJECT_CODE}/csv/${PROJECT_CODE}_sarek_samplesheet.csv \
    --outdir gs://ccg-bioinformatics/mnt/storage/centers/ccgd/Illumina/pipeline_sarek/${PROJECT_CODE} \
    -w gs://ccg-scratch/mnt/scratch/ccgd/global/nextflow_execution/${PROJECT_CODE}
```

## Run Somalier/Fingerprinting
1. Run somalier extract. Somalier has not been incorporated into the pipeline. Pass the recalibrated.csv samplesheet
into somalier. Manually creating recalibrated.csv is **NOT** required. It's produced by Sarek. It will take about an
1-2 hours to run. An example of the recalibrated.csv is as follows. The patient is the individual name. The sex can be 
NA if unknown (the column is not actually used by somalier - just output by sarek). The sample is the sample name.
The cram and crai point to the files that will be fingerprinted. These can also be genomic vcf files (created by 
GATK HaplotypeCaller). They probably can also bam and bai files (I haven't tried). The restrictions of what
cram and crai can be are set by
[somalier extract](https://github.com/brentp/somalier/blob/v0.2.19/README.md#somalier-extract), itself.
```
patient,sex,status,sample,cram,crai
A001,NA,0,CEPH1408-H1_A001,gs://ccg-bioinformatics/mnt/storage/centers/ccgd/nap24/pipeline_rnd/pipeline_sarek/20260114_CCG1200c_update_bed_hugo_and_filter_false/preprocessing/recalibrated/CEPH1408-H1_A001/CEPH1408-H1_A001.recal.cram,gs://ccg-bioinformatics/mnt/storage/centers/ccgd/nap24/pipeline_rnd/pipeline_sarek/20260114_CCG1200c_update_bed_hugo_and_filter_false/preprocessing/recalibrated/CEPH1408-H1_A001/CEPH1408-H1_A001.recal.cram.crai
A001,NA,1,A001_X101_49707-1-B1,gs://ccg-bioinformatics/mnt/storage/centers/ccgd/nap24/pipeline_rnd/pipeline_sarek/20260114_CCG1200c_update_bed_hugo_and_filter_false/preprocessing/recalibrated/A001_X101_49707-1-B1/A001_X101_49707-1-B1.recal.cram,gs://ccg-bioinformatics/mnt/storage/centers/ccgd/nap24/pipeline_rnd/pipeline_sarek/20260114_CCG1200c_update_bed_hugo_and_filter_false/preprocessing/recalibrated/A001_X101_49707-1-B1/A001_X101_49707-1-B1.recal.cram.crai
```
After you have a samplesheet.csv too run somalier. You can run it by the following command:
```bash
tmux
cd /data
PROJECT_CODE=[Your project code]
/opt/nextflow/25.10.2/nextflow-25.10.2-dist run /opt/ccg-nf-extensions/workflows/somalier/main.nf \
  -c /opt/ccg-nf-extensions/nextflow.config \
  -profile docker,google \
  --sample_sheet gs://ccg-bioinformatics/mnt/storage/centers/ccgd/Illumina/pipeline_sarek/${PROJECT_CODE}/csv/recalibrated.csv \
  --outdir gs://ccg-bioinformatics/mnt/storage/centers/ccgd/Illumina/pipeline_sarek/${PROJECT_CODE} \
  -w gs://ccg-scratch/mnt/scratch/ccgd/global/nextflow_execution/${PROJECT_CODE}
```
2. Create a groups_file.csv if necessary. The command to generate the groups file is the following:
```bash
/opt/ccg-rotation-suite/sarek/create_somalier_groups.py -p ${PROJECT_CODE} \
    -r /mnt/storage/centers/ccgd/Illumina/pipeline_sarek/${PROJECT_CODE}/csv/recalibrated.csv
cp somalier_groups.csv /mnt/storage/centers/ccgd/Illumina/pipeline_sarek/${PROJECT_CODE}/csv/
```
There is no need a groups_file.csv if all samples derive from different
organisms/patients/individuals. The structure of the groups_file is mentioned in the
[somalier relate documentation](https://github.com/brentp/somalier/blob/v0.2.19/README.md#somalier-relate).
Basically, comma seperated groups by line:
```
normal1,tumor1a,tumor1b
normal2,tumor2a
```
3. Run Somalier relate. The default coverage depth is 7. We are trying `-d 20` as our standard for high-coverage runs. It's possible that some project will benefit from adjusting this threshold.:
```bash
docker run \
    -v /mnt/storage/centers/ccgd/Illumina/pipeline_sarek/${PROJECT_CODE}:/mnt/storage/centers/ccgd/Illumina/pipeline_sarek/${PROJECT_CODE} \
    -w /mnt/storage/centers/ccgd/Illumina/pipeline_sarek/${PROJECT_CODE}/reports/somalier \
    -it quay.io/biocontainers/somalier:0.2.19--h0c29559_0 bash
somalier relate -d 20 --groups ../../csv/somalier_groups.csv $(find . -name "*.somalier")
```

## Run Picard HsMetrics
1. Convert the .bed to .interval_list if not already done.
See [step 5](./migrate_old_bait_and_pon.md#convert-riker-interval_list-to-sarek-bed-file--hsmetrics-interval_list)
2. Use the markduplicates.csv to run picard. After quality recalibration, all the bases outside of the bait regions
will be marked as low quality and unusable.
```bash
/opt/nextflow/25.10.2/nextflow-25.10.2-dist run /opt/ccg-nf-extensions/workflows/picard/collecthsmetrics/main.nf \
  -c /opt/ccg-nf-extensions/nextflow.config \
  -profile docker,google \
  --converted_csv gs://ccg-bioinformatics/mnt/storage/centers/ccgd/Illumina/pipeline_sarek/${PROJECT_CODE}/csv/markduplicates.csv \
  --bait_intervals  gs://ccg-bioinformatics/mnt/storage/centers/ccgd/reference_files/baits/GRCh37/TLP_v3_TE_98697476/TLP_v3_TE_98697476_targets_20260508.interval_list \
  --target_intervals gs://ccg-bioinformatics/mnt/storage/centers/ccgd/reference_files/baits/GRCh37/TLP_v3_TE_98697476/TLP_v3_TE_98697476_targets_20260508.interval_list \
  --outdir /mnt/storage/centers/ccgd/Illumina/pipeline_sarek/${PROJECT_CODE} \
  -w gs://ccg-scratch/mnt/scratch/ccgd/global/nextflow_execution/${PROJECT_CODE}
```
3. Copy the header. Then copy all the hsmetrics from all the samples into 1 file.
```bash
TODAY_DATE=$(date --iso-8601 | tr -d "-")
PROJECT_CODE_LOWERCASE=$(echo ${PROJECT_CODE} | tr "[:upper:]" "[:lower:]")
OUTPUT_FOLDER="../ccg_${PROJECT_CODE_LOWERCASE}_reports"
OUTPUT_FILE="${OUTPUT_FOLDER}/${TODAY_DATE}_${PROJECT_CODE_LOWERCASE}_hsmetrics_report.csv"

cd /mnt/storage/centers/ccgd/Illumina/pipeline_sarek/${PROJECT_CODE}/reports/picard
mkdir -p ${OUTPUT_FOLDER}
readarray -d '' hsmetrics_files < <(find . -name "*.CollectHsMetrics.coverage_metrics" -type f -print0)
# Copy header
head -n7 "${hsmetrics_files[0]}" | tail -n1 > ${OUTPUT_FILE}
for one_file in "${hsmetrics_files[@]}"; do
  SAMPLE_NAME=$(basename ${one_file} .CollectHsMetrics.coverage_metrics)
  echo $SAMPLE_NAME
  # Copy data from each file.
  head -n8 $one_file | tail -n1 | tr -d '\n' | sed "s/\t\t$//" >> ${OUTPUT_FILE}
  echo -e "${SAMPLE_NAME}" >> ${OUTPUT_FILE}
done
```

## Edge Cases
### PDX (Patient Derived Xenografts) Project
1. Run the project with a `-profile pdx`. For example:
```bash
/opt/nextflow/25.10.2/nextflow-25.10.2-dist run  /opt/sarek-dev-nap24/main.nf \
    -profile docker,google,ccg,POPv3_1_S3050142,pdx \
    --input gs://ccg-out/mnt/storage/centers/ccgd/Illumina/pipeline_sarek/SO31074a_no_xengstor/csv/20260905_so31074a_sarek_samplesheet.csv \
    --outdir gs://ccg-out/mnt/storage/centers/ccgd/Illumina/pipeline_sarek/SO31074a_bbsplit2 \
    -w gs://ccg-scratch/mnt/scratch/ccgd/global/nextflow_execution/SO31074a_bbsplit
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
feature, which just saves the last failure in the database, nextflow has to retrace its steps through the run's
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