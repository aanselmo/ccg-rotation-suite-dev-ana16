# Running Sarek
Directions on how to run Sarek

## Basic Kickoff
We use GATK.GRCh37 by default as the reference. You might have to change that in the future when CCG migrates
to GRCh38. I recommend using tmux instead of screen with nextflow. It looks much better in the
terminal (i.e. no screen tearing). For a basic tutorial of tmux: https://www.redhat.com/en/blog/introduction-tmux-linux
This typically takes 5-8 hours to run.
```bash
gcloud compute ssh --zone "us-central1-a" "nextflow4" --project "ccgd-174313"
cd /data  # I try to keep the nextflow run database here. When you run it, it appends to the hidden .nextflow folder.
tmux  # start new session
# Make sure to add the bait profile
/opt/nextflow/25.10.2/nextflow-25.10.2-dist run  /opt/sarek/main.nf \
    -profile docker,google,ccg,[bait_profile] \
    --input gs://ccg-bioinformatics/mnt/storage/centers/ccgd/nap24/pipeline_rnd/pipeline_sarek/SE81226a-ULPWGS/csv/SE81226a_sarek_samplesheet.csv \
    --outdir gs://ccg-bioinformatics/mnt/storage/centers/ccgd/nap24/pipeline_rnd/pipeline_sarek/SE81226a-ULPWGS \
    -w gs://ccg-bioinformatics/mnt/scratch/ccgd/global/nextflow_execution/SE81226a-ULPWGS
```

## Run Somalier/Fingerprinting
1. Run somalier extract. Somalier has not been incorporated into the pipeline. Pass the recalibrated.csv samplesheet
into somalier. It will take about an 1-2 hours to run.
```bash
cd /data
PROJECT_CODE=[Your project code]
/opt/nextflow/25.10.2/nextflow-25.10.2-dist run /opt/somalier/main.nf \
  -c /opt/somalier/nextflow.config \
  -profile docker,google \
  --sample_sheet gs://ccg-bioinformatics/mnt/storage/centers/ccgd/Illumina/pipeline_sarek/${PROJECT_CODE}/csv/recalibrated.csv \
  --outdir /mnt/storage/centers/ccgd/Illumina/pipeline_sarek/${PROJECT_CODE} \
  -w gs://ccg-scratch/mnt/scratch/ccgd/global/nextflow_execution/${PROJECT_CODE}
```
2. Run Somalier relate:
```bash
docker run \
    -v /mnt/storage/centers/ccgd/Illumina/pipeline_sarek/${PROJECT_CODE}/reports/somalier:/mnt/storage/centers/ccgd/Illumina/pipeline_sarek/${PROJECT_CODE}/reports/somalier \
    -w /mnt/storage/centers/ccgd/Illumina/pipeline_sarek/${PROJECT_CODE}/reports/somalier \
    -it quay.io/biocontainers/somalier:0.2.19--h0c29559_0 bash
# You may or may not need a groups file. It depends if you have matching samples from the same patient/individual
# or not.
somalier relate [--groups groups_file.csv] $(find . -name "*.somalier")
```

## Generate PON
1. Run Sarek. Pretty much the same except use the /opt/sarek-create-pon directory, create_pon profile, and make sure no
PON is set. Use the create_pon flag as the last profile to make sure the PON is not set. I set Sarek to use the Broaod
PON by default so that's why you have to do this.
```bash
/opt/nextflow/25.10.2/nextflow-25.10.2-dist run  /opt/sarek-create-pon/main.nf \
    -profile docker,google,ccg,create_pon \
    --input gs://ccg-bioinformatics/mnt/storage/centers/ccgd/nap24/pipeline_rnd/pipeline_sarek/SCB1046d_pon_popv3/csv/recalibrated.csv \
    --outdir gs://ccg-bioinformatics/mnt/storage/centers/ccgd/nap24/pipeline_rnd/pipeline_sarek/SCB1046d_pon_popv3 \
    -w gs://ccg-bioinformatics/mnt/scratch/ccgd/global/nextflow_execution/SCB1046d_pon_popv3
```
2. Create an args file for all of the vcf files. Head of the file looks like. It contains args to pass to
GenomicsDBImport.
```bash
--variant gs://ccg-bioinformatics/mnt/storage/centers/ccgd/nap24/pipeline_rnd/pipeline_sarek/SCB1046d_pon_popv3_2/variant_calling/mutect2/OP-P1-A6-19-136b/OP-P1-A6-19-136b.bcftools.vcf.gz
--variant gs://ccg-bioinformatics/mnt/storage/centers/ccgd/nap24/pipeline_rnd/pipeline_sarek/SCB1046d_pon_popv3_2/variant_calling/mutect2/OP-P1-B6-19-123b/OP-P1-B6-19-123b.bcftools.vcf.gz
--variant gs://ccg-bioinformatics/mnt/storage/centers/ccgd/nap24/pipeline_rnd/pipeline_sarek/SCB1046d_pon_popv3_2/variant_calling/mutect2/OP-P1-C6-19-138b/OP-P1-C6-19-138b.bcftools.vcf.gz
```
3. On your desktop, run the following to increase the memory of the nextflow4 server. Then, **reboot** nextflow4.
```
gcloud compute instances set-machine-type nextflow4 \
    --zone us-central1-a --machine-type n2-highmem-64
```
4. On the nextflow4 server, I would next remount the drives because they still don't automount well after reboot.
```
sudo umount -l /mnt/scratch
sudo umount -l /mnt/storage
sudo mount /mnt/scratch
sudo mount /mnt/storage
```
5. Run GenomicsDBImport. This typically takes 1 day to run the GenomicsDBImport.
See the Broad website for more info: https://gatk.broadinstitute.org/hc/en-us/articles/30332078729883-CreateSomaticPanelOfNormals-BETA
```bash
nohup /opt/jdk-17.0.12/bin/java -Xmx450G -Xms20G -jar /opt/gatk/4.6.1.0/gatk-package-4.6.1.0-local.jar \
  GenomicsDBImport \
  --arguments_file /mnt/storage/centers/ccgd/nap24/pipeline_rnd/pipeline_sarek/SO31074a_20251217_w_mnp_flag/SO31074a_pon_vcfs_2.args \
  --genomicsdb-workspace-path /data/pon_db_popv3_add_SO31074a_5 \
  --max-num-intervals-to-import-in-parallel 10 \
  -R /mnt/storage/centers/ccgd/reference_files/igenomes/Homo_sapiens/igenomes/GATK/GRCh37/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.fasta \
  -L /mnt/storage/centers/ccgd/reference_files/baits/GRCh37/POPv3.1_S3050142/POPv3_AllTargetedRegions.interval_list \
  &> /data/20251218_pon_db_popv3_add_SO31074a.log &
```
6. Run GATK CreateSomaticPanelOfNormals. 
```bash
nohup /opt/jdk-17.0.12/bin/java -jar /opt/gatk/4.6.1.0/gatk-package-4.6.1.0-local.jar \
  CreateSomaticPanelOfNormals \
  -R /mnt/storage/centers/ccgd/reference_files/igenomes/Homo_sapiens/igenomes/GATK/GRCh37/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.fasta \
  -V gendb://pon_db_popv3_add_SO31074a_5 \
  --create-output-variant-index \
  -O mutect2_pon_popv3.1_20251218.vcf.gz \
  &> /data/20251218_pon_create_popv3_add_SO31074a.log &
```
7. On your desktop, run the following to decrease the memory of the nextflow4 server. Then, **reboot** nextflow4.
```
gcloud compute instances set-machine-type nextflow4 \
    --zone us-central1-a --machine-type e2-custom-4-20480
```
8. Again refer to step 4 to remount the drives.
9. Create the CNVKit reference.cnn.
```bash
docker run \
    -v /mnt/storage/centers/ccgd/Illumina/pipeline_sarek/SE81226a-ULPWGS/variant_calling/cnvkit:/mnt/storage/centers/ccgd/Illumina/pipeline_sarek/SE81226a-ULPWGS/variant_calling/cnvkit \
    -v /mnt/storage/centers/ccgd/reference_files/igenomes/Homo_sapiens/igenomes/GATK/GRCh37/Sequence/WholeGenomeFasta:/mnt/GRCh37 \
    -it quay.io/biocontainers/cnvkit:0.9.10--pyhdfd78af_0   bash
cnvkit.py reference [*coverage.cnn] -f ucsc.hg19.fa -o Reference.cnn

cd /mnt/storage/centers/ccgd/Illumina/pipeline_sarek/SE81226a-ULPWGS/variant_calling/cnvkit
cnvkit.py reference \
    -f /mnt/GRCh37/human_g1k_v37_decoy.fasta \
    CEPH1408-P1-H11_S148/*coverage.cnn \
    CEPH1408-P1-H1_S73/*coverage.cnn \
    CEPH1408-P2-G9_S66/*coverage.cnn \
    CEPH1408-P2-H1_S7/*coverage.cnn \
    -o 20250119_ulpwgs_se81226a_cnvkit_reference.cnn
```
10. Copy the reference.cnn file, and mutect2_pon.vcf.gz to the bait folder.
```bash
cp [reference.cnn] /mnt/storage/centers/ccgd/reference_files/baits/GRCh37/[bait_set]/
cp [mutect2_pon.vcf.gz] /mnt/storage/centers/ccgd/reference_files/baits/GRCh37/[bait_set]/
```
11. Copy the sample fasta files, sample target/antitarget cnn files, and sample vcf files in case you want to add
   to the PON later or run on a different pipeline (you can run the fastas through the new pipeline).
```
gcloud storage cp "*.fastq.gz" /mnt/storage/centers/ccgd/reference_files/baits/GRCh37/[bait_set]/pon_members/fastqs/
gcloud storage cp $(find . -name "*.cnn") /mnt/storage/centers/ccgd/reference_files/baits/GRCh37/[bait_set]/pon_members/cnn/
gcloud storage cp $(find . -name "*.mutect2.vcf.gz") /mnt/storage/centers/ccgd/reference_files/baits/GRCh37/[bait_set]/pon_members/vcf/
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