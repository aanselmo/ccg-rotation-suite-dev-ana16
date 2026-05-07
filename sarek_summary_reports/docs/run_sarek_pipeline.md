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
PROJECT_CODE=[your_project_code]
mkdir /data/${PROJECT_CODE}
cd /data/${PROJECT_CODE}  # When you run nextflow, it adds a hidden .nextflow folder.
tmux  # start new session
# Make sure to add the bait profile
/opt/nextflow/25.10.2/nextflow-25.10.2-dist run  /opt/sarek/main.nf \
    -profile docker,google,ccg,[bait_profile] \
    --input gs://ccg-bioinformatics/mnt/storage/centers/ccgd/Illumina/pipeline_sarek/${PROJECT_CODE}/csv/${PROJECT_CODE}_sarek_samplesheet.csv \
    --outdir gs://ccg-bioinformatics/mnt/storage/centers/ccgd/Illumina/pipeline_sarek/${PROJECT_CODE} \
    -w gs://ccg-bioinformatics/mnt/scratch/ccgd/global/nextflow_execution/${PROJECT_CODE}
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
cd /data
PROJECT_CODE=[Your project code]
/opt/nextflow/25.10.2/nextflow-25.10.2-dist run /opt/somalier/main.nf \
  -c /opt/somalier/nextflow.config \
  -profile docker,google \
  --sample_sheet gs://ccg-bioinformatics/mnt/storage/centers/ccgd/Illumina/pipeline_sarek/${PROJECT_CODE}/csv/recalibrated.csv \
  --outdir /mnt/storage/centers/ccgd/Illumina/pipeline_sarek/${PROJECT_CODE} \
  -w gs://ccg-scratch/mnt/scratch/ccgd/global/nextflow_execution/${PROJECT_CODE}
```
2. Create a groups_file.csv if necessary. The command to generate the groups file is the following:
```bash
/opt/ccgd-rotation-suite/sarek_summary_reports/create_somalier_groups.py -p ${PROJECT_CODE} \
    -r /home/nap24/Documents/Pipeline/${PROJECT_CODE}/recalibrated.csv
```
There is no need a groups_file.csv if all samples derive from different
organisms/patients/individuals. The structure of the groups_file is mentioned in the
[somalier relate documentation](https://github.com/brentp/somalier/blob/v0.2.19/README.md#somalier-relate).
Basically, comma seperated groups by line:
```
normal1,tumor1a,tumor1b
normal2,tumor2a
```
3. Run Somalier relate:
```bash
docker run \
    -v /mnt/storage/centers/ccgd/Illumina/pipeline_sarek/${PROJECT_CODE}/reports/somalier:/mnt/storage/centers/ccgd/Illumina/pipeline_sarek/${PROJECT_CODE}/reports/somalier \
    -w /mnt/storage/centers/ccgd/Illumina/pipeline_sarek/${PROJECT_CODE}/reports/somalier \
    -it quay.io/biocontainers/somalier:0.2.19--h0c29559_0 bash
somalier relate [--groups groups_file.csv] $(find . -name "*.somalier")
```

## Generate PON
1. Run Sarek. Pretty much the same except use the /opt/sarek-create-pon directory, create_pon profile, and make sure no
PON is set. Use the create_pon flag as the last profile to make sure the PON is not set. I set Sarek to use the Broad
PON by default so that's why you have to do this. Also, go into the `[bait_profile]` and make sure no PON or
reference.cnn is set there.
```bash
/opt/nextflow/25.10.2/nextflow-25.10.2-dist run  /opt/sarek-create-pon/main.nf \
    -profile docker,google,ccg,[bait_profile],create_pon \
    --input gs://ccg-bioinformatics/mnt/storage/centers/ccgd/nap24/pipeline_rnd/pipeline_sarek/SCB1046d_pon_popv3/csv/recalibrated.csv \
    --outdir gs://ccg-bioinformatics/mnt/storage/centers/ccgd/nap24/pipeline_rnd/pipeline_sarek/SCB1046d_pon_popv3 \
    -w gs://ccg-bioinformatics/mnt/scratch/ccgd/global/nextflow_execution/SCB1046d_pon_popv3
```
2. Run bcftools norm to normalize variants into a single letter per line.
```bash
BUCKET_DIR="/mnt/storage/centers/ccgd/Illumina/pipeline_sarek/[your bucket folder]"
GENOME="/mnt/storage/centers/ccgd/reference_files/igenomes/Homo_sapiens/igenomes/GATK/GRCh37/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.fasta"
for file in $(find ${BUCKET_DIR} -name "*.mutect2.vcf.gz")
  do split_file=${file/mutect2.vcf.gz/split.vcf.gz}
  vt_file=${file/mutect2.vcf.gz/vt.vcf}
  dedup_file=${file/mutect2.vcf.gz/dedup.vcf.gz}
  echo ${split_file}
  /opt/bcftools/1.23/bin/bcftools norm --multiallelics -any --output-type z --write-index=tbi --output ${split_file} ${file}
  /opt/vt/vt decompose_blocksub ${split_file} | /opt/vt/vt normalize -r ${GENOME} - > ${vt_file}
  /opt/bcftools/1.23/bin/bcftools norm --rm-dup all --output-type z --output ${dedup_file} --write-index=tbi ${vt_file}
done
```
3. Start the server used to create the GATK PON and login:
```bash
gcloud compute instances create \
  --zone=us-central1-a \
  --source-machine-image=nextflow4image \
  --machine-type n2-highmem-8 \
  gatkpon3
gcloud compute instances start --zone "us-central1-a" --project "ccgd-174313" gatkpon3
gcloud compute ssh --zone "us-central1-a" --project "ccgd-174313" "gatkpon3"
# Remount the fstab entry because it often doesn't load correctly
sudo umount -l /mnt/storage
sudo mount /mnt/storage
```
```
4. Create an args file for all of the vcf files. Head of the file looks like. It contains args to pass to
GenomicsDBImport.
Exmaple .args file:
```bash
--variant gs://ccg-bioinformatics/mnt/storage/centers/ccgd/nap24/pipeline_rnd/pipeline_sarek/SCB1046d_pon_popv3_2/variant_calling/mutect2/OP-P1-A6-19-136b/OP-P1-A6-19-136b.bcftools.vcf.gz
--variant gs://ccg-bioinformatics/mnt/storage/centers/ccgd/nap24/pipeline_rnd/pipeline_sarek/SCB1046d_pon_popv3_2/variant_calling/mutect2/OP-P1-B6-19-123b/OP-P1-B6-19-123b.bcftools.vcf.gz
--variant gs://ccg-bioinformatics/mnt/storage/centers/ccgd/nap24/pipeline_rnd/pipeline_sarek/SCB1046d_pon_popv3_2/variant_calling/mutect2/OP-P1-C6-19-138b/OP-P1-C6-19-138b.bcftools.vcf.gz
```
5. Run GenomicsDBImport. This typically takes 1 day to run the GenomicsDBImport.
See the Broad website for more info: https://gatk.broadinstitute.org/hc/en-us/articles/30332078729883-CreateSomaticPanelOfNormals-BETA
An article about tuning GenomicsDBImport: https://gatk.broadinstitute.org/hc/en-us/articles/360056138571-GenomicsDBImport-usage-and-performance-guidelines
Useful params:
    * `--batch-size` number of samples to process at a time.
    * `--max-num-intervals-to-import-in-parallel` number of intervals in the bed file to process at a time.
    * `--merge-input-intervals` merges all the intervals on 1 contig (i.e. chromosome). Shouldn't matter if we merge it.
      I'd want to filter out the off-target germline SNVs variants anyway (i.e. keep them in the PON).
Unused params:
    * `--bypass-feature-reader` makes GenomicsDBImport faster. It skips the checks that GATK does, and reads with
      bcftools instead. Should be fine if you did the bcftools and vt normalization steps earlier. I don't use
      it because it requires copying google storage files to local disk. It's simpler to run it without it.
    * `--reader-threads` will only work for 1 interval. I have no idea why there would only be 1 interval in a bait bed
      file. Typically there are thousands. That's why I don't use it.
```bash
# Note: Max ram should be 80% of VM ram. Some should be leftover for c libraries (not just java).
nohup /opt/jdk-17.0.12/bin/java -Xmx50G -Xms3G -jar /opt/gatk/4.6.1.0/gatk-package-4.6.1.0-local.jar \
  GenomicsDBImport \
  --merge-input-intervals \
  --max-num-intervals-to-import-in-parallel 10 \
  --batch-size 50 \
  --tmp-dir /data/tmp \
  -R /mnt/storage/centers/ccgd/reference_files/igenomes/Homo_sapiens/igenomes/GATK/GRCh37/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.fasta \
  --arguments_file /mnt/storage/centers/ccgd/Illumina/pipeline_sarek/SureSelect_XTHS_Human_All_Exon_V8_S33266340_create_pon/[your args].args \
  --genomicsdb-workspace-path /data/pon_db_sureselect_v8_20260125_5 \
  -L /mnt/storage/centers/ccgd/reference_files/baits/GRCh37/SureSelect_XTHS_Human_All_Exon_V8_S33266340/SureSelect_XTHS_Human_All_Exon_V8_S33266340_Covered_20260122.bed \
  &> /data/20260125_pon_db_sureselect_v8_6.log &
```
6. Run GATK CreateSomaticPanelOfNormals. 
```bash
nohup /opt/jdk-17.0.12/bin/java -jar /opt/gatk/4.6.1.0/gatk-package-4.6.1.0-local.jar \
  CreateSomaticPanelOfNormals \
  -R /mnt/storage/centers/ccgd/reference_files/igenomes/Homo_sapiens/igenomes/GATK/GRCh37/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.fasta \
  --create-output-variant-index \
  -V gendb://pon_db_sureselect_v8_20260125_5 \
  -O /mnt/storage/centers/ccgd/Illumina/pipeline_sarek/[your_project]/mutect2_pon_popv3.1_20251218.vcf.gz \
  &> /data/20250125_pon_create_popv3_add_SO31074a.log &
```
7. Shut down the gatkpon compute node. Note that the GenomicsDB will be lost if this is run; it shouldn't matter
because it's a temporary database.
```bash
gcloud compute instances stop --zone "us-central1-a" --project "ccgd-174313" gatkpon3
gcloud compute instances delete --zone=us-central1-a --delete-disks=all gatkpon3
```
8. Again refer to step 4 to remount the drives.
9. Create the CNVKit reference.cnn.
```bash
docker run \
    -v /mnt/storage/centers/ccgd/Illumina/pipeline_sarek/SE81226a-ULPWGS/variant_calling/cnvkit:/mnt/storage/centers/ccgd/Illumina/pipeline_sarek/SE81226a-ULPWGS/variant_calling/cnvkit \
    -v /mnt/storage/centers/ccgd/reference_files/igenomes/Homo_sapiens/igenomes/GATK/GRCh37/Sequence/WholeGenomeFasta:/mnt/GRCh37 \
    -it quay.io/biocontainers/cnvkit:0.9.10--pyhdfd78af_0   bash
cd /mnt/storage/centers/ccgd/Illumina/pipeline_sarek/SE81226a-ULPWGS/variant_calling/cnvkit
# Format
cnvkit.py reference -f ucsc.hg19.fa -o Reference.cnn [*coverage.cnn]
# Example Command
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
PROJECT_CODE=[your project]
BAIT_FOLDER=[your bait folder]
cp [reference.cnn] /mnt/storage/centers/ccgd/reference_files/baits/GRCh37/${BAIT_FOLDER}/
cp [mutect2_pon.vcf.gz] /mnt/storage/centers/ccgd/reference_files/baits/GRCh37/${BAIT_FOLDER}/
```
11. Edit the sarek config file that you created for the bait set [earlier](./README.md#convert-riker-interval_list-to-sarek-bed-file)
Add the pon, pon_tbi, and cnvkit_reference.
12. Copy the sample fasta files, sample target/antitarget cnn files, and sample vcf files in case you want to add
   to the PON later or run on a different pipeline (you can run the fastas through the new pipeline).
```bash
gcloud storage cp "gs://ccg-bioinformatics/mnt/storage/centers/ccgd/Illumina/pipeline_sarek/${PROJECT_CODE}/demultiplex/*.fastq.gz" \
    gs://ccg-bioinformatics/mnt/storage/centers/ccgd/reference_files/baits/GRCh37/${BAIT_FOLDER}/pon_members/fastqs/
gcloud storage cp "gs://ccg-bioinformatics/mnt/storage/centers/ccgd/Illumina/pipeline_sarek/${PROJECT_CODE}/variant_calling/cnvkit/**coverage.cnn" \
    gs://ccg-bioinformatics/mnt/storage/centers/ccgd/reference_files/baits/GRCh37/${BAIT_FOLDER}/pon_members/cnvkit/
gcloud storage cp "gs://ccg-bioinformatics/mnt/storage/centers/ccgd/Illumina/pipeline_sarek/${PROJECT_CODE}/variant_calling/mutect2/**.dedup.vcf.gz" \
    gs://ccg-bioinformatics/mnt/storage/centers/ccgd/reference_files/baits/GRCh37/${BAIT_FOLDER}/pon_members/mutect2/
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