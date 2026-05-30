# Regenerate Old Bait and PON
If you're running a bait set for the first time through Sarek, it's a good idea to generate a PON before the actual run
you need to do or at least generate the bed file.

## Convert Riker .interval_list to Sarek .bed file
1. Sarek runs on .bed files (not GATK .interval_list files). Convert the interval list to a .bed file.
```bash
ssh argos-ccgd.dfci.harvard.edu
qrsh -pe pvm 3
less /mnt/storage/centers/ccgd/apps/Riker/latest/Riker.ccgd.clouddb.config
# pick the one that's used for full_targets in Riker
# (https://wiki.dfci.harvard.edu:8443/pages/viewpage.action?pageId=60738749)
```
2. Upload to gcloud.
```bash
gcloud storage cp [path from riker.config] gs://ccg-bioinformatics/[path from riker.config]
```
3. Run the convert_to_sarek_bed.py tool. It does the following:
    * Converts the .interval_list to .bed
    * Sorts the chromosomes numerically (1,2,3,etc.) instead of alphabetically (1,10,11,etc.).
    * Merges overlapping regions.
    * Annotates the regions in a way that CNVKit expects.
```bash
/opt/ccg-rotation-suite/sarek/convert_to_sarek_bed.py -i [input.interval_list]
```
4. Add the .bed to Sarek config.
```bash
cp /opt/sarek/conf/ccg/grch37/baits/POPv3.1_S3050142.config /opt/sarek/conf/ccg/grch37/baits/${BAIT_NAME}.config
vim /opt/sarek/conf/ccg/grch37/baits/${BAIT_NAME}.config  # Change POPv3 beds to your new bed
# Delete the pon, pon_tbi, and cnvkit_reference if you don't have a PON, yet.
vim /opt/sarek/conf/ccg/ccg.config
# Add the bait set as a new profile. Note: periods (.) are not allowed in profile names in nextflow. I would
# use the same bait name next to the supplier id for the profile name though.
cp -r /opt/sarek/conf/ccg /opt/sarek-create-pon/conf/
```

## Regenerate Old PON
1. Look on argos-ccgd for samples in the old PON.
```bash
ssh argos-ccgd.dfci.harvard.edu
cd /mnt/storage/centers/ccgd/bmw35/RobustCNV_PON
cd [YOUR PON]
cd coverage_and_calls
ls -tlh *.dedup.cleaned.txt
# Try to find files/samples in the last 5 years. Ideally, you find a instrument run that has a lot of normals so you
# don't have to do a lot of Sarek runs, either.
find . -name "*.dedup.cleaned.txt" -ctime -1825 \
    -exec sh -c "basename {} .dedup.cleaned.txt | sed -r 's/_L_[0-9]{6}//g' | sed 's/^\|$/\"/g'" \; | tr '\n' ', '
```
2. Query to find what instrument runs the samples were on.
```sql
SELECT rs.sample_collaborator_name, p.project_code, ir.run_name, pl.platform_location_description, bs.bait_set_name 
FROM received_sample rs
JOIN project p ON p.project_id = rs.project_id
JOIN illumina_lane_member ilm ON ilm.sample_id = rs.received_sample_id
JOIN capture cap ON ilm.capture_id = cap.capture_id
JOIN bait_set bs ON cap.bait_set_id = bs.bait_set_id
JOIN flowcell_content fc ON fc.flowcell_content_id = ilm.flowcell_content_id
JOIN instrument_run ir ON ir.instrument_run_id = fc.instrument_run_id
JOIN platform_location pl ON pl.platform_location_id  = fc.platform_location_id 
WHERE rs.sample_collaborator_name IN (
  [FROM FIND COMMAND]
) AND bs.bait_set_name = "Your bait set";
```
3. Restore the old Riker runs and proceed to [rerun old riker project directions](#rerun-old-riker-project)
4. Proceed to [Generate PON](./run_sarek_pipeline.md#generate-pon)


## Create Sarek Samplesheet After Demultiplexing
9. Run the following to get the fastqs to add to the samplesheet:
```bash
# For fastq_1 column
gcloud storage ls "gs://ccg-bioinformatics/mnt/storage/centers/ccgd/Illumina/pipeline_sarek/${PROJECT_CODE}/demultiplex/*R1*"
# For fastq_2 column
gcloud storage ls "gs://ccg-bioinformatics/mnt/storage/centers/ccgd/Illumina/pipeline_sarek/${PROJECT_CODE}/demultiplex/*R2*"
```
10. Grab the other columns with a query.
    * **Nota bene:** You might only need to grab one lane if you used `--no-lane-splitting` on the demultiplexing step.
```sql
SELECT i.individual_collaborator_name as patient,
IF(rec.genomic_source_id=2, 0, 1) AS status,
rec.sample_collaborator_name AS sample,
SUBSTRING_INDEX(pl.platform_location_description, "_", -1) AS lane
FROM project proj
JOIN received_sample rec ON proj.project_id  = rec.project_id
JOIN individual i ON i.individual_id = rec.individual_id 
JOIN illumina_lane_member ilm ON ilm.`sample_id` = rec.`received_sample_id`
JOIN flowcell_content fc ON fc.flowcell_content_id = ilm.flowcell_content_id 
JOIN instrument_run ir ON ir.instrument_run_id = fc.instrument_run_id 
JOIN platform_location pl ON fc.platform_location_id = pl.platform_location_id
WHERE proj.project_code = "[Your project code]"
/* Run specific PON samples to quicken Sarek's runtime later.
AND rec.sample_collaborator_name IN (
  [FROM FIND COMMAND]
)*/
ORDER BY rec.sample_collaborator_name;
```
11. Create a sarek_samplesheet.csv. See [above directions](#1-generate-the-sarek-samplesheet) for details.
12. Kickoff Sarek pipeline as in [run pipeline directions](./run_sarek_pipeline.md)

## Generate PON
1. Run Sarek. Pretty much the same except the following:
    * Use the create_pon flag as the last profile to make sure the PON is not set. I set Sarek to use the
      Broad PON by default so that's why you have to do this.
    * Go into the `[bait_profile]` and make sure no PON or reference.cnn is set there.
    * Status should be 1 for all normals if creating a PON. Normals should be treated like tumor
      samples when run through Sarek (when creating a PON). This is because HaplotypeCaller will run instead of Mutect2
      if the status is 0. You need Mutect2 to run in order to create a PON against Mutect2.
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
4. Create an args file for all of the vcf files. Head of the file looks like. It contains args to pass to
GenomicsDBImport.
Exmaple .args file:
```bash
--variant gs://ccg-bioinformatics/mnt/storage/centers/ccgd/nap24/pipeline_rnd/pipeline_sarek/SCB1046d_pon_popv3_2/variant_calling/mutect2/OP-P1-A6-19-136b/OP-P1-A6-19-136b.bcftools.vcf.gz
--variant gs://ccg-bioinformatics/mnt/storage/centers/ccgd/nap24/pipeline_rnd/pipeline_sarek/SCB1046d_pon_popv3_2/variant_calling/mutect2/OP-P1-B6-19-123b/OP-P1-B6-19-123b.bcftools.vcf.gz
--variant gs://ccg-bioinformatics/mnt/storage/centers/ccgd/nap24/pipeline_rnd/pipeline_sarek/SCB1046d_pon_popv3_2/variant_calling/mutect2/OP-P1-C6-19-138b/OP-P1-C6-19Process `NFCORE_SAREK:SAREK:FASTQ_PREPROCESS_GATK:FASTP` input file name collision -- There are multiple input files for each of the following file names: 4
9_10_KAPA3uM_14cyc_S29_L002_R1_001.fastq.gzt-usage-and-performance-guidelines
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
