# Sarek Directions
These are directions for running the Sarek pipeline and it's reporting scripts.

Steps of running sarek are:
1. [Generate the Sarek Samplesheet](#1-generate-the-sarek-samplesheet)
2. [Run the sarek pipeline](#2-run-the-sarek-pipeline)
3. [Run somalier/Fingerprinting](#3-run-somalier--fingerprinting)
4. [Run reporting scripts](#4-run-the-reporting-scripts)

## 1. Generate the Sarek Samplesheet
An example of a samplesheet is as follows
(as described in [Sarek docs](https://nf-co.re/sarek/3.7.1/docs/usage/)). The patient an id for the individual/organism.
The status is 0 for normal and 1 for tumor. The sample is the same name. The lane is the lane on the instrument. If
you have multiple lanes for 1 sample, there will need to be multiple lines in the samplesheet for 1 sample with the
sample name repeating (each line has unique fastqs and lanes, however). The fastq_1 column is the R1 fastq.gz file,
and the fastq_2 is the R2 fastq.gz file.
```
patient,status,sample,lane,fastq_1,fastq_2
A001,1,A001_X101_49707-1-B1,L005,gs://ccg-bioinformatics/mnt/storage/centers/ccgd/nap24/pipeline_rnd/pipeline_sarek/demult_revcomp_0089/A001_X101_49707-1-B1_S14_R1_001.fastq.gz,gs://ccg-bioinformatics/mnt/storage/centers/ccgd/nap24/pipeline_rnd/pipeline_sarek/demult_revcomp_0089/A001_X101_49707-1-B1_S14_R2_001.fastq.gz
A008,1,A008_X108_58797-1-C1,L005,gs://ccg-bioinformatics/mnt/storage/centers/ccgd/nap24/pipeline_rnd/pipeline_sarek/demult_revcomp_0089/A008_X108_58797-1-C1_S13_R1_001.fastq.gz,gs://ccg-bioinformatics/mnt/storage/centers/ccgd/nap24/pipeline_rnd/pipeline_sarek/demult_revcomp_0089/A008_X108_58797-1-C1_S13_R2_001.fastq.gz
```
```bash
gcloud storage cp ${PROJECT_CODE}_sarek_samplesheet.csv gs://ccg-bioinformatics/mnt/storage/centers/ccgd/Illumina/pipeline_sarek/${PROJECT_CODE}/csv/
```
🚧 **Under construction.** 🚧 I'm hoping to create something like the python script created by Edwin to run the
Scaffron pipeline. For now, you can use Edwin's [Quick Search](https://ccgd2.dfci.harvard.edu/minilims/html/quick_search.php)
to match up the index/barcode with the MBCF fastq.

![minilims_quick_search](./minilims_quick_search.png)

## 2. Run the Sarek pipeline
Please see the [run pipeline directions](./run_sarek_pipeline.md)

## 3. Run Somalier / Fingerprinting
Please see the [run somalier directions](./run_sarek_pipeline.md#run-somalierfingerprinting)

## 4. Run the reporting scripts
Please see the [run reporting scripts directions](./run_reporting.md)

## Edge Cases
### Regenerate Old Bait or Old PON
If you're running a bait set for the first time through Sarek, it's a good idea to generate a PON before the actual run
you need to do or at least generate the bed file.

#### Convert Riker .interval_list to Sarek .bed file
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
/opt/ccgd-rotation-suite/sarek_summary_reports/convert_to_sarek_bed.py -i [input.interval_list]
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

#### Regenerate Old PON
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

### Rerun old Riker Project
If the project has never been run on the bait before, you will probably want to regenerate a PON or at least the .bed
file. See the directions [above](#regenerate-old-bait-or-old-pon)
1. Restore the old run to argos.
```bash
ssh rc-stor12.dfci.harvard.edu
cd /data/backups/mnt/storage/centers/ccgd/Illumina/Basecalls/
nohup rsync -r 200127_SN218_1093_BHCCKWBCX3 /mnt/storage/centers/ccgd/Illumina/Basecalls/ &> ~/log/rsync_restore_200127_SN218_1093_BHCCKWBCX3.txt &
```
2. Add the following line at the top of the bcl2fastq_samplesheet.csv:
```
[Data],,,
```
3. Run this query:
```sql
SELECT SUBSTRING_INDEX(pl.platform_location_description, "_", -1) AS Lane,
rec.sample_collaborator_name AS Sample_ID,
first_barcode.barcode_sequence AS `index`, second_barcode.barcode_sequence AS `index2`,
plat.second_index_direction
FROM project proj
JOIN received_sample rec ON proj.project_id  = rec.project_id 
JOIN illumina_lane_member ilm ON ilm.`sample_id` = rec.`received_sample_id`
JOIN flowcell_content fc ON fc.flowcell_content_id = ilm.flowcell_content_id 
JOIN instrument_run ir ON ir.instrument_run_id = fc.instrument_run_id 
JOIN platform_location pl ON fc.platform_location_id = pl.platform_location_id
JOIN platform plat ON plat.platform_id = pl.platform_id 
JOIN library lib ON lib.sample_id  = rec.received_sample_id  
JOIN seq_barcode_set_member bsm ON lib.barcode_set_member_id = bsm.`seq_barcode_set_member_id`
JOIN seq_barcode first_barcode ON first_barcode.seq_barcode_id = bsm.first_seq_barcode_id
JOIN seq_barcode second_barcode ON second_barcode.seq_barcode_id = bsm.second_seq_barcode_id
WHERE proj.project_code = "[old riker project]"
/* Run specific PON samples to quicken Sarek's runtime later.
AND rec.sample_collaborator_name IN (
  [FROM FIND COMMAND]
)*/
ORDER BY pl.platform_location_description, rec.sample_collaborator_name ;
```
4. If the `second_index_direction` comes out as REVCOMP, please use https://arep.med.harvard.edu/labgc/adnan/projects/Utilities/revcomp.html
   to find the reverse complement of index2. You can just copy and paste the column out of DBeaver. Rows should be
   seperated by new lines. The site will figure out what you mean.
5. Make the output directory, and upload the samplesheet file to argos:
```bash
ssh argos-ccgd.dfci.harvard.edu
cd /mnt/storage/centers/ccgd/Illumina/pipeline_sarek/
PROJECT_CODE=[YOUR PROJECT CODE]
mkdir -p ${PROJECT_CODE}/csv ${PROJECT_CODE}/demultiplex
```
6. On your desktop run the following:
```bash
PROJECT_CODE=[YOUR PROJECT CODE]
touch ${PROJECT_CODE}_bcl2fastq_samplesheet.csv
libreoffice --calc ${PROJECT_CODE}_bcl2fastq_samplesheet.csv
# Right click -> Paste -> Paste Special -> Use Text Import Dialog to paste header:
[Data],,,
Lane,Sample_ID,index,index2
# Paste query data
# Save as -> Check mark Edit Filter Settings -> OK
ssh argos-ccgd.dfci.harvard.edu
PROJECT_CODE=[your project code]
mkdir -p /mnt/storage/centers/ccgd/Illumina/pipeline_sarek/${PROJECT_CODE}/csv/
# On your desktop
scp ${PROJECT_CODE}_bcl2fastq_samplesheet.csv argos-ccgd:/mnt/storage/centers/ccgd/Illumina/pipeline_sarek/${PROJECT_CODE}/csv/
```
Example:
```
[Data],,,
Lane,Sample_ID,index,index2
1,KM12-C2,AAGGACAC,TGTTGTGG
1,PON_Blood6-C3,TGCAGGTA,CATCGTGA
1,CEPH1408-H1,CGATAGAG,ATGGCGAA
1,PON_Blood1-G2,GTAGAGCA,CGGTCATA
```
7. Run bcl2fastq on argos. Example commands:
```bash
screen
qrsh -pe pvm 16
PROJECT_CODE=[YOUR PROJECT CODE]
# Add `--use-bases-mask Y98,I8n*,I8,Y98` if needed
/mnt/storage/centers/ccgd/apps/bcl2fastq/2.20.0.422/bin/bcl2fastq \
  --create-fastq-for-index-reads --no-lane-splitting --ignore-missing-positions --ignore-missing-controls \
  --ignore-missing-filter --ignore-missing-bcls -r 6 -w 6  \
  -R /mnt/storage/centers/ccgd/Illumina/Basecalls/[your illumina run] \
  --output-dir=/mnt/storage/centers/ccgd/Illumina/pipeline_sarek/${PROJECT_CODE}/demultiplex \
  --interop-dir=/tmp \
  --sample-sheet=/mnt/storage/centers/ccgd/Illumina/pipeline_sarek/${PROJECT_CODE}/csv/${PROJECT_CODE}_bcl2fastq_samplesheet.csv
```
8. Upload to gcloud after done demultiplexing. The fastq.gz files usually should be 700 MB to 1 GB each. If they aren't
something's wrong. Try switching to reverse complement of the 2nd barcode or just regular if you already tried the
reverse complement. The database isn't 100% accurate for the REVCOMP because chemistry on instruments can change.
```bash
screen ls
screen -d -r [old session for demultiplexing]
gcloud storage rsync -r --exclude=".*Undetermined.*" \
    /mnt/storage/centers/ccgd/Illumina/pipeline_sarek/${PROJECT_CODE} \
    gs://ccg-bioinformatics/mnt/storage/centers/ccgd/Illumina/pipeline_sarek/${PROJECT_CODE}
```
9. Run the following to get the fastqs to add to the samplesheet:
```bash
# For fastq_1 column
gcloud storage ls "gs://ccg-bioinformatics/mnt/storage/centers/ccgd/Illumina/pipeline_sarek/${PROJECT_CODE}/demultiplex/*R1*"
# For fastq_2 column
gcloud storage ls "gs://ccg-bioinformatics/mnt/storage/centers/ccgd/Illumina/pipeline_sarek/${PROJECT_CODE}/demultiplex/*R2*"
```
10. Grab the other columns with a query.
    * **Nota bene:** You might only need to grab one lane if you used `--no-lane-splitting` on the demultiplexing step.
    * **Nota bene bis:** Status should be 1 for all normals if creating a PON. Normals should be treated like tumor
      samples when run through Sarek (when creating a PON).
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

### Demultiplex Newer Project
If you get the following error from bcl2fastq, that may mean you have a newer NovaSeq run,
and you will need to use bcl-convert. Here's how you run it.
```
No cbcl files found for cycle 3. Attempt 3 of 3. Corrupt or missing .cbcl files were found for 3 consecutive cycles, verify the source of the data.
```
1. Run the [query in step 3](#rerun-old-riker-project) or find some way to get the
   data section of the samplesheet with the same headers `Lane,Sample_ID,index,index2`.
   **Nota bene** unlike for bcl2fastq index2 should always be in forward orientation.
   (not reverse compliment.)
2. Add the following sections to the samplesheet to follow the v2 specification.
   The BCLConvert_Settings section is optional. Other settings are described in the 
   [illumina docs](https://support-docs.illumina.com/SW/BCL_Convert/Content/SW/BCLConvert/SampleSheets_swBCL.htm)
   **Nota bene** I2 cycles to ignore should be specified in reverse order if the i2 is reverse (i.e. N8I8... not I8N8)
```
[Header],,,
FileFormatVersion,2,,
[BCLConvert_Settings],,,
OverrideCycles,R1:Y151;R2:Y151;I1:I8U9;I2:N8I8,,
BarcodeMismatchesIndex1,1,,
[BCLConvert_Data],,,
Lane,Sample_ID,index,index2

your query data...
```
3. Upload data to gcloud storage.
```bash
ssh argos-ccgd.dfci.harvard.edu
screen
qrsh -pe pvm 16
cd /mnt/storage/centers/ccgd/Illumina/Basecalls
gcloud storage cp -r [your run] gs://ccg-bioinformatics/mnt/storage/centers/ccgd/Illumina/Basecalls/
```
4. Start the server used to run bcl-convert and login:
```bash
# n2d-custom-48-97280 is 48 cores, 95 * 1024 MB RAM
gcloud compute instances create \
  --zone=us-central1-a \
  --source-snapshot=https://compute.googleapis.com/compute/v1/projects/ccgd-174313/global/snapshots/nextflow4-20260125 \
  --machine-type n2d-custom-48-97280 \
  --local-ssd=device-name=nvme1,interface=NVME \
  --local-ssd=device-name=nvme2,interface=NVME \
  --local-ssd=device-name=nvme3,interface=NVME \
  --local-ssd=device-name=nvme4,interface=NVME \
  --local-ssd=device-name=nvme5,interface=NVME \
  --local-ssd=device-name=nvme6,interface=NVME \
  --local-ssd=device-name=nvme7,interface=NVME \
  --local-ssd=device-name=nvme8,interface=NVME \
  bclconvert3
gcloud compute ssh --zone "us-central1-a" --project "ccgd-174313" "bclconvert3"
# Remount the fstab entry because it often doesn't load correctly
sudo umount -l /mnt/storage
sudo mount /mnt/storage
# Create a RAID disk setup for temporary saving
sudo apt install mdadm
sudo mdadm --create /dev/md0 \
  --level=0 \
  --raid-devices=8 \
  /dev/nvme0n1 /dev/nvme0n2 /dev/nvme0n3 /dev/nvme0n4 \
  /dev/nvme0n5 /dev/nvme0n6 /dev/nvme0n7 /dev/nvme0n8
sudo mkfs.ext4 /dev/md0
sudo mkdir /mnt/localssd
sudo mount /dev/md0 /mnt/localssd
sudo chown -R :google-sudoers /mnt/localssd
sudo chmod g+ws /mnt/localssd
cd /mnt/localssd
# copy files over
gcloud storage cp -r gs://ccg-bioinformatics/mnt/storage/centers/ccgd/Illumina/Basecalls/[your run] .
```
5. Run docker and then bcl-convert:
Docker image was made from: https://github.com/nf-core/demultiplex/blob/1.7.0/modules/nf-core/bclconvert/Dockerfile
```bash
tmux
docker run \
    -v /mnt/localssd:/mnt/localssd \
    -v /mnt/storage/centers/ccgd/Illumina/pipeline_sarek/SE81226a/csv:/mnt/storage/centers/ccgd/Illumina/pipeline_sarek/SE81226a/csv \
    -it us-central1-docker.pkg.dev/ccgd-174313/sarek/bclconvert:4.4.6 bash

# Don't run with --no-lane-splitting. It crashes if there are no reads from a sample in one of the lanes.
bcl-convert \
    --bcl-input-directory /mnt/localssd/20260122_LH00547_0196_B23GNT7LT3_JN13826-L001-L002-L003_AB13838-L007 \
    --sample-sheet /mnt/storage/centers/ccgd/Illumina/pipeline_sarek/SE81226a/csv/SE81226a_run196_bcl2fastq_samplesheet.csv \
    --output-directory /mnt/localssd/demultiplex_run196 \
    --fastq-gzip-compression-level 6
```
6. Upload fastqs back to gcloud storage.
```bash
gcloud storage rsync -r --exclude=".*Undetermined.*" /mnt/localssd/demultiplex_run196 gs://ccg-bioinformatics/mnt/storage/centers/ccgd/Illumina/pipeline_sarek/SE81226a/demultiplex_run196
```
7. Delete the old image:
```bash
gcloud compute instances stop --zone "us-central1-a" --project "ccgd-174313" --discard-local-ssd=true bclconvert3
gcloud compute instances delete --zone=us-central1-a --delete-disks=all bclconvert3
```