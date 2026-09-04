# Demultiplex

## Demultiplex Old Riker Project
If the project has never been run on the bait before, you will probably want to regenerate a PON or at least the .bed
file. See the directions [above](#regenerate-old-bait-or-old-pon)
1. Restore the old run to argos.
```bash
ssh rc-stor12.dfci.harvard.edu
cd /data/backups/mnt/storage/centers/ccgd/Illumina/Basecalls/
nohup rsync -r 200127_SN218_1093_BHCCKWBCX3 /mnt/storage/centers/ccgd/Illumina/Basecalls/ &> ~/log/rsync_restore_200127_SN218_1093_BHCCKWBCX3.txt &
```
2. On your desktop run the following:
```bash
PROJECT_CODE=[YOUR PROJECT CODE]
touch ${PROJECT_CODE}_bcl2fastq_samplesheet.csv
libreoffice --calc ${PROJECT_CODE}_bcl2fastq_samplesheet.csv
# Right click -> Paste -> Paste Special -> Use Text Import Dialog to paste header:
[Data],,,
Lane,Sample_ID,index,index2
```
3. Run this query and paste it below the header:
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
5. In Libreoffice: Save as -> Check mark Edit Filter Settings -> OK
Example:
```
[Data],,,
Lane,Sample_ID,index,index2
1,KM12-C2,AAGGACAC,TGTTGTGG
1,PON_Blood6-C3,TGCAGGTA,CATCGTGA
1,CEPH1408-H1,CGATAGAG,ATGGCGAA
1,PON_Blood1-G2,GTAGAGCA,CGGTCATA
```
5. Make the output directory, and upload the samplesheet file to argos:
```bash
PROJECT_CODE=[YOUR PROJECT CODE]
# {csv,demultiplex} is an example of brace expansion. It makes 2 dirs.
ssh argos-ccgd.dfci.harvard.edu mkdir -p /mnt/storage/centers/ccgd/Illumina/pipeline_sarek/${PROJECT_CODE}/{csv,demultiplex}
scp ${PROJECT_CODE}_bcl2fastq_samplesheet.csv argos-ccgd:/mnt/storage/centers/ccgd/Illumina/pipeline_sarek/${PROJECT_CODE}/csv/
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

## Demultiplex Newer Project
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
  --zone=us-central1-f \
  --project=ccgd-174313 \
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
gcloud compute ssh --zone "us-central1-f" --project "ccgd-174313" "bclconvert3"
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
PROJECT_CODE=[project code]
docker run \
    -v /mnt/localssd:/mnt/localssd \
    -v /mnt/storage/centers/ccgd/Illumina/pipeline_sarek/${PROJECT_CODE}$/csv:/mnt/storage/centers/ccgd/Illumina/pipeline_sarek/${PROJECT_CODE}/csv \
    -it us-central1-docker.pkg.dev/ccgd-174313/sarek/bclconvert:4.4.6 bash

# Don't run with --no-lane-splitting. It crashes if there are no reads from a sample in one of the lanes.
PROJECT_CODE=[project code]
bcl-convert \
    --bcl-input-directory /mnt/localssd/20260122_LH00547_0196_B23GNT7LT3_JN13826-L001-L002-L003_AB13838-L007 \
    --sample-sheet /mnt/storage/centers/ccgd/Illumina/pipeline_sarek/${PROJECT_CODE}/csv/${PROJECT_CODE}_bclconvert_samplesheet.csv \
    --output-directory /mnt/localssd/${PROJECT_CODE}_fastqs \
    --fastq-gzip-compression-level 6
```
6. Upload fastqs back to gcloud storage.
```bash
gcloud storage rsync -r --exclude=".*Undetermined.*" /mnt/localssd/${PROJECT_CODE}_fastqs gs://ccg-bioinformatics/mnt/storage/centers/ccgd/Illumina/pipeline_sarek/${PROJECT_CODE}/preprocessing/fastq
```
7. Delete the old image:
```bash
gcloud compute instances stop --zone "us-central1-f" --project "ccgd-174313" --discard-local-ssd=true bclconvert3
gcloud compute instances delete --zone=us-central1-f --delete-disks=all bclconvert3
```