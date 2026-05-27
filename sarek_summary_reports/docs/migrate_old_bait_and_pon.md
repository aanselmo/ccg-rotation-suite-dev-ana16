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

