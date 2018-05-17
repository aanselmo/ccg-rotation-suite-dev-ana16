#!/usr/bin/env Rscript
#version 2.0.0

args <- commandArgs(trailingOnly=TRUE)

usage <- "Usage: generate_report.R <project code> <project path>\n\n\tgenerates breakmer reports in Ryans format\n\n"

if(length(args)!=2){
	cat((prompt=usage), "\n")
	q(save="no")
}

project_code <- args[1]
project_path <- args[2]

require("XLConnect")
require("RMySQL")

############# global settings ##############################

#where to place the temporary files
targetPath <- "/media/isilon/bmw35/analysis/breakmer/breakmer_reports/"
#where the projects folder is located
projectPath <- "/media/isilon/Illumina/pipeline_riker/Project/"


#Connect to riker using project code to retreive relevant info
riker <- dbConnect(MySQL(), user="ccgd", password="cc9d_u53r", dbname="riker_db_july2013", host="ccgd2", client.flag=CLIENT_MULTI_STATEMENTS)
riker_ver <- "riker_db_july2013"

#profile database, when necessary
#riker <- dbConnect(MySQL(), user="camd2", password="c4nc3r", dbname="riker_db_profile_production_may_2013", host="staging-camd.partners.org", client.flag=CLIENT_MULTI_STATEMENTS)
#riker_ver <- "riker_db_profile_production_may_2013"
	
#Connect to local test database:
#riker <- dbConnect(MySQL(), user="root", password="damn1tjan3t", dbname="riker_local", host="127.0.0.1", client.flag=CLIENT_MULTI_STATEMENTS)
#riker_ver <- "riker_local"

project_description <- dbGetQuery(riker, statement= paste("SELECT project_name FROM project WHERE project_code = '",project_code,"';",sep=""))

parent_project_id <- project_description <- dbGetQuery(riker, statement= paste("SELECT parent_project_id FROM project WHERE project_code = '",project_code,"';",sep=""))

parent_project_code <- project_description <- dbGetQuery(riker, statement= paste("SELECT project_code FROM project WHERE project_id = '",parent_project_id,"';",sep=""))

project_name <- project_description <- dbGetQuery(riker, statement= paste("SELECT project_name FROM project WHERE project_code = '",project_code,"';",sep=""))

projectId <- paste(project_path,paste(project_code,project_name,sep="_"),sep="")

fileName <- paste(project_path,paste(format(Sys.Date(),"%y%m%d"),project_name,project_code,sep="_"),sep="")

#cat(paste(projectId, "breakmer.df.txt", sep='_'))

#system(paste("python link_breakmer_files.rcprime.py ",sep=""))

#######################################################################################
#-------------------- begin Ryan's plot_results.R section ----------------------------#
#######################################################################################

svs <- read.table( paste(projectId, "breakmer.df.txt", sep='_'), header=T, sep="\t",as.is=T) 
tFreqs <- table(svs$sampleid, svs$target_name)

print(tFreqs)
targetGeneSummary <- data.frame("target_region" = colnames(tFreqs), "Nevents" = colSums(tFreqs), "Nsamples_with_event"=colSums(tFreqs>0))
write.table(targetGeneSummary, quote=F, row.names=T, col.names=T, sep="\t", file=paste(projectId, "breakmer.target_summary.txt", sep="_"))

tFreqs2 <- table(svs$sampleid, svs$rearr_type)

print(tFreqs2)
print(dim(tFreqs2))
if(ncol(tFreqs2)==2 && c("rearrangement", "indel") %in% colnames(tFreqs2)){
	sampleSummary <- data.frame("sample_summary" = rownames(tFreqs), "Nevents"=rowSums(tFreqs), "Nindels"=tFreqs2[,1], "Nrearr"=tFreqs2[,2])
	sample_header <- c("Sample Name", "Number of events detected", "Number of indels", "Number of rearrangements") 
}else if(ncol(tFreqs2)==1 && colnames(tFreqs2) == "rearrangement"){
	sampleSummary <- data.frame("sample_summary" = rownames(tFreqs), "Nevents"=rowSums(tFreqs), "Nrearr"=tFreqs2[,1]) 
	sample_header <- c("Sample Name", "Number of events detected", "Number of rearrangements")
}else if(ncol(tFreqs2)==3){
	sampleSummary <- data.frame("sample_summary" = rownames(tFreqs), "Nevents"=rowSums(tFreqs), "Nindels"=tFreqs2[,1], "Nrearr"=tFreqs2[,2],  "Nrearrinv"=tFreqs2[,3])
	sample_header <- c("Sample Name", "Number of events detected", "Number of indels", "Number of rearrangements", "Number of rearrangement inversions") 
}else{
	print("ERROROROROROROROORORORO!!!!!!!")
}
#, "Ninversions"=tFreqs2[,3])
write.table(sampleSummary, quote=F, row.names=T, col.names=T, sep="\t", file=paste(projectId, "breakmer.sample_summary.txt", sep="_"))

#Rearranged genes	No. samples	No. of events
rearr_idx <- which(svs$rearr_type != 'indel')
rearr_pair_freq <- table(svs$sampleid[rearr_idx], svs$gene_list[rearr_idx])
rearr_pair_freq <- colSums(rearr_pair_freq > 0)
rearr_pair_df <- data.frame("rearr_genes" = factor(names(rearr_pair_freq), levels=names(rearr_pair_freq)[order(rearr_pair_freq, decreasing=F)], ordered=T), "Nsample" = rep(NA,length(rearr_pair_freq)), "Nevents" = rep(NA,length(rearr_pair_freq)))
rearr_pair_df$Nsample <- unlist(lapply(rearr_pair_df$rearr_genes,function(rearr_gene){length(unique(svs$sampleid[which(svs$gene_list == rearr_gene)]))}))
rearr_pair_df$Nevents <- unlist(lapply(rearr_pair_df$rearr_genes,function(rearr_gene){length(svs$sampleid[which(svs$gene_list == rearr_gene)])}))
write.table(rearr_pair_df, quote=F, row.names=F, col.names=T, sep="\t", file=paste(projectId, "breakmer.rear_pair_diff.txt", sep="_"))

library("ggplot2")
#------------------------------------------------------------------------------
results_by_target_bar <- function(data) {

  rearr_idx <- which(data$rearr_type != 'indel')
  indel_idx <- which(data$rearr_type == 'indel')

  rearr_target_freq <- table(data$sampleid[rearr_idx], data$target_name[rearr_idx])
  indel_target_freq <- table(data$sampleid[indel_idx], data$target_name[indel_idx])

  rearr_target_freq <- colSums(rearr_target_freq > 0)
  indel_target_freq <- colSums(indel_target_freq>0)

  rearr_df <- data.frame("target_gene" = factor(names(rearr_target_freq), levels=names(rearr_target_freq)[order(rearr_target_freq, decreasing=F)], ordered=T), "rearr_freq" = as.numeric(rearr_target_freq))
  xl <- "Target captured gene name"
  yl <- "Number of samples with rearr event"
  gp <- ggplot(rearr_df[which(rearr_df[,2] > 5),], aes(x=target_gene, y=rearr_freq) ) +
               geom_bar(stat="identity") +
               labs(x=xl, y=yl, title="Target capture gene rearrangement freq (> 5 samples)") +
               theme_bw() +
               guides(fill=guide_legend(title=NULL)) +
               theme(text = element_text(size=18)) +
               theme(legend.text=element_text(size=18)) +
               theme(strip.text=element_text(size=18)) +
               coord_flip() +
               theme(axis.text=element_text(size=14), axis.title=element_text(size=16))
  ggsave(filename="target_rearr_freq_bar.pdf", plot=gp, width=14, height=12)

  indel_df <- data.frame("target_gene" = factor(names(indel_target_freq), levels=names(indel_target_freq)[order(indel_target_freq, decreasing=F)], ordered=T), "indel_freq" = as.numeric(indel_target_freq))
  xl <- "Target captured gene name"
  yl <- "Number of samples with indel event"
  gp <- ggplot(indel_df[which(indel_df[,2] > 5),], aes(x=target_gene, y=indel_freq) ) +
               geom_bar(stat="identity") +
               labs(x=xl, y=yl, title="Target gene indel freq (> 5 samples)") +
               theme_bw() +
               guides(fill=guide_legend(title=NULL)) +
               theme(text = element_text(size=18)) +
               theme(legend.text=element_text(size=18)) +
               theme(strip.text=element_text(size=18)) +
               coord_flip() +
               theme(axis.text=element_text(size=16), axis.title=element_text(size=16))
  ggsave(filename="target_indel_freq_bar.pdf", plot=gp, width=14, height=10)

  rearr_pair_freq <- table(data$sampleid[rearr_idx], data$gene_list[rearr_idx])
  rearr_pair_freq <- colSums(rearr_pair_freq > 0)
  rearr_pair_df <- data.frame("rearr_genes" = factor(names(rearr_pair_freq), levels=names(rearr_pair_freq)[order(rearr_pair_freq, decreasing=F)], ordered=T), "rearr_freq" = as.numeric(rearr_pair_freq))
  xl <- "Rearranged genes"
  yl <- "Number of samples with rearr pairs"
  gp <- ggplot(rearr_pair_df[which(rearr_pair_df[,2] > 2),], aes(x=rearr_genes, y=rearr_freq) ) +
               geom_bar(stat="identity") +
               labs(x=xl, y=yl, title="Gene rearrangements freq (> 5 samples)") +
               theme_bw() +
               guides(fill=guide_legend(title=NULL)) +
               theme(text = element_text(size=18)) +
               theme(legend.text=element_text(size=18)) +
               theme(strip.text=element_text(size=18)) +
               coord_flip() +
               theme(axis.text=element_text(size=14), axis.title=element_text(size=16))
  ggsave(filename="rearr_pair_freq_bar.pdf", plot=gp, width=14, height=12)

  rearrs <- data[rearr_idx, ]
  indels <- data[indel_idx, ]

  xl <- "Split read counts"
  yl <- "Discordant read counts"
  gp <- ggplot(rearrs, aes(x=log2(max_splitread_counts+1), y=log2(max_discread_counts+1), color=max_repeat_overlap) ) +
               geom_point(aes(size=1/min_uniq_align), alpha=0.7) +
               labs(x=xl, y=yl, title="Rearrangement evidence") +
               theme_bw() +
               guides(fill=guide_legend(title=NULL)) +
               theme(text = element_text(size=18)) +
               theme(legend.text=element_text(size=18)) +
               theme(strip.text=element_text(size=18)) +
               coord_flip() +
               theme(axis.text=element_text(size=14), axis.title=element_text(size=16))
  ggsave(filename="rearr_counts_scatter.pdf", plot=gp, width=14, height=12)

  xl <- "log2(Read counts)"
  yl <- "Contig length"
  gp <- ggplot(rearrs, aes(x=log2(max_splitread_counts+max_discread_counts), y=len_contig, color=max_repeat_overlap) ) +
               geom_point(aes(size=1/min_uniq_align), alpha=0.7) +
               labs(x=xl, y=yl, title="Rearrangement evidence vs. contig length") +
               theme_bw() +
               guides(fill=guide_legend(title=NULL)) +
               theme(text = element_text(size=18)) +
               theme(legend.text=element_text(size=18)) +
               theme(strip.text=element_text(size=18)) +
               coord_flip() +
               theme(axis.text=element_text(size=14), axis.title=element_text(size=16))
  ggsave(filename="rearr_counts_vs_len_scatter.pdf", plot=gp, width=14, height=12)

  xl <- "log2(Read counts)"
  yl <- "Contig length"
  gp <- ggplot(indels, aes(x=log2(max_splitread_counts), y=len_contig, color=subprojectid) ) +
               geom_point(size=4, alpha=0.7) +
               labs(x=xl, y=yl, title="Indel evidence vs. contig length") +
               theme_bw() +
               guides(fill=guide_legend(title=NULL)) +
               theme(text = element_text(size=18)) +
               theme(legend.text=element_text(size=18)) +
               theme(strip.text=element_text(size=18)) +
               coord_flip() +
               theme(axis.text=element_text(size=14), axis.title=element_text(size=16))
  ggsave(filename="indel_counts_vs_len_scatter.pdf", plot=gp, width=14, height=12)

  xl <- "Target genes"
  yl <- "Read count support"
  gp <- ggplot(rearrs, aes(x=target_name, y=log2(max_splitread_counts+max_discread_counts)) ) +
               geom_boxplot() +
               labs(x=xl, y=yl, title="Indel evidence vs. contig length") +
               theme_bw() +
               guides(fill=guide_legend(title=NULL)) +
               theme(text = element_text(size=18)) +
               theme(legend.text=element_text(size=18)) +
               theme(strip.text=element_text(size=18)) +
               coord_flip() +
               theme(axis.text=element_text(size=7), axis.title=element_text(size=16))
  ggsave(filename="rearr_target_counts_box.pdf", plot=gp, width=14, height=20)
}

#######################################################################################
#---------------------- end Ryan's plot_results.R section ----------------------------#
#######################################################################################


wb <- loadWorkbook("utils/BreaKmer/BreaKmer_summary_tab.xlsx")
createSheet(wb, "Results")
createSheet(wb, "Filtered out")
createSheet(wb, "Sample summary")
createSheet(wb, "Target summary")
createSheet(wb, "Rearr pair counts")
createSheet(wb, "circos links all")
createSheet(wb, "circos annos all")

# function for setting a range of column widths
# wraps setColumnWidth from XLConnect and just loops over
# the column variable for however many columns we need
setColumnWidths <- function(workbook, sheet,  columns, currWidth=-1){
  for(currCol in columns){
    setColumnWidth(workbook, sheet, column=currCol, width = currWidth)
  }
}

#putting the outputted files together

# hs <- createStyle(textDecoration = "BOLD")

df <- read.table(paste(projectId, "breakmer.df.txt", sep="_"),sep="\t",header=FALSE,stringsAsFactors=FALSE)
filtered <- read.table(paste(projectId, "breakmer.filtered.txt", sep="_"),sep="\t",header=TRUE,stringsAsFactors=FALSE)
breakmer <- read.table(paste(projectId, "breakmer.txt", sep="_"),sep="\t",header=TRUE,stringsAsFactors=FALSE,fill=T)
sample_summary <- read.table(paste(projectId, "breakmer.sample_summary.txt", sep="_"),sep="\t",header=TRUE,stringsAsFactors=FALSE)
circos_annos <-  try(read.table(paste(projectId, "circos_annos.txt", sep="_"),sep="\t",header=FALSE,stringsAsFactors=FALSE))
if (class(circos_annos) == "try-error") {
	circos_annos <- NULL
}
circos_links <- try(read.table(paste(projectId, "circos_links.txt", sep="_"),sep="\t",header=FALSE,stringsAsFactors=FALSE))
if (class(circos_links) == "try-error") {
	circos_links <- NULL
}
rear_pair <- read.table(paste(projectId, "breakmer.rear_pair_diff.txt", sep="_"),sep="\t",header=TRUE,stringsAsFactors=FALSE)
target_summary <- read.table(paste(projectId, "breakmer.target_summary.txt", sep="_"),sep="\t",header=TRUE,stringsAsFactors=FALSE)


filtered_header <- c("subprojectid", "sampleid", "rearr_type", "event_description", "realign_strands", "target_name", "gene_list", "gene-gene_type", "breakpoint coverages", "splitread_counts", "discordantread_counts", "segment_values(realign_chr_cytoband:bp:refseq_annotated_gene:distance to annotated_gene:realign_strand)", "contig_sequence", "filtered")
colnames(filtered) <- filtered_header
writeWorksheet(wb,filtered,"Filtered out",startRow=1, header=T) #headerStyle=hs
setColumnWidths(wb, "Filtered out", column=1:length(filtered_header))


results_header <- c("Subproject ID", "Sample ID", "Rearrangement Type", "Event Description", "Realignment strands", "Target region name", "Gene List", "gene-gene type", "breakpoint coverages", "Splitread Counts", "Discordant Read Counts", "Segment Values (realign_chr_cytoband:bp:refseq_annotated_gene:distance to annotated_gene:realign_strand)", "Contig Sequence")
colnames(breakmer) <- results_header
writeWorksheet(wb,breakmer,"Results",startRow=1) # headerStyle=hs
setColumnWidths(wb, "Results",  columns=1:length(results_header))
setColumnWidths(wb, "Results",  columns=1, currWidth = 10000)


colnames(sample_summary) <- sample_header
writeWorksheet(wb,sample_summary,"Sample summary",startRow=1,header=T) 
setColumnWidths(wb, "Sample summary",  columns=1:length(sample_header))

if(!is.null(circos_annos)){
  writeWorksheet(wb,circos_annos,"circos annos all",startRow=1,header=F)
  setColumnWidths(wb, "circos annos all",  columns=1:ncol(circos_annos))
}
if(!is.null(circos_links)){
  writeWorksheet(wb,circos_links,"circos links all",startRow=1,header=F)
  setColumnWidths(wb, "circos links all",  columns=1:ncol(circos_links))
} 

rear_header <- c("Rearranged genes", "No. samples", "No. of events")
colnames(rear_pair) <- rear_header
writeWorksheet(wb,rear_pair,"Rearr pair counts",startRow=1 ) 
setColumnWidths(wb, "Rearr pair counts", columns=1:length(rear_header))

target_header <- c("Target region", "Nevents", "Nsamples_with_event")
colnames(target_summary) <- target_header
writeWorksheet(wb,target_summary,"Target summary",startRow=1)
setColumnWidths(wb, "Target summary",  columns=1:length(target_header))
saveWorkbook(wb,paste(fileName, "Breakmer_report.xlsx", sep="_"))

warnings()