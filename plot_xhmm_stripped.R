# Felix Richter
# 08/25/2016

# module load xhmm
# module load R
# R

library(xhmmScripts)
library(dplyr)
library(ggplot2)
library(tidyr)

basename = "/hpc/users/richtf01/chdiTrios/Felix/xhmm_test/hg19.all.xhmm"
JOB_TARGETS_TO_GENES = "/hpc/users/richtf01/chdiTrios/Felix/xhmm_test/hg19.annotated_targets.refseq.loci"

basename = "/sc/orga/projects/chdiTrios/ME_DEL/xhmm/xhmm_hg19/hg19.all.xhmm"
JOB_TARGETS_TO_GENES = "/sc/orga/projects/chdiTrios/ME_DEL/xhmm/xhmm_hg19/hg19.annotated_targets.refseq.loci"


###################
# load XHMM results
###################

# import PCA normalization scores
dataList = list()
PCA_NORM_Z_SCORES = readNamedMatrix(paste(basename, ".PCA_normalized.filtered.sample_zscores.RD.txt",
sep = ""))
dataList[["PCA_NORM_Z_SCORES"]] = PCA_NORM_Z_SCORES
remove(PCA_NORM_Z_SCORES)

# import targets (aka exons) and their locations
TARGETS_CHR_BP1_BP2 = targetsToChrBp1Bp2(colnames(dataList[["PCA_NORM_Z_SCORES"]]))
dataList[["CHROMOSOMES_START_WITH_CHR"]] = length(grep("^chr",
TARGETS_CHR_BP1_BP2[["chr"]], perl = TRUE)) > 0
dataList[["TARGETS_TO_GENES"]] = loadTargetsToGenes(JOB_TARGETS_TO_GENES,
CHROMOSOMES_START_WITH_CHR = TRUE)

XCNV_CALLS = loadXCNVcalls(paste(basename, ".xcnv", sep = ""))
dataList[["XCNV_CALLS"]] = XCNV_CALLS
NUM_ADD_TARGS = 2
SQ_THRESH = 60
PLOT_ONLY_PNG = TRUE
PLOT_DIR = "/sc/orga/projects/chdiTrios/ME_DEL/xhmm/xhmm_hg19/chr22_plots"

for (cnvInd in 1:nrow(XCNV_CALLS)) {
  cnv = XCNV_CALLS[cnvInd, ]
  sample = as.character(cnv["SAMPLE"])
  interval = as.character(cnv["INTERVAL"])
  SQ = as.numeric(cnv["Q_SOME"])
  type = as.character(cnv["CNV"])
  numTargets = as.numeric(cnv["NUM_TARG"])
  MARK_SAMPLES = sample
  plotName = paste(PLOT_DIR, "/sample_", sample,
    ".", interval, ".SQ_", SQ, sep = "")
  chrom_start_stop = targetsToChrBp1Bp2(interval)
  MARK_INTERVALS = list()
  MARK_INTERVALS[[paste("XHMM ", type, " (",
    numTargets, " targets)", sep = "")]] = list(as.numeric(chrom_start_stop$bp1),
    as.numeric(chrom_start_stop$bp2), col = "black",
    cex = 0.8)
  for (addTargs in c(0, NUM_ADD_TARGS)) {
    start_stop_inds = regionToStartStopInds(dataList,
      interval, NUM_ADD_TARGS = addTargs)
    usePlotName = plotName
    if (addTargs == 0) {
      usePlotName = paste(usePlotName, ".exact",
        sep = "")
    }
    plot_XHMM_targets(usePlotName, dataList,
      dataList[["TARGETS_TO_GENES"]],
      SAMPLE_FEATURES, SQ_THRESH, start_stop_inds[1],
      start_stop_inds[2], PLOT_ONLY_PNG = PLOT_ONLY_PNG,
      MARK_SAMPLES = MARK_SAMPLES, MARK_INTERVALS = MARK_INTERVALS,
      APPEND_REGION_NAME = FALSE)
  }
}











#
