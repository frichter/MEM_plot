module load xhmm
module load R
R

library(xhmmScripts)
library(dplyr)

# cat /sc/orga/projects/chdiTrios/Katie/MEM/WES/1-10188.putative_deletion_region_with_me_to_plot.txt

PLOT_PATH = "/hpc/users/richtf01/chdiTrios/Felix/xhmm_test/plots_test"
basename = "/hpc/users/richtf01/chdiTrios/Felix/xhmm_test/hg19.all.xhmm"
JOB_TARGETS_TO_GENES = "/hpc/users/richtf01/chdiTrios/Felix/xhmm_test/hg19.annotated_targets.refseq.loci"
# plot for all individuals on one figure:
usePlotName = "/hpc/users/richtf01/chdiTrios/Felix/xhmm_test/plot_fam/fam_1-10188"

# load ID, CHROM_ME START_ME STOP_ME
fam_file_list = list.files("/sc/orga/projects/chdiTrios/Katie/MEM/WES",
  pattern = "*.putative_deletion_region_with_me_to_plot.txt", full.names = TRUE)
sample_id_list = lapply(fam_file_list, function(x) gsub("/sc/orga/projects/chdiTrios/Katie/MEM/WES/(.*).putative_deletion_region_with_me_to_plot.txt", "\\1", x))
fam_table_list = lapply(fam_file_list, function(x) read.table(x, header = TRUE))
names(fam_table_list) = sample_id_list
sample_info_df = fam_table_list %>% bind_rows(.id = "SAMPLE") %>%
  group_by(SAMPLE) %>%
  summarise(CHROM = CHROM_ME[[1]], START = START_ME[[1]], END = END_ME[[n()]])


# set desired x axis
chr = "chr12"
start_bp = 21011310 - 50000
end_bp = 21377559 + 50000

# set z-score cut-offs for y-axis
# none for initial analysis, zoomed in for repeats
# ylim_upper_z = 10
# ylim_lower_z = -10

# Set defaults
SAMPLE_FEATURES = NULL
SQ_THRESH = 60
NUM_ADD_TARGS = 2
PLOT_READ_DEPTHS = TRUE
PLOT_PC_CORRS = TRUE
PLOT_ALL_CNVS = TRUE
USE_XCNV_TO_PLOT = NULL
INCLUDE_PEDIGREE_SAMPLES = NULL
PLOT_ONLY_PNG = TRUE
LIMIT_MEMORY = FALSE
type = "DEL"

dataList = list()
# import PCA normalization scores
PCA_NORM_Z_SCORES = readNamedMatrix(paste(basename, ".PCA_normalized.filtered.sample_zscores.RD.txt",
    sep = ""))
dataList[["PCA_NORM_Z_SCORES"]] = PCA_NORM_Z_SCORES
remove(PCA_NORM_Z_SCORES)

# Obtain target indices: base pair location and position in data of
# exons to plot
TARGETS_CHR_BP1_BP2 = targetsToChrBp1Bp2(colnames(dataList[["PCA_NORM_Z_SCORES"]]))
target_indices = which(TARGETS_CHR_BP1_BP2$bp1 >= start_bp &
    TARGETS_CHR_BP1_BP2$bp2 <= end_bp &
    TARGETS_CHR_BP1_BP2$chr == chr)

# confirm your targets (exons) cover the desired regions:
TARGETS_CHR_BP1_BP2$chr[target_indices]
TARGETS_CHR_BP1_BP2$bp1[target_indices]
TARGETS_CHR_BP1_BP2$bp2[target_indices]
# USE THIS TO CREATE THE ACTUAL INTERVAL

# import relationship between target (exon) names and coordinates
dataList[["CHROMOSOMES_START_WITH_CHR"]] = length(grep("^chr",
  TARGETS_CHR_BP1_BP2[["chr"]], perl = TRUE)) > 0
dataList[["TARGETS_TO_GENES"]] = loadTargetsToGenes(JOB_TARGETS_TO_GENES,
  CHROMOSOMES_START_WITH_CHR = dataList[["CHROMOSOMES_START_WITH_CHR"]])

# convert interval into format for plotting

addTargs = tail(target_indices, n = 1) - target_indices[1] + 1
interval = paste(TARGETS_CHR_BP1_BP2$chr[target_indices[[1]]], ":", TARGETS_CHR_BP1_BP2$bp1[target_indices[[1]]], "-", TARGETS_CHR_BP1_BP2$bp2[tail(target_indices, n = 1)], sep = "")
# interval = "chr1:861282-3397257"
chrom_start_stop = targetsToChrBp1Bp2(interval)
MARK_INTERVALS = list()
MARK_INTERVALS[[paste("XHMM ", type, " (",
  addTargs, " targets)", sep = "")]] = list(as.numeric(chrom_start_stop$bp1),
  as.numeric(chrom_start_stop$bp2), col = "black",
  cex = 0.8)

interval


start_stop_inds = c(target_indices[1], tail(target_indices, n = 1))
names(start_stop_inds) = c("startTarg", "stopTarg")
##### confirm that this can be used as input instead of regionToStartStopInds

# create y-axis cut-offs for graphing (do not use results for statistical analyses)
dataList_subset = dataList
# dataList_subset[["PCA_NORM_Z_SCORES"]][ dataList_subset[["PCA_NORM_Z_SCORES"]] > ylim_upper_z ] = ylim_upper_z
# dataList_subset[["PCA_NORM_Z_SCORES"]][ dataList_subset[["PCA_NORM_Z_SCORES"]] < ylim_lower_z ] = ylim_lower_z

# all individuals on one plot
# 1p36
# sample_list = c("1-02711", "1-04652", "1-05084", "1-06498", "1-07442")
sample = "1-10188"
MARK_SAMPLES = as.character(unlist(paste(sample, c("", "-01", "-02"), sep = "")))
# confirm samples are in PCA
lapply(MARK_SAMPLES, function(x) grep(x, rownames(dataList[["PCA_NORM_Z_SCORES"]])))


XCNV_CALLS = loadXCNVcalls(paste(basename, ".xcnv", sep = ""))
XCNV_CALLS_test = XCNV_CALLS[1:3, ]
XCNV_CALLS_test$SAMPLE = MARK_SAMPLES
XCNV_CALLS_test$INTERVAL = interval
XCNV_CALLS_test$KB = round((as.numeric(chrom_start_stop$bp2) - as.numeric(chrom_start_stop$bp1))/1000, digits = 2)
XCNV_CALLS_test$CHR = chrom_start_stop$chr
XCNV_CALLS_test$MID_BP = round((as.numeric(chrom_start_stop$bp2) + as.numeric(chrom_start_stop$bp1))/2, digits = 0)
XCNV_CALLS_test$TARGETS = paste(start_stop_inds[1], start_stop_inds[2], sep = "..")
XCNV_CALLS_test$NUM_TARG = start_stop_inds[2] - start_stop_inds[1]

dataList_subset[["XCNV_CALLS"]] = XCNV_CALLS_test
print(XCNV_CALLS_test)
plot_XHMM_targets(usePlotName, dataList_subset,
    dataList_subset[["TARGETS_TO_GENES"]],
    SAMPLE_FEATURES, SQ_THRESH, start_stop_inds[1],
    start_stop_inds[2], PLOT_ONLY_PNG = PLOT_ONLY_PNG,
    MARK_SAMPLES = MARK_SAMPLES, MARK_INTERVALS = MARK_INTERVALS,
    APPEND_REGION_NAME = FALSE)


# Trios on separate plots
sample_list = c("1-02711", "1-04652", "1-05084", "1-06498", "1-07442")
for(sample in sample_list) {
    MARK_SAMPLES = as.character(unlist(paste(sample, c("", "-01", "-02"), sep = "")))
    print(MARK_SAMPLES)
    XCNV_CALLS = loadXCNVcalls(paste(basename, sample, "xcnv", sep = "."))
    dataList[["XCNV_CALLS"]] = XCNV_CALLS
    dataList_temp[["XCNV_CALLS"]] = XCNV_CALLS
    print(XCNV_CALLS)
    # usePlotName = "/hpc/users/richtf01/chdiTrios/Felix/xhmm_test/test_1_avail_1p36_yZoom_xZoomSuper"
    usePlotName = paste("/hpc/users/richtf01/chdiTrios/Felix/xhmm_test/plot_fam/fam_", sample, sep = "")
    print(usePlotName)
    plot_XHMM_targets(usePlotName, dataList_temp,
        dataList_temp[["TARGETS_TO_GENES"]],
        SAMPLE_FEATURES, SQ_THRESH, start_stop_inds[1],
        start_stop_inds[2], PLOT_ONLY_PNG = PLOT_ONLY_PNG,
        MARK_SAMPLES = MARK_SAMPLES, MARK_INTERVALS = MARK_INTERVALS,
        APPEND_REGION_NAME = FALSE)
}
