# Felix Richter
# 08/25/2016
# MEM XHMM cross-over plotting

module load xhmm
module load R
R

library(xhmmScripts)
library(dplyr)

# sshfs richtf01@mothra.hpc.mssm.edu:/hpc/users/richtf01/chdiTrios/Felix/xhmm_test ~/Volumes/xhmm/

plot_prefix = "/hpc/users/richtf01/chdiTrios/Felix/xhmm_test/plot_fam/fam"
basename = "/hpc/users/richtf01/chdiTrios/Felix/xhmm_test/hg19.all.xhmm"
JOB_TARGETS_TO_GENES = "/hpc/users/richtf01/chdiTrios/Felix/xhmm_test/hg19.annotated_targets.refseq.loci"

################
# load MEM calls
################

# general file format = ID, CHROM_ME START_ME STOP_ME
# .chr15_deletions_me.txt
# .putative_deletion_region_with_me_to_plot.txt
fam_file_list = list.files("/sc/orga/projects/chdiTrios/Katie/MEM/WES",
  pattern = "*.chr15_deletions_me.txt", full.names = TRUE)
sample_id_list = lapply(fam_file_list, function(x) gsub("/sc/orga/projects/chdiTrios/Katie/MEM/WES/(.*).chr15_deletions_me.txt", "\\1", x))
fam_table_list = lapply(fam_file_list, function(x) read.table(x, header = TRUE))
names(fam_table_list) = sample_id_list
sample_info_df = fam_table_list %>% bind_rows(.id = "SAMPLE") %>%
  group_by(SAMPLE, CHROM_ME) %>%
  summarise(START = START_ME[[1]], END = END_ME[[n()]]) %>%
  ungroup() %>%
  rename(CHROM = CHROM_ME) %>%
  mutate(CHROM = paste("chr", CHROM, sep = ""))


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

########################
# what data can you use?
########################

# confirm that there are exons (aka targets) in the MEM loci
IsTargetAvail = function(CHROM, START, END, TARGETS_CHR_BP1_BP2) {
    # chr = paste("chr", CHROM, sep = "")
    chr = CHROM
    start_bp = START #- 50000
    end_bp = END #+ 50000
    total_targets = sum(TARGETS_CHR_BP1_BP2$bp1 >= start_bp &
        TARGETS_CHR_BP1_BP2$bp2 <= end_bp &
        TARGETS_CHR_BP1_BP2$chr == chr)
    return(total_targets > 0)
}

# add columns to determine which samples to process
sample_info_df = sample_info_df %>%
  rowwise() %>%
  # confirm that all samples are available in XHMM results:
  mutate(avail_in_xhmm = SAMPLE %in% rownames(dataList[["PCA_NORM_Z_SCORES"]])) %>%
  mutate(target_avail = IsTargetAvail(CHROM, START, END, TARGETS_CHR_BP1_BP2)) %>%
  ungroup()

######
# Plot
######

# wrapper function for plotting the family
PlotFamMEXHMM = function(SAMPLE, CHROM, START, END, dataList, TARGETS_CHR_BP1_BP2, plot_prefix) {
    print(paste("Starting with sample", SAMPLE))
    # set desired x axis
    # chr = paste("chr", CHROM, sep = "")
    chr = CHROM
    start_bp = START #- 50000
    end_bp = END #+ 50000
    type = "DEL"

    usePlotName = paste(plot_prefix, SAMPLE, chr, start_bp, end_bp, collapse = "_")

    # set z-score cut-offs for y-axis
    # none for initial analysis, possibly zoomed in for repeats
    ylim_upper_z = 5
    # ylim_lower_z = -10
    dataList[["PCA_NORM_Z_SCORES"]][ dataList[["PCA_NORM_Z_SCORES"]] > ylim_upper_z ] = ylim_upper_z
    # dataList[["PCA_NORM_Z_SCORES"]][ dataList[["PCA_NORM_Z_SCORES"]] < ylim_lower_z ] = ylim_lower_z


    # Obtain target indices: base pair location and position in data of
    # exons to plot
    target_indices = which(TARGETS_CHR_BP1_BP2$bp1 >= start_bp &
        TARGETS_CHR_BP1_BP2$bp2 <= end_bp &
        TARGETS_CHR_BP1_BP2$chr == chr)

    # convert interval into format for plotting
    addTargs = tail(target_indices, n = 1) - target_indices[1] + 1
    interval = paste(TARGETS_CHR_BP1_BP2$chr[target_indices[[1]]], ":", TARGETS_CHR_BP1_BP2$bp1[target_indices[[1]]], "-", TARGETS_CHR_BP1_BP2$bp2[tail(target_indices, n = 1)], sep = "")
    chrom_start_stop = targetsToChrBp1Bp2(interval)
    MARK_INTERVALS = list()
    MARK_INTERVALS[[paste("XHMM ", type, " (",
      addTargs, " targets)", sep = "")]] = list(as.numeric(chrom_start_stop$bp1),
      as.numeric(chrom_start_stop$bp2), col = "black",
      cex = 0.8)
    start_stop_inds = c(target_indices[1], tail(target_indices, n = 1))
    names(start_stop_inds) = c("startTarg", "stopTarg")
    print(MARK_INTERVALS)

    MARK_SAMPLES = as.character(unlist(paste(SAMPLE, c("", "-01", "-02"), sep = "")))

    XCNV_CALLS = loadXCNVcalls(paste(basename, ".xcnv", sep = ""))
    XCNV_CALLS_test = XCNV_CALLS[1:3, ]
    XCNV_CALLS_test$SAMPLE = MARK_SAMPLES
    XCNV_CALLS_test$INTERVAL = interval
    XCNV_CALLS_test$KB = round((as.numeric(chrom_start_stop$bp2) - as.numeric(chrom_start_stop$bp1))/1000, digits = 2)
    XCNV_CALLS_test$CHR = chrom_start_stop$chr
    XCNV_CALLS_test$MID_BP = round((as.numeric(chrom_start_stop$bp2) + as.numeric(chrom_start_stop$bp1))/2, digits = 0)
    XCNV_CALLS_test$TARGETS = paste(start_stop_inds[1], start_stop_inds[2], sep = "..")
    XCNV_CALLS_test$NUM_TARG = start_stop_inds[2] - start_stop_inds[1]

    dataList[["XCNV_CALLS"]] = XCNV_CALLS_test
    print(XCNV_CALLS_test)
    print(start_stop_inds)
    plot_XHMM_targets(usePlotName, dataList,
        dataList[["TARGETS_TO_GENES"]],
        binarySampleFeatures = NULL, SQ_THRESH = 60, start_stop_inds[1],
        start_stop_inds[2], PLOT_ONLY_PNG = TRUE,
        MARK_SAMPLES = MARK_SAMPLES, MARK_INTERVALS = MARK_INTERVALS,
        APPEND_REGION_NAME = FALSE)
    return(0)
}

t = sample_info_df %>% filter(SAMPLE == "1-04824", CHROM == 14)

plot_XHMM_targets

# test run a single individual
t = sample_info_df[1, ]
PlotFamMEXHMM(t$SAMPLE, t$CHROM, t$START, t$END, dataList, TARGETS_CHR_BP1_BP2)

# generate plots :)
sample_info_df[9:80, ] %>%
  filter(avail_in_xhmm) %>%
  filter(target_avail) %>%
  rowwise() %>%
  summarise(PlotFamMEXHMM(SAMPLE, CHROM, START, END, dataList, TARGETS_CHR_BP1_BP2))


######################
# adjustments for 1p36
######################

# 1p36 samples of interest:
sample_list = c("1-02711", "1-04652", "1-05084", "1-06498", "1-07442")

# chr15 samples of interest: load in first section
sample_list = sample_info_df %>% filter(avail_in_xhmm) %>% select(SAMPLE) %>% unlist
names(sample_list) = sample_list
# sample_re = paste(sample_list, collapse = "|")

# plotting windows:
# for 1p36:
chr_i = "chr1"
start_i = 1337395
end_i = 1431087

# chr15 group
chr_i = "chr15"
start_i = 22800000 # min(sample_info_df$START) - 1
end_i = max(sample_info_df$END) + 1
sum(TARGETS_CHR_BP1_BP2$bp1 >= start_i &
    TARGETS_CHR_BP1_BP2$bp2 <= end_i &
    TARGETS_CHR_BP1_BP2$chr == chr_i)

# plot by FAMILY
plot_prefix = "/hpc/users/richtf01/chdiTrios/Felix/xhmm_test/chr15/fam_"
# test on 1 individual
PlotFamMEXHMM(sample_lo_list_other[[1]], chr_i, start_i, end_i,
  dataList, TARGETS_CHR_BP1_BP2, plot_prefix)
# repeat for all other individuals
lapply(sample_lo_list_other[[2]], function(x)
  PlotFamMEXHMM(x, chr_i, start_i, end_i, dataList, TARGETS_CHR_BP1_BP2, plot_prefix))

# plotting all PROBANDS, PARENTS of interest, or both
PlotSampleList = function(sample_list, CHROM, START, END, dataList,
                          TARGETS_CHR_BP1_BP2, plot_prefix, basename,
                          id_index_to_keep) {
    # set desired x axis (ie location on genome)
    chr = CHROM # chr_i #
    start_bp = START # start_i # - 50000
    end_bp = END # end_i # + 50000
    type = "DEL"

    usePlotName = plot_prefix

    # set z-score cut-offs for y-axis
    ylim_upper_z = 5
    # ylim_lower_z = -10
    dataList[["PCA_NORM_Z_SCORES"]][ dataList[["PCA_NORM_Z_SCORES"]] > ylim_upper_z ] = ylim_upper_z
    # dataList[["PCA_NORM_Z_SCORES"]][ dataList[["PCA_NORM_Z_SCORES"]] < ylim_lower_z ] = ylim_lower_z

    dataList[["PCA_NORM_Z_SCORES"]] = dataList[["PCA_NORM_Z_SCORES"]][id_index_to_keep, ]

    # Obtain target indices: base pair location and position in data of
    # exons to plot
    target_indices = which(TARGETS_CHR_BP1_BP2$bp1 >= start_bp &
        TARGETS_CHR_BP1_BP2$bp2 <= end_bp &
        TARGETS_CHR_BP1_BP2$chr == chr)

    # convert interval into format for plotting
    numberOfTargets = tail(target_indices, n = 1) - target_indices[1] + 1
    interval = paste(TARGETS_CHR_BP1_BP2$chr[target_indices[[1]]], ":",
                     TARGETS_CHR_BP1_BP2$bp1[target_indices[[1]]], "-",
                     TARGETS_CHR_BP1_BP2$bp2[tail(target_indices, n = 1)],
                     sep = "")
    chrom_start_stop = targetsToChrBp1Bp2(interval)
    MARK_INTERVALS = list()
    MARK_INTERVALS[[paste("XHMM ", type, " (",
      numberOfTargets, " targets)", sep = "")]] = list(as.numeric(chrom_start_stop$bp1),
      as.numeric(chrom_start_stop$bp2), col = "black",
      cex = 0.8)
    start_stop_inds = c(target_indices[1], tail(target_indices, n = 1))
    names(start_stop_inds) = c("startTarg", "stopTarg")
    print(MARK_INTERVALS)

    MARK_SAMPLES = sample_list

    XCNV_CALLS = loadXCNVcalls(paste(basename, ".xcnv", sep = ""))
    XCNV_CALLS_test = XCNV_CALLS[1:length(sample_list), ]
    XCNV_CALLS_test$SAMPLE = MARK_SAMPLES
    XCNV_CALLS_test$INTERVAL = interval
    XCNV_CALLS_test$KB = round((as.numeric(chrom_start_stop$bp2) - as.numeric(chrom_start_stop$bp1))/1000, digits = 2)
    XCNV_CALLS_test$CHR = chrom_start_stop$chr
    XCNV_CALLS_test$MID_BP = round((as.numeric(chrom_start_stop$bp2) + as.numeric(chrom_start_stop$bp1))/2, digits = 0)
    XCNV_CALLS_test$TARGETS = paste(start_stop_inds[1], start_stop_inds[2], sep = "..")
    XCNV_CALLS_test$NUM_TARG = start_stop_inds[2] - start_stop_inds[1]

    dataList[["XCNV_CALLS"]] = XCNV_CALLS_test
    print(XCNV_CALLS_test)
    print(start_stop_inds)
    plot_XHMM_targets(usePlotName, dataList,
        dataList[["TARGETS_TO_GENES"]],
        binarySampleFeatures = NULL, SQ_THRESH = 60, start_stop_inds[1],
        start_stop_inds[2], PLOT_ONLY_PNG = TRUE,
        MARK_SAMPLES = MARK_SAMPLES, MARK_INTERVALS = MARK_INTERVALS,
        APPEND_REGION_NAME = FALSE)
    return(0)
}


# Choose which individuals to plot (Parents, Probands, or both)
# only probands:
id_index_to_keep = grep("-01$|-02$", row.names(dataList[["PCA_NORM_Z_SCORES"]]), invert = T)
# only parents
id_index_to_keep = grep("-01$|-02$", row.names(dataList[["PCA_NORM_Z_SCORES"]]))
# all/both
id_index_to_keep = grep("", row.names(dataList[["PCA_NORM_Z_SCORES"]]))


plot_prefix = "/hpc/users/richtf01/chdiTrios/Felix/xhmm_test/chr15/proband_only_zoomed"
PlotSampleList(sample_list, chr_i, start_i, end_i, dataList,
               TARGETS_CHR_BP1_BP2, plot_prefix, basename,
               id_index_to_keep)


##########################################################
# Identify other individuals with deletions in same region
##########################################################


# identify probands with data in specified target regions
# 1p35: targets 7:12, chr15 region: ???
# Extract read depth z-scores for these individuals
sample_df = dataList[["PCA_NORM_Z_SCORES"]][, target_indices]
sample_df = sample_df %>% as.data.frame %>% mutate(Blinded.ID = row.names(sample_df))

# What is the definition of "low"?
# Heuristic: use highest z-score in MEM probands + 0.5 per exon as cutoff
sample_list = c("1-02711", "1-04652", "1-05084", "1-06498", "1-07442") # for 1p36
sample_list = sample_info_df %>% filter(avail_in_xhmm) %>% select(SAMPLE) %>% unlist # for chr15

blinded_id_interest = sample_df$Blinded.ID %in% sample_list
perColCutOffs = apply(sample_df %>% filter(blinded_id_interest) %>% select(-Blinded.ID), 2, max)

# count the number of targets less than the cut-off for probands and parents
id_index_to_keep = grep("-01$|-02$", row.names(dataList[["PCA_NORM_Z_SCORES"]]), invert = T)
cutCount_vec = apply(sample_df %>% select(-Blinded.ID), 1, function(pca_values) sum(pca_values <= perColCutOffs + 0.5))
sample_df$cutCount_empirical = cutCount_vec

sample_lo_list = sample_df %>%
  filter(!grepl("-01$|-02$", Blinded.ID)) %>%
  filter(cutCount_empirical >= 60) %>% # for chr15
  # filter(cutCount_empirical >= 6 & `chr1:1387418-1387724` < -2) %>% # for 1p36
  select(Blinded.ID) %>% unlist
names(sample_lo_list) = NULL
sample_lo_list
# 6: 10, 5: 15, 4: 24

sample_lo_list_other = sample_lo_list[!(sample_lo_list %in% sample_list)]

# sample_lo_list = lapply(sample_lo_list, function(x) as.character(unlist(paste(x, c("-01", "-02"), sep = ""))))

plot_prefix = "/hpc/users/richtf01/chdiTrios/Felix/xhmm_test/chr15/parent_lo"
PlotSampleList(sample_lo_list, chr_i, start_i, end_i, dataList,
               TARGETS_CHR_BP1_BP2, plot_prefix, basename)


###############################
# plot the "de novos" by family
###############################

# plot the "de novos" by family
sample_lo_list_parent = gsub("-0[1-2]$", "", sample_lo_list)
parent_transmitted = lapply(sample_lo_list_parent, function(parent) parent %in%
                            sample_lo_list_proband) %>% unlist
names(sample_lo_list_proband) = sample_lo_list_proband
not_in_parents = lapply(sample_lo_list_proband, function(proband) !(proband %in%
                            sample_lo_list_parent)) %>% unlist
SAMPLE = sample_lo_list_proband[not_in_parents][[2]]
sample_list = as.character(unlist(paste(SAMPLE, c("", "-01", "-02"), sep = "")))
plot_prefix = paste("/hpc/users/richtf01/chdiTrios/Felix/xhmm_test/plot_fam/1p36/fam_", SAMPLE, sep = "")
PlotSampleList(sample_list, chr_i, start_i, end_i, dataList,
               TARGETS_CHR_BP1_BP2, plot_prefix, basename)




# What is the definition of "low"?

# NAIVE: use arbitrary cut-offs of 0, -1, and -2
# cutoff_list = c(0, -1, -2)
# names(cutoff_list) = c("Neg", "Neg1", "Neg2")
# cutCount_df = sapply(cutoff_list, function(cutoff)
#                      apply(sample_df, 1, function(pca_values) sum(pca_values < cutoff)),
#                      USE.NAMES = T)
# cutCount_df = cutCount_df %>% as.data.frame %>% mutate(Blinded.ID = row.names(cutCount_df))
# sample_df = sample_df %>% left_join(cutCount_df)
# sample_df$ID_means = rowMeans(sample_df[, 1:6])
# sample_df %>% filter(Blinded.ID %in% sample_list)
# sample_df %>% filter(Neg == 6) %>% dim
# # 62 samples with read depth all negative
# sample_df %>% filter(Neg1 >= 5) %>% dim
# # 31
# sample_df %>% filter(Neg1 == 6) %>% dim
# # 23 samples with all reads < -1
# sample_df %>% filter(Neg2 >= 4) %>% dim
# # 28
