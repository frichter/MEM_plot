# Felix Richter
# 08/25/2016
# XHMM exon read depth PCA z-score for each ME

# module load xhmm
# module load R
# R

library(xhmmScripts)
library(dplyr)
library(ggplot2)
library(tidyr)

basename = "/hpc/users/richtf01/chdiTrios/Felix/xhmm_test/hg19.all.xhmm"
JOB_TARGETS_TO_GENES = "/hpc/users/richtf01/chdiTrios/Felix/xhmm_test/hg19.annotated_targets.refseq.loci"

##############
# Load ME data
##############

# 3 sources of Mendelian errors

# 1 ME per line
me_loci = read.table("/sc/orga/projects/chdiTrios/Katie/MEM/WES/hg19.exons_me_intersect.bed", sep = "\t") %>%
  select(V4, V5, V6) %>%
  rename(SAMPLE = V4, CHROM = V5, ME_POS = V6) %>%
  mutate(SAMPLE = gsub("/sc/orga/projects/chdiTrios/Katie/MEM/WES/(.*).mendel.bed", "\\1", SAMPLE)) %>%
  mutate(CHROM = gsub("^", "chr", CHROM))

# 1 exon (and multiple MEs) per line
me_exon_loci = read.table("/sc/orga/projects/chdiTrios/Katie/MEM/WES/hg19.exons_me_list.txt", sep = "\t", header = T) %>%
  mutate(CHROM_EXON = gsub("^", "chr", CHROM_EXON))

# 1 exon (and multiple MEs) per line, overlapped with windows of interest (WOI)
me_exon_loci_woi = read.table("/hpc/users/richtf01/chdiTrios/Felix/xhmm_test/wes_mem_final_woi_by_exon.txt",
  sep = "\t", header = T, comment.char = "") %>%
  rename(CHROM_WOI = X.CHROM_WOI) %>%
  mutate(CHROM_EXON = gsub("^", "chr", CHROM_EXON),
         CHROM_WOI = gsub("^", "chr", CHROM_WOI),
         ID_WOI = as.character(ID_WOI),
         ID_ME_on_EXON = as.character(ID_ME_on_EXON))

# take only individuals with window of interest and a Mendelian error
me_exon_loci_woi = me_exon_loci_woi %>%
  rowwise() %>%
  mutate(woi_exon_ID_match = ID_WOI %in% strsplit((ID_ME_on_EXON), "[|]")[[1]])

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

############################################
# figure out which exons + IDs have z-scores
############################################

# label exons that have read depth z scores
me_exon_loci_woi = me_exon_loci_woi %>%
  rowwise() %>%
  # count how many exons each ME overlaps with (expecting 1)
  mutate(target_avail_count = sum(TARGETS_CHR_BP1_BP2$chr == CHROM_EXON &
      TARGETS_CHR_BP1_BP2$bp1 <= START_EXON &
      TARGETS_CHR_BP1_BP2$bp2 >= END_EXON)) %>%
      # note that XHMM final z-score targets are shorter than the
      # original exon reference targets
  ungroup()

# check that data looks good
me_exon_loci_woi %>% filter(target_avail_count != 0 & woi_exon_ID_match) %>% as.data.frame %>% head

# label the individuals that were processed by XHMM
me_exon_loci_woi = me_exon_loci_woi %>%
  rowwise() %>%
  mutate(avail_in_xhmm = ID_WOI %in% rownames(dataList[["PCA_NORM_Z_SCORES"]])) %>%
  ungroup()

###########################################
# Combine into 1 dataframe and plot results
###########################################

# annotate the ME loci with the corresponding z-score
me_loci_z = me_exon_loci_woi %>%
  filter(target_avail_count == 1, avail_in_xhmm, woi_exon_ID_match) %>%
  rowwise() %>%
  mutate(target = which(TARGETS_CHR_BP1_BP2$chr == CHROM_EXON &
                        TARGETS_CHR_BP1_BP2$bp1 <= START_EXON &
                        TARGETS_CHR_BP1_BP2$bp2 >= END_EXON)) %>%
  mutate(pca_norm_z = dataList[["PCA_NORM_Z_SCORES"]][ ID_WOI, target],
         exon_chr = TARGETS_CHR_BP1_BP2$chr[target],
         exon_start = TARGETS_CHR_BP1_BP2$bp1[target],
         exon_end = TARGETS_CHR_BP1_BP2$bp2[target],
         exon_target_name = colnames(dataList[["PCA_NORM_Z_SCORES"]])[target]) %>%
  ungroup()

head(me_loci_z) %>% as.data.frame()

# unique(me_loci_z$ID_WOI)
# unique(me_loci_z$exon_target_name)

# extract z-scores for non-ME individuals
non_me_id_indices = !(rownames(dataList[["PCA_NORM_Z_SCORES"]]) %in% me_loci_z$ID_WOI)
non_me_id_indices = non_me_id_indices & !grepl("-01$|-02$", rownames(dataList[["PCA_NORM_Z_SCORES"]]))
non_me_loci_z = dataList[["PCA_NORM_Z_SCORES"]][non_me_id_indices, unique(me_loci_z$exon_target_name)]
non_me_loci_z_long = non_me_loci_z %>%
  as.data.frame() %>%
  gather("exon_target_name", "pca_norm_z") #, `chr1:79403474-79403996`:`chr21:43523637-43524180`

# combine non-ME and ME individuals for graphing
loci_all_z = me_loci_z %>% select(exon_target_name, pca_norm_z) %>%
  bind_rows(non_me_loci_z_long, .id = "ID_w_ME") %>%
  mutate(ID_w_ME = ifelse(ID_w_ME == 1, "Yes", "No"))

p = ggplot(loci_all_z, aes(x = ID_w_ME, y = pca_norm_z)) +
  geom_boxplot() +
  xlab("Mendelian error present") +
  ylab("Read depth z-score (XHMM)") +
  ggtitle("IDs with MEs vs other\n(in exons with MEs)") +
  theme_bw()
ggsave("/hpc/users/richtf01/chdiTrios/Felix/xhmm_test/exons_w_ME.png", p, width = 3, height = 4)

p = loci_all_z %>%
  filter(ID_w_ME == "No") %>%
  ggplot(., aes(x = exon_target_name, y = pca_norm_z)) +
  geom_boxplot() +
  geom_boxplot(data = loci_all_z %>% filter(ID_w_ME == "Yes"), color = "red") +
  xlab("Exon target name") +
  ylab("Read depth z-score (XHMM)") +
  ggtitle("IDs with MEs (red) vs other") +
  theme_bw() +
  theme(legend.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
ggsave("/hpc/users/richtf01/chdiTrios/Felix/xhmm_test/exons_w_ME_by_target.png", p, width = 10, height = 4)

# t test for significance
id_w_me = loci_all_z %>% filter(ID_w_ME == "Yes") %>% select(pca_norm_z) %>% unlist %>% as.numeric
id_wo_me = loci_all_z %>% filter(ID_w_ME == "No") %>% select(pca_norm_z) %>% unlist %>% as.numeric
t.test(id_w_me, id_wo_me)$p.value
# 1.595864e-24

# long format of IDs + z-scores for merging
id_vec = rownames(dataList[["PCA_NORM_Z_SCORES"]])
z_score_id_df = dataList[["PCA_NORM_Z_SCORES"]][, unique(loci_all_z$exon_target_name)] %>%
  as.data.frame %>%
  mutate(Blinded.ID = id_vec) %>%
  gather("exon_target_name", "pca_norm_z", -Blinded.ID) 

# select individuals + targets with z-scores less than the hits with MEM scores
# potential MEM false negatives
non_mem_outliers = loci_all_z %>%
  group_by(exon_target_name) %>%
  filter(pca_norm_z == min(pca_norm_z) | pca_norm_z < -10) %>%
  filter(pca_norm_z < -5) %>%
  ungroup %>%
  filter(ID_w_ME == "No") %>%
  left_join(z_score_id_df)
write.table(non_mem_outliers, "/hpc/users/richtf01/chdiTrios/Felix/xhmm_test/non_mem_outliers.txt",
            sep = "\t", col.names = T, row.names = F, quote = F)

# select individuals + targets with z-scores less than the hits with MEM scores
# potential MEM false positives
mem_potential_fp_df = loci_all_z %>%
  filter(pca_norm_z > -1 & ID_w_ME == "Yes") %>%
  left_join(z_score_id_df)

write.table(mem_potential_fp_df, "/hpc/users/richtf01/chdiTrios/Felix/xhmm_test/mem_potential_fp_df.txt",
            sep = "\t", col.names = T, row.names = F, quote = F)

#####################
# possibly deprecated
#####################

me_loci = me_loci %>%
  rowwise() %>%
  mutate(target_avail_count = sum(TARGETS_CHR_BP1_BP2$chr == CHROM & TARGETS_CHR_BP1_BP2$bp1 <= ME_POS & TARGETS_CHR_BP1_BP2$bp2 >= ME_POS)) %>%
  ungroup()

# half of MEs did not overlap with any measured exon
table(me_loci$target_avail_count)
# 0: 6098, 1: 6174

# mark samples that were processed with XHMM
me_loci = me_loci %>%
  rowwise() %>%
  mutate(avail_in_xhmm = SAMPLE %in% rownames(dataList[["PCA_NORM_Z_SCORES"]])) %>%
  ungroup()

# identify continuous MEs
me_loci_z_cont = me_loci_z %>%
  group_by(SAMPLE, CHROM) %>%
  mutate(same_id_same_chrom_count = n()) %>%
  mutate(same_id_same_chrom = ifelse(n() > 2, TRUE, FALSE)) %>%
  ungroup() %>%
  mutate(ME_count = ifelse(same_id_same_chrom_count > 10, "> 10", "6-10")) %>%
  mutate(ME_count = ifelse(same_id_same_chrom_count < 6, "3-6", ME_count))

# plot a histogram of the z-scores
p = me_loci_z_cont %>%
  filter(same_id_same_chrom) %>%
  ggplot(., aes(pca_norm_z, fill = ME_count)) +
  geom_histogram(bins = 100) +
  # geom_vline(xintercept = mean(me_loci_z$pca_norm_z), col = "red") +
  ylab("Count of MEs for each exon z-score") +
  xlab("Exon z-score (binned)") +
  theme_bw()

ggsave("/hpc/users/richtf01/chdiTrios/Felix/xhmm_test/z_hist_cont.png", p, width = 6, height = 5)

write.table(me_loci_z, "/hpc/users/richtf01/chdiTrios/Felix/xhmm_test/z_loci_me.txt", sep = "\t", quote = F, row.names = F)
