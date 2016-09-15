# continuation of test.R

basename = JOB_PREFICES[job]
EXCLUDE_LARGE_MATRICES = c()

writeLines(paste("Loading data from '", basename, "' XHMM run...",
    sep = ""))
dataList = list()
GC = loadNamedVectorNoHeaderMayNotExist(paste(basename, ".locus_GC.txt",
    sep = ""))
dataList[["GC"]] = GC
Complexity = loadNamedVectorNoHeaderMayNotExist(paste(basename,
    ".locus_complexity.txt", sep = ""))
dataList[["Repeat-masked"]] = Complexity

#########
# Read RD
RD = readNamedMatrix(paste(basename, ".RD.txt", sep = ""))
#########


# indices in .xcnv file
# chr4:3445689-3450029 "1-02042", c(39556:39560)
(3445689+3450029)/2
3450029-3445689

# should be: chr5:140482544-140531117, 27.73
# problem is that indices of TARGETS is off
# indices in .xcnv file
RD["1-04044", c(51348:51353)]
# indices in RD dataframe
RD["1-04044", c(58422:58427)]
# chr4:159531708-159531925
mean(RD["1-04044", c(58422:58427)])
# 27.725
# these indices are probably for RD filtered

#  sort( sapply(ls(),function(x){object.size(get(x))}))
dataList[["RD"]] = RD
TARGETS = colnames(RD)
remove(RD)

TARGETS_CHR_BP1_BP2 = targetsToChrBp1Bp2(TARGETS)
TARGET_SIZES = targetsToSizes(TARGETS)
dataList[["TARGETS"]] = TARGETS
dataList[["TARGETS_CHR_BP1_BP2"]] = TARGETS_CHR_BP1_BP2
dataList[["Target size"]] = TARGET_SIZES
dataList[["CHROMOSOMES_START_WITH_CHR"]] = length(grep("^chr",
    TARGETS_CHR_BP1_BP2[["chr"]], perl = TRUE)) > 0

# "chr1:1-3400000"
# using only the chromXHMM generated targets, new interval is:
"chr1:762098-3397257"
# TARGETS = 1..529
# NUM_TARG = 529

# should be next chr:
# TARGETS_CHR_BP1_BP2$chr
# TARGETS_CHR_BP1_BP2$bp1
# TARGETS_CHR_BP1_BP2$bp2

target_indices = which(TARGETS_CHR_BP1_BP2$bp1 >= 140482544 &
    TARGETS_CHR_BP1_BP2$bp2 <= 140531117 &
    TARGETS_CHR_BP1_BP2$chr == "chr5")
# check if this does the same thing:

interval = "chr4:3445689-3450029" # 1-02042
interval = "chr5:140482544-140531117" # 1-04044
targetsToChrBp1Bp2(interval)

##################
# Read centered RD
filtered_centered_RD = readNamedMatrix(paste(basename,
    ".filtered_centered.RD.txt", sep = ""))
##################
# these don't match with anything so probably not necessary
TARGETS_CHR_BP1_BP2 = targetsToChrBp1Bp2(colnames(filtered_centered_RD))

# chr4:3445689-3450029
target_indices = which(TARGETS_CHR_BP1_BP2$bp1 >= 140482544 &
    TARGETS_CHR_BP1_BP2$bp2 <= 140531117 &
    TARGETS_CHR_BP1_BP2$chr == "chr5")
# indices in .xcnv file, should be -4.08/40.59, and for -02 -4.33/27.75
filtered_centered_RD["1-02042", c(39556:39560)]
mean(filtered_centered_RD["1-02042", c(39556:39560)])

# from this matrix
filtered_centered_RD["1-02042", c(39574:39578)]
mean(filtered_centered_RD["1-02042", c(39574:39578)])
# 15.34205

# indices from original RD matrix
filtered_centered_RD["1-02042", c(45036:45040)]
mean(filtered_centered_RD["1-02042", c(45036:45040)])


dataList[["filtered_centered_RD"]] = filtered_centered_RD

filtered_centered_RD.filtered_samples = scanVectorMayNotExist(paste(basename,
    ".filtered_centered.RD.txt.filtered_samples.txt", sep = ""))
dataList[["filtered_centered_RD.filtered_samples"]] = filtered_centered_RD.filtered_samples
filtered_centered_RD.filtered_targets = scanVectorMayNotExist(paste(basename,
    ".filtered_centered.RD.txt.filtered_targets.txt", sep = ""))
dataList[["filtered_centered_RD.filtered_targets"]] = filtered_centered_RD.filtered_targets

remove(filtered_centered_RD)

################
# Read PC matrix
PC = readNamedMatrix(paste(basename, ".RD_PCA.PC.txt",
    sep = ""))
################
dataList[["PC"]] = PC
remove(PC)

PC_SD = readNamedMatrix(paste(basename, ".RD_PCA.PC_SD.txt",
    sep = ""))
dataList[["PC_SD"]] = PC_SD
PC_LOADINGS = readNamedMatrix(paste(basename, ".RD_PCA.PC_LOADINGS.txt",
    sep = ""))
dataList[["PC_LOADINGS"]] = PC_LOADINGS
PCA_NORMALIZE_NUM_REMOVED = paste(basename, ".PCA_normalized.txt.num_removed_PC.txt",
    sep = "")
NUM_PC_REMOVED_VEC = scanVectorMayNotExist(PCA_NORMALIZE_NUM_REMOVED)
if (length(NUM_PC_REMOVED_VEC) == 1) {
    NUM_PC_REMOVED = as.numeric(NUM_PC_REMOVED_VEC[1])
}
else {
    PCA_NORMALIZE_OUT_FILE = paste(basename, ".PCA_normalized.txt.out",
        sep = "")
    if (!file.exists(PCA_NORMALIZE_OUT_FILE)) {
        stop(paste("Cannot find output log file '", PCA_NORMALIZE_OUT_FILE,
            "' either", sep = ""))
    }
    NUM_PC_REMOVED = as.numeric(system(paste("grep 'Removing first' ",
        PCA_NORMALIZE_OUT_FILE, " | awk '{print $3}'", sep = ""),
        intern = TRUE))
    if (length(NUM_PC_REMOVED) != 1) {
        stop(paste("Unable to determine number of PC removed from '",
            PCA_NORMALIZE_OUT_FILE, "'", sep = ""))
    }
}
dataList[["NUM_PC_REMOVED"]] = NUM_PC_REMOVED

##################
# reading PCA norm
PCA_NORMALIZED = readNamedMatrix(paste(basename, ".PCA_normalized.txt",
    sep = ""))
##################
dataList[["PCA_NORMALIZED"]] = PCA_NORMALIZED
remove(PCA_NORMALIZED)

###########################
# reading PCA norm z scores
PCA_NORM_Z_SCORES = readNamedMatrix(paste(basename, ".PCA_normalized.filtered.sample_zscores.RD.txt",
    sep = ""))
###########################
dataList[["PCA_NORM_Z_SCORES"]] = PCA_NORM_Z_SCORES
remove(PCA_NORM_Z_SCORES)

# confirm that colnames are location
colnames(dataList[["PCA_NORM_Z_SCORES"]])[1:5]
TARGETS_CHR_BP1_BP2 = targetsToChrBp1Bp2(colnames(dataList[["PCA_NORM_Z_SCORES"]]))

# chr4:3445689-3450029
target_indices = which(TARGETS_CHR_BP1_BP2$bp1 >= 1 &
    TARGETS_CHR_BP1_BP2$bp2 <= 3400000 &
    TARGETS_CHR_BP1_BP2$chr == "chr1")

TARGETS_CHR_BP1_BP2$chr[[target_indices[1]]]
TARGETS_CHR_BP1_BP2$bp1[[target_indices[1]]]
# TARGETS_CHR_BP1_BP2$bp2[[target_indices[1]]]
# TARGETS_CHR_BP1_BP2$bp1[[target_indices[421]]]
TARGETS_CHR_BP1_BP2$bp2[[target_indices[100]]]

TARGETS_CHR_BP1_BP2$chr[[target_indices[421]]]
TARGETS_CHR_BP1_BP2$bp1[[39556]]
TARGETS_CHR_BP1_BP2$bp2[[39556]]

PCA_NORM_Z_SCORES.filtered_samples = scanVectorMayNotExist(paste(basename,
    ".PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt",
    sep = ""))
dataList[["PCA_NORM_Z_SCORES.filtered_samples"]] = PCA_NORM_Z_SCORES.filtered_samples
PCA_NORM_Z_SCORES.filtered_targets = scanVectorMayNotExist(paste(basename,
    ".PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt",
    sep = ""))
dataList[["PCA_NORM_Z_SCORES.filtered_targets"]] = PCA_NORM_Z_SCORES.filtered_targets

grep("00050", rownames(dataList[["PCA_NORM_Z_SCORES"]]))

################################
# reading unkown variation of RD
RD_SAME_FILTERED = readNamedMatrix(paste(basename, ".same_filtered.RD.txt",
    sep = ""))
dataList[["RD_SAME_FILTERED"]] = RD_SAME_FILTERED
################################


#######################
# posterior prob of DIP
DIP_POST = readNamedMatrix(paste(basename, ".posteriors.DIP.txt",
    sep = ""))
#######################
dataList[["DIP_POST"]] = DIP_POST
remove(DIP_POST)


#######################
# posterior prob of DEL
DEL_POST = readNamedMatrix(paste(basename, ".posteriors.DEL.txt",
    sep = ""))
#######################
dataList[["DEL_POST"]] = DEL_POST
remove(DEL_POST)

#######################
# posterior prob of DUP
DUP_POST = readNamedMatrix(paste(basename, ".posteriors.DUP.txt",
    sep = ""))
#######################
dataList[["DUP_POST"]] = DUP_POST
remove(DUP_POST)

XCNV_CALLS = loadXCNVcalls(paste(basename, ".xcnv", sep = ""))
dataList[["XCNV_CALLS"]] = XCNV_CALLS
