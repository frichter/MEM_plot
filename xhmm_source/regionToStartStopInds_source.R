parseRegion <- function(chrom, start=NULL, stop=NULL) {
	if (is.null(start) || is.null(stop)) {
		parse1 = strsplit(chrom, ":")
		chrom = sapply(parse1, "[[", 1)

		parse2 = sapply(parse1, "[[", 2)
		bp1_bp2 = strsplit(parse2, "-")

		start = as.numeric(sapply(bp1_bp2, "[[", 1))
		stop = as.numeric(sapply(bp1_bp2, "[[", 2))
	}

	return(list(chrom=chrom, start=start, stop=stop))
}

regionToFinalInds <- function(xhmm_data, chrom, start=NULL, stop=NULL) {
	if (is.null(start) || is.null(stop)) {
		chrom_start_stop = parseRegion(chrom, start, stop)
		chrom = as.character(chrom_start_stop["chrom"])
		start = as.numeric(chrom_start_stop["start"])
		stop  = as.numeric(chrom_start_stop["stop"])
	}

	chrBp1Bp2 = targetsToChrBp1Bp2(colnames(xhmm_data[["PCA_NORM_Z_SCORES"]]))
	return(which(chrBp1Bp2[["chr"]] == chrom & chrBp1Bp2[["bp1"]] <= stop & start <= chrBp1Bp2[["bp2"]]))
}

regionToStartStopInds <- function(xhmm_data, chrom, start=NULL, stop=NULL, NUM_ADD_TARGS=2) {
	targFinalInds = regionToFinalInds(xhmm_data, chrom, start, stop)

	startTarg = min(targFinalInds)
	stopTarg = max(targFinalInds)

	chr = targetsToChrBp1Bp2(colnames(xhmm_data[["PCA_NORM_Z_SCORES"]]))[["chr"]]

	if (NUM_ADD_TARGS > 0) {
		potentialStartTargs = max(1, startTarg - NUM_ADD_TARGS):startTarg
		startTarg = min(potentialStartTargs[chr[potentialStartTargs] == chr[startTarg]])

		potentialStopTargs = stopTarg:min(length(chr), stopTarg + NUM_ADD_TARGS)
		stopTarg = max(potentialStopTargs[chr[potentialStopTargs] == chr[stopTarg]])
	}

	return(c(startTarg=startTarg, stopTarg=stopTarg))
}
