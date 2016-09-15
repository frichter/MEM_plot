# Felix Richter
# felix.richter@icahn.mssm.edu

# I started using dplyr in January thanks to Jason Homsy. It
# was a little tricky at first but has made my life A LOT easier
# %>% is a pipe, like "|" in bash
library(dplyr)

setwd("/Users/felixrichter/Volumes/minerva/")

Get_Min_DP_table = function(filename) {
  # The second line of this file has 3 columns- in R the number of columns has to 
  # be the same for every line in order to read it as a table. Here we skip
  # the second line. Alternatively pre-process with bash: 
  # sed '2d' 1-07196.mendel.vcf.txt > 1-07196.mendel.tbl.txt
  all_content = readLines(filename)
  skip_second = all_content[-2]
  mendel_table = read.table(textConnection(skip_second), sep = "\t", 
                            header = T, row.names = NULL, comment.char = "")
  mendel_min_dp = mendel_table %>%
    rowwise() %>%
    mutate(min_DP = min(child_read_depth, mom_read_depth, dad_read_depth)) %>% 
    select(X.CHROM, START, min_DP)
  return(mendel_min_dp)
}

# test out the function with a single file
filename = "1-07196.mendel.vcf.txt"
mendel_min_dp = Get_Min_DP_table(filename) 

# apply the function to all files
file_list = list.files(pattern = "*.mendel.vcf.txt")

# apply (and its derivatives eg lapply) is weird and only in R, but extremely useful
# the output is a list of tables, where 1 table 
mendel_min_dp_list = lapply(file_list, function(filename) Get_Min_DP_table(filename))

# here I combine the list of tables into one giant table
mendel_min_dp_all = mendel_min_dp_list %>% bind_rows

# this command provides the count matrix:
count(mendel_min_dp_all, min_DP)

# this command plots the histogram (number of counts per depth)
library(ggplot2)
ggplot(mendel_min_dp_all, aes(x = min_DP)) + geom_histogram(bins = 100)

