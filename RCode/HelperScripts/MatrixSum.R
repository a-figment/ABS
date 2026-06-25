#!/usr/bin/env Rscript
files = commandArgs(trailingOnly=TRUE) # Last argument is output location

cat("\nSumming matrices contained in the following files:\n")
print(data.frame(files=files[1:(length(files)-1)]), max=10)

cat("\nWriting output to:\n")
outfile <- files[length(files)]
print(data.table::as.data.table(outfile))

# Sum of multiple matrices
matrix_to_sum <- data.table::fread(files[1], colClasses=c("numeric"))
for (i in seq(2, length(files)-1)) {
    current_matrix <- data.table::fread(files[i], colClasses=c("numeric"))
    matrix_to_sum <- matrix_to_sum + current_matrix
}


data.table::fwrite(x=matrix_to_sum, file=outfile, append=FALSE, sep=" ", col.names=FALSE, row.names=FALSE)
cat("\nDone!\n")
