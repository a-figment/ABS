#!/usr/bin/env Rscript
files = commandArgs(trailingOnly=TRUE) # Last argument is output location

options(digits=18)

cat("\nAveraging matrices contained in the following files:\n")
print(data.frame(files=files[1:(length(files)-1)]), max=10)

# Online averaging
dims <- dim(data.table::fread(files[1], colClasses=c("numeric"))) 
matrix_to_average <- matrix(0., nrow=dims[1], ncol=dims[2])
for (i in seq(1, length(files)-1)) {
    current_matrix <- data.matrix(data.table::fread(files[i], colClasses=c("numeric")))
    matrix_to_average <- matrix_to_average + (current_matrix - matrix_to_average)/i
}
matrix_to_average <- data.matrix(matrix_to_average)

# Raw averaging
#matrix_to_average <- data.table::fread(files[1])
#for (i in seq(2, length(files)-1)) {
#    matrix_to_average <- matrix_to_average + data.table::fread(files[i])
#}
#matrix_to_average <- (matrix_to_average / (length(files) - 1))

cat("\nWriting output to:\n")
outfile <- files[length(files)]
print(data.table::as.data.table(outfile))
data.table::fwrite(x=matrix_to_average, file=outfile, sep=" ", col.names=FALSE)
