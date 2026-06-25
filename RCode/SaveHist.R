
args <- commandArgs(trailingOnly = TRUE)


domain <- function(dx, n) seq(min(dx,na.rm=TRUE), max(dx,na.rm=TRUE), length=n)
cat("Creating histogram from file ", args[1], " with ", args[2], " bins.\n")

fin <- args[1]
nBins <- as.numeric(args[2])
fout <- args[3]

# Expects column ordered matrix 
dat <- data.matrix(data.table::fread(fin))

cat("ABSOLUTE VALUES: ON")

# Create histogram from each column
hMids <- hDens <- c()
for (i in seq(ncol(dat))) {
    df <- abs(dat[,i])
    H <- hist(x=df, breaks=domain(df, nBins), plot=FALSE)
    hMids <- rbind(hMids, H$mids)
    hDens <- rbind(hDens, H$density)
}

cat("Writing histogram data to ", fout, ".\n")
data.table::fwrite(x=cbind(t(hMids), t(hDens)), file=fout, sep=" ", col.names=FALSE)

# Testing
#L <- data.matrix(data.table::fread("histout.dat"))
#K <- 1:(ncol(L)/2)
#mplot(L[,K], L[,-K], type='s')



