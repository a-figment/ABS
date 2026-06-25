GroupFiles <- function(fname, key) {
    Files <- lf(fname)
    for (i in seq(length(key))) 
        Files <- Files[grep(x=Files, pattern = key[i])]
    return(Files)
}

getSimulation <- function(File) {
    Time <- as.numeric(read.table(file=File, nrows=1, skip=1, header=FALSE))
    nICs <- read.table(file=File, nrows=1, skip=3, header=FALSE)$V1
    nCollisions <- read.table(file=File, nrows=1, skip=5, header=FALSE)$V1
    TotalTime <- read.table(file=File, nrows=1, skip=8, header=FALSE)$V1
    L <- as.numeric(read.table(file=File, nrows=1, skip=11, header=FALSE))
    Symbols <- read.table(file=File, nrows=1, skip=14, header=FALSE)$V1
    Sim <- list(Time=Time, nIC=nICs, nCollisions=nCollisions, TotalTime=TotalTime, L=L, Symbols=Symbols)
    return(Sim)
}


setwd("Data/")
DataOut <- "../DataSummary/PGFData/"


# Escape to absorbing boundary
Files <- GroupFiles(".dat", "s0")
Sim <- getSimulation("s0-b1-p1--Simulation.txt") # b0 == b1


MFPT <- ReadData(Files[grep(x=Files, pattern="MFPT")])[,1,]


FPT <- ReadData(Files[grep(x=Files, pattern="FP_Time")], Tpose=FALSE)


# w.r.t to 
ECDF <- apply(FPT, 2, 




# Forced escape to absorbing boundary
System1 <- lf("s1")


# Symbolic escapes
System2 <- lf("s2")



