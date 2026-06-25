#!/usr/bin/env Rscript

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

