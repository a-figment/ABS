
setwd("../Studies/S-JumpStatistics/Data/")

file.means <- gtools::mixedsort(lf("--ND-Means"))
file.vars <- gtools::mixedsort(lf("--ND-Vars"))
file.num <- gtools::mixedsort(lf("--ND.dat"))


# Detection statistics are written in terms of 

# DIM = (TIME [small to large], WIDTH [small to large], PARTICLE NUMBER) 
means <- ReadData(FileList = file.means, Tpose = FALSE)


num   <- ReadData(FileList = file.num, Tpose = FALSE)

rss <- rowSums(num, dims=2)
rms <- rowMeans(num, dims = 2)

rvs <- array(NA, dim = c(5, 5, dim(num)[3]))
for (i in seq(dim(num)[3])) {
    rvs[,,i] <- (num[,,i] - rms)^2
}
rvs <- rowMeans(rvs, dims=2)


