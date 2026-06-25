
####' 
##   PATH: S-JumpStatistics/Data/
##
####' 


img.path = "~/Dropbox/PHD/Papers/JumpStatistics/Images/"
setwd("../Studies/S-JumpStatistics/Data/")


t <- seq(0, 1e6, length=1e6)
x <- data.matrix(fread("p0-a1-d2--1-x.dat"))

xtp <- ReadDataList(lf("1-xtp"))
ytp <- ReadDataList(lf("1-ytp"))

spar(cx=2)
matplot(t, x, type='l', lwd=2)

lines(xtp[,2], ytp[,2], col='red', type='o', pch=19)

h <- 100
s1 <- 1:(nrow(x)-h+1)
s2 <- h:nrow(x)
pac <- x[s1] * x[s2]


