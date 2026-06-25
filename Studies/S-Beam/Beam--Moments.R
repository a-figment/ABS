

# Here we compare moments of displacement for discrete and continuous cases

###
#   CONTINUOUS MOMENTS
###

t <- seq(0, 1e6, length=1001)

setwd("Data/T0/")
# Top and bottom angles

at <-
ab <- 

f <- list.files(path = "Continuous/", pattern = "Moments")
msd <- tpose(apply(ReadData(f), 1, colMeans))
mplot(t, msd, log='xy')
LinFit(mat=msd, time=t, plot = TRUE, range = 2:501)

mu <- getPowerScaling(t = t, M = msd, z=0.9)
mplot(t, mu, log='xy')


###
#   DISPLACEMENT DISTRIBUTION
###

dx <- data.matrix(fread(lff("b0-t1", "DX")))

L <- getTSD_ecdf(abs(dx), n=1000, domain=exp.domain)
mplot(L$x, 1-L$y, type='s', log='xy')



###
#   DISCRETE MOMENTS
###


