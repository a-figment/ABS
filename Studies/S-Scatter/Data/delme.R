
source("~/Dropbox/Rlibs.R")

V <- data.matrix(fread("p0-a7-r0--V.dat"))

L <- getTSD_hist(dx=V, n=1000)

plot(L$x, L$y)
mplot(L$x, L$y, type='l')
