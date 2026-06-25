setwd("Data/")

r <- c(0.20000, 0.25000, 0.33333, 0.45000, 0.50000, 0.66667, 0.75000, 1.00000, 1.50000, 1.90000, 2.00000, 2.10000, 2.50000, 3.00000, 3.50000, 3.90000, 4.00000, 4.10000, 4.50000, 5.00000, 5.50000, 5.90000, 6.00000, 6.10000, 6.50000, 7.00000, 7.50000, 7.90000, 8.00000, 8.10000, 8.50000, 9.00000, 9.50000, 10.00000)

t <- c(16e3, 32e3)
t <- c(16e3, 32e3)/4
FH <- paste0("a7-r", c(10,13,14,15), "--DX.dat")
FH <- paste0("p0-a7-r", 2, "--DX.dat")
dx <- ReadDataList(gtools::mixedsort(lf("a1-")), idx=1)
dx <- ReadDataList(FileList=FH, idx=1)
dx <- data.matrix(fread(lf("a1-r2--")))

dim(dx)
L <- getTSD_hist(abs(dx), t=rep(t[1],length(FH)), n=500)#, domain=exp.domain)
mplot(L$x, L$y, type='s', log='y', lty='solid',)# , spa=FALSE, xaxt='n', xlim=c(0.35,0.6), )#main=paste0("r = ", r[i+7]))
spar(mfrow=c(3,3))
vbfh <- c(0.355, 0.325, 0.285, 0.305, 0.308, 0.333, 0.365)
for (i in 1:8) {
    mplot(L$x[,i], L$y[,i], type='s', log='y', lty='solid', spa=FALSE, xaxt='n', xlim=c(0.35,0.6))
    abline(v=vb.pi2.irr[i])
    #abline(v=vbfh[i])
    #axis(side=1, at=seq(0,0.5,by=0.05))
    axis(side=1, at=seq(0,1,by=0.05))
}

vb.pi3 <- c(0.4355, 0.421, 0.40, 0.373, 0.364, 0.334, 0.363, 0.4, 0.77, 0.962, 0.4, 0.446, 0.609, 0.77, 0.897, 0.982, 0.666, 0.688, 0.77, 0.857, 0.934, 0.988, 0.77, 0.784, 0.835, 0.896, 0.951, 0.991, 0.823, 0.833, 0.873, 0.919, 0.962, 0.857)
vb.pi3.err <- c(0.001, 0.001, 0.001, 0.001, 0.001, 0.002, 0.002, 0.0025, 0.003, 0.0035, 0.0035, 0.0035, 0.0032, 0.0031, 0.003, 0.003, 0.003, 0.003, 0.0034, 0.0035, 0.0041, 0.003, 0.0032, 0.003, 0.004, 0.003, 0.003, 0.0025, 0.0035, 0.0035, 0.0038, 0.004, 0.0035, 0.004)


vb.pi2 <- c(0.5,0.5001,0.5,0.5045,0.5,0.515,0.514, 0.4999, 0.6, 0.655,2/3,0.678, 0.714,0.75,0.778,0.796,0.80,0.804,0.819,0.833,0.847,0.855,0.857,.86,.867,0.876,0.883,0.8875,.89,0.891,0.895,0.90,0.905,0.908)

vb.pi3.irr.err <- c(0.001, 0.001, 0.001, 0.001, 0.001, 0.002, 0.002, 0.0025, 0.003, 0.0035, 0.0035, 0.0035, 0.0032, 0.0031, 0.003, 0.003, 0.003, 0.003, 0.0034, 0.0035, 0.0041, 0.003, 0.0032, 0.003, 0.004, 0.003, 0.003, 0.0025, 0.0035, 0.0035, 0.0038, 0.004, 0.0035, 0.004)

vb.pi3.irr     <- c(0.44, 0.428, 0.408, 0.38, 0.37, 0.34, 0.368, 0.402, 0.755, 0.94, 0.977, 0.425, 0.583, 0.74, 0.864, 0.946, 0.965, 0.983, 0.734, 0.82, 0.895, 0.948, 0.96, 0.972, 0.796, 0.856, 0.91, 0.95, 0.958, 0.97, 0.83, 0.878, 0.918, 0.958)
vb.pi2.irr <- c(0.512,0.51,0.509,0.515,0.51,0.519,0.52, 0.508, 0.598, 0.649,0.659,0.67,0.704,0.738,0.763,0.78,0.785,0.789,0.801,0.816,0.827,0.835,0.838,.84,.848,0.85,0.8595,0.862,.868,0.869,0.871,0.874,0.879,0.885)
#fwrite(x=cbind(r, vb.pi2, vb.pi3, vb.pi2.irr, vb.pi3.irr), file="VB--ALL.tkz", col.names=FALSE, sep=" ")
i <- i+1
print(paste0("i = ", i, " length(vb) = ", length(vb.pi2)))
{
    file <- lf(paste0("a6-r",i,"--"))
    dx <- data.matrix(fread(file))
    L <- getTSD_hist(abs(dx), t=t, n=750)#, domain=exp.domain)
    Lx <- L$x; Ly <- L$y
    mplot(Lx, Ly, type='s', log='y', xaxt='n',main=paste0("r = ", r[i+1]), spa=FALSE, xlim=c(0,1))
    axis(side=1, at=seq(0,1,by=0.05))
    abline(v=  vb.pi3.irr[i+1], col='red')
    #abline(v=c(vb.pi3.slow[i-7]+vb.pi3.err.slow[i-7]))
    #abline(v=c(vb.pi3.slow[i-7]-vb.pi3.err.slow[i-7]))
}
Lx <- Lx[,1]
    Ly <- Ly[,1]
    n <- which(Lx[] > 0.3)
    # FH: W = 5, LX > 0.2. 
    # IH: W = 7, LX > 0.3
    ids <- argmax(x=Lx[n], y=(log(Ly[n])), w=20)
    X <- Lx[n[1] + ids -1]
    Y <- (Ly)[n[1] + ids -1]
    abline(v=X)
}

mplot(r, vb.pi3, type='o', ylim=c(0.25, 1))
lines(r, vb.pi3-vb.pi3.err)
k <- seq(3)[1]
k <- 1


# For infinite horizon only! n_0 = 1
g <- function(r, a) {
    n_0 <- 1
    1/(1 + n_0 * r * (1 + cos(a)))
}

g.bar.old <- function(r, a) {
    g(r,a) * (r + (r - 2)*cos(a))
}


g.bar <- function(r, a) { 
    dy <- (1/2) * (1/tan(a/2))
    N <- floor((dy - dy*r*cos(a))/(dy*r*tan(a/2)*sin(a)))
    den1 <- 1 + r*(N + 1)*(1 + cos(a)) 
    den2 <- 1 + r*N*(1 + cos(a))
    num1 <- (1/r) - N + cos(a) * (N - 1)
    num2 <- 1 - (1/r) + N*(1 - cos(a))
    (num1/den1) + (num2/den2)
}

g.bar.1 <- function(r, a) {  # For a > pi/2
    #rr <- which(rr > 
    dy <- (1/2) * (1/tan(a/2))
    N <- floor(1/(r*(1-cos(a))))
    den2 <- 1 + r*N*(1 + cos(a))
    num2 <- 1 - (1/r) + N*(1 - cos(a))
    (num2/den2)
}

g.bar.2 <- function(r, a) {  # For a > pi/2
    dy <- (1/2) * (1/tan(a/2))
    N <- floor(1/(r*(1-cos(a))))
    den1 <- 1 + r*(N + 1)*(1 + cos(a)) 
    num1 <- (1/r) - N + cos(a) * (N + 1)
    (num1/den1) 
}

g.tilde <- function(r, n, a) { 
    g.bar(r,n+1,a) * ((1/r) - n + cos(a)*(n - 1)) + g.bar(r,n,a) * (1 - (1/r) + n*(1 - cos(a))) 
}
g.1 <-  function(r,n,a) {
    dy <- (1/2) * (1/tan(a/2))
    s <- 2*r*dy*sin(a) + 2*(n+1)*cos(a) + 1
    1/(2*sqrt(3)*r*dy + 4)
}

# alpha = pi/3 average speed finite horizon
phi <- (sqrt(5) + 1)/2
a <- (pi/2) * (phi-(3/5))
dy <- (1/2) * cot(a/2)

mplot(r[1:8], cbind(vb.pi2, vb.pi2.irr)[1:8,], type='o', col=c('red', 'blue'))
mplot(r, cbind(vb.pi2, vb.pi2.irr), type='o', col=c('red', 'blue'))
R <- seq(0.2, 10, length=10000)
plot(R, g.bar(r=R, a), type='l')
lines(R, g.bar(r=R, a), type='l')

lines(R, g.bar.1(r=R, a) + g.bar.2(r=R,a), type='l')
lines(R, g.bar(r=R, a), type='l')

lines(R, 1 - (1/ (1 + R*(1 + cos(a)))), type='l', col='red')

lines(R, g.bar(r=R, pi/2), type='l')
lines(R, g.tilde(R, 1, pi/3) + g.1(R,0,pi/3), type='l')

lines(R, g.bar.old(r=R, a), type='l', col='blue')


k <- seq(1, 10)
d <- 1/((1 - cos(a))*k)
abline(v=d)
abline(h=1/(1+(1+cos(a))/(1-cos(a))))

points(r, vb.pi2.irr)

eps <- 0.01
a <- runif((pi/2) - eps, pi/2 + eps, n=10000)
x <- sapply(1:length(R), function(x) g.bar(R[x], a))
dim(x)

L <- getTSD_ecdf(x, n=1000)
lf <- LinFit(mat=1-L$y, time=L$x, range=100:200 ,logs=FALSE)
plot(R, lf$ax)
mplot(L$x[,1], 1-L$y[,1], type='s')

# Density plots for pi/3
t <- c(16e3, 32e3)
id <- c(10,12,16)-1 + 8
r[id]
FH <- paste0("a2-r", id, "--DX.dat")
dx <- ReadDataList(FileList=FH, idx=2)
{

    L <- getTSD_hist(abs(dx), t=t[2], n=600)#, domain=exp.domain)
    Lx <- L$x; Ly <- L$y
    mplot(Lx, Ly, type='s', log='y', xaxt='n', spa=FALSE, xlim=c(0.2,1))
    axis(side=1, at=seq(0,1,by=0.05))
    abline(v=vb.pi2[id+1])
    abline(v=vb.pi3[id+1])
}
#fwrite(x=cbind(L$x, L$y), file="PdxPI2vbal.tkz", sep=" ", col.names=FALSE)
ys <- sapply(1:ncol(L$x), function(x) L$y[min(which(L$x[,x] > vb.pi3[id[x]+1]))-1,x])
ys <- sapply(1:ncol(L$x), function(x) L$y[min(which(L$x[,x] > vb.pi2[id[x]+1]))-1,x])
fwrite(file="temp", format(cbind(vb.pi2[id+1], ys) , digits = 4, trim = TRUE))
abline(h=ys)
    abline(v=vb.pi3[id+1])




