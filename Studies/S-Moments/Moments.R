
setwd("Data/")

# Length of trajectories and number of steps (see *Moments.cpp)
T <- 10^(seq(5,7))
nSteps <- c(101, 1001)

# Ratios
r <- seq(1/4, 10, length=391)[8:31]
r <- seq(1/5, 10, length=981)
d <- (1/2) * r

# Use data key to determine time
{
    a.file <- lf("Moments.dat")[1]
    t.idx <- grep(a.file, pattern="-t") + 1
    n.idx <- grep(a.file, pattern="-n") + 1

    t.idx <- 2
    n.idx <- 2

    times <- seq(0, T[t.idx], length=nSteps[n.idx])
    ttg <- c(times[2], times[11], times[101]) 
}


X <- tpose(fread(lf("R.dat")))
T <- tpose(fread(lf("L.dat")))
mplot(times, T, type='l', lwd=1.5, xlab="Times to observe the particle", ylab="Actual times observed")
max(as.vector(apply(T, 2, function(x) max(x - times))))
mplot(T[,87], X[,87])

Moments <- data.matrix(ReadData(FileList=lff("a7", "Moments"))[,,1])
Moments <- tpose(fread(lf("Moments")))[]
dim(Moments)

n <- 1000
tt <- 1e8
times <- c(seq(tt/n, tt, length=n))

# Plot moments as a function of time
mplot(times[-1], Moments[-1,], log='xy')

# Plot the power scaling as a function of time
mu <- getPowerScaling(t = times[-1], M = Moments[-1,], z = 0.9)
mplot(times[-1], mu[,], ylim=c(1,2))

# Scalings obtained from large mu
gamma <- mu[min(which(times > 8e4)),]

# Linear regression on log-log plot (blind)
LF <- LinFit(mat = Moments[-1,], time = times[-1], range = 200:999)
gamma <- LF$ax

# Non-linear regression
gamma <- c()
for (i in seq(ncol(Moments))) {
    df <- data.frame(x=times, y=Moments[,i])
    gamma <- c(gamma, coef(nlsLM(y ~ b*x^gamma, data=df, start=list(b=1, gamma=1)))[2])
}

# Parameter dependence of moment scalings
mplot(r[1:81], gamma, type='o', xaxt='n')#ylim=c(1.58, 1.85))
mplot(r[1:81], gamma, add=TRUE, col='blue')#ylim=c(1.58, 1.85))
mplot(r[1:81], gamma, add=TRUE, col='black')#ylim=c(1.58, 1.85))
abline(h=5/3)
abline(h=1.65)
axis(side=1, at=seq(0.2, 1, by=0.05))
mplot(r, gamma, type='p', col='blue', add=TRUE)


# Displacement distributions
DX <- ReadDataList(FileList=lff("a7", "DX"), idx = 3)
DX <- data.matrix(fread(lff("r30", "DX")))

L <- getTSD_hist(dx = DX, n=1000, t = ttg, eta=1.82/2)

mplot(L$x, L$y, type='s', log='y')



