
####' 
##   DATA PATH:  ~/Dropbox/ABS/Studes/S-JumpStatistics/Data/
##   IMAGE PATH: ~/Dropbox/PHD/Papers/JumpStatistics/Images/
##
##   -- Overview --
##   Analysis of jump events from a bottom up 'flow' definition and 
##   a top down 'coarse grained' definition. 
##
##   Date created: Wed 16 Feb 2022 11:17:52
##   Author: Jordan Orchard
##   
####' 

img.path = "~/Dropbox/PHD/Papers/JumpStatistics/Images/"
setwd("../Studies/S-JumpStatistics/Data/PaperPlots/")

alpha <- c(pi/seq(9, 2, by=-1), seq(2, 8)*pi/seq(3,9) )

alpha <- c(pi/seq(9, 2, by=-1), seq(2, 4)*pi/seq(3,5) )
phi <- (sqrt(5) + 1)/2
alpha.irr <- alpha * (phi - 3/5)

alpha.lab <- c(tex("\\pi/9"), 
	       tex("\\pi/8"), 
	       tex("\\pi/7"), 
	       tex("\\pi/6"), 
	       tex("\\pi/5"), 
	       tex("\\pi/4"), 
	       tex("\\pi/3"), 
	       tex("\\pi/2"), 
	       tex("2\\pi/3"), 
	       tex("3\\pi/4"), 
	       tex("4\\pi/5"), 
	       tex("5\\pi/6"), 
	       tex("6\\pi/7"), 
	       tex("7\\pi/8"), 
	       tex("8\\pi/9"))

alpha.lab.irr <- c(tex("(\\varphi - 3/5) \\pi/9"), 
	       tex("(\\varphi - 3/5) \\pi/8"), 
	       tex("(\\varphi - 3/5) \\pi/7"), 
	       tex("(\\varphi - 3/5) \\pi/6"), 
	       tex("(\\varphi - 3/5) \\pi/5"), 
	       tex("(\\varphi - 3/5) \\pi/4"), 
	       tex("(\\varphi - 3/5) \\pi/3"), 
	       tex("(\\varphi - 3/5) \\pi/2"), 
	       tex("(\\varphi - 3/5) 2\\pi/3"), 
	       tex("(\\varphi - 3/5) 3\\pi/4"), 
	       tex("(\\varphi - 3/5) 4\\pi/5"), 
	       tex("(\\varphi - 3/5) 5\\pi/6"), 
	       tex("(\\varphi - 3/5) 6\\pi/7"), 
	       tex("(\\varphi - 3/5) 7\\pi/8"), 
	       tex("(\\varphi - 3/5) 8\\pi/9"))


# ----------------- LAMINAR FLOW APPROACH -----------------

# Step times and step distances
{
    #nfront      <- data.matrix(fread(lf("nFront")))
    #jump        <- data.matrix(fread(lf("jump")))
    lambda.t    <- ReadDataList(lff("p0", "jump.dat"), idx=1)
    lambda.x    <- ReadDataList(lff("p0", "jump.dat"), idx=2)
    cc          <- ReadDataList(lff("p0", "jump.dat"), idx=3)
    nstep       <- ReadDataList(lff("p0", "jump.dat"), idx=4)
    th          <- ReadDataList(lff("p0", "jump.dat"), idx=5)
}

# Number of fronts as a function of N collisions
th <- c(seq(-pi/400, -1e-5, length=M/2), seq(1e-5, pi/400, length=M/2))
N <- seq(1, 2000, by=1)

#pdf("pi2--nfronts.pdf", width=12, height=9)
mplot(N, nfront, type='o', cex=1, col='blue', log='xy', xaxt='n', yaxt='n',
      xlab=tex("n"), ylab=tex("N(n)"))
axis(side=1, at=c(1,10,100,2000))
axis(side=2, at=c(1,10,100,1000,10000))
LinFit(mat = nfront, time=N, range = 100:990, plot=TRUE)
dev.off()


##
#  FIG. 1 --- IC dependence 2
## 
#pdf(paste0(img.path, "front-jump-time.pdf"), height=12, width=14)
id <- which(nstep[,1] == 1)
mplot(th[id,], lambda.t[id,], type='p', cex=0.4, xaxt='n',
	ylim=c(200,2000), #xlim=c(-pi/400, pi/400),
	xlab=tex("\\theta_0"), ylab=tex("\\tau"))
axis(side=1, at=c(-pi/600, 0, pi/600), labels=c(tex("-599\\pi/600"), "0", tex("599\\pi/600")))
mplot(th[-id,], lambda.t[-id,], type='p', cex=0.4, col='blue', add=TRUE)
dev.off()


# pi/3 no good coz two fronts
xi.bound <- 5000
x <- (lambda.x)
id <- which(nstep[,1] == 1)
id <- which(x < xi.bound[1])
#pdf(paste0(img.path, "icvbal.pdf"), height=12, width=14)
mplot(th[id,], (lambda.x[id,]/lambda.t[id,]), type='p', cex=0.1, xaxt='n',
      ylim=c(0.5, 0.515),
      xlim=c(-pi/600, pi/600),
      #ylim=c(0, 1),
      #xlim=c(-pi,pi),
      xlab=tex("\\theta_0"), ylab=tex("\\xi / \\tau"))
axis(side=1, at=c(-pi/600, 0, pi/600), labels=c(tex("-599\\pi/600"), "0", tex("599\\pi/600")))
mplot(th[-id,], (lambda.x[-id,])/lambda.t[-id,], type='p', cex=0.1, col='blue', add=TRUE)
id <- which(x < xi.bound[1] & nstep > 1)
mplot(th[id,], (lambda.x[id,])/lambda.t[id,], type='p', cex=0.1, col='green', add=TRUE)


xi.bound <- seq(10, 10000, by=10)
xi.bound <- exp(seq(0, log(10000), length=1000))
means <- sdvs <- c()
for (i in seq(length(xi.bound))) {
    id <- which(x > xi.bound[i])
    means <- c(means, mean(lambda.x[id,] / lambda.t[id,]))
    sdvs <- c(sdvs, sqrt(var(lambda.x[id,] / lambda.t[id,])))
}
mplot(xi.bound, means, type='o', log='x', ylim=c(0.4, 0.6))
ErrorBars.y(xi.bound, means, sdvs)


mean(x/lambda.t[id


y0 <- 0.5143
#abline(h=y0)
x <- seq(0, pi/500, length=1000)
lines(x, y0*(-1/(1350000*x+1) + 1))
x <- seq(-pi/500, -0.00001, length=1000)
lines(x, y0*(1/(500000*x+1) + 1))
dev.off()


## Checking ballistic speed as a function of collision number
jump <- data.matrix(fread(lf("jump.dat")))[1:500,]
hist(abs(jump[,2])/jump[,1], n=100)
mplot((jump[,2])/jump[,1], type='l')

##
#  FIG. --- Characteristic jump speed (CJS) as a function of aperture angle
#  We define the CJS as being the two most probable jump events, given by
#  their velocity |\xi| / \tau when 0 << \xi. i.e. a large jump limit 
## 
lambda.t     <- ReadDataList(lff("a3", "jump"), idx=1)
lambda.x     <- ReadDataList(lff("a3", "jump"), idx=2)

lambda.t     <- ReadDataList(lff("d24-i4-c1", "jump.dat"), idx=1)
lambda.x     <- ReadDataList(lff("d24-i4-c1", "jump.dat"), idx=2)

# Distribution of jump speeds
L <- getTSD_hist(dx = (lambda.x/lambda.t)[1:1000000,], n=300)
L <- getTSD_ecdf(dx = abs(lambda.x/lambda.t)[,3], n=2000)

L <- getTSD_hist(dx = (lambda.x/lambda.t)[,], n=1000)#, domain=exp.domain)
Lx <- L$x[,]
Ly <- L$y[,]

#fwrite(x=cbind(Lx, sweep(Ly,2,apply(Ly,2,max),'/')), file="JumpThreshDensity.tkz", sep=" ", col.names=FALSE)

#pdf("small-jv.pdf", width=14, height=10)
#mplot(Lx, sweep(Ly,2,apply(Ly,2,max),'/'), type='s', #xaxt='n', ylim=c(0, 4), xlim=c(0,0.3),
mplot(Lx, Ly, type='s', #xaxt='n', ylim=c(0, 4), xlim=c(0,0.3),
      xlab=TeX("$\\xi / \\tau$"), ylab="Density", lwd=1, log='')
spar(mfrow=c(4,5), cex=0.5)
Cspeed <- c()
vs <- c(0.2,0.2,0.2,0.2,0.25,0.31,0.38,0.4,0.4,0.65, 0.79)
for (i in seq(ncol(Lx))) {
    mplot(Lx[,i], Ly[,i], type='s', spa=FALSE,
    xlab=TeX("$\\xi / \\tau$"), ylab="Density", lwd=1.5, log='')

    # First estimates
    sly <- rev(sort(Ly[,i]))
    #id <- seq(nrow(Lx))
    #if (i > 6) {
    #    id <- which(Lx[,i] > 0.21)
    #} else {
    #    id <- which(Lx[,i] < 0.21)
    #}
    #id <- which(Lx[,i] > vs[i])
    #lxx <- Lx[id,i]
    #lyy <- Ly[id,i]
    #Cspeed <- c(Cspeed, lxx[which.max(lyy)])
    #abline(v=Cspeed[i], col='red')
    
    #abline(v = Cspeed[i], col='blue')
}

r <- seq(1/4, 1, length=76)
cs <- sapply(1:ncol(Ly), function(x) Lx[which.max(Ly[,x]), x])
mplot(r, cs, type='o')
abline(h=1/2)

# RATIONAL AND IRATIONALS ARE BASICALLY THE SAME HERE...
# d8 RATS
Cspeed[c(7,11)] <- c(0.31, 0.745)
#fwrite(file="CS--d8-rat.tkz", x=cbind(alpha,Cspeed), col.names=FALSE, sep=" ")
abline(v = Cspeed, col='blue')

# d8 IRAS
Cspeed[c(7,11)] <- c(0.31, 0.745)
#fwrite(file="CS--d8-ira.tkz", x=cbind(alpha,Cspeed), col.names=FALSE, sep=" ")
abline(v = Cspeed, col='blue')

# d24 RATS
#fwrite(file="CS--d24-rat.tkz", x=cbind(alpha,Cspeed), col.names=FALSE, sep=" ")

# Write rational data for d8 and d24 together
d8 <- fread("d8/rats/CS--d8-rat.tkz")
d24 <- fread("d24/rats/CS--d24-rat.tkz")
X <- cbind(d8, d24[,2])
#fwrite(file="CS--d8d24-rat.tkz", x=cbind(d8,d24[,2]), col.names=FALSE, sep=" ")

mplot(alpha*180/pi, X[,2:3], type='p', pch=19, ylim=c(0,1))
a <- sort(c(seq(a[1], a[length(a)],length=1000),a))
gd8 <- sapply(2:3, function(x) 1/(1 + x*(1/2)*(1 + cos(a))))
mplot(a*180/pi, gd8, add=TRUE)


a <- alpha.irr
a <- alpha
mplot(a, Cspeed, type='o',ylim=c(0,1))

# g(alpha)
dy <- (1/2) * (1/tan(a/2))
d <- dy
r <- d/dy
a <- sort(c(seq(a[1], a[length(a)],length=1000),a))
lines(a, 1/(1 + 1*r[1]*(1 + cos(a))), type='l')
lines(a, 1/(1 + 2*r[1]*(1 + cos(a))), type='l')
lines(a, 1/(1 + 3*r[1]*(1 + cos(a))), type='l')
lines(a, 1/(1 + 4*r[1]*(1 + cos(a))), type='l')
abline(v=1/(1 + seq(2,7)*r[1]*(1 + cos(a[i]))))
abline(h=1/2)

fwrite(x=cbind(alpha, Cspeed[rats], Cspeed[-rats]), file="CSpeeds-d8.tkz", col.names=FALSE, sep=" ")
LinFit(mat=Cspeed, time=alpha, range=1:6, logs=FALSE)

axis(side=1, at=seq(0, 0.3, by=0.1/2))
axis(side=1, at=seq(0, 0.3, by=0.1/2))
id <- argmax(x = L$x, y = L$y, w = 2)#[c(1:6, 12, 30)]
idd <- id[which(L$y[id,] > 0.02)][-c(1:3, 12, 13, 15:28, 30:50, 52:54, 56:80)]
lines(L$x[idd,], L$y[idd,])
lines(L$x[idd,], exp(5.5)*L$x[idd,]^2.1)
abline(v=L$x[idd,])
y <- diff(L$x[idd,])

legend('topleft', legend=c(tex("0 < \\xi < 1000"), tex("1 < \\xi < 1000"), tex("\\xi > 1000")), lty='solid', col=rbow(), lwd=2, cex=2)
axis(side=1, at=seq(-3/4, 3/4, by=1/4))
dev.off()


abline(v=0.5 + c(1/40, -1/40), lty='dashed')
legend('topleft', legend=c("Special velocities removed", "Ordinary distribution of velocities"), col=rbow(), lty='solid', lwd=2, cex=1.75)
abline(v=0.6 + c(25/1e3, -25/1e3))

# Distribution of jump time/distance
lambda.t     <- ReadDataList(lff("", "Pos"), idx=1)
lambda.t     <- ReadDataList(lff("", "jump"), idx=1)
lambda.t <- data.matrix(lambda.t[which(lambda.t > 10),])
L <- getTSD_ecdf(dx = lambda.t, n=1000, domain=exp.domain)
mplot((L$x), (1-L$y), type='s', log='xy')
mplot(L$x, 1-L$y, type='s', add=TRUE, col='blue')

LinFit(mat = 1-L$y, time = L$x, logs=TRUE, range = 10:400, plot=TRUE, show.range=TRUE)

######
lambda.t     <- ReadDataList(lff("", "Pos"), idx=1)
lambda.t <- data.matrix(lambda.t[which(lambda.t > 10),])
L <- getTSD_ecdf(dx = lambda.t, n=1000, domain=exp.domain)
mplot((L$x), (1-L$y), type='s', log='xy')
mplot((L$x), (1-L$y), type='s', col='black', add=TRUE)
LinFit(mat = 1-L$y, time = L$x, logs=TRUE, range = 1:800, plot=TRUE, show.range=TRUE)


lambda.t     <- ReadDataList(lff("", "jump"), idx=1)
L <- getTSD_ecdf(dx = lambda.t, n=1000, domain=exp.domain)
mplot(L$x, 1-L$y, type='s', col=rbow(), log='xy', ylab=tex("\\Psi(\\tau)"), xlab=tex("\\tau"))
mplot(L$x, 1-L$y, type='s', add=TRUE, col='blue')
LinFit(mat = 1-L$y, time = L$x, logs=TRUE, range = 400:800, plot=TRUE, show.range=TRUE)
legend('bottomleft', legend=paste0("N = ", 10^seq(0,4)), col=c('black', rbow()), lwd=2, cex=2)

curve(x^(-.85), from=10, to = 1e8, add=TRUE)

lambda.x     <- ReadDataList(lff("", "jump"), idx=2)
L <- getTSD_hist(dx = lambda.x/lambda.t, n=1000)
mplot(L$x, L$y, type='s')

abline(v=c(0.49, 0.52))
id <- which(L$x[,1] < 0.51 & L$x[,1] > 0.49)
sum(L$y[id,1])

mplot(log10(L$x), log10(1-L$y), type='s', log='', xlim=c(2, 5), ylim=c(-5, 2), xaxs='i', yaxs='i')



LinFit(mat = log10(1-L$y), time = log10(L$x), logs=FALSE, range = 200:600, plot=TRUE, show.range=TRUE)
LinFit(mat = log10(1-L$y), time = log10(L$x), logs=FALSE, range = 1:600, plot=TRUE, show.range=TRUE)
mplot(log10(L$x), log10(1-L$y), type='s', log='',)# xlim=c(2, 5), ylim=c(-5, 2), xaxs='i', yaxs='i')
LinFit(mat = log10(1-L$y), time = log10(L$x), logs=FALSE, range = 100:300, plot=TRUE, show.range=TRUE)
LinFit(mat = log10(1-L$y), time = log10(L$x), logs=FALSE, range = 400:800, plot=TRUE, show.range=TRUE)
LinFit(mat = log10(1-L$y), time = log10(L$x), logs=FALSE, range = 1:500, plot=TRUE, show.range=TRUE)


# Speed corresponding to the two greatest probabilities
CJS <- sapply(1:ncol(L$y), function(i) {
		   nna <- length(which(!is.na(L$y[,i])))
		   m1 <- sort(L$y[,i])[nna]
		   L$x[which(L$y[,i] == m1),i]})
CJS <- tpose(CJS)

#pdf(paste0(img.path, "characteristic_speeds4.pdf"), width=15, height=12)
# r <- 1/2 -- red and blue circles, rational and irrational
mplot(alpha[1:15], abs(CJS[1:15]), type='p', xaxt='n', xlab=tex("\\alpha"), ylab=tex("\\nu"), cex=2, col=rbow()[1], ylim=c(0,1))
mplot(alpha[1:15], abs(CJS[16:30]), cex=2, type='p', col=rbow()[5], add=TRUE)
# r = 3/4
mplot(alpha[1:15], abs(CJS[1:15]), type='p', xaxt='n', xlab=tex("\\alpha"), ylab=tex("\\nu"), cex=2, col=rbow()[1], pch=0, add=TRUE)
mplot(alpha[1:15], abs(CJS[16:30]), cex=2, type='p', col=rbow()[5], pch=0, add=TRUE)

# ballistic speed for 4-bounce orbits
x <- seq(0, pi, length=100)
dy <- (1/2)*1/tan(alpha/2)
dy <- (1/2)*1/tan(x/2)

d1 <- (1/2)*dy
d2 <- (3/4)*dy
lines(x, 1/(1 + 4*d1*sin(x)), lwd=2, lty='dashed', col='green')
lines(x, 1/(1 + 4*d2*sin(x)), lwd=2, lty='dashed', col='black')

# x-axis
id <- c(1,5,7,8,9,11,15)
axis(side=1, at=alpha[id], labels=alpha.lab[id])#, line = 2) 
abline(h=1/2, lty='dashed')

# Inset
par(fig = c(0.58, 0.95, 0.12, 0.58), new=TRUE, cex=.75)
#pdf(paste0(img.path, "pi3pi22pi3_jspeed.pdf"), width=15, height=12)
L <- getTSD_hist(dx = lambda.x/lambda.t, n=50, breaks=(seq(-10,log(1), length=1000)))#, domain=exp.domain) 
L <- getTSD_hist(dx = lambda.x/lambda.t, n=200)
K <- 7:9
matplot(L$x[,K], L$y[,K], type='s', xlab=tex("\\xi / \\tau"), ylab="Density", col=rbow(), lwd=3, lty='solid', log='', xlim=c(-.7, .7))#, ylim=c(0,15))
K <- K + 15
matplot(L$x[,K], L$y[,K], type='s', xlab=tex("\\xi / \\tau"), ylab="Density", col=rbow(), lwd=3, lty='dashed', add=TRUE)
#legend('topleft', legend=c(alpha.lab[7:9], alpha.lab.irr[7:9]), col=rep(rbow()[1:3], 2), lty=c(rep('solid', 3), rep('dashed', 3)), lwd=3, cex=2)
dev.off()

# ---------------------------------------------------------

##
#  FIG. MOMENTS --- 
## 

img.path = "~/Dropbox/PHD/Papers/JumpStatistics/Images/"
setwd("../Studies/S-JumpStatistics/Data/")

# Moment scalings based on time dep. scaling function
getMu <- function(time, moment.files, z=0.95) { 
    mom <- vector(mode='list', length=length(moment.files))
    for (i in seq(length(moment.files))) {
	mom[[i]] <- colMeans(fread(moment.files[i]))
    }
    mom <- data.matrix(list.cbind.na(mom))
    print(dim(mom))
    mu <- getPowerScaling(t = time, M = mom, z = z)
    return(mu)
}

#Scalings based on non-linear fit
getMom <- function(moment.files) {
    mom <- vector(mode='list', length=length(moment.files))
    for (i in seq(length(moment.files))) 
	mom[[i]] <- colMeans(fread(moment.files[i]))
    mom <- data.matrix(list.cbind.na(mom))
    return(mom)
}

#Scalings based on non-linear fit
getScalings <- function(time, moment.files, p) {
    mom <- vector(mode='list', length=length(moment.files))
    for (i in seq(length(moment.files))) 
	mom[[i]] <- colMeans(fread(moment.files[i]))
    mom <- data.matrix(list.cbind.na(mom))
    coefs <- sapply(1:ncol(mom), function(x) {
	df <- data.frame(t=t, x=mom[,x])
	coef(nlsLM(x ~ b*t^a, data=df, start=list(a=p[x], b=1), control = nls.lm.control(maxiter=200)))
    })
}

# Initial velocity angles
Vx <- ReadDataList(lf("P0.dat"), idx=3)
Vy <- ReadDataList(lf("P0.dat"), idx=4)
th <- data.matrix(atan2(Vy, Vx))
L <- getTSD_ecdf(dx = th, n = 1000)
mplot(L$x, L$y, type='s')
x <- seq(-pi,pi, length=100)
lines(x, seq(length(x))/length(x))

# Jump velocities after filtering
lambda.t     <- ReadDataList(lff("p", "jump.dat"), idx=1)
lambda.x     <- ReadDataList(lff("p", "jump.dat"), idx=2)
L <- getTSD_hist(dx = lambda.x / lambda.t, n = 200)
mplot(L$x, L$y, type='s')

# Scratch - checking individual moments
mu <- getMu(time=t, gtools::mixedsort(lff("a23", "Moment_1")), z=0.95)
mom <- getMom(gtools::mixedsort(lff("a23", "Moment_1")))
g2 <- getScalings(t, gtools::mixedsort(lff("", "Moment_1")), p=rep(2, 27))

tau0 <- c(0,1,10,100)
mplot(tau0, g2[1,1:4], type='p', cex=2, ylim=c(1.8, 2))
mplot(tau0, g2[1,5:8], type='p', cex=2, col='blue', add=TRUE)
ErrorBars.y(x = tau0, g2[1,], sdv=g2[1,])

mplot(t, mu[,], ylim=c(1.5,2))#, add=TRUE)
mplot(t, mom[,1:2], log='xy')
LinFit(mat = mom[,1:4], time = t, range = seq(10, 100), plot=TRUE)
lines(t, 0.1*t^2)
abline(h=2)

T <- 1e6
t <- seq(0, T, length=T/1000+1)

# Fig-MSD

tau0 <- xi0 <- c(0,1,10,100,100)
g2.reg <- getScalings(t, gtools::mixedsort(lff("c3", "Moment_1")), p=rep(2, 3))[1,]
spar(mfrow=c(1,2))
a21.g2.tau <- getScalings(t, gtools::mixedsort(lff("c0", "a21")), p=rep(2, 4))[1,]
a22.g2.tau <- getScalings(t, gtools::mixedsort(lff("c0", "a22")), p=rep(2, 4))[1,]
a23.g2.tau <- getScalings(t, gtools::mixedsort(lff("c0", "a23")), p=rep(2, 4))[1,]

x <- cbind(tau0[1:4], a21.g2.tau, a22.g2.tau, a23.g2.tau)
x <- rbind(x, c(tau0[5], g2.reg))
print(xtable::xtable(x), file="tau0_msd_scalings.txt")

a21.g2.xi <- getScalings(t, gtools::mixedsort(lff("c1", "a21")), p=rep(2, 4))[1,]
a22.g2.xi <- getScalings(t, gtools::mixedsort(lff("c1", "a22")), p=rep(2, 4))[1,]
a23.g2.xi <- getScalings(t, gtools::mixedsort(lff("c1", "a23")), p=rep(2, 4))[1,]

x <- cbind(xi0[1:4], a21.g2.xi, a22.g2.xi, a23.g2.xi)
x <- rbind(x, c(tau0[5], g2.reg))
print(xtable::xtable(x), file="xi0_msd_scalings.txt")


p <- seq(1/2, 10, by=1/2)
pp <- seq(-10, 20, by=1/2)

# NEW IC DISTRIBUTION /T8
mu.rat <- getScalings(time=t, gtools::mixedsort(lff("a7", "Moment")), p=p)
mu.irr <- getScalings(time=t, gtools::mixedsort(lff("a22", "Moment")), p=p)
p0.irr <- data.matrix(fread(lff("a23", "P0")))
p0.rat <- data.matrix(fread(lff("a7", "P0")))

# NORMAL IC DISTRIBUTION /T10
reg.mu.rat <- getScalings(time=t, gtools::mixedsort(lff("a6", "Moment")), p)
reg.mu.irr <- getScalings(time=t, gtools::mixedsort(lff("a22", "Moment")), p)

spar(cx=2.5)
par(mfcol=c(1,2))
# RATIONAL -- Regular and new IC's
matplot(p, mu.rat[1,], pch=19, type='p', col='blue', cex=2, xaxs='i', yaxs='i', xlim=c(0,11), ylim=c(0,10), xlab=tex("p"), ylab=tex("\\gamma_p"))
matplot(p, reg.mu.rat[1,], pch=1, type='p', col='red', cex=2, lwd=2)#, add=TRUE)

FitMomentSpectrum(p, mu.rat[1,], force.ballistic=FALSE, offs=0, PLOT=TRUE)
FitMomentSpectrum(p, mu.rat[1,], force.ballistic=TRUE, offs=+0.2, PLOT=TRUE)

    df <- data.frame(x=p, y=reg.mu.rat[1,])
    LF1 <- lm(y ~ 0 + ., data=df, subset=1:5)
    LF2 <- LinFit(mat = reg.mu.rat[1,], time = p, range=5:20, logs = FALSE)
    lines(pp, coef(LF1)*pp, lty='dashed', col='black', lwd=2)
    lines(pp, LF2$ax * pp + LF2$b, lty='dashed', col='black', lwd=2)
    
    lines(pp, pp - 0.65)
    lines(pp, (mu.rat[1,4]/2)*pp)

    df <- data.frame(x=p, y=mu.rat[1,])
    LF2 <- LinFit(mat = mu.rat[1,], time = p, range=6:20, logs = FALSE)
    lines(pp, LF2$ax * pp + LF2$b, lty='dashed', col='black', lwd=2)
    lines(c(2,2), c(0,2), lty='dotted', lwd=2)
    lines(c(0,2), c(2,2), lty='dotted', lwd=2)

## Inset -- Initial distribution of angles - RATIONAL
{
    par(mfcol=c(1,2))
    par(fig=c(0.05, 0.3, 0.450, 0.955), new=TRUE)#, mgp=c(2,1,0))
    par(cex=0.5, mgp=c(2,.75,0))
    theta0 <- cbind(runif(n=nrow(p0.rat), -pi, pi), atan2(p0.rat[,4], p0.rat[,3]))
    L <- getTSD_ecdf(dx = theta0, n = 1000)
    matplot(L$x, L$y, type='s', xlab=tex("\\theta_0"), ylab="", col=c('red', 'blue'), lty='solid', lwd=2, xaxt='n')
    axis(side=1, at=c(-3, -2, -1, 1, 2, 3))
    abline(h=c(0,1))
}

# IRRATIONAL -- REGULAR AND NEW IC's
matplot(p, mu.irr[1,], pch=19, type='p', col='blue', cex=2, xaxs='i', yaxs='i', xlim=c(0,11), ylim=c(0,11), xlab=tex("p"), ylab=tex("\\gamma_p"))
matplot(p, reg.mu.irr[1,], pch=1, type='p', col='red', cex=2, lwd=2)#, add=TRUE)

FitMomentSpectrum(p, mu.irr[1,], force.ballistic=FALSE, offs=0, PLOT=TRUE)
FitMomentSpectrum(p, mu.irr[1,], force.ballistic=TRUE, offs=0.1, PLOT=TRUE)

    df <- data.frame(x=p, y=reg.mu.irr[1,])
    LF1 <- lm(y ~ 0 + ., data=df, subset=1:14)
    lines(pp, coef(LF1)*pp, lty='dashed', col='black', lwd=2)
    LF2 <- LinFit(mat = reg.mu.irr[1,], time = p, range=14:20, logs = FALSE)
    lines(pp, LF2$ax * pp + LF2$b, lty='dashed', col='black', lwd=2)

    LF2 <- LinFit(mat = mu.irr[1,], time = p, range=1:20, logs = FALSE)
    lines(pp, LF2$ax * pp + LF2$b, lty='dashed', col='black', lwd=2)
    lines(c(2,2), c(0,2), lty='dotted', lwd=2)
    lines(c(0,2), c(2,2), lty='dotted', lwd=2)

## Inset -- Initial distribution of angles
{
    par(mfcol=c(1,2))
    par(fig=c(0.55, 0.8, 0.450, 0.955), new=TRUE)#, mgp=c(2,1,0))
    par(cex=0.5, mgp=c(2,.75,0))
    theta0 <- cbind(runif(n=nrow(p0.rat), -pi, pi), atan2(p0.irr[,4], p0.irr[,3]))
    L <- getTSD_ecdf(dx = theta0, n = 1000)
    matplot(L$x, L$y, type='s', xlab=tex("\\theta_0"), ylab="", col=c('red', 'blue'), lty='solid', lwd=2, xaxt="n")
    axis(side=1, at=c(-3, -2, -1, 1, 2, 3), line=0)#, line = 0, pos = -0.1)
    abline(h=c(0,1))
}

## Inset -- mu(t)
par(fig=c(0.05, 0.28, 0.450, 0.955), new=TRUE, cex=0.75, mgp=c(3,1,0))
matplot(t, mu.rat[,4], type='l', col='blue', lwd=2, xlab=tex("t"), ylab=tex("\\mu(t)"), ylim=c(1.4,1.8))
matplot(t, reg.mu.rat[,4], type='l', col='red', cex=2, lwd=2, add=TRUE)
abline(h=1.66)

# Inset
lambda.t     <- ReadDataList(lff("p", "jump.dat"), idx=1)
lambda.x     <- ReadDataList(lff("p", "jump.dat"), idx=2)
L <- getTSD_ecdf(dx = lambda.x[,]/lambda.t[,], n=200)#, domain=exp.domain)

par(fig=c(0.08, 0.48, 0.450, 0.955), new=TRUE, cex=0.75, mgp=c(3,1,0))
matplot(L$x[,], 1-L$y[,], col=c('red', 'blue'),
	type='s', 
	lty='solid', lwd=3, 
	ylab="", 
	xlab=tex("\\xi / \\tau"))

# ---------------- CHANGING IC MSD ----------------
tx <- c(0, 1, 10, 100, 120)
T <- 1e6
t <- seq(0, T, length=T/1000+1)

# PI/2
f <- c(list.files(path = "PI2_IRR/xi0_tau0_1", "Moment_3", full.names = TRUE), 
       list.files(path = "PI2_IRR/xi0_tau0_10", "Moment_3", full.names=TRUE),
       list.files(path = "PI2_IRR/xi0_tau0_100", "Moment_3", full.names=TRUE))
mu <- getMu(time=t, moment.files=f[-c(1,3,5)])
mplot(t, mu, add=TRUE)
gamma.2 <- getScalings(time=t, moment.files=f, p=rep(2,length(f)))
mplot(tx, c(1.999, gamma.2[1,c(1,3,5)], 1.81), 
      xlab=tex("\\xi_0 (red), \\tau_0 (blue)"), ylab=tex("\\gamma_2"),
      cex=2, type='p', xaxt='n',
      ylim=c(1.8, 2), xlim=c(0,120)) # tau_0 = tx, xi_0 = free -> e0
axis(side=1, at=c(1,10,100,120), labels=c(1,10,100,TeX(paste0("\\infty"))))
mplot(tx, c(1.999, gamma.2[1,-c(1,3,5)],1.81), 
      cex=2, type='p', col='blue', pch=1, add=TRUE) # tau_0 = free, xi_0 = tx -> e1
abline(h=1.81)

# 2PI/3
f <- c(list.files(path = "2PI3_IRR/xi0_tau0_1", "Moment_3", full.names = TRUE), 
       list.files(path = "2PI3_IRR/xi0_tau0_10", "Moment_3", full.names=TRUE),
       list.files(path = "2PI3_IRR/xi0_tau0_100", "Moment_3", full.names=TRUE))
mu <- getMom(moment.files=f)
mplot(t, mu, log='xy')
gamma.2 <- getScalings(time=t, moment.files=f, p=rep(2,length(f)))
mplot(tx, c(gamma.2[1,c(1,3,5)], 1.75), 
      xlab=tex("\\xi_0 (red), \\tau_0 (blue)"), ylab=tex("\\gamma_2"),
      cex=2, type='p', xaxt='n',
      ylim=c(1.8, 2), xlim=c(0,120)) # tau_0 = tx, xi_0 = free -> e0
axis(side=1, at=c(1,10,100,120), labels=c(1,10,100,TeX(paste0("\\infty"))))
mplot(tx, c(gamma.2[1,-c(1,3,5)],1.81), 
      cex=2, type='p', col='blue', pch=1, add=TRUE) # tau_0 = free, xi_0 = tx -> e1
abline(h=1.81)



# Initial angle distribution
f <- c(list.files(path = "xi0_tau0_1",  "P0", full.names=TRUE), 
       list.files(path = "xi0_tau0_10", "P0", full.names=TRUE),
       list.files(path = "xi0_tau0_100","P0", full.names=TRUE))
Vx <- ReadDataList(f, idx=3)
Vy <- ReadDataList(f, idx=4)
th <- data.matrix(atan2(Vy, Vx))
L <- getTSD_ecdf(dx = th[,c(1,3,5)], n = 100)
L <- getTSD_hist(dx = th[,-c(1,3,5)], n = 1000)
mplot(L$x, L$y, type='s')

# Jump velocity distribution
f <- c(list.files(path = "xi0_tau0_1",  "jump", full.names=TRUE), 
       list.files(path = "xi0_tau0_10", "jump", full.names=TRUE),
       list.files(path = "xi0_tau0_100","jump", full.names=TRUE))
lambda.t <- ReadDataList(f, idx=1)
lambda.x <- ReadDataList(f, idx=2)
L <- getTSD_hist(dx = lambda.x/lambda.t, n = 200)
mplot(L$x, L$y, type='s')
# ---------------------------------------------------------


# ---------------- COARSE GRAINED APPROACH ----------------
img.path = "~/Dropbox/PHD/Papers/JumpStatistics/Images/"
setwd("../Studies/S-JumpStatistics/Data/")

cnt     <- ReadDataList(lff("", "_jump.dat"), idx=1)
t0      <- ReadDataList(lff("", "w50_jump.dat"), idx=2)
x0      <- ReadDataList(lff("", "w50_jump.dat"), idx=3)
jt      <- ReadDataList(lff("", "w10_jump.dat"), idx=4)
jt      <- ReadDataList(lff("", "_jump.dat"), idx=4)
jd      <- ReadDataList(lff("", "_jump.dat"), idx=5)
id <- apply(cnt, 2, function(x) which(x == 0))

num.ics <- apply(cnt, 2, function(x) length(which(x == 0)))
ranges <- lapply(id, function(x) sapply(1:(length(x)-1), function(y) seq(x[y], x[y+1]-1)))
nJumpPerIC <- lapply(ranges, function(x) lapply(x, function(y) length(y))) # Number of jumps per IC for each window

W <- 10^seq(1,4)
W <- c(rep(5,3), rep(50,3), rep(500,3), rep(5000,2))

# Distributions
L <- getTSD_hist(dx = jt, n=200, domain=exp.domain)
mplot(L$x, L$y, type='s', log='xy',)# ylim=c(0,10))

id <- apply(abs(jd), 2, function(x) which(x > 100))
jd0 <- list.cbind.na(sapply(1:length(id), function(x) jd[-id[[x]],x]))
jt0 <- list.cbind.na(sapply(1:length(id), function(x) jt[-id[[x]],x]))
L <- getTSD_ecdf(dx = jt0[,], n=1000, domain=exp.domain, t=1)#W[-4])

L <- getTSD_ecdf(dx = jt, n=1000, domain=exp.domain)

ids <- paste0("a", c(6:8, 21:23))
datx <- daty <- list()
for (i in seq(length(ids))) {
    jt <- ReadDataList(lff(ids[i], "_jump.dat"), idx=4)

    L <- getTSD_ecdf(dx = jt, n=1000, domain=exp.domain)
    mplot(L$x, 1-L$y, type='s', add=TRUE)
    LinFit(time=L$x, mat=1-L$y, range=40:370, plot=TRUE)
    #datx[[i]] <- L$x
    #daty[[i]] <- L$y
}
datx <- list.cbind.na(datx)
daty <- list.cbind.na(daty)
#fwrite(file="JTD--WeakBall.tkz", x=cbind(datx, daty), sep=" ", col.names=FALSE)


mplot(L$x[,], 1-L$y[,], type='s', log='xy', )
      #xlim=c(1000,3e4),
      #ylim=c(0.01, 1))
LinFit(time=L$x, mat=1-L$y, range=40:370, plot=TRUE, show.range=TRUE)
LinFit(time=L$x, mat=L$y, range=10:150, plot=TRUE, show.range=TRUE)

# Reconstruction
t    <- seq(0,1e6, length=1e6/100 + 1)
t    <- data.matrix(fread(lf("w500_1--t.dat")))
x    <- data.matrix(fread(lf("w500_1--x.dat")))
T    <- apply(rbind(t0[1,], jt), 2, function(x) cumsum(x))
X    <- apply(rbind(x0[1,], jd), 2, function(x) cumsum(x))

s <- seq(100000)
mplot(t[s], x[s], lwd=0.5)#, xlim=c(1,1e6))

mm <- getMM(x = t, y = x, w = 5)
mm$df
ShowJTw(df=mm$df, eps = 0)
ShowJT(t0,x0,jt,jd)

j <- 1
lines(T[,j], X[,j])
points(T[,j], X[,j], pch=19, col='green', cex=0.5)

##
#  FIG. MOTIVATING EXAMPLE -- 
## 
img.path = "~/Dropbox/PHD/Papers/JumpStatistics/Images/"
setwd("../Studies/S-JumpStatistics/Data/PaperPlots/Fig-CGMotivation")

    x <- data.matrix(fread(lf("x.dat")))
    t <- data.matrix(fread(lf("t.dat")))
    jump <- data.matrix(fread(lf("N0_jump.dat")))
    #pdf(paste0(img.path, "cg_motivation.pdf"), width=15, height=10)
    spar(mfrow = c(1,1))
    mplot(t, x, spa=FALSE, xlab="t", ylab="x(t)", col='black', lwd=2)
    mm <- getMM(x = t, y = x, w = 1000)
    ShowJTw(df = mm$df, pt.cx = 0.1)
    
    s <- seq(420000, 440000)
    rect(xleft = t[s[1]], xright=t[s[length(s)]], ybottom=min(x[s]), ytop=max(x[s]), lty='solid', lwd=1, border='red')
    
    spar(fig = c(0.1, 0.485, 0.45, 0.95), cex=0.75, new=TRUE)
    mplot(t[s], x[s], spa=FALSE, xlab="t", ylab="x(t)", col='black', lwd=1, xaxs='i', yaxs='i')
    
    jt <- cumsum(c(0,jump[,1]))
    jd <- cumsum(c(x[1,],jump[,2]))
    odd.pnts <- seq(1, length(jt), by=2)
    even.pnts <- seq(2, length(jt)-1, by=2)
    
    points(jt[odd.pnts], jd[odd.pnts], cex=0.2, col='red', pch=19)
    points(jt[even.pnts], jd[even.pnts], cex=0.2, col='blue', pch=19)
    
    dev.off()


# ---------------- COMPARISON WITH FLOW METHOD ----------------
img.path = "~/Dropbox/PHD/Papers/JumpStatistics/Images/"
setwd("../Studies/S-JumpStatistics/Data/PaperPlots/Fig-CGFlowCompare/")

cnt     <- ReadDataList(lff("", "w5_jump.dat"), idx=1)
t0      <- ReadDataList(lff("", "_jump.dat"), idx=2)
x0      <- ReadDataList(lff("", "_jump.dat"), idx=3)
jd      <- ReadDataList(lff("c0", "w5"), idx=5)
jt      <- ReadDataList(lff("c0", "w"), idx=4)
jt.flow <- ReadDataList(lff("a", "--jump"), idx=1)
jd.flow <- ReadDataList(lff("a", "--jump"), idx=2)

id <- apply(cnt, 2, function(x) which(x == 0))
num.ics <- apply(cnt, 2, function(x) length(which(x == 0)))
ranges <- lapply(id, function(x) sapply(1:(length(x)-1), function(y) seq(x[y], x[y+1]-1)))
nJumpPerIC <- lapply(ranges, function(x) lapply(x, function(y) length(y))) # Number of jumps per IC for each window

# Filter small distances / times
l <- xi0.filter(jt.flow, jd.flow, 100000)
jt0 <- l$jt
jd0 <- l$jd
apply(jt0, 2, function(x) length(which(is.na(x))) == length(x))

# Jump Velocity distribution
L <- getTSD_hist(dx = (jd.flow / jt.flow), n=1000)
mplot(L$x, L$y, type='s', log='', ylim=c(0,10))# ylim=c(0,100), xlim=c(0,1))
mplot(L$x, L$y, type='s', log='', lty='dashed', add=TRUE)# ylim=c(0,100), xlim=c(0,1))

getSV <- function(x, y) {
    v <- c()
    for (i in seq(ncol(y))) 
	v <- c(v, x[which.max(y[,i]),i])
    return(v)	
}

v <- getSV(abs(L$x), L$y)
mplot(v[-1], type='o', log='',)# ylim=c(0.499, 0.5))

# Jump time distribution 
L <- getTSD_ecdf(dx = abs(jt.flow)[,], n=1000, domain=exp.domain)
mplot(L$x[,], (1-L$y[,]), #add=TRUE,
      type='s', log='xy', lwd=3, #xaxt='n',
      #xlim=c(1e0,1e7), ylim=c(1e-6, 1),
      xlab=tex("\\tau"), ylab=tex("\\Psi(\\tau)"))

# Multi scaling power law
X <- L$x[,]
Y <- 1-L$y[,]
LF1 <- LinFit(mat = Y, time=X,  plot=TRUE, range=300:850,   show=TRUE)
LF1 <- LinFit(mat = Y, time=X,  plot=TRUE, range=1:450,   show=TRUE)
LF2 <- LinFit(mat = Y, time=X,  plot=TRUE, range=650:780, show=TRUE)

Lf <- getTSD_hist(dx = (jt.flow), n=1000, domain=exp.domain)
Lf <- getTSD_ecdf(dx = abs(jt.flow), n=1000, domain=exp.domain)
mplot(Lf$x[,], 1-Lf$y[,], type='s', col='blue', add=TRUE)
mplot(Lf$x[,], 1-Lf$y[,], type='s', log='xy', )

X <- Lf$x[,]
Y <- 1-Lf$y[,]
LFf1 <- LinFit(mat = Y, time=X,  plot=TRUE, range=1:150,   show=FALSE)
LFf2 <- LinFit(mat = Y, time=X,  plot=TRUE, range=350:600, show=FALSE)


clf <- getCrossoverPoints(LF1, LF2)
clff <- getCrossoverPoints(LFf1, LFf2)
clf$inty  <- as.numeric(format(clf$inty, digits=1))
clff$inty <- as.numeric(format(clff$inty, digits=1))
tab <- cbind(LF1$ax, LF2$ax, clf$intx, LFf1$ax, LFf2$ax, clff$intx)
print(xtable::xtable(tab, digits = 2), file=paste0(img.path, "../scales-flow-cg-jt.txt"))


# Plot sequence for flow / cg comparison (jump times)
fn <- c("pi3", "pi2", "2pi3", "ipi3", "ipi2", "i2pi3") # Individual plots
fn <- c("irpi3", "irpi2", "ir2pi3") # Side by side panel of rational/irrational
for (j in seq(6)) {
    pdf(paste0(img.path, fn[j], "-flow_cg_jt.pdf"), height=10, width=20)
    spar(mfrow=c(1,2))
    j <- 2
    {
	# CG method
	L <- getTSD_ecdf(dx = abs(jt)[,j], n=1000, domain=exp.domain)
	mplot(L$x[,], (1-L$y[,]),spa = FALSE,  
	      type='s', log='xy', lwd=4, xaxt='n',
	      xlim=c(1e3,1e7), ylim=c(1e-5, 1),
	      xlab=tex("\\tau"), ylab=tex("\\Psi(\\tau)"))

	# CG Fits
	LF1 <- LinFit(mat = 1-L$y, time=L$x,  plot=FALSE, range=1:450,   show=FALSE, lty='dashed', col='blue')
	LF2 <- LinFit(mat = 1-L$y, time=L$x,  plot=FALSE, range=450:750, show=FALSE, lty='dotdash', col='blue')
	clf <- ConnectLinFits(100, 1e7, LF1, LF2, extend=1000)
	mplot(clf$lx, clf$ly, lty=c('dashed', 'dotdash'), lwd=3, col='blue', add=TRUE, spa=FALSE)
	dvline(clf$intx, y0=1e-7, y1=clf$inty, lty='dashed', lwd=1) 

	# Flow method
	Lf <- getTSD_ecdf(dx = data.matrix(jt.flow[,j]), n=1000, domain=exp.domain)
	mplot(Lf$x[,], 1-Lf$y[,], type='s', col='blue', lwd=4, add=TRUE, spa=FALSE)

	# Flow fits
	LFf1 <- LinFit(mat = 1-Lf$y, time=Lf$x,  plot=FALSE, range=1:150,   show=FALSE)
	LFf2 <- LinFit(mat = 1-Lf$y, time=Lf$x,  plot=FALSE, range=350:600, show=FALSE)
	clf1 <- ConnectLinFits(100, 1e8, LFf1, LFf2, extend=5000)
	mplot(clf1$lx, clf1$ly, lty=c('dashed', 'dotdash'), lwd=2, col='red', add=TRUE, spa=FALSE)
	dvline(clf1$intx, y0=1e-7, y1=clf1$inty, lty='dashed', lwd=1) 
	xtlab <- format(sort(c(1e3, 1e5, 1e6, 1e7, min(clf$intx, clf1$intx))), digits=1)
	axis(side=1, at=xtlab)


	rates <- format(c(LF1$ax, LF2$ax, LFf1$ax, LFf2$ax), digits=2)
	legend('topright', legend=TeX(paste0("$\\tau^{", rates, "}")), col=c(rep('blue',2), rep('red',2)), lty=c('dashed', 'dotdash', 'dashed', 'dotdash'), cex=2, lwd=3)
    }
    dev.off()
}




mplot(Lf$x[,], 1-Lf$y[,], type='s', log='y', add=TRUE, lty='dashed')
f <- function(x, a, b) { exp(-(x/a)^b) }
df <- data.frame(x=L$x[,], y=1-L$y[,])
NLF <- nlsLM(y ~ 1-f(x, a, b), data=df, start=list(a=200, b=2.0), subset=46:1000,
	     control = nls.lm.control(maxiter = 200))
X <- L$x[46:1000,]
cf <- coef(NLF)
print(cf)
lines(X, (1-f(X, a=cf[1], b=cf[2])), lwd=3, lty='dashed')#, c=cf[3])))
lines(X, (1-f(X, a=cf[1], b=cf[2])), lwd=3, lty='dashed')#, c=cf[3])))
# ---------------------------------------------------------


### 
#   JUMP TIME DISTRIBUTION AS A FUNCTION OF XI0 -- Flow method
### 
jt.flow <- ReadDataList(lff("a", "--jump"), idx=1)
jd.flow <- ReadDataList(lff("a", "--jump"), idx=2)

xi0 <- c(seq(1,1000,by=100), seq(1000, 5000, by=250))
xi0 <- seq(1, 500, by=10)
xi0.filter <- function(jt, jd, xi0) {
    jtf <- jdf <- list()
    for (i in seq(ncol(jt))) {
	id <- which(abs(jd[,i]) < xi0)
    	if (length(id) > 0) {
    	    jtf[[i]] <- jt[-id,i]
    	    jdf[[i]] <- jd[-id,i]
    	} else {
    	    jtf[[i]] <- jt[,i]
    	    jdf[[i]] <- jd[,i]
	}
    }
    jtf <- list.cbind.na(jtf)
    jdf <- list.cbind.na(jdf)
    return(list(jt=jtf, jd=jdf))
}

j <- 3
njt <- njd <- list()
for (i in seq(length(xi0))) {
    l <- xi0.filter(jt.flow[,j], jd.flow[,j], xi0[i])
    njt[[i]] <- l$jt
    njd[[i]] <- l$jd
}
njt <- list.cbind.na(njt)
njd <- list.cbind.na(njd)

Lf <- getTSD_ecdf(dx = data.matrix(njt), n=5000, domain=exp.domain,)
mplot(Lf$x[,], 1-Lf$y[,], type='s', log='xy', 
      ylim=c(1e-1, 1), 
      xlim=c(1e3, 5e3))
LF <- LinFit(mat = 1-Lf$y[,], time=Lf$x[,],  plot=TRUE, range=300:600, show=TRUE)
LF <- LinFit(mat = 1-Lf$y[,], time=Lf$x[,],  plot=TRUE, range=800:2000, show=TRUE)
mplot(xi0, 3+LF$ax, type='o')
