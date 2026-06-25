
img.path = "~/Dropbox/PHD/Papers/JumpStatistics/Images/"
setwd("../Studies/S-JumpStatistics/Data/")

####' 
##   PATH: S-JumpStatistics/Data/
##   Distribution of jump times and jump distances
##    
####' 
f <- gtools::mixedsort(lf("-d2--w10_jump.dat"))[1]
jump <- data.matrix(fread(f[1]))
jump <- ReadDataList(FileList = lf("d2--w10_jump"), idx=2)

id <- c(1, which(diff(jump[,1]) < 0))
ranges <- sapply(1:(length(id)-1), function(x) seq(id[x], id[x+1]))

times <- lapply(ranges, function(x) cumsum(c(jump[x[1],2], jump[x,2])))
posit <- lapply(ranges, function(x) cumsum(c(jump[x[1],3], jump[x,3])))

{
    s <- sample(1:length(times), size=5)
    #pdf(paste0(img.path, "x_w.pdf"), width=12, heigh=12)
    spar(cx=2.15, col='white')
    par(mgp=c(3,1,0), cex.main=2.15)
    tt <- list.cbind.na(sapply(s, function(x) times[[x]]))
    xx <- list.cbind.na(sapply(s, function(x) posit[[x]]))
    matplot(tt, xx, type='o', pch=19, ylab=tex("x_w(t')"), xlab=tex("t'"), lty='solid', col=rbow(), xaxt='n', xlim=c(0,8e5), lwd=2, cex=1.25)
    axis(side=1, at=seq(0, 1e6, by=1e5))
    abline(v=0, h=0, lwd=1.5)
    dev.off()
}

# Collect jump distances within a certain range of times
mid <- seq(1e4, 1e6, length=50)
eps <- rep(25000, length(mid))
eps <- c(mid[1]/2, diff(mid)/2)
X <- T <- X.std <- list()
for (j in seq(length(mid))) {
    x <- t <- x.std <- c()
    for (i in seq(length(times))) {
        ids <- which(times[[i]] > mid[j] - eps[j] & times[[i]] < mid[j] + eps[j])

        if (length(ids) > 0) {
	   x <- c(x, mean(posit[[i]][ids]))
	   t <- c(t, mean(times[[i]][ids]))

	   if (length(ids) > 1) {
		x.std <- c(x.std, sqrt(var(posit[[i]][ids])))
	   } else {
		x.std <- c(x.std, 0)
	   }
	}
    }
    X[[j]] <- x
    X.std[[j]] <- x.std
    T[[j]] <- t
}
X <- list.cbind.na(X)
T <- list.cbind.na(T)
X.std <- list.cbind.na(X.std)

# Distribution of compressed trajectories
{
    s <- unique(round(exp(seq(log(1), log(nrow(T)), length=400))))
    pdf(paste0(img.path, "xw_dist.pdf"), width=18, height=13)
    spar(cx=2.5, col='white')
    cols <- rainbow(10)
    matplot(T[s,], X[s,], pch=19, cex=0.7, log='', xlab="t", ylab=tex("x_w(t)"), 
	    ylim=c(-12e5, 5e5), col=cols)

    hid <- c(2, 4, 8)
    L <- getTSD_hist(dx = X[,hid], n = 120, plot=FALSE, logs='y')
    
    par(fig=c(0.055, 0.525, 0.075, 0.62), new=TRUE, cex=0.75)
    matplot(L$x[,], L$y[,], type='s', col=cols[hid], lty='solid', lwd=2, log='y', ylim=c(3e-7, 1.5e-5),
	    xlab=tex("x_w"), ylab=tex("P(x_w, t')"))

    dev.off()
}


# MSD
{
    apply(T, 2, function(x) sqrt(var(x, na.rm=TRUE)))
    apply(X, 2, function(x) sqrt(var(x, na.rm=TRUE)))
    Tm <- colMeans(T, na.rm=TRUE)
    msd <- colMeans(X^2, na.rm=TRUE)
    msdp <- colMeans((X + X.std)^2, na.rm=TRUE)
    msdm <- colMeans((X - X.std)^2, na.rm=TRUE)
    
    matplot(Tm, msd, log='xy', type='o', pch=19)
    lines(Tm, exp(1.5)*Tm^(1.66), type='l', lty='dashed', col='red')
    lines(T, exp(1.5)*T^(1.59), type='l', lty='dashed', col='red')
    
    spar(cx=2, col='white')
    mu <- getPowerScaling(t = Tm, M = msd, z = 0.8)
    mu <- getPowerScaling(t = Tm, M = msdp, z = 0.8)

    plot(Tm, mu, type='o', pch=19, ylim=c(0.5, 2), ylab=tex("\\gamma(t'; z = 0.8)"), xlab="t")
    abline(h=5/3)
}


sc <- sapply(mid, function(x) mean(seq(x - eps, x + eps)^(5/6)))
L <- getTSD_hist(dx = X, n=200, plot=FALSE, t = mid, eta = 0.9)

matplot(L$x[,], L$y[,], type='s', col=rbow(), lty='solid', lwd=2, log='y', ylim=c(5e-8, 5e-5))


LF <- LinFit(time=L$x[,1], mat=log(L$y[,1]), range = seq(110, 180), show.range=TRUE, logs=FALSE)
lines(L$x[,1], exp(LF$b[1])*L$x[,1]^LF$ax[1], lwd=3, col='magenta', lty='dashed')


Ljt <- getTSD_hist(dx = jump, n=200, plot=TRUE, domain=exp.domain)


spar(cx=2)
Ljt <- getTSD_ecdf(dx = jump, n=200, domain=exp.domain)
Ljt <- getTSD_hist(dx = jump[,2], n=200, plot=TRUE, domain=exp.domain)
Ljt <- getTSD_ecdf(dx = data.matrix(jump[,2]), n=200, domain=exp.domain)
Ljd <- getTSD_hist(dx = abs(jump[,3]), n=200, plot=TRUE, domain=exp.domain)
L <- Ljt


# ---- Jump time plot with insets ----
jump <- ReadDataList(FileList = lf("d2--w10_jump"), idx=2)

matplot(L$x[,], 1-L$y[,], type='s', col=rbow(), lty='solid', lwd=2, log='xy',
	xlab=tex("\\xi / \\tau"),
    	ylab="Density")


LF <- LinFit(time=Ljt$x[,1], mat=1-Ljt$y[,1], range = seq(1,70), show.range=TRUE, logs=TRUE)
lines(Ljt$x[,1], exp(LF$b[1])*Ljt$x[,1]^LF$ax[1], lwd=3, col='magenta', lty='dashed')


mpd <- data.matrix(fread("a1-d2--w10_MaxPD.dat", fill=TRUE))
mpd <- ReadDataList(lf("d2--w10_MaxPD.dat"), idx=1, written.as.row = TRUE)

ttc <- c(5e4, 250000, 500000)

L <- getTSD_ecdf(dx = mpd, n=500, domain=exp.domain)

spar(cx=2.15, col='white')
par(mgp=c(3.5,1.25,0), cex.main=2.15)
matplot(L$x[,], 1-L$y[,], type='s', lwd=2, lty='solid', 
	xlab=TeX("$max \\{ D_{p} \\}$"), 
	ylab="Prob.",log='x', col=rbow())

legend('topright', legend=TeX(paste0("$\\tau \\in \\[", ttc-1e4, ", ", ttc+1e4, "\\]$")), 
       col=rbow(), lwd=2, cex=2)
	#main=TeX(paste0("$\\tau \\in \\[", tau.to.check[slice]-25, ", ", tau.to.check[slice]+25, "\\]$")), 
	#log='x', col=rainbow)



# ---------------------------------------



pl <- power.law(alpha = 2, x0=1000, n = 100000)
L <- getTSD_hist(dx = 2*pl, n=200, plot=FALSE, domain=exp.domain)

x <- seq(1000, 5e5, by=100)
lines(x, exp(15)*x^(-3))
lines(L$x, L$y, type='s', col='red')




{
    #pdf(paste0(img.path, "jump_speed.pdf"), width=14, height=12)
    spar(cx=2.5, col='white', mfrow=c(2,2))
    par(mar = c(8, 8, 1.3, 1.3), oma = c(2, 2, 2, 2), mgp=c(5,2,0))
    
    matplot(L$x[,], L$y[,], type='s', xlim=c(-.5, .5), ylim=c(0,15),
    	col=rbow(), lty='solid', lwd=2, 
    	xlab=tex("\\xi / \\tau"),
    	ylab="Density", 
    	xaxt='n')
    axis(side=1, at=c(-0.5, -0.25, 0, 0.25, 0.5))
    legend('topright', legend=TeX(paste0("$\\alpha = ", c("\\pi /2", "\\pi /3")[2])),
           col=rbow()[1:2], lty='solid', cex=2, lwd=2.5)
    
    matplot(L$x[,2], L$y[,2], type='s', ylim=c(0,15),
    	col=rbow()[2], lty='solid', lwd=2, 
    	xlab=tex("\\xi / \\tau"),
    	ylab="Density", 
    	xaxt='n')
    axis(side=1, at=c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75))
    legend('topright', legend=TeX(paste0("$\\alpha = ", c("\\pi /2", "\\pi /3")[1])),
           col=rbow()[2], lty='solid', cex=2, lwd=2.5)
    
    trm <- which(jt[,1] > 10000) 
    matplot(jt[trm,1], jd[trm,1], pch=19, cex=0.5, col=rbow()[1], 
	    xlab=tex("\\tau"), ylab=tex("\\xi"), ylim=c(-5e4, 5e4), xlim=c(0,5e5))
    trm <- which(jt[,2] > 10000) 
    matplot(jt[trm,2], jd[trm,2], pch=19, cex=0.5, col=rbow()[2], 
	    xlab=tex("\\tau"), ylab=tex("\\xi"), ylim=c(-5e4, 5e4), xlim=c(0,5e5))

    dev.off()
}

####' 
##   PATH: S-JumpStatistics/Data/window_compare/
##   Distribution of jump times 
##    
####' 
f <- gtools::mixedsort(lf("_jump"))
f <- gtools::mixedsort(lf("_tau"))
f <- gtools::mixedsort(lf("_jd"))
f <- gtools::mixedsort(lf("-counts-"))
tau <- list()
for (i in seq(length(f)))
    tau[[i]] <- data.matrix(fread(f[i]))
tau[[i+1]] <- data.matrix(fread("p0-a1-d20---t3.dat"))#[1:1500000,]
tau <- list.cbind.na(tau)

jd <- 
jd <- tau

sjd <- apply(jd, 2, sort, na.last = TRUE)
sjt <- apply(tau, 2, sort, na.last=TRUE)

plot(sjt[1:100000,], sjd[1:100000,]^2, type='l')



L <- getTSD_hist(dx=jd/tau, n = 100, prune=FALSE)
L <- getTSD_hist(dx=abs(tau), n = 100, prune=FALSE, domain=exp.domain)
L <- getTSD_hist(dx=tau, n = 300, prune=FALSE, domain=function(dt, t, n) seq(101, 10000, by=20))#, domain=exp.domain)
L <- getTSD_hist(dx=tau[-which(tau == 0, arr.ind=TRUE),], n = 1000, prune=FALSE)#, domain=exp.domain)

#pdf(paste0(img.path, "jd_dist.pdf"), width=14, height=12)
spar(cx=2.25, col='white')
par(mgp=c(3.25,1,0))
matplot(L$x[,], L$y[,], type='s', lty='solid', log='', lwd=2, xlab=tex("\\tau"), ylab=tex("\\Psi(\\tau)"), col=rainbow(10)[-5])

matplot(L$x[,], L$y[,], type='s', lty='solid', log='xy', lwd=2, xlab=TeX("$| \\xi |$"), ylab=TeX("$\\psi(| \\xi |)$"), col=rainbow(10)[-5], yaxt="n")
axis(side=2, at=c(1e-11, 1e-9, 1e-7, 1e-5, 1e-3, 1e-1)) 


dev.off()


#
s <- 1:10000
matplot(tau[s,1], jd[s,1] / tau[s,1], type='l')



B <- exp.domain(abs(tau[,3]), 1, n=200)
h <- hist(abs(tau[,3]), breaks=B, plot=FALSE)
mids <- h$mids
d <- h$density
lines(mids, d, type='s', log='xy')
X <- seq(5e6)
f <- ecdf(tau[,4])(X)
plot(X, 1-f, type='s', log='xy')


LF <- LinFit(time=L$x[,1], mat=L$y[,1], range = seq(40,95), show.range=TRUE, logs=TRUE)
lines(L$x[,1], exp(LF$b[1])*L$x[,1]^LF$ax[1], lwd=3, col='magenta', lty='dashed')
LF <- LinFit(time=L$x[,2], mat=L$y[,2], range = seq(60,99), show.range=TRUE, logs=TRUE)
lines(L$x[,2], exp(LF$b[1])*L$x[,2]^LF$ax[1], lwd=3, col='magenta', lty='dashed')
LF <- LinFit(time=L$x[,3], mat=L$y[,3], range = seq(60,95), show.range=TRUE, logs=TRUE)
lines(L$x[,3], exp(LF$b[1])*L$x[,3]^LF$ax[1], lwd=3, col='magenta', lty='dashed')

x <- seq(1, 1e7)
lines(x, exp(13.2)*x^(-3), lty='dashed', lwd=3)
lines(x, exp(9.6)*x^(-2.7), lty='dashed', lwd=3)

legend('topright', legend=c(paste0("w = ", c(10, 50, 200, 400)), tex("\\tau = s"), tex("\\tau^{-3}")), 
			    col=c(rainbow(10)[-5][1:5], 'black'), 
			    lty=c(rep('solid' ,5), 'dashed'), 
			    cex=3, lwd=3)

B <- seq(1, nrow(tau))
Q <- apply(tau, 2, function(x) x/sum(x,na.rm=TRUE))
matplot(B, Q, type='s', log='xy', )#xlim=c(0,1e5), ylim=c(0,0.005))


####' 
##   PATH: S-JumpStatistics/Data/window_compare/
##   Distribution of jump distances
##    
####' 
f <- gtools::mixedsort(lf("_jd"))[-1]
jd <- list()
for (i in seq(length(f)))
    jd[[i]] <- data.matrix(fread(f[i]))
jd <- list.cbind.na(jd)

L <- getTSD_hist(dx=abs(jd), n = 100, prune=FALSE, domain=exp.domain)

spar(cx=3, col='white')
matplot(L$x[,], L$y[,], type='s', lty='solid', log='xy', lwd=2, xlab=tex("\\tau"), ylab=tex("\\Psi(\\tau)"), col=rainbow(10)[-5])

LF <- LinFit(time=L$x[,1], mat=L$y[,1], range = seq(40,95), show.range=TRUE, logs=TRUE)
lines(L$x[,1], exp(LF$b[1])*L$x[,1]^LF$ax[1], lwd=3, col='magenta', lty='dashed')
LF <- LinFit(time=L$x[,2], mat=L$y[,2], range = seq(60,99), show.range=TRUE, logs=TRUE)
lines(L$x[,2], exp(LF$b[1])*L$x[,2]^LF$ax[1], lwd=3, col='magenta', lty='dashed')
LF <- LinFit(time=L$x[,3], mat=L$y[,3], range = seq(60,95), show.range=TRUE, logs=TRUE)
lines(L$x[,3], exp(LF$b[1])*L$x[,3]^LF$ax[1], lwd=3, col='magenta', lty='dashed')


####' 
##   PATH: S-JumpStatistics/Data/window_compare/
##   Reconstructing positions x(t) using jump time/distance
##    
####' 
t <- seq(0, 1e6, length=1e6+1)
x <- data.matrix(fread(lf("0-x.dat")))

tau <- data.matrix(fread(lf("w10000_tau.dat")))
xi  <- data.matrix(fread(lf("w10000_jd.dat")))

s <- seq(100000)
matplot(t, x, type='l')
mm <- getMM(x=t, y=x, w=10000)
ShowJT(mm$df, eps=0)

#' Reconstructed x'(t')
t0 <- mm$df$xl[1]
xi0 <- mm$df$yl[1]

tp <- c(t0, t0 + cumsum(mm$tau))
xi <- mm$df$yr - mm$df$yl
m <-  xi / mm$tau

X[1] <- xi0
X[2] <- X[1] + m[1] * mm$tau[1]
X[3] <- X[2] + m[2] * mm$tau[2]
X[4] <- X[3] + m[3] * mm$tau[3]

points(tp, X, pch=19, col='red', cex=2)
lines(tp, X)

####' 
##   PATH: S-JumpStatistics/Data/window_compare/
##   Role of the window width 
##    
####' 


T <- 1e6
fact <- 1
t <- seq(0, T, length=(T+1)/fact)
x <- ReadDataList(lf("p0-a0-d2"))
#pdf(paste0(img.path, ""), width=12, height=16)
spar(cx=2.5, col='white')
i <- 15
w <- 2500
s <- seq(nrow(x))
s <- seq(1, 280000)
s <- seq(1, 10000)

matplot(t[s], x[s,i], type='l', lty='solid', ylab=tex("x(t_k)"), xlab=tex("t_k"), lwd=1)
mm <- getMM(x = t, y = x[,i], w = w[j])
ShowJT(df = mm$df, text=FALSE, eps=0, pt.cx=0.65)



####' 
##   PATH: S-JumpStatistics/Data/window_compare/
##   Zooming in on time series
##    
####' 

T <- 1e6
fact <- 1
t <- seq(0, T, length=(T+1)/fact)
x <- ReadDataList(lf("p0-a0-d2"))
#pdf(paste0(img.path, ""), width=12, height=16)
spar(cx=2.5, col='white')
i <- 15
w <- 2500

s <- seq(1, 10000)
s <- seq(1, 100000)
s <- seq(nrow(x))
matplot(t[s], x[s,i], type='l', lty='solid', ylab=tex("x(t_k)"), xlab=tex("t_k"), lwd=1)

mm <- getMM(x = t[s], y = x[s,i], w = 10)
mm <- getMM(x = t[s], y = x[s,i], w = 500)
mm <- getMM(x = t[s], y = x[s,i], w = 5000)
ShowJT(df = mm$df, text=FALSE, eps=0, pt.cx=0.65)

m.star <- list()
m.star[[3]] <- mm$tau #abs(mm$df$yr - mm$df$yl) / mm$tau

plot(sort(m.star[[1]]),  1-seq(length(m.star[[1]])) / length(m.star[[1]]), type='s', log='x', xlim=c(1, max(m.star[[3]])))
lines(sort(m.star[[2]]), 1-seq(length(m.star[[2]])) / length(m.star[[2]]), type='s', col='red')
lines(sort(m.star[[3]]), 1-seq(length(m.star[[3]])) / length(m.star[[3]]), type='s', col='blue')

# Shows different jump events for different windows
T <- 1e6
fact <- 1
t <- seq(0, T, length=(T+1)/fact)
x <- ReadDataList(lf("-x.dat"))
pdf(paste0(img.path, "window_compare.pdf"), width=12, height=16)
spar(cx=2.5, col='white', mfrow=c(3,1))
par(mar = c(8, 8, 1.3, 1.3), oma = c(2, 2, 2, 2), mgp=c(5,2,0))
i <- 15
w <- rev(c(1000, 2500, 5000))
for (j in seq(length(w))) {
    s <- seq(1, ifelse(j == 5, yes=10000, no=nrow(x)), by=100)
    matplot(t[s], x[s,i], type='l', lty='solid', ylab=tex("x(t_k)"), xlab=tex("t_k"), lwd=1)
    dy <- max(x[s,i]) - min(x[s,i])
    text(x = t[1]+10, y = max(x[s,i]) - 0.1*dy , labels=j, cex=4)
    mm <- getMM(x = t, y = x[,i], w = w[j])
    ShowJT(df = mm$df, text=FALSE, eps=0, pt.cx=0.65)
}; dev.off()


####' 
##    Understanding jumps which are much smaller than window size
##    Conclusion: It seems regardless of how times are defined, there will
##    will always be cases in which the jump time << window size. This is
##    true in particular for when the chosen window size is too large.
####' 

T <- 1e6
fact <- 1
t <- seq(0, T, length=(T+1)/fact)
x <- data.matrix(fread(lf("-x.dat")))
x <- ReadDataList(lf("-x.dat"))#, written.as.row=FALSE)
tau <- data.matrix(fread(lf("--su_tau.dat")))
i <- 1

pdf(paste0(img.path, "-ex.pdf"), width=22, height=18)
s <- seq(20000, 25000)
s <- seq(nrow(x))
spar(cx=2, col='white')
matplot(t[s], x[s,i], type='l', lty='solid', pch=19, ylab=tex("x(t_k)"), xlab=tex("t_k"))

w <- c(50, 500, 5000)
m <- list()
for (k in seq(length(w))) {
    mm <- getMM(x = t, y = x[,i], w = w[k])

    m[[k]] <- ((mm$df$yr - mm$df$yl) / mm$tau)
}
m  <- list.cbind.na(m)
ShowJT(df = mm$df, text=FALSE, eps=0)

idx.max <- am(w,  x[,i])
idx.min <- am(w, -x[,i])
points(t[idx.min], x[idx.min,i], pch=19, cex=2.5,col='green')
points(t[idx.max], x[idx.max,i], pch=19, cex=2.5,col='red')


####' 
##    PATH: S-JumpStatistics/Data/minmin_possible/
##    Showing that a min-min or max-max detection is possible
##    but that it is not always irrelevant. i.e. jumps may occur between 
##    two detected local min./max.
####' 
am <- function(w, y) {
    y.max <- rollmax(zoo(y), k=2*w, align='center')
    delta <- y.max - y[-c(1:(w-1), length(y)+1-1:w)]
    which(delta <= 0) + w - 1
}

T <- 1e7
t <- seq(0, T, length=(T+1)/10)
x <- ReadDataList(lf("-x.dat"))

pdf(paste0(img.path, "minmin-ex.pdf"), height=18, width=12)
spar(cx=2.5, col='white', mfrow=c(3,1))
par(mar = c(8, 8, 1.3, 1.3), oma = c(2, 2, 2, 2), mgp=c(5,2,0))
w <- 5000
s <- seq(1, length(t)-1, by=100)
k <- 1
for (i in c(9, 1, 4)) {
    matplot(t[s], x[s,i], type='l', lty='solid', pch=19, ylab=tex("x(t_k)"), xlab=tex("t_k"))

    dy <- max(x[s,i]) - min(x[s,i])
    text(x = t[s[length(s)]]-10, y = max(x[s,i]) - 0.1*dy , labels=k, cex=4)
    k <- k + 1

    idx.max <- am(w,  x[,i])
    idx.min <- am(w, -x[,i])
    points(t[idx.min], x[idx.min,i], pch=19, cex=2,col='green')
    points(t[idx.max], x[idx.max,i], pch=19, cex=2,col='red')

    legend('topleft', legend=c("Local minimum", "Local maximum"), col=c('green', 'red'), pch=19, cex=2)
}
dev.off()



####' 
##    Visualising the window for whatever data is considered
####' 

K <- nrow(x)
am <- function(w, y) {
    y.max <- rollmax(zoo(y), k=2*w, align='center')
    delta <- y.max - y[-c(1:(w-1), length(y)+1-1:w)]
    which(delta <= 0) + w - 1
}
y.max.compare <- y.min.compare <- c()
idx.max <- am(w, x)
idx.min <- am(w, -x)
for (i in seq(1, K - 2*w + 1)) {

    spar(cx=2, col='white')
    plot(t[s], x[s], type='l', cex=2, xlab=tex("x_k, k = 1, 2, ..., K"))#, ylim=c(-100,100))
    ShowJT(df = mm$df, text=FALSE, eps=1000)

    #rect(xleft=t[i], xright=t[i + 2*w - 1], ytop=1e6, ybot=-1e6, col=rgb(0,0,1,alpha=0.2))
    abline(v=t[i], lwd=2, col='blue')
    abline(v=t[i + 2*w - 1], lwd=2, col='blue')

    j <- seq(i, i+2*w - 1)
    wm <- which.max(-x[j]) + i - 1
    wM <- which.max(x[j]) + i - 1
    y.min.compare <- c(y.min.compare, max(-x[j]))
    y.max.compare <- c(y.max.compare, max(x[j]))

    #text(x=t[wm]+2.6, y=x[wm]-0.28, labels=TeX(paste0("$m_i + x(t_{i+w-1}) = m_{", i, "} + x(t_{", i+w-1, "}) = $",  format(y.min.compare[i] + x[i+w-1],digits=3))), cex=2, col='blue')
    #text(x=t[wM]-2.6, y=x[wM]-0.28, labels=TeX(paste0("$M_i - x(t_{i+w-1}) = M_{", i, "} - x(t_{", i+w-1, "}) = $",  format(y.max.compare[i] - x[i+w-1],digits=3))), cex=2, col='blue')

    xoff <- w
    yoff <- -1
    text(x=t[i]+xoff, y=max(x[s])+yoff, labels=TeX(paste0("$m_i + y(x_{i+w-1}) = m_{", i, "} + y(x_{", i+w-1, "}) = $", format(y.min.compare[i] + x[i+w-1], digits=3))), cex=2, col='black')
    text(x=t[i]+xoff, y=max(x[s])-1+yoff, labels=TeX(paste0("$M_i - y(x_{i+w-1}) = M_{", i, "} 0 y(x_{", i+w-1, "}) = $", format(y.max.compare[i] - x[i+w-1], digits=3))), cex=2, col='black')

    points(t[wm], x[wm], pch=19, col='green')
    points(t[wM], x[wM], pch=19, col='red')


    points(t[idx.min], x[idx.min], pch=19, col='green')
    points(t[idx.max], x[idx.max], pch=19, col='red')

    Sys.sleep(0.1)

}




####' 
##    Coarse graining of some IH data 
####' 

setwd("../Studies/S-JumpStatistics/Data/fp_gt_100_final/")

f <- gtools::mixedsort(lf("_tau"))
tau <- list()
for (i in seq(length(f)))
    tau[[i]] <- data.matrix(fread(f[i]))
tau[[i+1]] <- data.matrix(fread("p0-a1-d20---t3.dat"))
tau <- list.cbind.na(tau)

L <- getTSD_hist(dx=tau[,], n = 200, prune=FALSE, domain=exp.domain)

{
    img.path = "~/Dropbox/PHD/Papers/JumpStatistics/Images/"
    pdf(paste0(img.path, "ih_coarseg.pdf"), width=22, height=18)
    spar(cx=3, col='white')
    par(mar = c(8, 8, 1.3, 1.3), oma = c(2, 2, 2, 2), mgp=c(5,2,0))
    matplot(L$x[,-5], L$y[,-5], type='s', lty='solid', log='xy', lwd=4, xlab=tex("\\tau"), ylab=tex("\\Psi(\\tau)"), col=rainbow(10)[-5])
    lines(L$x[1:190,5], L$y[1:190,5], type='s', lty='solid', lwd=4, col=rainbow(10)[-5][5])
    
    x <- seq(1, 1e4)
    lines(x, exp(13.2)*x^(-3), lty='dashed', lwd=3)
    lines(x, exp(9.6)*x^(-3), lty='dashed', lwd=3)
    legend('topright', legend=c(paste0("w = ", c(10, 50, 200, 400)), tex("\\tau = s"), tex("\\tau^{-3}")), 
    			    col=c(rainbow(10)[-5][1:5], 'black'), 
    			    lty=c(rep('solid' ,5), 'dashed'), 
    			    cex=4, lwd=4)
    dev.off()
} 



####' 
##    R^2 as an acceptance criterion
####' 

setwd("../Studies/S-JumpStatistics/Data/R2_acceptance/")
{ 
    T <- 1e4
    t <- seq(0, T, length=T+1)
    x <- data.matrix(fread("p0-a1-d2--0-x2.dat"))
    s <- seq(length(t))
    s <- seq(1000)
    
    spar(cx=2, col='white')
    matplot(t[s], x[s], type='l', lty='solid', lwd=1, col='black', xlab=tex("t_k"), ylab=tex("x(t_k)"), 
    	pch=19, cex=0.75)
    
    mm <- getMM(x = t, y = x, w = 200)
    ShowJT(df = mm$df, wid = 1, text=FALSE, eps=1000)
    
    r2 <- c()
    for (i in seq(nrow(mm$df))) {
        X <- seq(mm$df$xl[i]-1, mm$df$xr[i]-1)
        Y <- x[X]
        #lines(X, Y, col='red')
        L <- lm(Y~X)
        r2 <- c(r2, summary(L)$r.square)
        cof <- coef(lm(Y ~ X))
        lines(X, cof[2]*X + cof[1], col=rainbow(10)[-5][i], lwd=5)
    }
    legend('topright', legend=TeX(paste0("$R^2 = $", format(r2, digits=3))), col=rainbow(10)[-5][1:nrow(mm$df)], lwd=3, cex=2)
}


####' 
##    Maximum distance criterion
####' 

setwd("../Studies/S-JumpStatistics/Data/maxdist_acceptance/")

perp_dist <- function(P1, P2, x0, y0) {
    num <- abs( (P2[1] - P1[1]) * (P1[2] - y0) - (P1[1] - x0) * (P2[2] - P1[2]))
    den <- sqrt((P1[1] - P2[1])^2 + (P1[2] - P2[2])^2)
    return(num / den)
}

{ 
    T <- 1e4
    t <- seq(0, T, length=T+1)
    x <- data.matrix(fread("p0-a1-d2--0-x2.dat"))
    s <- seq(length(t))
    s <- seq(50, 350)
    
    spar(cx=2, col='white')
    matplot(t[s], x[s], type='l', lty='solid', lwd=1, col='black', xlab=tex("t_k"), ylab=tex("x(t_k)"), 
    	pch=19, cex=0.75)
    
    mm <- getMM(x = t, y = x, w = 20)
    ShowJT(df = mm$df[4:5,], wid = 1, text=FALSE, eps=1000)
    
    id <- seq(4,5)
    for (j in seq(length(id))) {
	X <- seq(mm$df$xl[id[j]], mm$df$xr[id[j]])
    	Y <- x[X+1]

    	df <- data.matrix(mm$df)
    	lines(df[4,1:2], df[4,3:4],  col='red',   lwd=3.5)
    	lines(df[5,1:2], df[5,3:4],  col='magenta', lwd=3.5)
    	points(df[4:5,1:2], df[4:5,3:4], col='green', pch=19)

	P1 = as.vector(c(df[id[j],1], df[id[j],3]))
	P2 = as.vector(c(df[id[j],2], df[id[j],4]))

	#points(P1, P2, col='green', pch=19)
	tt <- seq(0,1,length=length(X))
	xx <- P1[1] + tt*(P2[1] - P1[1])
	yy <- P1[2] + tt*(P2[2] - P1[2])
	#points(xx, yy)

	pd <- c()
    	for (k in seq(1, length(X), by=3)) {
    	    x0 <- X[k]
	    y0 <- Y[k]

    	    pd <- c(pd, perp_dist(P1, P2, x0, y0))
	    lines(c(xx[k], x0), c(yy[k], y0), col='blue')

	    #print(paste0(pd, " ",  sqrt((x0 - xx[k])^2 + (y0 - yy[k])^2)))
	    #lines(c(x0, x0 + pd*
    	}
	print(max(pd))
    }
    legend('topright', legend=paste0("Max. deviation = ", c(1.2, 0.70)), col=c('red', 'magenta'), lwd=3, cex=2)
}


###
##
###







