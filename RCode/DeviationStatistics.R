
####' 
##   PATH: S-JumpStatistics/Data/
##   We consider the mean deviation from the line connecting two turning points
##   s.t. that line segment spans a particular jump time.
##
##   The idea is that small windows should produce a mean deviation which is
##   small when the jump time is large. This implies that the tails of the
##   jump time distribution obtained for small windows, should accurately represent 
##   large jumps, but also only large jumps which are "fast". That is, the slope
##   of the line connecting two turning points is large. This is problematic for 
##   an asymptotic picture, which does not distinguish between the speed of jumps. 
##
##   Thu 27 Jan 2022 16:02:27. UPDATE. 
##   Small windows leave us with 
####' 


img.path = "~/Dropbox/PHD/Papers/JumpStatistics/Images/"
setwd("../Studies/S-JumpStatistics/Data/")

cnt     <- ReadDataList(lff("a1", "_jump.dat"), idx=1)
t0      <- ReadDataList(lff("", "_jump.dat"), idx=2)
x0      <- ReadDataList(lff("", "_jump.dat"), idx=3)
jt      <- ReadDataList(lff("a", "w500_jump.dat"), idx=4)
jd      <- ReadDataList(lff("a", "w500_jump.dat"), idx=5)
mepd	<- ReadDataList(lff("a", "_jump.dat"), idx=6)
x       <- ReadDataList(lff("a", "xt.dat"), idx=4)
vel.x   <- ReadDataList(lff("a", "vel.dat"), idx=1)
vel.y   <- ReadDataList(lff("a", "vel.dat"), idx=2)
th	<- atan(vel.y / vel.x)


id <- apply(cnt, 2, function(x) which(x == 0))
num.ics <- apply(cnt, 2, function(x) length(which(x == 0)))
ranges <- lapply(id, function(x) sapply(1:(length(x)-1), function(y) seq(x[y], x[y+1]-1)))
nJumpPerIC <- lapply(ranges, function(x) lapply(x, function(y) length(y))) # Number of jumps per IC for each window

# Replicates initial conditions for each jump time
thetas <- list()
for (i in seq(length(nJumpPerIC))) {
    non.na <- length(which(!is.na(th[,i]))) - 1
    tmp <- list()
    for (j in seq(non.na)) 
	tmp[[j]] <- rep(th[j,i], nJumpPerIC[[i]][[j]])
    thetas[[i]] <- tmp
}

eps <- 3 * pi/180
j <- 1
val <- jd[ranges[[j]][[1]],1] 
val <- jt[ranges[[j]][[1]],1] 
matplot(thetas[[j]][[1]], val, pch=19, cex=.2, log='', xlim=c(-eps, eps), ylim=c(-300,300)) #ylim=c(50, 10000))
matplot(thetas[[j]][[1]], val, pch=19, cex=.2, log='y', xlim=c(-eps, eps), ylim=c(50, 10000))
for (i in seq(2,length(ranges[[j]]))) {
    val <- jd[ranges[[j]][[i]],1] 
    #val <- jt[ranges[[j]][[i]],1] 
    matplot(thetas[[j]][[i]], val, pch=19, cex=.2, add=TRUE)
}
abline(h=0, v=0)
x <- seq(0, eps, length=10000) - 0.035
a <- c(8, 2, 2/4, 2/8)
sapply(a, function(b) lines(x + 0.035, 1/(2*b*x), col='blue'))


spar(cx=2, col='white')
matplot((mepd*jt)/sqrt(jt + jd), jt, pch=19, col=rbow(), cex=0.2, log='xy')
matplot(jd/jt, mepd[,], pch=19, cex=0.2, log='')
matplot(atan(vel.y/vel.x), jd/jt, pch=19, cex=0.2, log='')

s <- seq(2, nrow(jt)-1, by=10)
n <- 1
F <- n:(nrow(jt)-n)
La <- (n+1):nrow(jt)
#pdf(paste0(img.path, "correl_jtjd.pdf"), height=9, width=18)
par(fig=c(0.07, 0.315, 0.095, 0.52), new=TRUE, cex=0.75)
plot(jt[F,1][s], jt[La,1][s], cex=0.2, pch=19, log='xy', xlab=tex("\\tau_i"), ylab=tex("\\tau_{i+1}"))
par(fig=c(0.57, 0.81, 0.095, 0.52), new=TRUE, cex=0.75) # jd
plot(abs(jd[F,1][s]), abs(jd[La,1][s]), cex=0.2, pch=19, log='', xlab=tex("\\xi_i"), ylab=tex("\\xi_{i+1}"), xlim=c(1e4,1e5), ylim=c(1e4,1e5))
dev.off()


###
#   Distribution of jump times 
###
legend('topright', cex=3, col=rbow(), lty='solid', lwd=3, legend=c(TeX(paste0("$\\alpha = ", c("\\pi/3", "\\pi/2", "2\\pi/3")))))
# Probabilities
#pdf(paste0(img.path, "cg_jt_fits.pdf"), width=14, height=10)
id <- apply(abs(jd), 2, function(x) which(x < 1e6))
L <- getTSD_hist(dx = (jd[id[[2]],]/jt[id[[2]],]), n= 100)#, domain=exp.domain)

spar(cx=3, col='white', mfrow=c(1,2))
L <- getTSD_ecdf(dx = jt, n= 1000, domain=exp.domain)
mplot(L$x[,1:3], 1-L$y[,1:3], type='s', spa=FALSE, lwd=5, log='xy', ylim=c(0.0001, 1), xlim=c(2e4,5e6), xlab=tex("\\tau"), ylab="Prob.")
L <- getTSD_ecdf(dx = abs(jd), n= 1000, domain=exp.domain)
mplot(L$x[,1:3], 1-L$y[,1:3], type='s', spa=FALSE, lwd=5, log='xy', ylim=c(0.0001, 1), xlim=c(1e4,1e6), xlab=tex("\\xi"), ylab="Prob.")
x <- seq(1e3, 1e7, by=100)
lines(x, exp(16.9)*x^(-1.8), lwd=3, lty='dashed')

LF <- LinFit(time = L$x[,3], mat = L$y[,3], range = 150:199, logs=TRUE, show.range=FALSE, plot=TRUE)
LF <- LinFit(time = L$x[,3], mat = L$y[,3], range = 100:150, logs=TRUE, show.range=FALSE, plot=TRUE)
LF <- LinFit(time = L$x[,3], mat = L$y[,3], range = 50:100, logs=TRUE, show.range=FALSE, plot=TRUE)
LF <- LinFit(time = L$x[,1:3], mat = 1-L$y[,1:3], range = 1:400, logs=TRUE, show.range=TRUE, plot=TRUE)
LF <- LinFit(time = L$x[,1], mat = 1-L$y[,1], range = 100:400, logs=TRUE, show.range=TRUE, plot=TRUE)
LF <- LinFit(time = L$x[,2], mat = 1-L$y[,2], range = 500:650, logs=TRUE, show.range=TRUE, plot=TRUE)
LF <- LinFit(time = L$x[,3], mat = 1-L$y[,3], range = 500:650, logs=TRUE, show.range=TRUE, plot=TRUE)
LF <- LinFit(time = exp(L$x[s,4]), mat = log(1-L$y[s,4]), range = 800:899, show.range=TRUE, plot=TRUE)


#for (j in seq(ncol(L$y))) {
for (j in seq(1,3)) {
#for (j in seq(1,3)) {
    df <- data.frame(x=L$x[,j], y=1-L$y[,j])
    f <- function(x, a, b) { exp(-(x/a)^b) }
    #NLF <- nlsLM(y ~ 1-f(x,a,b), data=df, start=list(a=20, b=-0.9),control = nls.lm.control(maxiter = 200))#,subset=s)
    NLF <- nlsLM(y ~ 1-f(x, a, b), data=df, start=list(a=2000, b=0.4),control = nls.lm.control(maxiter = 200))#,subset=s)
    X <- L$x[,j]
    cf <- coef(NLF)
    print(cf)
    lines(X, (1-f(X, a=cf[1], b=cf[2])), lwd=3, lty='dashed')#, c=cf[3])))
    #lines(X, (f(X, a=cf[1], b=cf[2])), lwd=2, lty='dashed')#, c=cf[3])))
}


LM <- lm(formula = log(1-L$y[,4]) ~ L$x[,4], subset=800:899)
co <- LM$coefficients
lines(L$x[,4], L$x[,4] * co[2] + co[1], lty='dashed', lwd=1)
lines(x, exp(16.5)*x^(-1.442))
lines(x, exp(16.9)*x^(-1.8), lwd=3, lty='dashed')
x <- seq(-1, 1, length=100)
lines(L$x[,2], 3*(f(L$x[,2], a=2, b=4.2) / (L$x[,2]^(-.6))))

x <- seq(0,1, by=0.01)
X <- inv.strexp(x,a=90000, b=-1.5)
plot(x, X)
stres <- strexp.sample(x, a=91000, b=1.44, n=10000)
L <- getTSD_ecdf(dx = stres, dx

j <- 3
df <- data.frame(x=L$x[,j], y=L$y[,j])
f <- function(x, a, b) { ((a^b)/(gamma(b))) * (exp(-a/x))/(x^(1+b)) }



legend('bottomleft', legend=paste0("w = ", 10^seq(1,4)), col=rbow(), lty='solid', cex=2, lwd=2)
coefs <- list()
for (j in seq(1, ncol(L$y))) {
    df <- data.frame(x=L$x[,j], y=1-L$y[,j])
    s <- c(800,920)
    LM <- lm(formula = log(y) ~ x, subset=s[1]:s[2], data = df)
    coefs[[j]] <- LM$coefficients
    co <- LM$coefficients
    lines(df$x, df$x * co[2] + co[1], lty='dashed', lwd=1)
    #abline(v=c(L$x[s[1],j], L$x[s[2],j]))
}
coefs <- tpose(list.cbind.na(coefs))
diff(coefs[,1])
mean(coefs[,2])
sqrt(var(coefs[,2]))

par(fig=c(0.62, 0.97, 0.475, 0.95), new=TRUE, cex=0.75, mgp=c(3,1.5,0))
s <- seq(2e4, 2e6, by=100000)
matplot(s, s, pch=19, cex=0.2, log='xy', col='white', xlab=tex("\\tau_i"), ylab=tex("\\tau_{i+1}"), xaxt='n', yaxt='n')
axis(side=1, at=c(2e4, 2e5, 2e6))
axis(side=2, at=c(2e4, 2e5, 2e6))
id <- which(lapply(ranges[[2]], function(x) length(x)) > 4)
for (i in seq(1, length(id), by=10)) {
    rg <- ranges[[2]][[id[i]]]
    lrg <- length(ranges[[2]][[id[i]]])
    f <- rg[1:(lrg-1)]
    l <- rg[2:(lrg)]
    points(jt[f,2], jt[l,2], pch=19, cex=0.2, col='black', )
    
}
dev.off()

LF <- LinFit(time = L$x, mat = 1-L$y, range = 110:300, logs=TRUE, show.range=TRUE, plot=TRUE)
LF <- LinFit(time = L$x, mat = 1-L$y, range = 400:700, logs=TRUE, show.range=TRUE, plot=TRUE)
LF <- LinFit(time = L$x, mat = 1-L$y, range = 850:950, logs=TRUE, show.range=TRUE, plot=TRUE)


lines(df$x, df$x * co[2] + co[1], lty='dashed', lwd=1)

tau <- lambda.t
omega <- list()
for (j in seq(ncol(tau))) {
    nas <- which(is.na(tau[,j]))
    print(length(nas))
    temp.tau <- tau[,j]
    if (length(nas) > 0)
	temp.tau <- tau[-nas,j]
    om <- c()
    for (i in seq(1,length(temp.tau)-1)) {
        om[i] <- temp.tau[i] / (temp.tau[i] + temp.tau[i+1])
    }
    omega[[j]] <- om
}
omega <- list.cbind.na(omega)
Lw <- getTSD_hist(dx=data.matrix(omega), n = 100, prune=FALSE)#, domain=exp.domain)

par(fig=c(0.055, 0.40, 0.075, 0.6), new=TRUE, cex=0.75, mgp=c(3,1.5,0))
par(fig=c(0.67, 0.97, 0.475, 0.95), new=TRUE, cex=0.75, mgp=c(3,1.5,0))
matplot(Lw$x, Lw$y, type='l', pch=19, xlim=c(0,1), col=rbow(), lty='solid', log='', lwd=2,
	ylab=tex("P(\\omega_2)"), xlab=tex("\\omega_2"), xaxt="n")
axis(side=1, at=c(0, 0.25, 0.5, 0.75, 1))


L <- getTSD_ecdf(dx = abs(jd), n= 1000, domain=exp.domain)
L <- getTSD_hist(dx = (jd)/jt, n= 200)#, domain=exp.domain)
matplot(L$x, L$y, type='s', lwd=2, lty='solid', log='', col=rbow(), #xlim=c(-1e6, 1e6), 
	xlab=tex("|\\xi |"),
	ylab="Density")
x <- seq(0, 5e5)
lines(x, 100000/((x^3)*log(x)))
LF <- LinFit(time = L$x, mat = 1-L$y, range = 10:600, logs=TRUE, show.range=TRUE, plot=TRUE)

lines(L$x[,1], exp(17)*L$x[,1]^(-2.85), lty='dashed', col='black', lwd=2)
legend('topright', legend=c(TeX(paste0("$\\alpha = ", c("\\pi/3", 
							"\\pi/2", 
							"2\\pi/3", 
						        "(\\varphi - 5/3)\\pi/3", 
						        "(\\varphi - 5/3)\\pi/2", 
						        "(\\varphi - 5/3)2\\pi/3"), "$")), TeX("$\\tau^{-2.85}$")),
		   col=c(rbow()[1:6], 'black'), lty=c(rep('solid',6), 'dashed'), lwd=2, cex=2)

# Densities
L <- getTSD_hist(dx = (jd/jt)[,-4], n= 200)#, domain=exp.domain)
L <- getTSD_hist(dx = (x), n= 250, )
#pdf(paste0(img.path, "speed_w_dep.pdf"), width=15, height=10)
spar(cx=2, col='white')
matplot(L$x, L$y, type='s', lwd=2, lty='solid', log='', col=rbow(), #ylim=c(0, 4),
	xlab=tex("\\xi / \\tau"),
	ylab="Density")
matplot(L$x, L$y, type='s', lwd=2, lty='solid', log='y', col=rbow(), add=TRUE)
legend('topleft', legend=paste0("w = ", c(10, 100)), col=rbow(), lty='solid', cex=2, lwd=2)
legend('topleft', legend=paste0("w = ", 10*c(1000, 1e4, 1e5)), col=rbow(), lty='solid', cex=2, lwd=2)
dev.off()
LF <- LinFit(time = L$x[,4], mat = L$y[,4], range = 140:195, logs=TRUE, show.range=TRUE, plot=TRUE)
matplot(jd / jt, jt, pch=19, cex=0.2, col=rbow(), log='y')


spar(cx=2, col='white')
W <- c(10, 100, 1000, 1e4, 1e5)
L <- getTSD_ecdf(dx = mepd, n= 1000, domain=exp.domain)
matplot(L$x, 1-L$y, type='s', lwd=2, lty='solid', log='x', col=rbow(), xlim=c(.01,1e3),
	xlab=tex("\\tau"),
	ylab="Prob.")
hid <- apply(1-L$y, 2, function(x) min(which(x < 0.5)))
abline(v=L$x[hid[2],2] - L$x)
sapply(1:(length(hid)-1), function(x) log(L$x[hid[x+1],x+1]) - log(L$x[hid[x],x]))
legend('topleft', legend=paste0("w = ", W), col=rbow(), lty='solid', cex=2, lwd=2)


# Visualising flow model
{
    K <- 1
    CC           <- ReadData(lff("a1", "CC"))[,1,1]
    vx           <- ReadData(lff("a1", "-Vx.dat"))
    vy           <- ReadData(lff("a1", "-Vy.dat"))
    y            <- ReadData(lff("a1", "-y.dat"))
    t            <- ReadData(lff("a1", "-t.dat"))
    x            <- ReadData(lff("a1", "-x.dat"))
    tt           <- ReadData(lff("a1", "-tt.dat"))
    xx           <- ReadData(lff("a1", "-xx.dat"))
    lambda.t     <- ReadDataList(lff("a7", "_jump.dat"), idx=1)
    lambda.x     <- ReadDataList(lff("a7", "_jump.dat"), idx=2)

    s <- seq(1,1000000)
    mplot((t)[s], x[s], type='l', pch=19)
    mplot(cumsum(tt)[s], xx[s], type='l', pch=19, add=TRUE, col='blue')
    mplot(cumsum(tt)[s], xx[s], type='l', pch=19, col='blue')
    mplot(cumsum(lambda.t), cumsum(lambda.x), lwd=3)

    points(cumsum(lambda.t), cumsum(lambda.x), pch=19, col='blue', cex=1)
}


# DUMB CUNTERY ENSUES
cc <- c(1,cumsum(c(CC)))
lt <- lx <- c()
lt[1] <- 0
lx[1] <- x[1]
for (i in seq(2,length(cc)-1)) {
    if (i == 1)
        S <- seq(cc[i], cc[i+1])
    else
	S <- seq(cc[i]-1, cc[i+1])
    lt[i] <- sum(t[S])
    lx[i] <- x[cc[i+1]] - x[cc[i]]
}
lt <- c(0, lt)
lx <- c(x[1], lx)#, x[cc[length(cc)]+1] - x[cc[length(cc)-1]+1])
points(c(0,cumsum(lt)), c(x[1], cumsum(lx)), pch=19, col='red')
lines(cumsum(lt), cumsum(lx), pch=19, col='red', lwd=2)

cc <- c(1,cumsum(c(CC)))
lt[1] <- 0; lx[1] <- x[1]
lt <- lx <- c()
for (i in seq(1,length(cc)-1)) {
    if (i == 1)
	S <- seq(cc[i], cc[i+1])
    else 
	S <- seq(cc[i]+1, cc[i+1])
    lt[i] <- sum(t[-1][S])
    if (i == 1)
	lx[i] <- x[1]
    else if (i == 2)
	lx[i] <- x[cc[i]+1] - x[cc[i-1]]
    else
	lx[i] <- x[cc[i]+1] - x[cc[i-1]+1]
}
lt <- c(0, lt)
lx <- c(lx, x[cc[length(cc)]+1] - x[cc[length(cc)-1]+1])
points(cumsum(lt), cumsum(lx), pch=19, col='red')
lines(cumsum(lt), cumsum(lx), pch=19, col='red')

### Time averages

nsteps <- nrow(lambda.x)
X.t <- function(i) sum(lambda.x[1:i,1])
X.t <- cumsum(lambda.x)

Delta <- seq(1, 30000, by=100)
msd <- c()
for (j in seq(length(Delta))) {
    temp <- 0
    for (i in seq(1, nsteps - Delta[j])) {
        #temp <- temp +  (X.t(i+Delta[j]) - X.t(i))^2
        temp <- temp +  (X.t[i+Delta[j]] - X.t[i])^2
    }
    msd[j] <- (1/(nsteps - Delta[j])) * temp
}

mplot(Delta, msd, type='o')

###
##  Flow definition of a step
###

P            <- data.matrix(fread("a7-d8-e1--P0.dat"))

L <- getTSD_hist(dx = atan(P[,4] / P[,3]), n=75)
mplot(L$x, L$y, type='s')
plot((P[,1]), (P[,2]), cex=0.2, pch=19, asp=1)


th           <- ReadDataList(lff("", "theta.dat"))
vel          <- ReadDataList(lff("p1-a1", "_jump.dat"), idx=4)
tt           <- ReadDataList(lff("p", "TL"), idx=1)
lambda.t     <- ReadDataList(lff("p", "jump.dat"), idx=1)
lambda.x     <- ReadDataList(lff("p", "jump.dat"), idx=2)
CC           <- ReadDataList(lff("p", "_jump.dat"), idx=3)
apply(lambda.t, 2, function(x) length(which(!is.na(x))))

id <- which(lambda.t == 0)
L <- getTSD_hist(dx = th, n=75)
L <- getTSD_hist(dx = atan(lambda.x[,]/lambda.t[,]), n=200)#, domain=exp.domain)
L <- getTSD_hist(dx = lambda.x[,]/lambda.t[,], n=200)#, domain=exp.domain)
L <- getTSD_ecdf(dx = abs(lambda.x), n=1000, domain=exp.domain)
L <- getTSD_ecdf(dx = tt, n=1000, domain=exp.domain)

#pdf(paste0(img.path, "th0_msd.pdf"), height=9, width=18)
spar(cx=2.25, col='white', mfrow=c(1,1))
K <- c(1)
#matplot(L$x[,], 1-L$y[,], log='',
matplot(L$x[,], L$y[,], log='', xlim=c(-.1, .1),
	type='s', col=rbow(),
	lty='solid', lwd=3, 
	ylab="Density", 
	xlab=tex("\\theta_0 = \\tan^{-1} \\[ \\xi / \\tau \\], \\tau_0 > 1000"))
legend('topright', lwd=4, cex=1.5, col=rbow(), legend=TeX(paste0("$n = ", c(0,1,2,3,4,6,8)[K], "$")))
LF <- LinFit(mat = 1-L$y[,1], time=L$x[,1], range = 10:200, logs = TRUE, plot = TRUE, show.range=TRUE)

legend('topright', lwd=4, cex=1.5, col=rbow(), legend=TeX(paste0("$\\tau_0 = ", c(50, 100), "$")))
legend('topright', lwd=4, cex=1.5, col=rbow(), legend=TeX(paste0("$\\xi_0 = ", c(50, 100), "$")))


par(fig=c(0.06, 0.45, 0.125, 0.705), new=TRUE, cex=0.75, mgp=c(3,1.5,0))
matplot(seq(1,4e5,by=1000), seq(1,4e5,by=1000), lty='solid', type='l', col='black', xlab=tex("\\tau_i"), ylab=tex("\\tau_{i+1}"))
K <- 7
for (i in seq(length(K))) {
    len <- length(which(!is.na(lambda.t[,K[i]])))
    s <- seq(1, len, by=10)
    l <- (2:len)
    f <- (1:(len-1))
    matplot(lambda.t[f[s],K[i]], 
	    lambda.t[l[s],K[i]], pch=19, cex=0.2, col=rbow()[K[i]], add=TRUE)
}
dev.off()

rect(xleft=-0.28, xright=0.28, ybottom=0, ytop=10, col=rgb(1,0,0,alpha=0.1))
rect(xleft=pi/2, xright=pi/2-0.28, ybottom=0, ytop=10, col=rgb(1,0,0,alpha=0.1))
rect(xleft=-pi/2, xright=-pi/2+0.28, ybottom=0, ytop=10, col=rgb(1,0,0,alpha=0.1))
matplot(L$x[,-(1:3)], L$y[,-(1:3)], type='s', col=rbow(), lty='dotted', lwd=4,add=TRUE)

legend('topright', lwd=4, cex=1.5, col=rbow(), legend=TeX(paste0("$n = ", c(0,2^seq(1,6)), "$")))
legend('topright', lwd=4, cex=2, col=rbow()[1:3], legend=TeX(paste0("$\\xi > ", c(5, 50, 250), "$")))
legend('topright', lwd=4, cex=2, col=rep(rbow()[1:3],2), lty=c(rep('solid', 3), rep('dotted', 3)), legend=c(TeX(paste0("$\\alpha = ", c("\\pi/3", "\\pi/2", "2\\pi/3", 
						        "(\\varphi - 5/3)\\pi/3", 
						        "(\\varphi - 5/3)\\pi/2", 
						        "(\\varphi - 5/3)2\\pi/3$")))))
dev.off()

#pdf(paste0(img.path, "ic_jump_vel.pdf"), height=12, width=18)
spar(cx=2.25, col='white', mfrow=c(2,2))
par(cex.main=2.25)
spar(cx=2.5, col='white')
s <- seq(1, nrow(lambda.x), by=10)

for (i in seq(99)) {
    vel <- data.matrix(fread(paste0("p0-a3-d8--v", i, "_vel.dat")))[,1]
    matplot(vel[s], lambda.x[s,1]/lambda.t[s,1], ,pch=19, cex=0.2, xlab=tex("\\theta_0"), ylab=tex("\\xi / \\tau"), xlim=c(-pi/2, pi/2))
}
dev.off()

v <- c(tex("\\alpha = \\pi/4"),
       tex("\\alpha = \\pi/3"),
       tex("\\alpha = \\pi/2"),
       tex("\\alpha = 2\\pi/3"))


L <- getTSD_hist(dx = lambda.t, n = 200, domain=exp.domain)
L <- getTSD_ecdf(dx = lambda.t/tt, n = 1000, domain=exp.domain)
matplot(L$x[,], 1-L$y[,], type='s', lty='solid', lwd=2, col=rbow(), log='xy',)# xlim=c(2000,8e4))
LF <- LinFit(time = L$x[,1:2], mat = 1-L$y[,1:2], range = 300:875, logs=TRUE, show.range=TRUE, plot=TRUE)
LF <- LinFit(time = L$x[,1:2], mat = 1-L$y[,1:2], range = 790:875, logs=TRUE, show.range=TRUE, plot=TRUE)
LF <- LinFit(time = L$x, mat = 1-L$y, range = 500:820, logs=TRUE, show.range=TRUE, plot=TRUE)
LF <- LinFit(time = L$x, mat = 1-L$y, range = 10:100, logs=TRUE, show.range=TRUE, plot=TRUE)
LF <- LinFit(time = L$x, mat = 1-L$y, range = 932:999, logs=TRUE, show.range=TRUE, plot=TRUE)

id <- apply(L$x[,-(1:3)], 2, function(x) min(which(x > 35000)))

plot(diag(sapply(id, function(x) 1-L$y[x,])), type='o', pch=19, log='xy')

###
#   Illustration of the flow definition
###
{
    #vel      <- data.matrix(fread("p0-a1-d2--n0_vel.dat"))
    t <- cumsum(data.matrix(fread("p0-a0-d0-e1--n0_0--t.dat")))
    x        <- data.matrix(fread("p0-a0-d0-e1--n0_0--x.dat"))
    y        <- data.matrix(fread("p0-a0-d0-e1--n0_0--y.dat"))
    lt       <- data.matrix(fread("p0-a2-d8-e0--n0_jump.dat")[,1])
    lx       <- data.matrix(fread("p0-a2-d8-e0--n0_jump.dat")[,2])
    #x        <- ReadDataList(lf("Traj"1, idx=1)
    #y        <- ReadDataList(lf("Traj"1, idx=2)
}

#pdf(paste0(img.path, "nbb_ex.pdf"), width=20, height=10)
s <- seq(1, 1000, by=1)
spar(cx=3, col='white', mfrow=c(1,1))
par(mgp=c(4.2,1.80,0), oma=c(2,0,2,2), mar=c(5,7,2,2))
matplot(t[s], x[s,], type='o', pch=19, cex=0.8, col='white', xlab="t", ylab="x(t)")
old.lb <- 1
mxlt <- max(which(cumsum(lt[,1]) < t[s[length(s)]]))
mxlt <- nrow(lt)
for (i in seq(1,mxlt)) {
    lb <- old.lb
    ub <- max(which(t <= sum(lt[1:i,])))+1
    old.lb <- ub
    ss <- seq(lb, ub, by=1)
    colo <- ifelse(i %% 2 == 0, 'red', 'blue')
    #lwd <- ifelse(colo == 'red', 5, 2)
    lwd <- 3
    lines(t[ss], x[ss,], col=colo, lwd=lwd)
}

clt <- cumsum(lt)
clx <- cumsum(lx)
plot(t, x, ylim=c(min(clx, x), max(clx, x)), type='l') 
lines(clt, clx, col='red')


### 
#   Illustration of flow in trajectory
###
x <- data.matrix(x[,,1])
y <- data.matrix(y[,,1])
PlotChannelBoundary(from=-1, to=2, xVert = c(-1, 0, 1, 1, 0, -1), yVert = c(0, .3, 0, .1, .4, .1), mult=2)
PlotChannelBoundary(from=2, to=40, xVert = c(-1/2, 0, 1/2, 1/2, 0, -1/2), yVert = c(0, .28868, 0, 0.28868, 0.57735, 0.28868))
PlotChannelBoundary(from=0.5, to=4, xVert = c(-1/2, 0, 1/2, 1/2, 0, -1/2), yVert = c(0, .28868, 0, 0.28868, 0.57735, 0.28868))
matplot(x[1:120,1,1], y[1:120,1,1], type='l', add=TRUE, col=rbow(), lty='solid')

PlotChannelBoundary(from=-1, to=0, xVert = c(-1, 0, 1, 1, 0, -1), yVert = c(0, 10, 0, .51, 10.51, .51), mult=2)
PlotChannelBoundary(from=14, to=16, xVert = c(-1/2, 0, 1/2, 1/2, 0, -1/2), yVert = c(0, 1/2, 0, 1/4, 3/4, 1/4))
old.lb <- 1
points(x[1,,1], y[1,,1], col='green', pch=19)
for (i in seq(1,10)) {
    lb <- old.lb
    ub <- max(which(t <= sum(lt[1:i,])))+1
    old.lb <- ub
    print(paste0(lb, " " , ub))
    ss <- seq(lb, ub)
    colo <- ifelse(i %% 2 == 0, 'red', 'blue')
    #colo <- c('red', 'blue', 'magenta')[i]
    lines(x[ss,], y[ss,], col=colo, lwd=2)

    #for (j in seq(length(ss))) {
    #    rx <- round(x[ss[j]])
    #    if (x[ss[j]]-rx < 0) {
    #        arrows(x0=x[ss[j]], y0=y[ss[j]], 
    #    	   x1=x[ss[j]] + 0.035*lhs[1],
    #               y1=y[ss[j]] + 0.035*lhs[2], code=2, length=0.075, lwd=2, col='black')
    #    } else { 
    #        arrows(x0=x[ss[j]], y0=y[ss[j]], 
    #    	   x1=x[ss[j]] + 0.035*lhs[1],
    #               y1=y[ss[j]] - 0.035*lhs[2], code=2, length=0.075, lwd=2, col='black')
    #    }
    #}
}
dev.off()

lhs <- c(1/2, 0.28868)  / sqrt(1/4 + 0.28868^2)
lhs <- c(1/2, 0.5) 
th <- atan2(diff(y) , diff(x))
v <- cbind(cos(th), sin(th))

arrows(x0=-1/2, x1=0, y0=0, y1=0.28868, code=2, length=0.1, lwd=2)

j <- 1:5
v[j,1] * lhs[1] + v[j,2] * lhs[2]

###
#   Evolution of jumps
###
O <- order(jt[,2])
sjt <- jt[O,2]
sjd <- jd[O,2]
v <- jd / jt
usjt <- unique(sjt)[seq(1, length(O), by=500)]
mids <- seq(10^3, 10^5, by=1000)
mids <- seq(-0.6, 0.6, by=0.025)
ub <- mids + 0.025/xlim=c(1e4, 1e6),2
lb <- mids - 0.025/2
result <- std <- c()
j <- 3
for (i in seq(length(mids))) {
    #id <- which(usjt[i] == sjt)

    id <- which(v[,j] < ub[i] & v[,j] > lb[i])
    #v <- c(v, mean(sjd[id]^2))
    result <- c(result, mean(jt[id]))
    std <- c(result, sqrt(var(jt[id])))
}
s <- seq(1, nrow(v), by=1)
plot(v[s,j], jt[s,j], pch=19, col='black', log='y', cex=0.2, #ylim=c(0,1e4),
     xlab=tex("\\xi / \\tau"), 
     ylab=tex("\\tau"))
x <- seq(0, 1, length=1000)
lines(c(x,-x), 1e3/c(x,x), col='blue')
legend('topright', legend=tex("\\sim 1/\\tau"), col='blue', lty='solid')
lines(mids, result, type='o', pch=19, col='blue', lwd=2)
ErrorBars.x(x = mids, y = result, sdv=0.025/2)
ErrorBars.y(x = mids, y = result, sdv=std)

L <- getTSD_hist(dx = sjd[id], n=50, plot=TRUE, log='')
plot(usjt[1:length(v)], v, type='l', log='',) #ylim=c(1e5, 1e11))
LF <- LinFit(time = usjt[1:length(v)], mat = v, range = 1:100, logs=TRUE, show.range=TRUE, plot=TRUE)
plot(mids[1:length(v)], v, type='o', log='xy', pch=19) #ylim=c(1e5, 1e11))
LF <- LinFit(time = mids[1:length(v)], mat = v, range = 20:100, logs=TRUE, show.range=TRUE, plot=TRUE)


###
#   Distribution of mean perpendicular distances > \tau
#   We show that large windows => greater deviations
###
L <- getTSD_ecdf(dx = mepd, n= 1000, domain=exp.domain)

w <- c(10, 100, 1000)
getIDBetween <- function(x, lb, ub) { which(x > lb & x < ub) }
taus.gt <- c(10000, 100000, 500000)
taus.gt <- c(0.1, 0.2, 0.3, 0.4, 0.5)
ub <- taus.gt + 0.1*taus.gt
lb <- taus.gt - 0.1*taus.gt

mean.pds <- list()
for (j in seq(ncol(jt))) {
    meps <- list()
    for (i in seq(length(taus.gt))) {
	#id <- getIDBetween(jt[,j], lb[i], ub[i])
	#id <- getIDBetween(jd[,j], lb[i], ub[i])
	id <- getIDBetween(abs(jd[,j])/jt[,j], lb[i], ub[i])

	print(id)
	if (length(id) > 0)
	    meps[[i]] <- mepd[id,j]
    }
    mean.pds[[j]] <- list.cbind.na(meps)
}

# All window sizes, for a range of jump times
mpds <- list.cbind.na(lapply(mean.pds, function(x) x[,2]))
L <- getTSD_ecdf(dx = mpds, n= 1000, domain=exp.domain)

# All range of jump times for fixd window size
L <- getTSD_ecdf(dx = mean.pds[[3]], n= 1000, domain=exp.domain)


matplot(L$x, 1-L$y, type='s', lwd=2, lty='solid', log='x', col=rbow(), 
   	xlab=tex("<D_p>"),
    	ylab="Prob.",
    	main=TeX("$\\tau \\in \\[9 \\cdot 10^4, 11 \\cdot 10^4\\]$"))

{
    pdf(paste0(img.path, "mpd.pdf"), width=12, height=10)
    spar(cx=2, col='white')
    par(cex.main=2)
    matplot(L$x, 1-L$y, type='s', lwd=2, lty='solid', log='x', col=rbow(), 
    	xlab=tex("<D_p>"),
    	ylab="Prob.",
    	main=TeX("$\\tau \\in \\[9 \\cdot 10^4, 11 \\cdot 10^4\\]$"))
    abline(v=10, lwd=2)
    arrows(x0 = 10, x1 = 6, y0 = c(.25, .5, 0.75), y1 = c(.25, .5, 0.75), lwd=2, angle = 45, code = 2, length = 0.15)
    legend('topright', legend=TeX(paste0("$w = ", w)), lty='solid', lwd=4, col=rbow(), cex=2)
    dev.off()
}

matplot(L$x, 1-L$y, type='s', lwd=3, lty='dashed', col=rbow(), add=TRUE)
matplot(L$x, 1-L$y, type='s', lwd=3, lty='dotted', col=rbow(), add=TRUE)

legend('topright', legend=TeX(paste0("$\\tau \\in \\[", lb, ", ", ub, "\\]$")),
       lty='solid', lwd=4, col=rbow(), cex=3)


t <- seq(0, 1e7, length=(1e7+1)/10)
x <- ReadDataList(lf("w10_0--x.dat"))
mm <- getMM(x = t, y = x, w = 10)
xi <- mm$df$yr - mm$df$yl
ShowJT(df=mm$df[c(2,3),], eps=0)
ShowJT(df=mm$df, eps=0)
mmdf <- data.matrix(mm$df)
lines(mmdf[c(2,3),1:2], mmdf[c(2,3),3:4], col='magenta')
mm$tau


small.dev  <- apply(mpd, 2, function(x) length(which(x < 5)))
num.sample <- apply(mpd, 2, function(x) length(which(!is.na(x))))
small.dev / num.sample

matplot(jump[,2], jump[,4:5], pch=19, cex=1)
matplot(jt[,-2], mpd[,-2], pch=19, cex=0.5)
L <- getTSD_ecdf(dx = mpd, n = 1000, domain=exp.domain)
matplot(L$x, 1-L$y, type='s', log='x', col=rbow(), lwd=2)
O <- order(jt[,1])
jt[O,]

mdp <- ReadData(lf("MeanPDPP.dat"), Tpose = FALSE)
vdp <- ReadData(lf("VarPDPP.dat"), Tpose = FALSE)

mdp <- vdp

tau.to.check <- c(100, 1000, 5000, 10000)
slice <- 2
w <- 10^seq(1, 4)
wx <- mdp[,,slice]
wx <- mdp[,slice,]

l <- list()
for (i in seq(ncol(wx))) 
    l[[i]] <- wx[-which(wx[,i] == 0),i]
l <- list.cbind.na(l)

L <- getTSD_ecdf(dx = l[,1:4], n=750, domain=exp.domain)

pdf(paste0(img.path, "Dp_deviation_t1k.pdf"), width=14, height=12)
spar(cx=2.15, col='white')
par(mgp=c(3.5,1.25,0), cex.main=2.15)
matplot(L$x[,], 1-L$y[,], type='s', lwd=2, lty='solid', 
	xlab=TeX("$< D_{p} >$"), 
	ylab="Probability",
	main=TeX(paste0("$\\tau \\in \\[", tau.to.check[slice]-25, ", ", tau.to.check[slice]+25, "\\]$")), 
	log='x', col=rainbow)

abline(h=0, lwd=2)
legend('topright', legend=c(paste0("w = ", w[1:ncol(L$x)]), "Mean"), lty=c(rep('solid', ncol(L$x)), 'dashed'), lwd=4, col=c(rainbow[1:ncol(L$x)], 'black'), cex=3)
#legend('topright', legend=paste0("t = ", tau.to.check[1:ncol(L$x)]), lty='solid', lwd=4, col=rainbow, cex=3)

apply(l, 2, function(x) length(which(!is.na(x))))

means <- apply(l[,1:4], 2, mean, na.rm=TRUE)
for (i in seq(length(means))) {
    wmeans <- min(which(L$x[,i] > means[i]))
    yc <- 1 - L$y[wmeans,i]
    xc <- L$x[wmeans,i]
    lines(c(xc, xc), c(0, yc), 
          col='black', lty='dashed', lwd=2)
    lines(c(0.001, xc), c(yc, yc), 
          col='black', lty='dashed', lwd=2)
}
dev.off()


###
##  Instead of looking at mean deviation, consider max deviation
##  PATH: S-JumpStatistics/Data/small_large_deviations/
###

img.path = "~/Dropbox/PHD/Papers/JumpStatistics/Images/"
setwd("../Studies/S-JumpStatistics/Data/small_large_deviations")


mpd <- data.matrix(fread("a1-d2--w10_MaxPD.dat", fill=TRUE))
mpd <- data.matrix(fread("a0-d2--w10_MaxPD.dat", fill=TRUE))
mpd <- ReadDataList(lf("MaxPD"), idx = 2)

ttc <- c(5e4, 250000, 500000)

L <- getTSD_ecdf(dx = mpd, n=750, domain=exp.domain)
L <- getTSD_ecdf(dx = data.matrix(mpd[,1]), n=750, domain=exp.domain)

spar(cx=2.15, col='white')
par(mgp=c(3.5,1.25,0), cex.main=2.15)
matplot(L$x[,], 1-L$y[,], type='s', lwd=2, lty='solid', 
	xlab=TeX("$< D_{p} >$"), 
	ylab="Probability",log='x', col=rbow())

legend('topright', legend=TeX(paste0("$\\tau \\in \\[", ttc-1e4, ", ", ttc+1e4, "\\]$")), col=rbow(), lwd=2)
	#main=TeX(paste0("$\\tau \\in \\[", tau.to.check[slice]-25, ", ", tau.to.check[slice]+25, "\\]$")), 
	#log='x', col=rainbow)


matplot(L$x[,], 1-L$y[,], type='s', lwd=2, lty='dashed', add=TRUE,
	ylab="Probability",log='x', col=rbow)

###
##  Comparison of jump detections for small/large deviations of the MEAN perp. dist
##  PATH: S-JumpStatistics/Data/small_large_deviations/
###
t <- seq(0, 1e6, length=1e6)
x <- data.matrix(fread("p0-a1-d2--333-sd-x.dat"))
x <- data.matrix(fread("p0-a1-d2--w10_3536-sd-x.dat"))
x <- data.matrix(fread("p0-a1-d2--w100_3602-ld-x.dat"))

x <- ReadData(lf(".dat"), Tpose = FALSE)[,1,]
x[,1:4] <- x[,c(1,4,2,3)]

s <- seq(1e6)
s <- seq(1000, 4000)
matplot(t[s], x[s,], type='l')


img.path = "~/Dropbox/PHD/Papers/JumpStatistics/Images/"
#pdf(paste0(img.path, "small_large_deviations.pdf"), width=18, height=16)
spar(cx=3, col='white', mfrow=c(2,2))
par(mar = c(8, 8, 1.3, 1.3), oma = c(2, 2, 2, 2), mgp=c(6,2,0))
w <- c(10, 100, 100, 100)
ww <- c(10, 10, 100, 100)
Dp <- c(0.191, 0.163, 5.84, 6.32)
for (i in seq(ncol(x))) {

    mm <- getMM(x=t, y=x[,i], w=w[i])
    
    id <- min(which(mm$tau > 975 & mm$tau < 1025))

    s <- seq(mm$df$xl[id]-200, mm$df$xr[id]+200)
    print(s)
    matplot(t[s], x[s,i], type='l', col='black', lwd=2, lty='solid', xlab='t', ylab='x')

    dff <- data.matrix(mm$df)
    ss <- seq(id, id)
    lines(dff[ss,1:2], dff[ss,3:4], col='magenta', lwd=1.5)
    ss <- seq(id, id)
    ShowJT(df=mm$df[ss,], text=FALSE, eps=0, lwd = 1.5, pt.cx=1.5)


    dy <- abs(max(x[s,i]) - min(x[s,i]))
    dx <- abs(max(t[s]) - min(t[s]))
    text(x=min(t[s] + 0.01*dx), y=max(x[s,i]) - 0.01*dy, labels=i, cex=2.5)

    text(x=min(t[s] + 500 + 0.3*dx), y=max(x[s,i]) - 0.6*dy, labels=paste0("w = ", ww[i]), cex=2.5)

    text(x=min(t[s] + 500 + 0.3*dx), y=max(x[s,i]) - 0.7*dy, labels=TeX(paste0("$<D_p> = $", Dp[i])), cex=2.5)
}
dev.off()

# MOMENTS
Mom <- ReadData(lff("a", "Moment"))
Mom <- ReadData(lf(c("Moment_3")))

T <- 1e5
t <- seq(0, T, length=T/10 + 1)
msd <- apply(Mom, 3, rowMeans)
coefs <- sapply(1:ncol(msd), function(x) {
    df <- data.frame(t=t, x=msd[,x])
    coef(nlsLM(x ~ b*t^a, data=df, start=list(a=1.5, b=1)))
})
matplot(t, msd, type='l', col=rbow(), log='xy', lwd=2, lty='solid')

mu <- getPowerScaling(t = t, M = msd, z = 0.9)
matplot(t, mu[,], type='l', ylim=c(1, 2), col=rbow(), ylab=tex("\\mu(t)"), xlab="t", lwd=3, lty='solid')#, add=TRUE)
abline(h=c(1.66,1.7), lty='dashed', lwd=2)
text(x=50000, y=1.63, labels=tex("\\mu = 1.66"), cex=3)
text(x=80000, y=1.715, labels=tex("\\mu = 1.7"), cex=3)

p <- seq(1/2, 10, by=1/2)
mplot(p, mu[800,], pch=19, type='o', asp=1 )
points(p, mu[800,], pch=19, type='o', col='blue')

df <- data.frame(x=p, y=mu[800,])
LF1 <- lm(y ~ 0 + ., data=df, subset=1:6)
LF2 <- LinFit(mat = mu[800,], time = p, range=6:20, logs = FALSE)
pp <- seq(-10, 20, by=1/2)
lines(pp, coef(LF1)*pp, lty='dashed', col='black', lwd=2)
lines(pp, LF2$ax * pp + LF2$b, lty='dashed', col='black', lwd=2)

LinFit(mat = mu[8,], time = p, range=1:6, logs = FALSE, plot=TRUE)
LinFit(mat = mu[8,], time = p, range=6:20, logs = FALSE, plot=TRUE)
lines(p, 0.92*p - 0.14)
lines(p, 0.85*p - 0.08)

### 
#  Alternating velocity model
###


### Get base piece ###
getBasePiece <- function(v, x0, n, alpha) {
    times <- power.law(alpha = alpha, x0 =x0, n = n)
    times <- runif(n=n, min=0, max=2)
    vs <- rep(v, length(times)) * (-1)^(seq(1,length(times)))
    t <- c(0, cumsum(times))
    x <- c(0, cumsum(vs * times))
    return(list(x=t, y=x))
}

getRW <- function(T=1e5, b) {
    res <- list(x=0, y=0)
    curr.time <- b$x[length(b$x)]
    base.time <- curr.time
    i <- 1
    while(res$x[length(res$x)] < T) {
	jump.time <- power.law(alpha = alpha, x0=base.time, n = 1) # Duration of time from current time to jump
	while (jump.time > T)
	    jump.time <- power.law(alpha = alpha, x0=base.time, n = 1) # Duration of time from current time to jump
	#jump.time <- runif(n = 1, min = x0, max = 100*x0)
	curr.time <- jump.time + res$x[length(res$x)]
	print(paste0("Current time = ", curr.time, ", jump.time = ", jump.time)) 
	bsf <- list(x=b$x, y= b$y)  # Base building block
	bsb <- list(x=b$x, y=-b$y) # Base building block
	#bsf <- getBasePiece(v, x0, n, alpha) 
	#bsb <- getBasePiece(v, x0, n, alpha) 

	r <- runif(1)
	#if (i %% 2 == 0) {
	if (r < .5) {
	    #st <- bsf$x[length(bsf$x)]
	    st <- res$x[length(res$x)] 
	    while (st < curr.time) {
		bsf <- continue(bsf$x, bsf$y)
		#st <- bsf$x[length(bsf$x)]
		st <- res$x[length(res$x)] + bsf$x[length(bsf$x)]
	    }
	    res <- concat(x0=res$x, y0=res$y, x1=bsf$x, y1=bsf$y)
	} else {
	    #st <- bsb$x[length(bsb$x)]
	    st <- res$x[length(res$x)] 
	    while (st < curr.time) {
		bsb <- continue(bsb$x, bsb$y)
		#st <- bsb$x[length(bsb$x)]
		st <- res$x[length(res$x)] + bsb$x[length(bsb$x)]
	    }
	    res <- concat(x0=res$x, y0=res$y, x1=bsb$x, y1=bsb$y)
	}
	i <- i + 1
    }
    return(res)
}

{ # Manual
    v <- 1/2
    x0 <- 10
    n <- 8
    alpha <- 0.9
    b <- getBasePiece(v, x0, n, alpha)
    pdf(paste0(img.path , "steps_rw.pdf"), height=15, width=18)
    spar(mfrow=c(2,2), cx=2.5)
    plot(b$x, b$y, type='o', pch=19, lwd=2, xlab=tex("t"), ylab=tex("x_t"))
    #dev.off()
}

T <- 1e2
nwalks <- 10
X <- matrix(0, nrow=nwalks, ncol=101)
for (i in seq_along(1:nwalks)) {

    b <- getBasePiece(v, x0, n, alpha)
    rw <- getRW(T, b)
    rw <- getRW(1e2+rw$x[length(rw$x)], rw)
    rw <- getRW(1e3+rw$x[length(rw$x)], rw)
    rw <- getRW(1e4+rw$x[length(rw$x)], rw)
    rw <- getRW(1e5+rw$x[length(rw$x)], rw)

    s <- seq(length(rw$x))
    plot(rw$x[s], rw$y[s], type='l', pch=19, lwd=2, xlab=tex("t"), ylab=tex("x_t"))
    
    time <- seq(0, T, by=100)
    f <- approxfun(x=rw$x, y=rw$y)(time)

    X[i,] <- f
    #X[[i]] <- f
}
#X <- list.cbind.na(X)

concat <- function(x0, x1, y0, y1) {
    x.concat <- c(x0, x0[length(x0)] + x1)
    y.concat <- c(y0, y0[length(x0)] + y1)
    return(list(x=x.concat, y=y.concat))
}
continue <- function(x, y, s=1) {
    x.cont <- c(x, x[length(x)] + x)
    y.cont <- c(y, y[length(x)] + s*y)
    return(list(x=x.cont, y=y.cont))
}


spar()
matplot(time, tpose(X), type='l', col=rbow(), lwd=2, lty='solid', xlab="t", ylab=tex("x_t"))
L <- getTSD_hist(dx = (X[,c(100, 500, 1000)/10]), n=100, t = time[c(100,500,1000)/10], eta=1)#, domain=exp.domain)
matplot(L$x, L$y, type='s', lwd=2, log='y')

p <- seq(1/2, 7, by=1/2)
mom <- sapply(p, function(x) colMeans(abs(X)^x))
matplot(time, mom, log='xy', type='l')

mu <- getPowerScaling(t = time, M = mom, z = .8)
matplot(time, mu, type='l', lty='solid', lwd=2, log='x')
abline(v=time[80])
plot(p, mu[80,], pch=19, type='o')


##
#  FIG. JUMP STATISTICS -- \tau > \tau_0
## 

jt <- ReadDataList(lf("a8"), idx = 4)
jd <- ReadDataList(lf("a8"), idx = 5)

W <- c(5, 50, 500, 5000)
jt <- ReadDataList(lf("--w"), idx = 4)
jd <- ReadDataList(lf("--w"), idx = 5)

mm.jt <- apply(jt, 2, function(x) c(min(x, na.rm=TRUE), max(x, na.rm=TRUE)))
mm.jd <- apply(abs(jd), 2, function(x) c(min(x, na.rm=TRUE), max(x, na.rm=TRUE)))
mm.jd / mm.jt

id <- apply(abs(jd), 2, function(x) which(x < 1000))
jd0 <- list.cbind.na(sapply(1:length(id), function(x) jd[-id[[x]],x]))
jt0 <- list.cbind.na(sapply(1:length(id), function(x) jt[-id[[x]],x]))

L <- getTSD_ecdf(dx = abs(jt), n=1000, domain=exp.domain, t = 1)
L <- getTSD_ecdf(dx = abs(jt0), n=1000, domain=exp.domain, t = 1)
mplot(L$x, 1-L$y, type='s', log='xy', xlim=c(1e3, 1e5))

LF <- LinFit(mat = 1-L$y[,], time=L$x[,],  plot=TRUE,	range=2:400, show=TRUE)
LF <- LinFit(mat = 1-L$y[,], time=L$x[,],  plot=TRUE,	range=400:700, show=TRUE)

L <- getTSD_hist(dx = (jd)/jt, n=100)#, domain=exp.domain, t = 1)
mplot(L$x, L$y, type='s')

##
#  FIG. TIME AVERAGES --- We consider the time averaged MSD
#  of trajectories determined from a long jump event
## 

# Reconstructed trajectories
j <- 3 # Window size
X.t <- t <- list()
for (i in seq(nrow(id)-1)) {
    dom <- c()
    if (i > 1)
	dom <- seq(id[i,j]+1, id[i+1,j])
    else
	dom <- seq(id[i,j], id[i+1,j])
    X.t[[i]] <- cumsum(c(x0[1,j], jd[dom,j]))
    t[[i]]   <- cumsum(c(t0[1,j], jt[dom,j]))
}
X.t <- list.cbind.na(X.t); X.t <- X.t[-nrow(X.t),]
t <- list.cbind.na(t); t <- t[-nrow(t),]
X.t <- cumsum(c(x0[1,j], jd[,j]))
t   <- cumsum(c(t0[1,j], jt[,j]))

tr <- seq(0, 1e8, by=100)
x <- data.matrix(fread(lf("0_x.dat")))
mplot(tr, x)
mplot(t[,1], X.t[,1], type='o', col='blue', add=TRUE)
mplot(t[], X.t[], type='o', col='blue', cex=0.5, add=TRUE)

tamsd <- deltas <- list()
for (k in seq(ncol(X.t))) {
    nonna <- length(which(!is.na(X.t[,k])))
    Delta <- seq(2, nonna-1, by=1) # Lag time
    msd <- c()
    for (j in seq(length(Delta))) {
	temp <- 0
	s <- seq(1, nonna - Delta[j])
	msd[j] <- (1/(nonna - Delta[j])) * sum((X.t[s+Delta[j],k] - X.t[s,k])^2)

	#for (i in seq(1, nonna - Delta[j])) 
	#    temp <- temp + (X.t[i+Delta[j],k] - X.t[i,k])^2

	#msd[j] <- (1/(nonna - Delta[j])) * temp
    }
    tamsd[[k]] <- msd
    deltas[[k]] <- Delta
}
tamsd <- list.cbind.na(tamsd)
deltas <- list.cbind.na(deltas)

#pdf(paste0(img.path, "tamsd_ex.pdf"), height=12, width=12)
mplot(deltas, tamsd, lwd=3, type='l', log='xy', xlab=tex("\\Delta"), ylab=tex("\\bar{\\delta}^2"))
mplot(deltas, tamsd, lwd=3, type='l', lty='dashed', add=TRUE)
mplot(deltas, tamsd, lwd=3, type='l', lty='dotted', add=TRUE)
dev.off()
