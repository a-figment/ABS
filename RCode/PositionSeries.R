setwd("~/Dropbox/ABS/Studies/S-HPO/Data/Positions/")

T <- function(x) seq(1, 1e7+1, length=x)
i <- 5
L <- 10^seq(2,5) + 1
t <- list.cbind.na(sapply(L, function(x) T(x)))
Xc <- tpose(fread("p0-a5-r1-t0--Xc.dat", sep=" "))[1:10000,1]
Yc <- tpose(fread("p0-a5-r1-t0--Yc.dat", sep=" "))[1:10000,1]
FP <- tpose(fread("p0-a5-r1-t0--FP.dat", sep=" "))[1:10000,1]
P <- ReadDataList(FileList=lf("-Xt.dat"), rowidx = i)

spar(cx=2)
matplot(t[,rev(1:3)], P[,rev(1:3)], type='l', 
	#xlim=c(0,1e7),
	lwd=3, lty='solid',
	xlab=tex("t_k"),
	ylab=tex("x(t_k)"),
	col=rainbow(10)[c(1,3,5)])
y.smooth <- loess(P[,3] ~ t[,3], span=0.1)$fitted
lines(T(L[3]), y.smooth, type='l', col='black')

legend('topleft', legend=c(tex(paste0("\\Delta t = ", rev(L[1:3]-1))),
			   "Smoothed data"), 
	lwd=2, col=c(rainbow(10)[c(1,3,5)], 'black'), lty='solid', cex=2)


par(fig=c(0.4,0.95,0.075,0.5), new=TRUE)
matplot(t[1:1001,3], P[1:1001,3], type='l', 
	lwd=1.5, lty='solid',
	xlab=tex("t_k"),
	ylab=tex("x(t_k)"),
	col=rainbow(10)[5])
matplot(t[1:101,2], P[1:101,2], type='l', add=TRUE,
	lwd=1.5, lty='solid',
	xlab=tex("t_k"),
	ylab=tex("x(t_k)"),
	col=rainbow(10)[3])
matplot(t[1:11,1], P[1:11,1], type='l', add=TRUE,
	lwd=1.5, lty='solid',
	xlab=tex("t_k"),
	ylab=tex("x(t_k)"),
	col=rainbow(10)[1])


tc <- c(0, cumsum(sqrt(diff(Xc)^2 + diff(Yc)^2)))
spar(cx=2)
matplot(tc, Xc, type='l', 
	#xlim=c(0,1e7),
	lwd=2, lty='solid',
	xlab=tex("t_k"),
	ylab=tex("x(t_k)"),
	col=rainbow(10)[c(1,3,5)])
