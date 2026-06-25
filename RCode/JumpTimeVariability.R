setwd("../Studies/S-JumpStatistics/Data/")


####' 
##   PATH: S-JumpStatistics/Data/
##   We consider the variability measure (similarity index) used in  
##   Orchard_2021 on escape times, to probe varibility of jump times  
####' 
f <- gtools::mixedsort(lf("_tau"))[-1]
f <- gtools::mixedsort(lf("_jd"))[-1]
tau <- list()
for (i in seq(length(f)))
    tau[[i]] <- data.matrix(fread(f[i]))
tau <- list.cbind.na(tau)

f <- gtools::mixedsort(lf("-d2--w10_counts.dat"))[2]
counts <- fread(f)
counts <- rowSums(counts)
jump <- counts

f <- gtools::mixedsort(lf("-d2--w10_jump.dat"))
jump <- ReadDataList(f, idx = 3)

tau <- data.matrix(jump)

L <- getTSD_hist(dx=abs(tau), n = 150, prune=FALSE, domain=exp.domain, plot=TRUE)

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


L <- getTSD_hist(dx=data.matrix(omega), n = 100, prune=FALSE)#, domain=exp.domain)
spar(cx=3, col='white')
matplot(L$x, L$y, type='o', pch=19, xlim=c(0,1), col=rainbow(10)[-5], lty='solid', log='y', 
	ylab=tex("P(\\omega_2)"), xlab=tex("\\omega_2"))
W <- 10^seq(2, 5)
legend('topright', paste0("w = ", W), col=rainbow(10)[-5], pch=19, lty='solid')
       


##
#  Plotting jump time distribution and similary index
##

img.path = "~/Dropbox/PHD/Papers/JumpStatistics/Images/"
pdf(paste0(img.path, "jt_variability.pdf"), width=14, height=12)

tau.comb <- data.matrix(as.vector(tau))
L.comb <- getTSD_hist(dx=tau.comb, n = 150, prune=FALSE, domain=exp.domain)
L <- getTSD_hist(dx=tau, n = 150, prune=FALSE, domain=exp.domain, plot=TRUE)
spar(cx=2.5, col='white')
matplot(L$x[,], L$y[,], type='s', lty='solid', log='xy', lwd=2.5, xlab=tex("\\tau"), ylab=tex("\\Psi(\\tau)"), col=rainbow(10)[-5])
matplot(L.comb$x[,], L.comb$y[,], type='s', lty='solid', lwd=2.5, col='black', add=TRUE)

x <- seq(1, 1e7, by=100)
lines(x, exp(10)*x^(-3), lty='dashed', lwd=2)
lines(x, exp(21.6)*x^(-2.7), lty='dashed', lwd=2)

legend('topright', legend=c(paste0("w = ", W[-1]), tex("\\tau^{-2.7}")), 
			    col=c(rainbow(10)[-5][1:3], 'black'), 
			    lty=c(rep('solid' ,3), 'dashed'), 
			    cex=2, lwd=3)

legend('topright', legend=c(paste0("w = ", 10), tex("\\tau^{-2.7}")), 
			    col=c(rainbow(10)[-5][1:3], 'black'), 
			    lty=c(rep('solid' ,3), 'dashed'), 
			    cex=2, lwd=3)

Lw <- getTSD_hist(dx=data.matrix(omega), n = 100, prune=FALSE)#, domain=exp.domain)

par(fig=c(0.075, 0.5, 0.075, 0.55), new=TRUE, cex=0.75)
matplot(Lw$x, Lw$y, type='l', pch=19, xlim=c(0,1), col=rainbow(10)[-5], lty='solid', log='y', 
	ylab=tex("P(\\omega_2)"), xlab=tex("\\omega_2"), xaxt="n")
axis(side=1, at=c(0, 0.25, 0.5, 0.75, 1))
W <- 10^seq(2, 5)
#legend('topright', paste0("w = ", W), col=rainbow(10)[-5], pch=19, lty='solid')
dev.off()
