img.path = "~/Dropbox/PHD/Papers/JumpStatistics/Images/AdditionalResults/"
setwd("../Studies/S-JumpStatistics/Data/")

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

jt      <- ReadDataList(lff("c0", "jump.dat"), idx=1)
jd      <- ReadDataList(lff("c0", "jump.dat"), idx=2)
jv	<- data.matrix(abs(jd) / jt)


# @brief wgt = @x which-greater-than. Becomes which less then if @x -> -@x
# @x[vector] 
# @eps[numeric] Th
xwgt <- function(x, eps) {
    id <- which(x > eps)
    if (length(id) > 0) {
	return(x[id,])
    } else {
	cat("Threshold condition failed. Returning input.\n")
	return(x)
    }
}

# @return Head and tail of a vector
# @x[vector] 
# @n Number of elements from head/tail
headtail <- function(x, n) {
    return(c(head(x, n), tail(x,n)))
}

matrix_convert <- function(x) { if (!is.matrix(x)) { return(data.matrix(x)) } }

getFastSlowJumpVelocity <- function(xi0, tau0, jt, jd, n=3) {
    if (length(xi0) != length(tau0)) { stop("xi0 must equal tau0 in length") }
    jt <- matrix_convert(jt)
    jd <- matrix_convert(jd)
    jv <- jd / jt

    Vf.xi0 <- Vf.tau0 <- vector(mode='list', length=length(xi0))
    for (i in seq(length(xi0))) {
	id <- which(jd > xi0[i]) # Jump displacements > xi0
        Vf.xi0[[i]]  <- headtail(sort(jv[id,]), n)
	id <- which(jt > tau0[i])
        Vf.tau0[[i]] <- headtail(sort(jv[id,]), n)
    }
    Vf.xi0  <- list.cbind.na(Vf.xi0)
    Vf.tau0 <- list.cbind.na(Vf.tau0)
    cat("Fast/slow jumps conditioned by a smallest DISPLACEMENT xi0.\n")
    print(noquote(format(Vf.xi0, digits=2, scientific=FALSE)))
    cat("\nFast/slow jumps conditioned by a smallest TIME tau0.\n")
    print(noquote(format(Vf.tau0, digits=2, scientific=FALSE)))

    result <- cbind(slow.tau0=Vf.tau0[1,],  fast.tau0=Vf.tau0[2*n,], 
		   slow.xi0=Vf.xi0[1,], fast.xi0=Vf.tau0[2*n,])
    print(result)
    return(result)
}

xi0 <- tau0 <- 1000
xi0 <- tau0 <- c(0, 1, 10, 100, 1000)

# A single ensemble of jumps with no conditioning. Shows limiting behaviour of fast/slow jumps 
# in terms of small jump displacements or times.
{
    #pdf(paste0(img.path, "fsjv-tau0-xi0.pdf"), width=18, height=9)
    xi0 <- tau0 <- c(seq(0, 10), seq(11, 100, by=10), seq(100, 1000, by=100))
    fsjv <- getFastSlowJumpVelocity(xi0, tau0, jt, abs(jd), n=3)
    
        spar(mfrow=c(1,2))
        mplot(xi0, fsjv[,c(1,3)], type='p', log='y', cex=2, spa=FALSE,
              xlab=tex("\\xi_0, \\tau_0"),
              ylab=tex("V_s"))
        legend('bottomright', legend=c(tex("\\tau > \\tau_0"), tex("| \\xi | > \\xi_0")), pch=19, col=rbow(), cex=2)
        mplot(xi0, fsjv[,-c(1,3)], type='p', log='', cex=1.5, spa=FALSE,
              xlab=tex("\\xi_0, \\tau_0"),
              ylab=tex("V_f"))
        legend('topright', legend=c(tex("\\tau > \\tau_0"), tex("| \\xi | > \\xi_0")), pch=19, col=rbow(), cex=2)
    #dev.off()
}

## Several ensembles of jumps conditioned for multiple xi0/tau0 values. 
## Checking if results match with the above 
{
    fsjv <- c()
    for (i in seq(ncol(jt))) 
        fsjv <- rbind(fsjv, getFastSlowJumpVelocity(xi0, tau0, jt[,i], abs(jd[,i]), n=3))
    mplot(c(1, 10, 100, 1000), fsjv, type='o')
    mplot(fsjv, type='o', add)
}
Vs <- ReadDataList(lff("", "Vf.dat"), idx=1)
Vf <- ReadDataList(lff("", "Vs.dat"), idx=1)
mm <- c(min(Vs), max(Vf))
L <- getTSD_hist(dx = cbind(Vs, Vf), n = 1000)
mplot(L$x, L$y, type='s')



### 
#    xi0 = tau0 = 1000, N = M = 10^5, all rational and irrationals and 15 pore sizes
###
setwd("PaperPlots/Fig-FSJV/")

r <- seq(1/4, 1, length=25)
alpha     <- c(pi/seq(9, 2, by=-1), seq(2, 8)*pi/seq(3,9)) #
alpha.irr <- c(pi/seq(9, 2, by=-1), seq(2, 8)*pi/seq(3,9)) * ((sqrt(5)+1)/2 - 3/5) # x

did <- c(0, 1, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24) + 1
d <- sapply(1/(2*tan(alpha/2)), function(x) x*r[did])
R <- r[did]

a.vf <- a.vs <- list()
#for (i in seq(length(alpha), 2*length(alpha)-1)) { # IRRATIONALS
for (i in seq(0, length(alpha)-1)) { # RATIONALS
    Vf <- ReadDataList(gtools::mixedsort(lff(paste0("-a", i, "-"), "Vf.dat")), idx=1)
    Vs <- ReadDataList(gtools::mixedsort(lff(paste0("-a", i, "-"), "Vs.dat")), idx=1)
    a.vf[[i+1]] <- apply(Vf, 2, max, na.rm=TRUE)
    a.vs[[i+1]] <- apply(Vs, 2, max, na.rm=TRUE)
}
a.vf <- list.cbind.na(a.vf)#[,-(1:15)]
a.vs <- list.cbind.na(a.vs)#[,-(1:15)]

{
    #pdf(paste0(img.path, "RAT--VfVs.pdf"), width=18, height=9)
    spar(mfrow=c(1,2))
    rid <- c(1,6,14)
    id <- c(1,7,8,9,15)

    mplot(alpha[], tpose(a.vf[,])[,rid], spa=FALSE, 
          xlab=tex("\\alpha"), ylab=tex("V_f"),
          cex=2, type='p', xaxt='n', 
          ylim=c(0,1), pch=19)
    axis(side=1, at=alpha[id], labels=alpha.lab[id])#, line = 2) 
    legend('bottomright', legend=paste0("r = ", format(R[rid], digits=2)), pch=c(19), col=rbow(), cex=2)
    
    mplot(alpha.irr[], tpose(a.vs[,])[,rid], spa=FALSE, 
          xlab=tex("\\alpha"), ylab=tex("V_s"),
          cex=2, type='p', xaxt='n', )
          #ylim=c(0,.2))
    axis(side=1, at=alpha[id], labels=alpha.lab[id])#, line = 2) 
    legend('topleft', legend=paste0("r = ", format(R[rid], digits=2)), pch=c(19), col=rbow(), cex=2)
    #dev.off()
}

vb <- list()
for (i in seq(length(did))) { # RATIONALS
    Vf <- ReadDataList(gtools::mixedsort(lff(paste0("d", did[i]-1, "-"), "i4-c1--Vf.dat")), idx=1)
    Vs <- ReadDataList(gtools::mixedsort(lff(paste0("d", did[i]-1, "-"), "i4-c1--Vs.dat")), idx=1)
    L <- getTSD_hist(dx = Vf, n = 1000)
    #mplot(L$x, L$y, type='s', ylim=c(0,10))
    wm <- apply(L$y, 2, function(x) which.max(x))
    vb[[i]] <- sapply(1:length(wm), function(x) L$x[wm[x],x])
}
vb <- list.cbind.na(vb)[1:15,]


mplot(alpha, (vb)[,9], type='o', xaxt='n', ylim=c(0, .8))
axis(side=1, at=alpha[id], labels=alpha.lab[id])#, line = 2) 
abline(h=1/2)
legend('topleft', legend=paste0("r = ", format(R, digits=2)), pch=c(19), col=rbow(), cex=2)
mplot(R, tpose(vb), type='o', xaxt='n', ylim=c(0, 1))

# Compare with ordinary
lambda.t     <- ReadDataList(lf(c(".dat")), idx=1)
lambda.x     <- ReadDataList(lf(c(".dat")), idx=2)

# Distribution of jump speeds
L <- getTSD_hist(dx = abs(lambda.x/lambda.t), n=500)#, domain=exp.domain)
mplot(L$x, L$y, type='s')

CJS <- sapply(1:ncol(L$y), function(i) {
		   nna <- length(which(!is.na(L$y[,i])))
		   m1 <- sort(L$y[,i])[nna]
		   L$x[which(L$y[,i] == m1),i]})
CJS <- tpose(CJS)

#pdf(paste0(img.path, "characteristic_speeds4.pdf"), width=15, height=12)
# r <- 1/2 -- red and blue circles, rational and irrational
mplot(alpha[1:15], abs(CJS[1:15]), type='p', xaxt='n', xlab=tex("\\alpha"), ylab=tex("\\nu"), cex=2, col=rbow()[1], ylim=c(0,1))
mplot(alpha, vb[,9], add=TRUEq)


# Write data in 3D format for gnuplot
#for (i in seq(ncol(vb))) {
#    fwrite(file="vb3d.dat", sep = " ", x=cbind(alpha, R[i], vb[,i]), append = TRUE, col.names=FALSE)
#    write("",file="vb3d.dat",append=TRUE)
#}







