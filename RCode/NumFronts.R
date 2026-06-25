img.path = "~/Dropbox/PHD/Papers/JumpStatistics/Images/AdditionalResults/"
setwd("../Studies/S-JumpStatistics/Data/")

alpha <- c(pi/seq(9, 2, by=-1), seq(2, 8)*pi/seq(3,9) )
alpha.lab <- c(tex("\\pi/9"), 
	       tex("\\pi/8"), 
	       tex("\\pi/7"), 
	       tex("\\pi/6"), 
	       tex("\\pi/5"), 
	       tex("\\pi/4"), 
	       tex("\\pi/3"), # 7
	       tex("\\pi/2"), # 8
	       tex("2\\pi/3"),# 9 
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
	       tex("(\\varphi - 3/5) \\pi/3"), # 7
	       tex("(\\varphi - 3/5) \\pi/2"), # 8
	       tex("(\\varphi - 3/5) 2\\pi/3"),# 9
	       tex("(\\varphi - 3/5) 3\\pi/4"), 
	       tex("(\\varphi - 3/5) 4\\pi/5"), 
	       tex("(\\varphi - 3/5) 5\\pi/6"), 
	       tex("(\\varphi - 3/5) 6\\pi/7"), 
	       tex("(\\varphi - 3/5) 7\\pi/8"), 
	       tex("(\\varphi - 3/5) 8\\pi/9"))

alpha.ir <- (((sqrt(5)+1)/2) - 3/5) * pi/2
r <- seq(1/4, 1, length=50) # = d / dy

# Data 
nf      <- ReadDataList(lff("a22", "NPT.dat"), idx=1)
jv      <- ReadDataList(lff("a7", "NPT.dat"), idx=2)
jvstd   <- ReadDataList(lff("a7", "NPT.dat"), idx=3)

M <- 1e5
n <- seq(1, 10000, length=100)
t <- seq(0, 1e4, length=101)
t <- exp(seq(0, log(1e4), length=10))-1
nf <- data.matrix(rowSums(nf))


# Example of average ballistic speed as a function of n
{
    #pdf(paste0(img.path, "vbalofn.pdf"), width=9, height=9)
    n <- n[-1]
    jv <- jv[-1,]
    id <- sapply(c(1/4, 1/2, 3/4), function(x) min(which(r >= x)))
    mplot(n, jv[,id], type='o', ylim=c(0.35, 0.45), xlab=tex("n"), ylab=tex("\\nu"))
    s <- c(2, 24, 49, 74, 98)
    for (i in seq(length(s))) {
	ss <- c(s[i] - 1, s[i], s[i] + 1)
	sapply(1:length(ss), function(x) ErrorBars.y(n[ss[x]], jv[ss[x],id[x]], sqrt(jvstd / nf)[ss[x],id[x]]))
    }
    legend('topleft', legend=paste0("r = ", format(r[id], digits=2)), pch=c(19), col=rbow(), cex=2)
    dev.off()
}

# Example of num propagating decaying with n
{
    #pdf(paste0(img.path, "Nofn-example.pdf"), width=9, height=9)
    mplot(t[1:nrow(nf)], nf[,]/M, type='o', log='xy', #ylim=c(0,0.1),
    	xlab=tex("n"), ylab="N(n)")
    LinFit(time=t, nf/M, range=3:60, plot=TRUE, lty='dashed')
    
    df <- data.frame(x=t, y=nf/M)
    f <- function(x, b, gamma) b*x^(-gamma)
    f <- function(x, b, gamma) b*exp(-x^(gamma))
    NLF <- nlsLM(y ~ f(x,b,gamma), start=list(b=0, gamma=0), data=df, subset=1:100)
    lines(df$x, f(df$x, coef(NLF)[1], coef(NLF)[2]))

    legend('topright', pch=19, lty='solid', lwd=2,
           legend=alpha.lab[7:9], col=rbow(), cex=2)
    #dev.off()
}


# Analysis of the mean jump velocity obtained from N(n) propagating
{
    # The jump velocities jv are themselves means. Compute the mean of the means
    # for each width d. This is an average over the number of collisions n which
    # technically corresponds to a different experiment. 

    #pdf(paste0(img.path, "nubar.pdf"), width=9, height=9)
    jv      <- ReadDataList(lff("a7", "TB.dat"), idx=2)
    mu   <- apply(jv, 2, mean) # Mean of means
    std  <- sqrt(apply(jv, 2, var)) 
    sem  <- std / length(mu) # SEM
    mplot(r, mu, type='o', ylim=c(0.3, 0.425), xlab=tex("r"), ylab=tex("\\bar{\\nu}"))
    ErrorBars.y(r, mu, sem, col=rbow()[1])

    jv      <- ReadDataList(lff("a22", "TB.dat"), idx=2)
    mu   <- apply(jv, 2, mean) # Mean of means
    std  <- sqrt(apply(jv, 2, var)) 
    sem  <- std / length(mu) # SEM
    mplot(r, mu, col=rbow()[2], type='o', add=TRUE)
    ErrorBars.y(r, mu, sem, col=rbow()[2])
    dev.off()
}


{
    nf      <- ReadDataList(lff("a", "TB.dat"), idx=1)
    jv      <- ReadDataList(lff("a", "TB.dat"), idx=2)
    jvstd   <- ReadDataList(lff("a", "TB.dat"), idx=3)
    #pdf(paste0(img.path, "nofn--ddep.pdf"), width=18, height=9)
    LF <- LinFit(time=n, nf/M, range=2:100, plot=TRUE,show.range = TRUE, lty='dashed', logs=TRUE)

    fwrite(x=cbind(r, LF$ax[1:50], LF$ax[-(1:50)]), file="nPropScalings.tkz", sep=" ", col.names=FALSE)
    mplot(n,nf[,1], log='xy')
    spar(mfrow=c(1,2))
    mplot(r, LF$ax[1:50], ylim=c(-0.9, -0.75), type='o', spa=FALSE, xaxt='n', xlab=tex("r"), ylab=tex("\\gamma"))
    mplot(r, LF$ax[-(1:50)], col=rbow()[2], type='o', add=TRUE, spa=FALSE)
    axis(side=1, at=c(0.25, 0.5, 0.75, 1))

    mplot(r, LF$b[1:50],  type='o', ylim=c(0.9, 2.5), spa=FALSE, xaxt='n', xlab=tex("r"), ylab=tex("b"))
    mplot(r, LF$b[-(1:50)],  type='o', col=rbow()[2], spa=FALSE, add=TRUE)
    axis(side=1, at=c(0.25, 0.5, 0.75, 1))
    dev.off()
}



mplot(diff(mu), diff(LF$ax[1:50]), type='p')
cor.test(diff(mu), diff(LF$ax[1:50]))

#print(xtable::xtable(cbind(LF$ax, LF$b)), file=paste0(img.path, "nofn--scalings.dat"))


# Image of propagating trajectories

{
    dx <- 1/2
    varphi <- (sqrt(5) + 1)/2
    alpha <- c(pi/seq(9, 2, by=-1), seq(2,8)*pi/seq(3,9))
    alpha <- c(alpha, alpha * (varphi - 3/5))
    ratio <- seq(1/4, 10, length=40)
}

KeyID <- c(7,1,0)	# Input id 
{
    a <- alpha[KeyID[1]+1] 
    r <- ratio[KeyID[2]+1]
    dy <- 1/(2*tan(a/2))
    d <- dy * r
    PlotChannelBoundary(xVert=c(-dx, 0, dx, dx, 0, -dx), yVert=c(0,  dy, 0,  d, d+dy, d), from=0, to=5)
}


t <- ReadDataList(lf("Pos"), idx = 1)[,-1,]
x <- ReadDataList(lf("Pos"), idx = 2)
y <- ReadDataList(lf("Pos"), idx = 3)


t <- ReadDataList(lf("Pos"), idx=1)[,-1]
x <- ReadDataList(lf("Pos"), idx=2)[,-1]
L <- getTSD_hist(dx = (x/t)[1:900,], n=300)#, t = 10^seq(2,4))
mplot(L$x, L$y, type='s')#, ylim=c(0,100/2))



L <- getTSD_ecdf(dx = t[,1:10], n=1000, domain=exp.domain)#, t = c(1,10, 100, 1000))
mplot(L$x, 1-L$y, type='s', log='xy')#, ylim=c(0,100/2))
LF <- LinFit(mat = 1-L$y, time=L$x, plot=TRUE, range = 10:300, show.range=TRUE)
