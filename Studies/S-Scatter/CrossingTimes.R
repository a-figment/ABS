

ArrowAtPeak <- function(peak, lx, ly, eps=0.2) {
    id <- min(which(lx >= peak))
    ypeak <- ly[id]
    fact <- log(ly[1]) - log(tail(ly[which(!is.infinite(log(ly)))])[1])
    print(fact)
    arrows(x0=peak, x1=peak, y0=ypeak+fact*eps*2, y1=ypeak+0.001*fact, length=0.1, code = 2, col='red', lwd=2)
}


# Produce a heat map out of the crossing times

# Crossing times
d <- 1/2
dx <- dy <- 1/2
{
    kappa0 <- function(d, h0, t0) { ceiling((d - h0 - tan(t0))/(d*(tan(t0) -1))) }
    gamma0 <- function(d, h0, t0) { ceiling((h0 - 1 - d)/(d*(1 - tan(t0)))) }
    alpha0 <- function(d, h0, t0) { ceiling((-h0 - tan(t0))/(d*(1 + tan(t0)))) }

    t0 <- function(h,a) -h*cot(a)
    t1 <- function(h,a) d * (alpha0(d,h,a) + 3) * sec(a)
    t2 <- function(h,a) -(h + d)*csc(a)
    t3 <- t7 <- function(h,a) (3*dx + 2*d)*sec(a)
    t4 <- t13 <- t11 <- function(h,a) (d + 2*dx)*sec(a)
    t5 <- function(h,a) -h*csc(a)
    t6 <- t8 <- t9 <- function(h,a) (2*dx + 2*d)*sec(a)
    t10 <- function(h,a) (2*dx - h)*csc(a)
    t12 <- function(h,a) (2*dx + kappa0(d,h,a)*d)*sec(a)
    t14 <- function(h,a) (d + 2*dx - h) * csc(a)
    t15 <- t16 <- function(h,a) 2*dx*sec(a)
    t17 <- function(h,a) (2*dy + 2*d) * csc(a)
    t18 <- function(h,a) (2*dx + d*(1 + gamma0(d,h,a)))*csc(a)
    t19 <- function(h,a) 2*dx*sec(a)

    times <- list(t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19)
}

# ========================================================================================
# ============================= HEAT MAP OF CROSSING TIMES ===============================

#png("CrossingTimes--Heatmap.png", width=1940, height=1250)
spar(mfrow=c(1,2))
{
    SV <- data.matrix(fread(lf("ST.dat")))
    d <- 1/2
    #mplot(SV[,1], SV[,2]*180/pi, pch=19, type='p', cex=0.2)
    mplot(0, 0, ylim=c(-90, 90), xlim=c(0,d), col='white', xaxs='i', yaxs='i', spa=FALSE,
	  ylab=tex("\\theta_0"), 
	  xlab=tex("h_0"),yaxt='n')
    axis(side=2, at=seq(-90,90,by=30))
    pal <- leaflet::colorBin("inferno", domain = NULL, bins=10000, reverse=TRUE)
    points(SV[,1], SV[,2], col=pal(SV[,3]), pch=15, cex=.8)

    pts <- sort(SV[,3])[10^seq(1,log10(nrow(SV)))]
    legend('topright', pch=15, cex=2, horiz = TRUE,
	   legend=round(pts, digits=4),
	   col=pal(pts))
}
# Write data with colors

X <- cbind(SV[,1], SV[,2], tpose(col2rgb(SV[,3])))
X <- cbind(SV[,1], SV[,2], SV[,3])
colnames(X) <- c("x", "y", "R", "G", "B")
fwrite(file="CTHeat.tkz", sep=" ", x=X)
fwrite(file="pal.pal", sep=" ", x=tpose(col2rgb(SV[,3])), col.names=FALSE)
dev.off()

# Example 
plot(NA, NA, asp=1, xlim=c(0, 1), ylim=c(0,pi))
x <- seq(0.01, pi-0.01, length=500)
y <- c()
for (i in seq(length(x))) 
    y <- c(y, rep(x[i], length=length(x)))
X <- rep(x, length(x))
Y <- y
Z <- X^2 + Y^2

pal <- leaflet::colorBin("inferno", domain = NULL, bins=1000, reverse=TRUE)
points(X, Y, col=pal(Z), pch=15, cex=.8)

fwrite(file="~/Dropbox/PHD/Thesis/Tikz/example.tkz", sep=" ", x=cbind(X, Y, Z), col.names=FALSE)
fwrite(file="example.tkz", sep=" ", x=cbind(X, Y, Z), col.names=FALSE)

plot(X, Y, pch=19, cex=0.5, asp=1)
y <- seq(-1, 1, length=100)
plot(0,0,xlim=c(-1,1),ylim=c(-1,1))
for (i in seq(length(x))) {
    for (j in seq(length(x))) {
	points(x[i], x[j], )
    }
}

g <- sapply(x, function(x) y)
gp <- grid::grid.points(x=x, y=y)
f <- function(x, y) { x/y } 

plot(g[,1], g[1,], pch=19, cex=0.1, ylim=c(-1, 1), xlim=c(-1,1))
sapply(2:length(x), function(x) points(g[,x], g[x,], pch=19, cex=0.1, add=TRUE))
pal <- leaflet::colorBin("inferno", domain = NULL, bins=1000, reverse=TRUE)
points(SV[,1], SV[,2], col=pal(SV[,3]), pch=15, cex=.8)


SV <- data.matrix(fread(lf("ST.dat")))
rgl::plot3d(x=SV[,1], y=SV[,2], z=SV[,4], zlim=c(0,10))

# ========================================================================================
# Visualising the changes in region (sort of directed graph)
#pdf("RegionVisits.pdf", width=18, height=18)
spar(mfrow=c(2,2))
{
    SV <- data.matrix(fread(lf("ST.dat")))
    theta <- seq(0, 2*pi, length=20)
    a <- 1.5

    mplot(a*cos(theta), a*sin(theta), asp=1, type='p', xaxt='n', yaxt='n', ylab="", xlab="",
	  ylim=c(-1.6, 1.6), spa=FALSE)
    text(x=(a+.15)*cos(theta), y=(a+0.15)*sin(theta), labels=seq(0,18), cex=2)
    lines(a*cos(theta[SV[-1,5]+1]), a*sin(theta[SV[-1,5]+1]), lwd=2.5)

    bid <- SV[-1,5]
    X <- a*cos(theta[bid+1])
    Y <- a*sin(theta[bid+1])
    n <- length(X)
    myArrow <- function(x0, y0, x1, y1, cut = 1, ...){
	x.new <- (1 - cut) * x0 + cut * x1
	y.new <- (1 - cut) * y0 + cut * y1
	# segments(x0, y0, x1, y1, ...)
	arrows(x0=x0, y0=y0, x1=x.new, y1=y.new, ...)
    }
    myArrow(x0=X[1:(n-1)], x1=X[2:n], y0=Y[1:(n-1)], y1=Y[2:n], 
	    cut=0.5, angle=45, length=0.2)
    points(a*cos(theta[SV[2,5]+1]), a*sin(theta[SV[2,5]+1]), col='green', pch=19, cex=2.5)
    points(a*cos(theta[SV[nrow(SV),5]+1]), a*sin(theta[SV[nrow(SV),5]+1]), col='red', pch=19, cex=2.5)

}
dev.off()
# ========================================================================================

source("../TrajectoryCompare.R")
PlotTrajectory(1/2)

# ========================================================================================
#pdf("RegionVisitations.pdf", width=14, height=14)

HeatMap <- function(x, y, z, coll="RdYlBu") {
}

#png("RegionCounts--Heatmap.png", width=1940, height=1250)
spar(mfrow=c(1,2))
{
    SV <- data.matrix(fread(lf("--RGC.dat")))
    d <- 1
    mplot(0, 0, ylim=c(-90, 90), xlim=c(0,d), col='white', xaxs='i', yaxs='i', spa=FALSE,
	  ylab=tex("\\theta_0"), 
	  xlab=tex("h_0"),yaxt='n')
    axis(side=2, at=seq(-90,90,by=30))

    cols = colorRampPalette(c("yellow", "red"))(nUniq)
    pal <- leaflet::colorNumeric("RdYlBu", domain = NULL, reverse=TRUE)
    points(SV[,1], SV[,2]*180/pi, col=pal(SV[,3]), pch=15, cex=.4)
    
    pts <- sort(unique(SV[,3]))
    legend('topright', pch=15, cex=2, horiz = FALSE,
	   legend=round(pts, digits=4),
	   col=pal(pts))
}
dev.off()
# ========================================================================================
# ========================================================================================
#' @brief Fraction of particles propagating and not propagating over time

TID <- function(d) {
    if (d < 1) {
	c(1,2,3,4,5,6,7,12,18) 
    } else {
	c(4,5,6,12,18,19)
    }
}

B <- data.matrix(fread(lf("r1--BID.dat")))
whichTransmitting <- function(i) sapply(TID, function(x) which(x == (B[,i])))

# Number of transmitting particles at time n
nTransmitting <- sapply(2:ncol(B), function(i) sum(unlist(lapply(whichTransmitting(i), function(x) length(x)))))

#pdf("TransRatio.pdf", width=12, height=12)
{
    Y <- (nrow(B)-nTransmitting)/nTransmitting
    mplot(1:(ncol(B)-1), Y, ylim=c(1.2, 2),
	  xlab=tex("n"), ylab=tex("N_L / N_T"),
	  log='', type='o', col='black')
    lines(1:(ncol(B)-1), Y, type='o', col='blue')
    abline(h=mean(Y[60:100]), lty='dashed')
}
#dev.off()

# ========================================================================================
#' @brief More resolved displacement distributions

L <- data.matrix(fread(lf("L.dat")))
mL <- apply(L, 2, mean)
abs(rbind(t,t) - apply(L, 2, minmax))

R <- data.matrix(fread(lf("r1--V.dat")))
t <- c(36,60) * 1e3

L <- getTSD_hist(dx=abs(R/L), n=1000)
L <- getTSD_hist(dx=abs(R), n=1000)
mplot(L$x, L$y, type='s', log='y',) #)ylim=c(0.00005, 1), xlim=c(0,0.6), spa=FALSE)
mplot(L$x, L$y, type='s', log='y', add=TRUE, col='red')


Vbs <- list.cbind.na(cbind(fread("../VBAL--FINE/Vbal--FH-Peaks.tkz")[,1], fread("../VBAL--FINE/Vbal--IH-Peaks.tkz")[,1]))

# ============= LOADING HISTOGRAMS FROM CLUSTER =============
L <- data.matrix(fread("DX--FINE/a7-r0--Vhist.dat"))[,-c(1,4)]
N <- 1:(ncol(L)/2)
Lx <- L[,N]; Ly <- L[,-N]
#spar(mfrow=c(3,1), cx=1.5)
id <- which(is.infinite(log(L[,3])))
Ly <- data.matrix(Ly[-id,1]); Lx <- data.matrix(Lx[-id,1])
mplot(Lx, Ly, type='s', log='y', xlim=c(0.2,0.75),
#mplot(L[,N], L[,-N], type='s', log='y', xlim=c(0.2,0.75),
      xlab=tex("x_{n_t} / s_{n_t}"), ylab=tex("s_{n_t}^{-1} P(x_{n_t} / s_{n_t})"),
      lwd=1, spa=FALSE)

# Rescale according to a fixed time 
Lx <- sweep(L[,N], 2, t, '/')
Ly <- sweep(L[,-N], 2, t, '*')
mplot(Lx, Ly, type='s', log='y',
      xlab=tex("x_{n_t} / s_{n_t}"), ylab=tex("s_{n_t}^{-1} P(x_{n_t} / s_{n_t})"),
      lwd=1, spa=FALSE)

# Bounding curve above which we identify large fluctuations / fronts 
{
    id <- which(Lx[,2] > 0.075 & Lx[,2] < 0.1)
    x0 <- Lx[id[1],2]
    x1 <- Lx[id[length(id)],2]
    F0 <- Ly[id[1],2]
    F1 <- Ly[id[length(id)],2]
    m <- log(F1/F0) / (x1 - x0)
    if (max(Lx) < 0.8){
	curve(expr=exp(m*x + 1.3*F0), lty='dashed', from=0.001, to=1, add=TRUE)
	#curve(expr=exp(m*x + 1.65*F0), lty='dashed', from=0.001, to=1, add=TRUE)
    } else {
	curve(expr=exp(0.88*m*x + 0.9*F0), lty='dashed', from=0.001, to=1, add=TRUE)
    }
    #abline(v=Lx[minmax(id),3])
}

# Extract peak location algorithmically
{
    n <- which(Lx[,1] > 0.2)
    # FH: W = 5, LX > 0.2. 
    # IH: W = 7, LX > 0.3
    ids <- argmax(x=Lx[n,1], y=(log(Ly[n,1])), w=5)
    X <- Lx[n[1] + ids -1,1]
    Y <- (Ly)[n[1] + ids -1,1]
    FH.peak.loc <- cbind(X,Y)
    eps <-.05
    #arrows(x0=X, x1=X, y0=Y+0.05*max(log(Ly)), y1=Y+eps, length=0.1, code = 2, col='blue', lwd=2)

    # ----- FINITE AND INFINITE HORIZON PEAK LOCATIONS ----- 
    #
    if (max(Lx) < 0.8){
	idsToKeep <- c(5,8,12,19,20,22,25,33,34,35,37,39)
    }else{
	idsToKeep <- c(2, 3, 4, 5, 7, 8, 9, 17, 18, 19, 20, 23)
    }

    # Finding the peaks
    Peaks <- Lx[n[1] + ids -1,1]
    sapply(1:length(Peaks), function(x) ArrowAtPeak(Peaks[x], lx=Lx[,1], ly=(Ly[,1]), eps=0.5))

    # Exact locations of the peaks
    Peaks.x <- Lx[n[1] + ids[idsToKeep] - 1,1]
    Peaks.y <- Ly[n[1] + ids[idsToKeep] - 1,1]
    sapply(1:length(Peaks.x), function(x) ArrowAtPeak(Peaks.x[x], lx=Lx[,1], ly=(Ly[,1]), eps=0.5))
    fwrite(x=cbind(Peaks.x, Peaks.y), file="Vbal_FH_PeakLocs.tkz", col.names=FALSE, sep=" ")

    # Defining the neighbourhoods of peaks (for finding trajectories of a particular speed)
    Peaks <- Lx[n[1] + ids[idsToKeep] - 2,1]
    Peaks.var <- Lx[n[1] + ids[idsToKeep]+1,1]
    abline(v=c(Peaks, Peaks.var))

    Peaks.LB <- Lx[n[1] + ids[idsToKeep] - 3,1]
    Peaks.UB <- Lx[n[1] + ids[idsToKeep] + 2,1]
    abline(v=c(Peaks.LB, Peaks.UB))
    fwrite(x=fm(cbind(Peaks.LB, Peaks.UB),6), file="Vbal_IH.tkz", col.names=FALSE, sep=" ")

    # ----- INFINITE HORIZON d = 1 PEAK LOCATIONS ----- 
}

# ========================================================================================

# ===== Plot the initial conditions =====
IC <- data.matrix(fread(lf("ST.dat")))
spar(mfrow=c(1,2))
IC <- cbind(tpose(fread(lf("BranchHeights.dat"))), tpose(fread(lf("BranchAngles.dat"))))
mplot(IC[,1], IC[,2]*180/pi, type='p', pch=19, cex=0.1, xlim=c(0, 1/2), ylim=c(-90, 90), spa=FALSE)
SV <- data.matrix(fread(lf("r0--ST.dat")))
mplot(SV[,1], SV[,2]*180/pi, cex=0.1, col='red', type='p', xlim=c(0, 1), ylim=c(-90, 90), add=TRUE, spa=FALSE)
abline(h=0)

# ====== Load jump distances and jump times ======
d <- 1/2
xi <- ReadDataList(lf("p0-a7-r0--ST"), idx=3)
tau <- ReadDataList(lf("p0-a7-r0--ST"), idx=4)
ST <- data.matrix(fread(lf("--ST")))

nubar <- apply(abs(xi/tau), 2, mean)
mean(nubar)
mean(abs(xi/tau))

#pdf("Vbal--FHIH.pdf", width=14, height=18)
spar(mfrow=c(2,1), cx=1.5)
{
    L <- getTSD_hist(dx=abs(xi/tau), n=1000)#, normalise.by.range=TRUE)
    { 
	Ly <- sweep(L$y, 2, apply(L$y,2,max,na.rm=TRUE), '/') # Normalise each hist by their height
	Ly <- L$y
	sids <- c(3,4,5,6,7) 
	sids <- c(10)
	mplot(L$x[,-sids], Ly[,-sids], type='s', ylim=c(0,10),
	      xaxt='n',  xlim=c(0,1),
	      xlab=tex("\\xi/\\tau"), ylab="Density", spa=FALSE)
	xt <- seq(0, 1, by=.1)
	axis(side=1, at=xt, labels=xt)
	#legend('topleft', legend=1:ncol(Ly), col=rbow(), lty=LTY(rbow()), cex=2, lwd=2)
	#abline(v=fread(lf("VB_FH_PeakLocs"))$V1)
    }
    dev.off()


    # ===== Trajectories which were found to have a particular ballistic speed =====

    STn <- data.matrix(fread(lf("r0--VB5.dat")))
    SCD <- SplitColumnData(STn, 50000)
    maxTime <- 60000 # Truncate data if time exceeds this
    t <- SCD$ByObservation[[4]]

    # Assume time scales are roughly the same - truncate based on the smallest minIDs
    minIDs <- apply(t, 2, function(x) min(which(x > maxTime)))
    OBS <- lapply(SCD$ByObservation, function(x) x[1:min(minIDs),])

    t <- OBS[[4]]
    vb <- OBS[[3]] / t
    mid <- mean(abs(vb[,1]),na.rm=TRUE)
    eps <- 0.03
    mplot(t, abs(vb[]), type='l', ylim=c(mid-eps, mid+eps))


    # Writing the unique regions visited by several trajectories having a particular speed  
    files <- lf("--VB")
    unids <- list()
    for (i in seq(length(files))) {
	STn <- data.matrix(fread(lf(files[i])))
	SCD <- SplitColumnData(STn, 50000)
	maxTime <- 60000 # Truncate data if time exceeds this
	t <- SCD$ByObservation[[4]]
	minIDs <- apply(t, 2, function(x) min(which(x > maxTime)))
	OBS <- lapply(SCD$ByObservation, function(x) x[1:min(minIDs),])

	unid <- sort(unique(as.vector(unlist(apply(OBS[[5]][ltk,], 2, unique)))))
	unids[[i]] <- unid
	#cat(file="UniqueRegions--FH.dat", append = TRUE, format(FH.peaks[i],digits=6), sort(unid), "\n")
    }
    
    ltk <- seq(nrow(OBS[[5]]) - 1e3, nrow(OBS[[5]]))
    for (i in seq(ncol(OBS[[5]]))) {
	print(table(OBS[[5]][ltk,i]))
    }
    sapply(1:ncol(OBS[[5]]), function(x) table(OBS[[5]][,x]))
    apply(OBS[[5]], 2, table)
    apply(OBS[[5]], 2, unique)
    which(OBS[[5]][,3] == 13)


    mplot(STn[1:maxTimes[1],4], vb[1:maxTimes[1]], type='l', col=rbow()[1])
    sapply(2:length(maxTimes), function(x) lines(STn[maxTimes[


    hist(STn[,5])$counts
    L <- getTSD_hist(dx=abs(vb), n=1000)#, normalise.by.range=TRUE)
    mplot(L$x[,], L$y[,], type='s', log='y', xlim=c(0.3,0.55)) #ylim=c(0,50),

}
dev.off()


# ====== Distribution of escape times (escaping from a particular region) ======

#pdf("ptau--FHIH.pdf", width=20, height=12)
spar(mfrow=c(1,2))
{
    L <- getTSD_ecdf(dx=abs(tau), n=1000, domain=exp.domain)
    mplot(L$x[,], 1-L$y[,], type='s', spa=FALSE,
	  log='xy', xlab=tex("\\tau"), ylab=tex("p_j(\\tau)"))
    legend('topright', legend=paste0("Region ", TID(d)), col=rbow(), lty=LTY(rbow()), lwd=2, cex=2)
    LinFit(mat=1-L$y[,], t=L$x[,], range=550:750)
    LinFit(mat=1-L$y[,], t=L$x[,], range=350:450)
    curve(1/x, from=0.01, to=1e7, lty='dotted', add=TRUE)
}
dev.off()

# ======= Distribution of ballistic speeds for several different regions =======
d <- 1
xi <- ReadDataList(lf("r1--ST"), idx=3)
tau <- ReadDataList(lf("r1--ST"), idx=4)
pdf("Vbal--IH.pdf", width=18, height=14)
{
    spar(mfrow=c(2,3), cx=1.75)
    L <- getTSD_hist(dx=xi/tau, n=3000)
    for (i in seq(1,ncol(xi))) {
	mplot(L$x[,i], L$y[,i], type='s', xlab=tex("\\xi/\\tau"), ylab="Density", spa=FALSE)
	text(x=0, y=0.9*max(L$y[,i]), labels=paste0("Region ", TID(d)[i]), cex=2, col='blue')
    }
}
dev.off()

# ======= Heat map of the ballistic speeds (a bit meaningless) =======
#png("Vbal--Heatmap.png", width=1940, height=1050)
spar(mfrow=c(1,2))
{
    SV <- data.matrix(fread(lf("ST.dat")))
    d <- 1
    mplot(0, 0, ylim=c(-90, 90), xlim=c(0,d), col='white', xaxs='i', yaxs='i', spa=FALSE,
	  ylab=tex("\\theta_0"), 
	  xlab=tex("h_0"),yaxt='n')
    axis(side=2, at=seq(-90,90,by=30))
    vb <- SV[,3]/SV[,4]
    id <- which(is.nan(vb))
    if (length(id) > 0) {
	SV <- SV[-id,]
	vb <- SV[,3]/SV[,4]
    }
    pal <- leaflet::colorBin("plasma", domain = NULL, bins=10000, reverse=TRUE)
    pal <- leaflet::colorNumeric("inferno", domain = NULL, reverse=TRUE, na.color = 'gray')
    points(SV[,1], SV[,2]*180/pi, col=pal(vb), pch=15, cex=1)
    pts <- c(min(vb), 0.5, max(vb))
    legend('topright', pch=15, cex=2, horiz = TRUE,
	   legend=round(pts, digits=4),
	   col=pal(pts))
}
dev.off()

# ======= Initial conditions of the ballistic fronts =======
Peaks <- ReadDataList(lf("PeakLocs"), idx=1)

foundPeaks.FH <- c(3,4,5,6,7,8)+1
foundPeaks.IH <- c(3,6,7,8)+1
Peaks.FH <- Peaks[foundPeaks.FH,1]
Peaks.IH <- Peaks[foundPeaks.IH,2]

r <- 0
d <- list(1/2, 1)[[r+1]]
K <- list(foundPeaks.FH, foundPeaks.IH)[[r+1]] 
h0 <- ReadDataList(lf(paste0("r", r)), idx=1)[,K]
t0 <- ReadDataList(lf(paste0("r", r)), idx=2)[,K]
vb <- ReadDataList(lf(paste0("r", r)), idx=5)[,K]


## 25/03/2024 ##
# h_0 | theta_0 | branch id | speed at t_0 = 32k | speed at t_1 = 64k
dat <- data.matrix(fread(lf(paste0("VB-FH.dat"))))

table(dat[,3])
table(dat[which(dat[,2] > -0.12 & dat[,2] < 0.12),3])

id <- which(dat[,2] > pi/2)
dat[id,2] <- pi - dat[id,2]
dat <- dat[-id,]


# Initial conditions causing the sub-leading front v ~ 1/2 
summary(dat[,2] * 180/pi)
plot(dat[,1], dat[,2]*180/pi, pch=19, type='p', cex=0.2, ylim=c(-90, 270), col=dcols[dat[,3]+3], ylab="theta", xlab="h")
abline(h=c(-90,90))
plot(dat[,1], dat[,2]*180/pi, pch=19, type='p', cex=0.3, ylim=c(269.99, 270), col='red')
plot(dat[,1], dat[,2]*180/pi, pch=19, type='p', cex=0.3, ylim=c(89.95, 90.001), col='red')
plot(dat[,1], dat[,2]*180/pi, pch=19, type='p', cex=0.1, ylim=c(-89.995, -89.88), col='red', xlab="h", ylab="theta")
plot(dat[,1], dat[,2]*180/pi, pch=19, type='p', cex=0.3, ylim=c(-6, 6), col=dcols[dat[,3]+3])
abline(v=c(0.179, 0.25))

plot(dat[,1], dat[,2]*180/pi, pch=19, type='p', cex=0.3, ylim=c(-4.09, -4.06), col=dcols[dat[,3]+3])
plot(dat[,1], dat[,2]*180/pi, pch=19, type='p', cex=0.3, ylim=c(4.06, 4.09), col=dcols[dat[,3]+3])
matplot(dat[,1], -dat[,2]*180/pi, pch=19, type='p', cex=0.3, ylim=c(4.06, 4.09), col='red', add=TRUE)


{
    # Examining the initial conditions near -pi/2
    dat <- data.matrix(fread(lf(paste0("tst.dat"))))
    # Remove some points right near -pi/2
    a <-  -1.57074
    abline(h=a)
    fil <- which(dat[,2] > a)
    dat <- dat[fil,]

    # 
    dim(dat)
    fil <- which(dat[,2] < -1.5706)
    sam <- sample(fil, size=120000)
    dat <- dat[-sam,]
    plot(dat[,1], dat[,2], pch=20, type='p', cex=0.1, , col='red', xlab="h", ylab="theta")
    abline(v=c(0, 1/3, 2/5, 3/7,4/9, 5/11), lty='dashed', lwd=0.7)
}

# With region 6 as a boundary
d <- 1/2
plot(dat[,1], dat[,2], pch=19, type='p', cex=0.5,)
plot(dat[,1], dat[,2], pch=19, type='p', cex=0.2, col=dcols[dat[,3]+3], ylab="theta", xlab="h")


# Plot region lines for extra clarity
plot(dat[,1], dat[,2], pch=19, type='p', cex=0.3,ylim=c(-0.1, 0.1), col=dcols[dat[,3]+3], ylab="theta", xlab="h")
mean((dat[which(dat[,2] < -0.065 & dat[,2] > -0.08),2]))
mean((dat[which(dat[,2] > -0.065 & dat[,2] < -0.052),2]))
minmax(dat[which(dat[,2] < -0.065 & dat[,2] > -0.08 & dat[,1] >0.15 & dat[,1] < 0.25),1])
minmax(dat[which(dat[,2] < -0.065 & dat[,2] > -0.08 & dat[,1] >0.25 & dat[,1] < 0.41),1])
abline(v=c(0.1793, 0.24677), lty='dashed')
abline(v=c(0.03515, 0.10386), lty='dashed')
abline(v=c(0.32256, 0.3923), lty='dashed')
abline(h=c(0.0713, -0.0713), lty='dashed')
abline(h=c(0.0555, -0.0555), lty='dashed')
abline(v=0.15)

abline(h=c(-0.065, -0.052), lty='dashed')
# theta = asec(2m/7) = acos(7/2m)
abline(h=-0.0713)
abline(h=acos(7))

abline(v=c(0.5-0.03515, 0.5-0.10386), lty='dashed')

h0 <- .178
ts <- -0.0713
h0 - 2*tan(ts)
H <- function(h) { h + 2*tan(ts) }

# First iteration of endpoints [0.178, 0.2467]
points(x=H(0.178), y=ts)
points(x=H(0.2467), y=ts)
c(H(0.178), H(0.2467))
points((H(0.178) + H(0.2467) ) / 2, y=ts)
m <- (H(0.178) + H(0.2467) ) / 2
H(H(h0))

# Second iteration of endpoints [0.178, 0.2467]
H <- function(h) { h + 0.5* (1 + 5*tan(ts)) }
points(x=H(0.035157), y=ts, col='red')
points(x=H(m), y=ts, col='red')

H <- function(h) { h + 0.5* (3*tan(ts) + 1) }
points(x=H(0.10385), y=ts, col='red')
points(x=H(m), y=ts, col='red')



plot(dat[,1], dat[,2], pch=19, type='p', cex=0.3,ylim=c(-0.6, 0.6), col=dcols[dat[,3]+3], ylab="theta", xlab="h")
legend(col=unique(dcols[c(3,4,6,7,12)+3]), x='topleft', legend=paste0("Region ", c(3,4,6,7,12)), pch=19)
d <- 1
{
    # Region 3
    curve(-atan((x)/(d)), from=0, to=1/3, add=TRUE) 
    # Region 4
    curve(-atan(x), from=0, to=0.5, add=TRUE) 
    # Region 5
    curve(-atan(x/(d+1)), from=0, to=0.5, add=TRUE) 

    #curve(-atan(x/(d+2)), from=0, to=0.5, add=TRUE) 
    #curve(-atan(x/(d+2.5)), from=0, to=0.5, add=TRUE) 
    #curve(-atan(x/(d+3)), from=0, to=0.5, add=TRUE) 
    #curve(-atan(x/(d+4)), from=0, to=0.5, add=TRUE) 
    #curve(-atan(x/(d+5)), from=0, to=0.5, add=TRUE) 

    # Region 6
    curve(-atan(x/(2*d+1)), from=0, to=0.5, add=TRUE) 
    curve(atan((d-x)/(2*d+1)), from=0, to=0.5, add=TRUE) 
    # Region 7
    curve(atan((0.5-x)/(2*d+0.5)), from=0, to=0.5, add=TRUE) 
    # Region 8
    curve(atan((0.5-x)/(d+0.5)), from=0, to=0.5, add=TRUE) 
    curve(atan((1-x)/(2*d+1)), from=0, to=0.5, add=TRUE) 
    # Region 9
}

abline(v=0.2)
abline(h=-atan(0.2))
abline(h=-atan(0.66*0.2))
-atan(0.66*0.2) * 180/pi
abline(v=0.4)
h <- 0.2
h1 <- 0.5 * (1 - 0.2)
1.5 * sec(-atan(0.2))

# Distribution of ballistic speeds 
# | speed at t_0 = 32k | speed at t_1 = 64k
V <-  data.matrix(fread(lf("--V.dat")))
t <- c(36,60) * 1e3
L <- getTSD_hist(dx=(V), n=500)
mplot(L$x, L$y, type='s', log='y',) #)ylim=c(0.00005, 1), xlim=c(0,0.6), spa=FALSE)

## 25/03/2024 ##

## 05/04/2024 ##
# Analysing the initial condition data that is used in the thesis and paper
dat <- data.matrix(fread(lf(paste0("FHh0t0.tkz"))))
dat <- data.matrix(fread(lf(paste0("IHh0t0.tkz"))))
dat[,2] <- dat[,2] * pi/180.0
n <- -which(dat[,2] > -0.0001 & dat[,2] < 0.0001)
plot(dat[,1], dat[,2], pch=19, type='p', cex=0.2, ylim=c(-.1, .1), ylab="theta", xlab="h", yaxt="n")
plot(dat[n,1], dat[n,2], pch=19, type='p', cex=0.2, ylim=c(-.1, .1), ylab="theta", xlab="h", yaxt="n")
axis(side=2, labels=seq(-0.1, 0.1, by=0.01), at=seq(-0.1, 0.1, by=0.01))


### Prepping and rewriting the old data ###
dat <- data.matrix(fread(lf(paste0("FHh0t0.tkz"))))
dat[,2] <- dat[,2] * pi/180 # Radians
n <- -which(dat[,2] > -0.0001 & dat[,2] < 0.0001) # Prune
dat <- dat[n,]

dat <- data.matrix(fread(lf(paste0("IHh0t0.tkz"))))
dat[,2] <- dat[,2] * pi/180 # Radians
plot(dat[,1], dat[,2], pch=19, type='p', cex=0.2, ylim=c(-.1, .1), ylab="theta", xlab="h", yaxt="n")

# Assign a number to distinguish the different regions for coloring FINITE HORIZONT
colcol <- rep(-1, nrow(dat))
for (i in seq(nrow(dat))) {
    h <- dat[i,1]
    a <- dat[i,2]
    if (a < -atan(h/(d+0.5)) & a > -atan(h/d)) { # Region 1
	colcol[i] <- 3
    } else if (a < -atan(h/(d+1)) & a > -atan(h) ) { # R4
	colcol[i] <- 4
    } else if (a < atan(d-h) & a > atan((d-h)/(2*d+1))) {
	colcol[i] <- 7
    } else if (a > atan((d-h)/d) & a < atan((d+0.5-h)/(d+0.5))) {
	colcol[i] <- 12
    } else if (a > atan((0.5-h)/(0.5)) & a < atan(d - h)) {
	colcol[i] <- 19
    } else {
	colcol[i] <- 6
    }
}

dat <- data.matrix(fread(lf(paste0("IHh0t0.tkz"))))
dat[,2] <- dat[,2] * pi/180 # Radians
plot(dat[,1], dat[,2], pch=19, type='p', cex=0.2, ylim=c(-.1, .1), ylab="theta", xlab="h", yaxt="n")

# Assign a number to distinguish the different regions for coloring INFINITE HORIZON
colcol <- rep(-1, nrow(dat))
for (i in seq(nrow(dat))) {
    h <- dat[i,1]
    a <- dat[i,2]
    if (a > -atan(h) & a < -atan(h/(1 + d))) {
	colcol[i] <- 4
    } else if (a > atan((0.5-h)/(0.5)) & a < atan(d - h)) {
	colcol[i] <- 19
    } else if (a > -atan(h/(1 + 2*d)) & a < atan((0.5-h)/(1 + 2*d - 0.5))) {
	colcol[i] <- 6
    } else {
	colcol[i] <- 12
    }
}
new.dat <- cbind(dat, colcol)
# Check
plot(new.dat[,1], new.dat[,2], pch=19, type='p', cex=0.2, ylim=c(-.1, .1), ylab="theta", xlab="h", yaxt="n", col=dcols[new.dat[,3]])

# Write
fwrite(x=new.dat, file="VBIC--redo/IHh0t0.tkz", sep=" ", col.names=FALSE)

###########################################



horizontals <- c(mean((dat[which(dat[,2] < -0.096 & dat[,2] > -0.101),2])),
		 mean((dat[which(dat[,2] < -0.065 & dat[,2] > -0.080),2])),
                 mean((dat[which(dat[,2] > -0.065 & dat[,2] < -0.052),2])),
                 mean((dat[which(dat[,2] > -0.051 & dat[,2] < -0.044),2])),
                 mean((dat[which(dat[,2] > -0.051+0.008 & dat[,2] < -0.044+0.008),2])),
                 mean((dat[which(dat[,2] > -0.051+2*0.008 & dat[,2] < -0.048+2*0.008),2])),
                 mean((dat[which(dat[,2] > -0.051+2.5*0.008 & dat[,2] < -0.048+2.5*0.008),2])),
                 mean((dat[which(dat[,2] > -0.051+3*0.008 & dat[,2] < -0.049+3*0.008),2])))

diff(horizontals[2:(length(horizontals))] / horizontals[1:(length(horizontals)-1)])

abline(h=c(horizontals,-horizontals))
abline(h=asin(horizontals))
bounds <- c(0.498558, 0.501723)
(bounds[1] + bounds[2])/2
mean(bounds)

abline(h=c(-0.096, -0.101), col=c('red', 'blue'))
abline(h=c(-0.051, -0.048)+2.5*0.008, col=c('red', 'blue'))
abline(h=c(-0.051, -0.049)+3*0.008, col=c('red', 'blue'))
## 05/04/2024 ##


######### 05/04/2024 ############

#pdf("FrontIC--FH-h0t0.pdf", height=12, width=16)
spar(mfrow=c(2,3))
for (i in seq(ncol(h0))) {
#    mplot(vb[,i], h0[,i], ylab=tex("h_0"), xlab=tex("\\nu"), ylim=c(0,d),
#	  main=tex(paste0("$\\nu = ", format(Peaks.FH[i],digits=5), "$")),
#	  spa=FALSE, pch='.', type='p', cex=.5)
#    mplot(vb[,i], t0[,i]*180/pi, yaxt='n', ylab=tex("\\theta_0"), xlab=tex("\\nu"), ylim=c(-90, 270),
    mplot(h0[,i], t0[,i]*180/pi,  ylab=tex("\\theta_0"), xlab=tex("h_0"), 
	  ylim=c(-7, 7), #yaxt='n',
	  main=tex(paste0("$\\nu = ", format(Peaks.FH[i],digits=5), "$")),
	  spa=FALSE, pch='.', type='p', cex=2)
     axis(side=2, at=seq(-90, 270, by=30))
}
dev.off()



#pdf("FrontIC--IH-05zoom.pdf", height=10, width=18)
spar(mfrow=c(1,2))
i <- 6
mplot(vb[,i], t0[,i]*180/pi,  ylab=tex("\\theta_0"), xlab=tex("h_0"), xaxs='i', yaxs='i', 
mplot(h0[,i], t0[,i]*180/pi,  ylab=tex("\\theta_0"), xlab=tex("h_0"), xaxs='i', yaxs='i', 
      ylim=c(-6, 6), 
      main=tex(paste0("$\\nu = ", format(Peaks.IH[i],digits=5), "$")),
      spa=FALSE, pch='.', type='p', cex=2)

# Tikz data
id <- which((t0[,i]*180/pi > -6 & t0[,i]*180/pi < 6))# & (t0[,i]*180/pi > 0.1 & t0[,i]*180/pi < 6))
id <- id[seq(1, length(id), by=1)]
plot(h0[id,i] ,t0[id,i]*180/pi)

fwrite(x=cbind(h0[id,i], t0[id,i]*180/pi), file="FHh0t0.tkz", sep=" ", col.names=FALSE)


source("../../RegionLines-FH.R")
DrawRegionLines(d=d,h=seq(0,d,length=1000))
abline(v=c(0,d))
source("../../RegionLines-IH.R")
DrawRegionLines(d=d,h=seq(0,d,length=1000))
dev.off()

#' @brief Draws a polygon between two curves
#' @f,g Vector satisfying f < g
#' @x Arugment to @f,g
#' @... Arguments to ?polygon
ShadeBetween <- function(x1, x2, f1, f2, ...) {
    X <- c(x1, rev(x2))
    Y <- c(f1, rev(f2))
	polygon(X, Y, ...)
}

# Consts
#dcols <- rev(dcols)
oh <- dy <- dx <- 1/2 
dx2 <- 2*dx
rtd <- 180/pi
d <- 1/2
H <- seq(0,d, length=1000)
p.lwd <- 1
l.lwd <- 1.5

# PLOT SEQUENCE
#pdf("Regions-FH-270.pdf", width=16, height=20)
source("../../RegionLines-FH.R")

mplot(H, a1(d, H), ylim=c(-90, 270), xlim=c(0,d), col='white', xaxs='i', yaxs='i', #spa=FALSE,
      ylab=tex("\\theta_0"), 
      xlab=tex("h_0"),yaxt='n')
axis(side=2, at=seq(-90,270,by=90))

dcols <- rep('orange', length=100)
# -90 to 90
{
    # Branch 0
    ShadeBetween(H, H, a0(d,H), a1(d,H),  col=dcols[1], lwd=p.lwd,)

    # Branches 1, 2, 3
    p0id <- PHI0.id(d,H); p1id <- PHI1.id(d,H)
    ShadeBetween(H, H, a1(d,H), a2(d,H), col=dcols[3], lwd=p.lwd,)
    ShadeBetween(H[p0id], H[p0id], a1(d,H[p0id]),   phi0(d,H)[p0id], col=dcols[2], lwd=p.lwd, )
    ShadeBetween(H[p0id], H[p1id], phi0(d,H)[p0id], phi1(d,H)[p1id], col=dcols[3], lwd=p.lwd,)
    ShadeBetween(H[p1id], H, phi1(d,H)[p1id], a2(d,H),               col=dcols[4], lwd=p.lwd,)
    ShadeBetween(H, H, a2(d,H), a3(d,H), col=dcols[2], lwd=p.lwd,) 
    # Branches 3 - 8
    A1 <- list(a3, a4, a5, a6, a7, a8)
    for (i in seq(length(A1)-1)) 
	ShadeBetween(H, H, A1[[i]](d,H), A1[[i+1]](d,H), col=dcols[i+4], , lwd=p.lwd,) 

    # Branches 9, 10, 11
    p2id <- PHI2.id(d,H); p3id <- PHI3.id(d,H)
    ShadeBetween(H,       H, a7(d,H), a8(d,H),  col=dcols[11],  , lwd=p.lwd,)
    ShadeBetween(H,       H[p2id], a7(d,H), phi2(d,H)[p2id],  col=dcols[10], lwd=p.lwd,)
    ShadeBetween(H[p2id], H[p3id], phi2(d,H)[p2id], phi3(d,H)[p3id], col=dcols[11], lwd=p.lwd,)
    ShadeBetween(H[p3id], H[p3id], phi3(d,H)[p3id], a8(d,H[p3id]), col=dcols[12], lwd=p.lwd, )

    A2 <- list(a12, a13, a14, a15, a16, a17, a18, a19, a20)
    for (i in seq(length(A2)-1)) 
	ShadeBetween(H, H, A2[[i]](d,H), A2[[i+1]](d,H), col=dcols[i+14], lwd=p.lwd,) 

    # Lines separating each region
    HP <- seq(0,d, length=10000)
    sapply(1:100, function(x) 
	   lines(HP[PROP1.id(d,HP,x)], PROP1(d,HP,x)[PROP1.id(d,HP,x)], lty='dashed', lwd=1))
    lines(H[p0id], phi0(d,H)[p0id], lwd=l.lwd)
    lines(H[p1id], phi1(d,H)[p1id], lwd=l.lwd)
    lines(H[p2id], phi2(d,H)[p2id], lwd=l.lwd)
    lines(H[p3id], phi3(d,H)[p3id], lwd=l.lwd)
    lines(H, a1(d,H), lwd=l.lwd)
    lines(H, a2(d,H), lwd=l.lwd)

    sapply(1:100, function(x) 
	   lines(HP[PROP2.id(d,HP,x)], PROP2(d,HP,x)[PROP2.id(d,HP,x)], lty='dashed', lwd=1))
    sapply(1:2, function(x)   
	   lines(HP[PROP3.id(d,HP,x)], PROP3(d,HP,x)[PROP3.id(d,HP,x)], lty='dashed', lwd=1))

    for (i in seq(length(A1))) {
	lines(H, A1[[i]](d,H), lwd=l.lwd)
    }
    for (i in seq(length(A2))) {
	lines(H, A2[[i]](d,H), lwd=l.lwd)
    }
}

# -90 to 270
{
    # Branch 0
    ShadeBetween(H, H, 180-a1(d,H), 180-a0(d,H),  col=dcols[1], lwd=p.lwd,)

    # Branches 1, 2, 3
    p0id <- PHI0.id(d,H); p1id <- PHI1.id(d,H)
    ShadeBetween(H, H, 180-a2(d,H), 180-a1(d,H), col=dcols[3], lwd=p.lwd,)
    ShadeBetween(H[p0id], H[p0id], 180-phi0(d,H[p0id]),   180-a1(d,H)[p0id], col=dcols[2], lwd=p.lwd, )
    ShadeBetween(H[p1id], H[p0id], 180-phi1(d,H)[p1id], 180-phi0(d,H)[p0id], col=dcols[3], lwd=p.lwd,)
    ShadeBetween(H, H[p1id], 180-a2(d,H), 180-phi1(d,H)[p1id], col=dcols[4], lwd=p.lwd,)
    ShadeBetween(H, H, 180-a3(d,H), 180-a2(d,H), col=dcols[2], lwd=p.lwd,) 
    # Branches 3 - 8
    A1 <- list(a3, a4, a5, a6, a7, a8)
    for (i in seq(length(A1)-1))
	ShadeBetween(H, H, 180-A1[[i+1]](d,H), 180-A1[[i]](d,H), col=dcols[i+4], , lwd=p.lwd,) 

    # Branches 9, 10, 11
    p2id <- PHI2.id(d,H); p3id <- PHI3.id(d,H)
    ShadeBetween(H,       H, 180-a8(d,H), 180-a7(d,H),  col=dcols[11],  , lwd=p.lwd,)
    ShadeBetween(H[p2id], H, 180-phi2(d,H)[p2id],180-a7(d,H), col=dcols[10], lwd=p.lwd,)
    ShadeBetween(H[p3id], H[p2id], 180-phi3(d,H)[p2id], 180-phi2(d,H)[p3id], col=dcols[11], lwd=p.lwd,)
    ShadeBetween(H[p3id], H[p3id], 180-a8(d,H[p3id]), 180-phi3(d,H)[p3id], col=dcols[12], lwd=p.lwd, )

    A2 <- list(a12, a13, a14, a15, a16, a17, a18, a19, a20)
    for (i in seq(length(A2)-1)) 
	ShadeBetween(H, H, 180-A2[[i+1]](d,H), 180-A2[[i]](d,H), col=dcols[i+14], lwd=p.lwd,) 

    # Lines separating each region
    HP <- seq(0,d, length=10000)
    sapply(1:100, function(x) 
	   lines(HP[PROP1.id(d,HP,x)], 180-PROP1(d,HP,x)[PROP1.id(d,HP,x)], lty='dashed', lwd=1))
    lines(H[p0id], 180-phi0(d,H)[p0id], lwd=l.lwd)
    lines(H[p1id], 180-phi1(d,H)[p1id], lwd=l.lwd)
    lines(H[p2id], 180-phi2(d,H)[p2id], lwd=l.lwd)
    lines(H[p3id], 180-phi3(d,H)[p3id], lwd=l.lwd)
    lines(H, 180-a1(d,H), lwd=l.lwd)
    lines(H, 180-a2(d,H), lwd=l.lwd)

    sapply(1:100, function(x) 
	   lines(HP[PROP2.id(d,HP,x)], 180-PROP2(d,HP,x)[PROP2.id(d,HP,x)], lty='dashed', lwd=1))
    sapply(1:2, function(x)   
	   lines(HP[PROP3.id(d,HP,x)], 180-PROP3(d,HP,x)[PROP3.id(d,HP,x)], lty='dashed', lwd=1))

    for (i in seq(length(A1))) {
	lines(H, 180-A1[[i]](d,H), lwd=l.lwd)
    }
    for (i in seq(length(A2))) {
	lines(H, 180-A2[[i]](d,H), lwd=l.lwd)
    }

    # Borders
    lines(c(0,d), c(-90,-90), lwd=l.lwd)
    lines(c(0,0), c(-90,270), lwd=l.lwd)
    lines(c(d,d), c(-90,270), lwd=l.lwd)
    lines(c(0,d), c(270,270), lwd=l.lwd)
}


