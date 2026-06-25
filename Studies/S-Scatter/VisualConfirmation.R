args <- commandArgs(trailingOnly = TRUE)

source("~/Dropbox/Rlibs.R")

setwd("../Studies/S-Scatter/Data/")
file <- paste0("Visuals-", as.numeric(args[1]), ".pdf")
print(file)
pdf(file, width=18, height=18)
spar(mfrow = c(5,4))
spar(mfrow = c(4,4))
par(mar=c(0,0,0,0))
d <- dx <- dy <- 1/2
#for (i in seq(0,18)) {
for (i in c(0,4,5,6,8,9,10,11,12,13,14,15,16,17,18,19)) {
{
    i <- 7
    x <- data.matrix(fread(lff(paste0("p", i, "-"), "--x.dat")))
    y <- data.matrix(fread(lff(paste0("p", i, "-"), "--y.dat")))

    # Height/position from the billiard simulation trajectory 
    BilliardH <- tpose(fread(lff(paste0("p", i, "-"), "--RawH.dat")))
    BilliardX <- tpose(fread(lff(paste0("p", i, "-"), "--RawX.dat")))
    
    # Height and position from scattering map calculation
    SMH <- tpose(fread(lff(paste0("p", i, "-"), "--H.dat"))[,-1])
    SMX <- tpose(fread(lff(paste0("p", i, "-"), "--R.dat"))[,-1])
    SMA <- tpose(fread(lff(paste0("p", i, "-"), "--A.dat"))[,-1])
    SML <- tpose(fread(lff(paste0("p", i, "-"), "--L.dat")))
    
    #SMH <- tpose(fread(paste0("p", i, "-a7-r0--H.dat")))
    #print(tpose(fread(lff(paste0("p", i, "-"), "--A"))) * 180/pi)
    #print(tpose(fread(lff(paste0("p", i, "-"), "--H.dat")))) 
    
    PlotChannelBoundary(xVert=c(0, dx, 2*dx, 2*dx, dx, 0), yVert=c(0, dy, 0, d, dy+d,d), from=-1,to=2)
    plot(c(0,1), c(0,d+dy+.1), col='white', asp=1, xaxt='n', yaxt='n', xlab='', ylab='')
    polygon(x=c(0, dx, 2*dx, 2*dx, dx, 0), y=c(0, dy, 0, d, dy+d,d), lwd=2)

    # Truncate trajectory to show only the first passing
    sub <- 1:(min(which(x[,1] > 1 | x[,1] < 0)) -1)
    pm1 <- ifelse(x[sub[length(sub)]+1] < 0, 0, 1)
    x <- c(x[sub], pm1)
    y <- c(y[sub], SMH[1])
#
    # Show trajectory lengths
    bTime <- sum(sqrt(diff(x)^2 + diff(y)^2))
    #text(x=-0.05, y=1.05, labels = paste0("(", LETTERS[i+1], ")"), cex=2)
    #text(x=-0.05, y=d+dy+.05, labels = paste0('<B>', i+1, '</B>'), cex=2)
    text(x=-0.05, y=d+dy-0.05, labels = paste0("(", i, ")"), cex=2.5, font=5)
    text(x=1, y=d+dy-0.1, labels = paste0("BT = ", format(bTime, digits=6)), cex=1.5)
    text(x=1, y=d+dy-0.25, labels = paste0("CT = ", format(SML[2,1], digits=6)), cex=1.5)

    # Plot trajectory (billiard simulation)
    lines(x, y, lwd=2)
    idx <- 1:min(nrow(SMH), nrow(BilliardH))

    if (pm1 == 0) {
	curve(tan(SMA[1])*x + SMH[1], from=-0.1, to=pm1, add=TRUE, lty='dotted', lwd=2)
    } else {
	curve(tan(SMA[1])*(x-pm1) + SMH[1], from=pm1, to=1.1, add=TRUE, lty='dotted', lwd=2)
    }
    #points(BilliardX[idx,], BilliardH[idx,], cex=2, lwd=2)
    points(SMX[1,], SMH[1,], pch=19, col='red')
    points(BilliardX[1,], BilliardH[1,], pch=19, col='green')
}
}
dev.off()


{
    i <- 0
    d <- 1
    x <- data.matrix(fread(lff(paste0("p", i, "-"), "--x.dat")))
    y <- data.matrix(fread(lff(paste0("p", i, "-"), "--y.dat")))


    BilliardH <- tpose(fread(lff(paste0("p", i, "-"), "--RawH.dat")))
    BilliardX <- tpose(fread(lff(paste0("p", i, "-"), "--RawX.dat")))

    # Height and position from scattering map calculation
    SMH <- tpose(fread(lff(paste0("p", i, "-"), "--H.dat")))
    SMX <- tpose(fread(lff(paste0("p", i, "-"), "--R.dat")))
    SMA <- tpose(fread(lff(paste0("p", i, "-"), "--A.dat")))
    SML <- tpose(fread(lff(paste0("p", i, "-"), "--L.dat")))

    mplot(x[1:1000])
    mplot(SMX, col='red')
    lines(BilliardX)
    
PlotChannelBoundary(xVert=c(0, dx, 2*dx, 2*dx, dx, 0), yVert=c(0, dy, 0, d, dy+d,d), from=-2,to=2)
lines(x[1:15],y[1:15], lwd=2)
points(SMX[1:8], SMH[1:8], col='red', cex=2, pch=19)
points(BilliardX[1:8], BilliardH[1:8], pch=19, cex=1)
}

{
# Height/position from the billiard simulation trajectory 
idx <- 1:min(nrow(SMH), nrow(BilliardH)) 
idx <- 1:7
plot(BilliardH[idx,], SMH[idx,], type='l', )
plot(BilliardX[idx,], SMX[idx,], type='l', asp=1)
}
	curve(tan(SMA[5])*(x-pm1) + SMH[5], from=pm1, to=1.1, add=TRUE, lty='dotted', lwd=2)




