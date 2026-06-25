
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
H <- seq(0,d, length=1000)

# PLOT SEQUENCE
#pdf("Regions-Transmit.pdf", width=18, height=12)
source("RegionLines-FH.R")
p.lwd <- 0.1
l.lwd <- 2.5
spar(mfrow=c(1,2))
d <- 1/2
mplot(H, a1(d, H), ylim=c(-90, 90), xlim=c(0,d), col='white', xaxs='i', yaxs='i', spa=FALSE,
      ylab=tex("\\theta_0 \\in \\Pi_{j}  (h_0)"), 
      xlab=tex("h_0"),yaxt='n')
axis(side=2, at=seq(-90,90,by=30))

# Checking
SV <- data.matrix(fread(lf("SV.dat")))
points(SV[,1], SV[,2]*180/pi, col=ifelse(SV[,4]==0, 'blue', 'red'), pch=19, cex=1.5)
# BLUE = NULL MOVE, RED = TRANSMIT
{


    ShadeBetween(H, H, a0(d,H), a1(d,H),  col='blue',  , lwd=p.lwd,)
    ShadeBetween(H, H, a1(d,H), a6(d,H), col='red', lwd=p.lwd,) # 4
    ShadeBetween(H,       H, a6(d,H), a8(d,H),  col='blue',  , lwd=p.lwd,)
    ShadeBetween(H, H, a12(d,H), a13(d,H), col='red', lwd=p.lwd,) 
    ShadeBetween(H, H, a13(d,H), a19(d,H), col='blue', lwd=p.lwd,) 
    ShadeBetween(H, H, a19(d,H), a20(d,H), col='red', lwd=p.lwd,) 

    # Lines separating each region
    # Borders
    lines(c(0,d), c(-90,-90), lwd=l.lwd)
    lines(c(0,0), c(-90,90), lwd=l.lwd)
    lines(c(d,d), c(-90,90), lwd=l.lwd)
    lines(c(0,d), c(90,90), lwd=l.lwd)
}

d <- 1
H <- seq(0,d, length=1000)
source("RegionLines-IH.R")
mplot(H, a1(d, H), ylim=c(-90, 90), xlim=c(0,d), col='white', xaxs='i', yaxs='i', spa=FALSE,
      ylab=tex("\\theta_0 \\in \\Pi_{j}  (h_0)"), 
      xlab=tex("h_0"),yaxt='n')
axis(side=2, at=seq(-90,90,by=30))

{
    ShadeBetween(H, H, a6(d,H), a19(d,H),  col='blue', border='blue',lwd=p.lwd,)
    ShadeBetween(H, H, a0(d,H), a1(d,H),  col='blue', lwd=p.lwd,)
    ShadeBetween(H, H, a12(d,H), a13(d,H),col='red',  border='red',lwd=p.lwd,)
    lines(H, a1(d,H))
    lines(H, a6(d,H))
    ShadeBetween(H, H, a1(d,H), a6(d,H),  col='red',  border='red', lwd=p.lwd,)
    ShadeBetween(H, H, a21(d,H), a22(d,H),col='red',  border='red', lwd=p.lwd,)
    ShadeBetween(H, H, a19(d,H), a20(d,H),col='red',  lwd=p.lwd,)


    lines(H, a21(d,H), col='white')

    lines(H, a6(d,H), col='white')
    lines(H, a22(d,H))
    lines(H, a13(d,H), col='white')
    lines(H, a19(d,H), col='black')
    #lines(H, a22(d,H))
    #ShadeBetween(H, H, a6(d,H), a5(d,H),  col='blue', border='blue',lwd=p.lwd,)
    #ShadeBetween(H, H, a7(d,H), a19(d,H), col='blue', lwd=p.lwd,)

    # Borders
    lines(c(0,d), c(-90,-90), lwd=l.lwd)
    lines(c(0,0), c(-90,90), lwd=l.lwd)
    lines(c(d,d), c(-90,90), lwd=l.lwd)
    lines(c(0,d), c(90,90), lwd=l.lwd)
}
dev.off()



# FINITE HORIZON d = 1/2
At <- 0.30127 + 0.099160 + 0.19227
Al <- 0.56599 + 0.13756 + 0.27456
At/Al
Al/At



