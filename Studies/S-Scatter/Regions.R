
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

# PLOT SEQUENCE
#pdf("Regions-FH.pdf", width=16, height=16)
source("RegionLines-FH.R")
#spar(mfrow=c(1,2))
mplot(H, a1(d, H), ylim=c(-90, 90), xlim=c(0,d), col='white', xaxs='i', yaxs='i', #spa=FALSE,
      ylab=tex("\\theta_0 \\in \\Pi_{j}  (h_0)"), 
      xlab=tex("h_0"),yaxt='n')
axis(side=2, at=seq(-90,90,by=30))

p.lwd <- 0.1
l.lwd <- 2.5
{
    # Branch 0
    ShadeBetween(H, H, a0(d,H), a1(d,H),  col=dcols[1],  , lwd=p.lwd,)

    # Branches 1, 2, 3
    p0id <- PHI0.id(d,H); p1id <- PHI1.id(d,H)
    ShadeBetween(H, H, a1(d,H), a2(d,H), col=dcols[3],  , lwd=p.lwd,)
    ShadeBetween(H[p0id], H[p0id], a1(d,H[p0id]),   phi0(d,H)[p0id], col=dcols[2],  , lwd=p.lwd, )
    ShadeBetween(H[p0id], H[p1id], phi0(d,H)[p0id], phi1(d,H)[p1id], col=dcols[3],  , lwd=p.lwd,)
    ShadeBetween(H[p1id], H, phi1(d,H)[p1id], a2(d,H),               col=dcols[4],  , lwd=p.lwd,)

    ShadeBetween(H, H, a2(d,H), a3(d,H), col=dcols[2], , lwd=p.lwd,) 
    # Branches 3 - 8
    A1 <- list(a3, a4, a5, a6, a7, a8)
    for (i in seq(length(A1)-1)) 
	ShadeBetween(H, H, A1[[i]](d,H), A1[[i+1]](d,H), col=dcols[i+4], , lwd=p.lwd,) 

    # Branches 9, 10, 11
    p2id <- PHI2.id(d,H); p3id <- PHI3.id(d,H)
    #lines(H, a7(d,H))
    #lines(H[p2id], phi2(d,H)[p2id])
    #lines(H, a8(d,H))
    ShadeBetween(H,       H, a7(d,H), a8(d,H),  col=dcols[11],  , lwd=p.lwd,)
    ShadeBetween(H,       H[p2id], a7(d,H), phi2(d,H)[p2id],  col=dcols[10],  , lwd=p.lwd,)
    ShadeBetween(H[p2id], H[p3id], phi2(d,H)[p2id], phi3(d,H)[p3id], col=dcols[11],  , lwd=p.lwd,)
    ShadeBetween(H[p3id], H[p3id], phi3(d,H)[p3id], a8(d,H[p3id]), col=dcols[12],  , lwd=p.lwd, )

    A2 <- list(a12, a13, a14, a15, a16, a17, a18, a19, a20)
    for (i in seq(length(A2)-1)) 
	ShadeBetween(H, H, A2[[i]](d,H), A2[[i+1]](d,H), col=dcols[i+14], , lwd=p.lwd,) 


    # Lines separating each region
    HP <- seq(0,d, length=10000)
    sapply(1:100, function(x) lines(HP[PROP1.id(d,HP,x)], PROP1(d,HP,x)[PROP1.id(d,HP,x)], lty='dashed', lwd=1))
    lines(H[p0id], phi0(d,H)[p0id], lwd=l.lwd)
    lines(H[p1id], phi1(d,H)[p1id], lwd=l.lwd)
    lines(H[p2id], phi2(d,H)[p2id], lwd=l.lwd)
    lines(H[p3id], phi3(d,H)[p3id], lwd=l.lwd)
    lines(H, a1(d,H), lwd=l.lwd)
    lines(H, a2(d,H), lwd=l.lwd)

    sapply(1:100, function(x) lines(HP[PROP2.id(d,HP,x)], PROP2(d,HP,x)[PROP2.id(d,HP,x)], lty='dashed', lwd=1))
    sapply(1:2, function(x)   lines(HP[PROP3.id(d,HP,x)], PROP3(d,HP,x)[PROP3.id(d,HP,x)], lty='dashed', lwd=1))

    for (i in seq(length(A1))) {
	lines(H, A1[[i]](d,H), lwd=l.lwd)
    }
    for (i in seq(length(A2))) {
	lines(H, A2[[i]](d,H), lwd=l.lwd)
    }

    # Borders
    lines(c(0,d), c(-90,-90), lwd=l.lwd)
    lines(c(0,0), c(-90,90), lwd=l.lwd)
    lines(c(d,d), c(-90,90), lwd=l.lwd)
    lines(c(0,d), c(90,90), lwd=l.lwd)
}
dev.off()

# =============================================================================

d <- 1
#eps <- 0.0001
#H <- c(seq(0,eps,length=10), seq(0,d, length=100), seq(d-eps, d, length=10))
H <- seq(0,d, length=1000)
source("RegionLines-IH.R")
# PLOT SEQUENCE
#pdf("Regions-IH.pdf", width=16, height=16)
#spar(mfrow=c(1,2))
mplot(H, a1(d, H), ylim=c(-90, 270), xlim=c(0,d), col='white', xaxs='i', yaxs='i', #spa=FALSE,
      ylab=tex("\\theta_0 \\in \\Pi_{j} \\  (h_0)"), 
      xlab=tex("h_0"),yaxt='n')
axis(side=2, at=seq(-90,270,by=30))

p.lwd <- 2
l.lwd <- 2.5
{

    As <- list(a0, a1, a3, a4, a6, a7, a5)
    colidx <- c(1,2,5,6,7,10)
    for (i in seq(length(As)-1)) {
	ShadeBetween(H, H, As[[i]](d,H), As[[i+1]](d,H),  col=dcols[colidx[i]], border='black', lwd=p.lwd,)
    }

    ShadeBetween(H, H, a5(d,H), phi3(d,H),  col=dcols[11], border='black', lwd=p.lwd,) #phi2 == a5
    id <- which(phi3(d,H) < a21(d,H))
    ShadeBetween(H, H[id], phi3(d,H), a21(d,H)[id],  col=dcols[12], border='black', lwd=p.lwd,)
    A1 <- list(a12,a13,a14,a15,a17,a18,a19,a20)
    for (i in seq(length(A1)-1)) 
	ShadeBetween(H, H, A1[[i]](d,H), A1[[i+1]](d,H), col=dcols[ifelse(i>=5,i+1,i)+14], lwd=p.lwd)  


    ShadeBetween(H, H, a21(d,H), a22(d,H), col=dcols[28], lwd=p.lwd,)
    # Borders
    lines(c(0,d), c(-90,-90), lwd=l.lwd)
    lines(c(0,0), c(-90,90), lwd=l.lwd)
    lines(c(d,d), c(-90,90), lwd=l.lwd)
    lines(c(0,d), c(90,90), lwd=l.lwd)
}
dev.off()

#Heights <- tpose(fread(lf("BranchHeights")))
#Angles <- tpose(fread(lf("BranchAngles")))
#mplot(Heights, Angles*180/pi, col=dcols, type='p', add=TRUE)
#mplot(Heights, Angles*180/pi, col=dcols, type='p', cex=0.5)

# ========================= FINITE HORIZON D = 1 ============================
# ===========================================================================




