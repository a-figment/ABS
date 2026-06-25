library(viridis)
library(RColorBrewer)

setwd("../Studies/S-Videos/Data/")
TemplateName <- "t0--Positions_"
TemplateName <- "--Beam"

x <- ReadDataList(lf(TemplateName), REV=TRUE, colidx=1)
y <- ReadDataList(lf(TemplateName), REV=TRUE, colidx=2)


n <- round(ReadData(lf("-Beam.dat"))[5,1,] * 150) # Particles per beam

V <- getChannelVertices(0.5, pi/3, NA, -1/2)
PolygonBaseplate(V[,1], V[,2], -1, 1)
PlotChannelBoundary(V[,1], V[,2], -1, 1)
points(x[1,], y[1,], col='blue', cex=0.3, pch=19)
for (i in seq(nrow(x))) {
    #PlotChannelBoundary(V[,1], V[,2], -3, 3)
    #points(x[1:(i-1),], y[1:(i-1),], col='white', cex=0.1, pch=19)
    points(x[i,], y[i,],     col='blue', cex=0.1, pch=19)
}


##### Save frames
# Experimenting with different beam colors
cols <- rainbow(length(n))
cols <- viridis_pal(option="D")(length(n))
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
cols <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# Vertices
V <- getChannelVertices(dx=1/2, alpha=NA, d=.38, m=NA, dy=0.56)
for (i in seq(500)) {


    png(paste0(i, "-frame.png"), width=1440, height=1024)
    #par(mar=c(0,0,0,0), bg='lightgray')
    #PolygonBaseplate(V[,1], V[,2], -1, 1)
    #PlotChannelBoundary(V[,1], V[,2], -1, 1)
    
    
    # Lorentz
    par(mar=c(0,0,0,0), bg='lightgray')
    r <- 3/10
    th <- seq(0, 2*pi, length=200)
    cx <- rep(c(-1/2, 1/2, 3/2), 3)
    cy <- c(rep(3/2,3), rep(1/2, 3), rep(-1/2, 3))
    plot(x=NA, y=NA, xlim=c(-1/2, 3/2), ylim=c(-1/2, 3/2), asp=1, xaxt="n", yaxt="n", xlab="", ylab=" ")
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white")
    for (j in seq(length(cx))) {
	polygon(r*cos(th) + cx[j], r*sin(th) + cy[j], col='grey', lty='solid')
    }
    


    # Regular plotting
    points(x[i,], y[i,], col='blue', cex=0.3, pch=19)
    dev.off()
     
    # For beam colors 
    #s <- seq(n[1]-1)
    #points(x[i,s], y[i,s], col=cols[1], cex=0.4, pch=19)
    #for (j in seq(2,length(n)-2)) {
    #    s <- seq(sum(n[1:j]) - j + 1, sum(n[1:(j+1)])-j)
    #    points(x[i,s], y[i,s], col=cols[j+1], cex=0.4, pch=19)
    #}
    #points(x[i,(ncol(x)-n[length(n)]):ncol(x)], y[i,(ncol(x)-n[length(n)]):ncol(x)], col=cols[j+1], cex=0.4, pch=19)
    #dev.off()
}


pdf("~/Dropbox/PHD/Papers/SlicingRates/Images/InitialBeams.pdf", width=1440, height=1024)

png("~/Dropbox/PHD/Papers/SlicingRates/Images/InitialBeams.png", width=1440, height=1024)
par(mfrow=c(2,1), oma=rep(0,4), mar=rep(4.5, 4))
PolygonBaseplate(V[,1], V[,2],    -2,2)
PlotChannelBoundary(V[,1], V[,2], -2,2)
i <- 5000
s <- seq(n[1]-1)
points(x[i,s], y[i,s], col=cols[1], cex=0.2, pch=19)
for (j in seq(2,length(n)-2)) {
    s <- seq(sum(n[1:j]) - j + 1, sum(n[1:(j+1)])-j)
    points(x[i,s], y[i,s], col=cols[j+1], cex=0.2, pch=19)
}
points(x[i,(ncol(x)-n[length(n)]):ncol(x)], y[i,(ncol(x)-n[length(n)]):ncol(x)], col=cols[j+1], cex=0.4, pch=19)
dev.off()




################## POLYGONAL BILLIARDS ####################
###########################################################
#' @brief Parallel channel vertices from internal angle
getChannelVertices <- function(dx, alpha, d, m, dy) {

    options <- c("\n1: (dx, alpha, d)\n2: (dx, alpha, m)\n3: (dx, dy, d)\nInput: ")

    input <- readline(prompt=options) 

    x  <- c(-dx, 0, dx, dx, 0, -dx)
    y  <- c()
    if (input == 1) {
        dy <- 1/(2*tan(alpha/2)) 
        y  <- c(0, dy, 0, d, d + dy, d)
    } else if (input == 2) {
        dy <- 1/(2*tan(alpha/2)) 
        d  <- dy + m * dy
        y  <- c(0, dy, 0, d, d + dy, d)
    } else if (input == 3) {
        y  <- c(0, dy, 0, d, d + dy, d)
    } else {
        print("Invalid options.")
    }

    return(cbind(x,y))
}

PolygonBaseplate <- function(xVert, yVert, from, to) {
    plot(x=NA, y=NA, xlim=c(from, to), ylim=c(min(yVert), max(yVert)), asp=1, xaxt="n", yaxt="n", xlab="", ylab=" ")
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")
}

PlotChannelBoundary <- function(xVert, yVert, from, to) {
    #PolygonBaseplate(xVert, yVert, from, to)
    for (i in seq(floor(from), floor(to))) {
        polygon(xVert + i, yVert, col='white', border="white")
        lines(xVert[1:3] + i, yVert[1:3], col='black', lwd=2.5)
        lines(xVert[4:6] + i, yVert[4:6], col='black', lwd=2.5)
    }
}
###########################################################
