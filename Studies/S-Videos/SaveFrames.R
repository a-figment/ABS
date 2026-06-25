

TemplateName <- "t0--Positions_"

x <- ReadDataList(lf(TemplateName), colidx=1)
y <- ReadDataList(lf(TemplateName), colidx=2)

################## POLYGONAL BILLIARDS ####################
###########################################################
#' @brief Parallel channel vertices from internal angle
getChannelVertices <- function(dx, alpha, d) {

    dy <- 1/(2*tan(alpha/2)) 
    x  <- c(-dx, 0, dx, dx, 0, -dx)
    y  <- c(0, dy, 0, d, d + dy, d)

    return(cbind(x,y))
}
#' @brief Parallel channel vertices from dy
getChannelVertices <- function(dx, dy, d) {
    x  <- c(-dx, 0, dx, dx, 0, -dx)
    y  <- c(0, dy, 0, d, d + dy, d)
    return(cbind(x,y))
}
###########################################################


PolygonBaseplate <- function(xVert, yVert, from, to) {
    plot(x=NA, y=NA, xlim=c(from, to), ylim=c(min(yVert), max(yVert)), asp=1, xaxt="n", yaxt="n", xlab="", ylab=" ")
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray")

}

PlotChannelBoundary <- function(xVert, yVert, from, to) {
    for (i in seq(floor(from), floor(to))) {
        polygon(xVert + i, yVert, col='white', border="white")
        lines(xVert[1:3] + i, yVert[1:3], col='black', lwd=2.5)
        lines(xVert[4:6] + i, yVert[4:6], col='black', lwd=2.5)
    }
}

V  <- getChannelVertices(1/2, 0.56, 0.5)


png(
PolygonBaseplate(V[,1], V[,2], -1, 1)
PlotChannelBoundary(V[,1], V[,2], -1, 1)
points(x[1,], y[1,], col='blue', cex=0.1, pch=19)
for (i in seq(nrow(x))) {
    #PlotChannelBoundary(V[,1], V[,2], -3, 3)
    #points(x[1:(i-1),], y[1:(i-1),], col='white', cex=0.1, pch=19)
    points(x[i,], y[i,],     col='blue', cex=0.1, pch=19)
}


##### Save frames
for (i in seq(1000)) {
    #png(
    #  paste0(i, "-frame.png"), 
    #  width     = 3.25,
    #  height    = 3.25,
    #  units     = "in",
    #  res       = 300,
    #)


    png(paste0(i, "-frame.png"), width=2560, height=1440)#, res=1200)#, units='cm'Colorspace)
    par(mar=c(0,0,0,0), bg='lightgray')
    #png(paste0(i, "-frame.png"))
    PolygonBaseplate(V[,1], V[,2], -1, 1)
    PlotChannelBoundary(V[,1], V[,2], -1, 1)
    points(x[i,], y[i,],     col='blue', cex=0.1, pch=19)
    dev.off()
}




space <- 10
ics <- 23
for (i in seq(1, nrow(x) - space, by=1)) {
    s <- seq(i, i + space)

    lines(c(x[i,ics], x[i+space,ics]), c(y[i,ics], y[i+space,ics]), col='blue', cex=0.1, pch=19)
}
dev.off()





