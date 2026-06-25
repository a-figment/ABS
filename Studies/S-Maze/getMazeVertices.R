

M <- 3 # M x M cell

getRandomCell <- function(M) { round(c(runif(1,M,n=1), runif(1,M,n=1))) }

# 0,1,2,3 in clockwise order starting from the bottom edge of a cell
# For left and right 
#' @i,j Index of the current cell
getRandomNeighbour <- function(M, i, j) {
    if (i == 1 && j == 1) { # Top left
	return (sample(c(0,1), 1))
    } else if (i == 1 && j == M) {
	return (sample(c(0,3), 1))
    } else if (i == M && j == 1) {
	return (sample(c(1,2), 1))
    } else if (i == M && j == M) {
	return (sample(c(2,3), 1))
    } else if (i == 1) { # Top row
	return (sample(c(0,1,3), 1)) 
    } else if (i == M) { # Bottom row
	return (sample(c(1,2,3), 1)) 
    } else if (j == 1) { # Left column
	return (sample(c(0,1,2), 1)) 
    } else if (j == M) { # Right column
	return (sample(c(0,2,3), 1)) 
    } else {
	return (sample(c(0,1,2,3), 1)) 
    }
}

# @i,j index of the current cell
# @edge number between 0-3 
getLocationFromEdge <- function(i, j, edge) {
    if (edge == 0) {
	return (c(i, j+1))
    } else if (edge == 1) {
	return (c(i+1, j))
    } else if (edge == 2) {
	return (c(i, j-1))
    } else if (edge == 3) {
	return (c(i-1, j)) 
    } else {
	print("Invalid edge")
    }
}

StartingCell <- SC <- getRandomCell(M)
CellsVisited <- matrix(0, nrow=M, ncol=M)
CellsVisited[SC[1], SC[2]] <- 1

nUnvisited <- M*M - 1
while (nUnvisited > 0) {

    # Random (valid) edge from current cell
    RN <- getRandomNeighbour(M, SC[1], SC[2])
    # Location of the random neighbour
    RNLoc <- getLocationFromEdge(SC[1], SC[2], RN)

    if (CellsVisited[RNLoc[1],RNLoc[2]] == 0) { # Random neighbour not visited




    }
    nUnvisited <- nrow(which(CellsVisited == 0, arr.ind=TRUE))
}


library(grid)
linesCount1 <- 0
lineAdd1 <- function(n1, x1, y1, x2, y2) {
    if ((x1 != x2) || (y1 != y2)) {
	linesCount1 <<- linesCount1 + 1
	grid.lines(c(x1,x2),c(y1,y2),default.units="native",gp=gpar(col="#000000",lwd=2))
	#cat(sprintf("%d. ( %d , %d ) to ( %d , %d )\n",linesCount1,x1,y1,x2,y2))

	cat(file="../lines.dat", append = TRUE,
	    sprintf("%d %d %d %d\n",x1,y1,x2,y2))
	    #sprintf("%f %f %f %f\n",x1,y1,x2,y2))
    }
}
mazeMake2 <- function(n1, rlen1, clen1, side1, x1, y1) {
    if ((rlen1 <= 1) || (clen1 <= 1)) {
	invisible(return())
    }
    x3 <- x1+side1*clen1
    y3 <- y1+side1*rlen1
    n2 <- n1+1
    r1 <- runif(1)
    c1 <- runif(1)
    if (rlen1 > clen1) {
	cpos1 <- floor(c1*clen1)
	rpos1 <- floor(r1*(rlen1-1))+1
	x2 <- x1+side1*cpos1
	y2 <- y1+side1*rpos1
	lineAdd1(n2,x1,y2,x2,y2)
	lineAdd1(n2,x2+side1,y2,x3,y2)
	mazeMake2(n2,rpos1,clen1,side1,x1,y1)
	mazeMake2(n2,rlen1-rpos1,clen1,side1,x1,y2)
    }
    else {
	cpos1 <- floor(c1*(clen1-1))+1
	rpos1 <- floor(r1*rlen1)
	x2 <- x1+side1*cpos1
	y2 <- y1+side1*rpos1
	lineAdd1(n2,x2,y1,x2,y2)
	lineAdd1(n2,x2,y2+side1,x2,y3)
	mazeMake2(n2,rlen1,cpos1,side1,x1,y1)
	mazeMake2(n2,rlen1,clen1-cpos1,side1,x2,y1)
    }
    invisible()
}
mazeMake1 <- function(rlen1, clen1, side1, gap1) {
    x1 <- side1*gap1
    y1 <- side1*gap1
    x3 <- x1+side1*clen1
    y3 <- y1+side1*rlen1
    lineAdd1(0,x1,y3,x3,y3)
    lineAdd1(0,x1,y1,x3,y1)
    lineAdd1(0,x1,y1,x1,y3-side1)
    lineAdd1(0,x3,y1+side1,x3,y3)
    mazeMake2(1,rlen1,clen1,side1,x1,y1)
    invisible()
}
mazeMake0 <- function(side1=1, rlen1=3, clen1=3) {
    gap1  <- 0
    linesCount1 <<- 0
    mazeMake1(rlen1,clen1,side1,gap1)
}


# Drawing the maze (remove lines.dat on each run)
PlotMaze <- function(side, len) {
    Lines <- data.matrix(fread("../lines.dat"))
    mplot(Lines[1,c(1,3)], Lines[1,c(2,4)], asp=1, xaxt='n', yaxt='n', ylab="", xlab="",
	  xlim=c(min(Lines[,c(1,3)]), max(Lines[,c(1,3)])),
	  ylim=c(min(Lines[,-c(1,3)]), max(Lines[,-c(1,3)])))
    for (i in seq(nrow(Lines))) {
	lines(Lines[i,c(1,3)], Lines[i, c(2,4)], lwd=3)
	#points(Lines[i,c(1,3)], Lines[i, c(2,4)], lwd=2, pch=10)
    }

    lines(c(0,0), c(len-side, len), col='green', lwd=3)
    lines(c(len,len), c(0, side), col='red', lwd=3)
}
len <- 3
mazeMake0(rlen1=len, clen1=len)
PlotMaze(side=1, len=len)

lines(x,y,lwd=0.2, col='blue')

sum(sqrt(diff(x)^2 + diff(y)^2))

x <- tpose(fread(lf("x.dat")))
y <- tpose(fread(lf("y.dat")))

x <- ReadData(lf("-x.dat"))[,1,]
y <- ReadData(lf("-y.dat"))[,1,]

exit.id <- apply(x, 2, function(x) min(c(which(x < 10)-1, which(x > 110)-1))) # Exit it
exit.id <- apply(x, 2, function(x) min(which(x > 110)-1)) # Exit it
#for (i in seq(ncol(x))) {
#    ids <- seq(1, ifelse(is.infinite(exit.id[i]), nrow(x), exit.id[i]))
#    mplot(x[,i],y[,i], add=TRUE)
#    points(x[,i],y[,i], add=TRUE)
#}

X <- list.cbind.na(sapply(1:length(exit.id), function(i) {
	    ids <- 1:ifelse(is.infinite(exit.id[i]), nrow(x), exit.id[i]);
	    x[ids,i] }))

Y <- list.cbind.na(sapply(1:length(exit.id), function(i) {
	    ids <- 1:ifelse(is.infinite(exit.id[i]), nrow(x), exit.id[i]);
	    y[ids,i] }))
X <- x
Y <- y

escapers <- which(!is.infinite(exit.id))
for (i in seq(nrow(x))) {
for (i in seq(1375)) {
    png(paste0(i, "-frame.png"), width=1440, height=1440)
    PlotMaze()
    points(X[i,-escapers], Y[i,-escapers], pch=19, cex=0.5, col='blue')
    points(X[i,escapers], Y[i,escapers], pch=19, cex=0.5, col='red')
    dev.off()
    #Sys.sleep(1)
}

points(x,y)


ET <- tpose(fread(lf("ET.dat")))
L <- getTSD_ecdf(dx=ET, n=1000, domain = exp.domain)
mplot(L$x, 1-L$y, log='xy')
LinFit(mat=1-L$y, time=L$x, plot=TRUE, range = 200:700, show.range=TRUE)

A <- tpose(fread(lf("Angle.dat")))

mplot(A*180/pi, ET, pch=19, cex=0.1, type='p', log='y', xlim=c(0,90))

######### Barrier billiard (testing) #########
xv <- c(0,1,1,0,0)
yv <- c(0,0,1,1,0)
mplot(xv, yv, lwd=2, asp=1)
lines(c(1/4, 3/4), c(1/2, 1/2), lwd=2)
x <- tpose(fread(lf("x.dat")))
y <- tpose(fread(lf("y.dat")))
points(x,y)
lines(x,y)
###############################################

