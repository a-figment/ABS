

#' @brief Given the fraction of particles in each cell, determine
#'	  the location at which they are sliced (in order to produce that area)
#' @L[vector] Fractions of particles at a given time
#' @w Width of elementary cells (also used for plotting)
#'
getSliceLocations <- function(w, L, print=TRUE) {
    S <- c()
    for (i in seq(length(L)-1)) 
	S[i] <- -w + sum(L[1:i])
    if (print)
	print(data.frame(S=diff(c(-w, S, w)), L=L, isSame=abs(diff(c(-w,S,w))-L) < 1e-15))
    return(S)

}

#' @brief A plotting function which resets the grid
#' 
reset <- function(xmm=5, ym=10, w=1/2, N=1000, grid=TRUE, clines=FALSE) {
    xlim <- c(-xmm, xmm)
    ylim <- c(0, ym)
    mplot(0,0, type='p', col='white', cex=.01, xlab="Distance (m)", ylab="Time (n)", xaxt='n', yaxt='n', xlim=xlim, ylim=ylim)

    #
    axis(side=1, at=seq(-N, N), labels=seq(-N,N)) 
    axis(side=2, at=seq(0,N)+1/2, labels = seq(0,N))

    # Partition cells
    if (grid)
	abline(h=seq(0, N, by=1), v=seq(-N,N, by=1/2)[-seq(1,4*N+1,by=2)], lty='solid') # Grid

    # Draw the first mass
    rect(xleft=-w, xright=w, ybottom=0, ytop=1, col = 'green')

    # Bounding cone
    curve(x-w, from = 0,  to = N, add=TRUE)
    curve(-x-w, from = 0, to = -N, add=TRUE)

    if (clines) 
	abline(v=seq(-N, N), lty='dashed')
}


#' @brief Draws the sliced particles 
#' @n Number of time units to draw
#'
Draw <- function(fracs, w=1/2, n=10, hfcol='green', sfcol='red', arrows=TRUE, print=FALSE) {

    # Draw the first couple of slices
    S <- getSliceLocations(w, fracs[[1]])
    rect(xleft=-w-1, xright=S[1]-1, ytop=1, ybottom=1+1, col=hfcol)
    rect(xleft=S[1]+1, xright=w+1, ytop=1, ybottom=1+1,  col=hfcol)

    S <- getSliceLocations(w, fracs[[2]])
    rect(xleft=-w-2, xright=S[1]-2, ytop=2, ybottom=2+1, col=hfcol)
    rect(xleft=S[1], xright=S[2], ytop=2,   ybottom=2+1, col=sfcol)
    rect(xleft=S[2]+2, xright=w+2, ytop=2,  ybottom=2+1, col=hfcol)

    # Draw the rest of them
    for (j in seq(3,n)) {

	S <- getSliceLocations(w, fracs[[j]], print=FALSE)

	# Heading fronts
	rect(xleft=-w-j, xright=S[1]-j, ytop=j, ybottom=j+1, col=hfcol)
	rect(xleft=S[j]+j, xright=w+j, ytop=j, ybottom=j+1,  col=hfcol)

	K <- seq(-j, j, by=2)
	for (i in seq(2, length(K)-1)) {
	    rect(xleft=K[i]+S[i-1], xright=S[i]+K[i],     ytop=j, ybottom=j+1, col=sfcol)
	    #rect(xleft=K[i+1]+S[i], xright=S[i+1]+K[i+1], ytop=j, ybottom=j+1, col='red')

	    if (arrows && j > 2) {
		# Arrows and lines pointing to slice location
		J <- j - 2
		s <- getSliceLocations(w, fracs[[J]], print=FALSE)
		arrows(x0=K[i]+s[i-1]+1,   x1=K[i]+s[i-1]+1, y1=J, y0=J+.5, col='black', len=0.1, lwd=2)
		dvline(x=K[i]+s[i-1]+1, y0=J-1, y1=J, lty='dashed', col='black')
	    }

	}
    }
}




# Obtain the fractions of particles in each cell at a given time (counts)
fracs <- data.matrix(fread(lff("HeadingCounts")))
fracs <- ReadData(lf("HeadingCounts"))
M <- sum(fracs[1,])

# Remove zeros
raw <- sapply(1:nrow(fracs), function(x) as.vector(fracs[x,-which(fracs[x,]==0)])/M)

# Plotting info
reset(xmm=5, ym=10, w=1/2, grid=TRUE, clines=FALSE)
Draw(raw, n=20, arrows=TRUE, print=FALSE)


R <- list.cbind.na(raw)

mplot(R[1:101,100], type='o', pch=0.5)

# Location of slicer as a function of 
#pdf("../Images/Lmn.pdf", width=14, height=12)
{
    n <- c(1,10,50,250,500)
    n <- c(1,10)
    mplot(1:10, tpose(R[n,]), type='o', log='xy', xlab="Time (m)", ylab=tex("L(m,n)"))
    legend('bottomleft', legend=paste0("n = ", n), pch=19, lty='solid', col=rbow(), cex=2)
}
#dev.off()

LinFit(mat=R[1,], time=1:1000, range=100:1000, plot=TRUE)
LinFit(mat=R[10,], time=1:1000, range=300:1000, plot=TRUE)

mplot(diag(R), type='o', log='xy')
mplot(R[,800], type='o', log='y')

L <- getTSD_ecdf(dx = R[,100] ,n=1000, domain=exp.domain)
mplot(L$x, 1-L$y, type='s', log='xy')

LinFit(mat=diag(R), time=1:1000, range = 10:1000, plot=TRUE)


###
#  The mathematical definition of a slicer
###

reset(xmm=5, ym=10, w=1/2, grid=TRUE, clines=FALSE)
Draw(raw, n=20, arrows=TRUE, print=FALSE)
A <- raw

#' Gets the location of a slicer at (n,m)
L <- function(A, n, m, w=1/2) { 

    A <- A[[n]]

    if (m < 1) { # m < -n
	-w
    } else if (m > length(A)) {
	w
    } else {
	-w + sum(A[1:m])
    }
}

n <- 1
abline(v=L(A,n,0), col='blue')
abline(v=L(A,n,1), col='blue')
abline(v=L(A,n,2), col='blue')

for (i in seq(2,8)) {
    abline(v=L(A,i,2), col=rbow()[i], lwd=2)
}

n <- 2
abline(v=L(A,n,0)-n,     col='blue')
abline(v=L(A,n,1)-n,     col='blue')
abline(v=L(A,n,2)+n,     col='blue')
abline(v=L(A,n,n+1)+n,   col='blue')

diff(c(L(A,n,0),  
L(A,n,1), 
L(A,n,2), 
L(A,n,n+1)))

n <- 3

abline(v=L(A,n,0), col='blue')

abline(v=L(A,n,1), col='blue')
abline(v=L(A,n,n), col='blue')

abline(v=L(A,n,2)+2, col='blue')
abline(v=L(A,n,n+1), col='blue')


raw[[1]]
raw[[2]]





