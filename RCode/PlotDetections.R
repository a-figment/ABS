
img.path = "~/Dropbox/PHD/Papers/JumpStatistics/Images/"
setwd("../Studies/S-JumpStatistics/Data/")


###
#   Intuition for mean perp. distance deviation
#   Plots a single jump event (the first detected) along with the window
#   size and the mean perp. distance
###
W <- c(5, 50, 500, 5000) # Must have 4 window sizes for panel plot (ideal case)
W <- c(10000, 50000, 100000, 250000) /10# Must have 4 window sizes for panel plot (ideal case)
W <- c(10000, 25000, 50000, 75000) # Must have 4 window sizes for panel plot (ideal case)
W <- c(250, 500, 2500, 5000) # Must have 4 window sizes for panel plot (ideal case)
W <- c(50, 500, 2500, 25000) # Must have 4 window sizes for panel plot (ideal case)
W <- c(1000, 5000, 10000)
) # Must have 4 window sizes for panel plot (ideal case)
t <- seq(0, 1e7, by=10)
t <- seq(0, 1e6, by=1)
t <- seq(0, 1e5, by=1)

# Single traj
j <- 3
i <- 5
jump.file <- lf(paste0("_jump.dat"))
jump.file <- lf(paste0("w", W[j], "_jump.dat"))
jump <- data.matrix(fread(jump.file))
x.file <- gtools::mixedsort(lff(paste0("w", W[j], "_"), "x.dat"))
x.file <- rev(gtools::mixedsort(lf("--x.dat")))
x <- data.matrix(fread(x.file[i]))
s <- seq(1e6)
matplot(t[s], x[s,], lwd=2, lty='solid', type='l', xlab="", ylab="")

ShowJT(t0=jump[,2], x0=jump[,3], tau=jump[,4], jd=jump[,5], lwd=2, cex=1, jump.line=TRUE)

# System parameters
alpha <- "a1"
d <- "d2"
dxs <- 100
mpd <- 0

show.all <- TRUE
file.name.key <- paste0(alpha, "-", d, "-dxs", dxs, "-mpd", mpd, "--") 
AS <- ifelse(show.all, "A", "S") 

for (j in seq(length(W))) {
    jump.file <- lf(paste0("w", W[j], "_jump.dat"))
    jump <- data.matrix(fread(jump.file))
    id <- which(jump[,1] == 0)
    ranges <- sapply(1:(length(id)-1), function(x) seq(id[x], id[x+1]-1))
    #ranges <- lapply(ranges, function(x) { if (x[1] != 1) { x + 1 } else { x }})
    x.file <- gtools::mixedsort(lff(paste0("w", W[j], "_"), "x.dat"))

    length(id) == length(x.file) # Required

    file.group <- as.list(split(x.file, ceiling(seq_along(x.file)/4)))
    id.group <- as.list(split(id, ceiling(seq_along(id)/4)))
    ranges.group <- as.list(split(ranges, ceiling(seq_along(ranges)/4)))
    
    pdf(file=paste0(img.path, "/jump_examples/", file.name.key, "w", W[j], "_", AS, "jump.pdf"), width=18, height=14)
    #pdf(file=paste0(img.path, "delme.pdf"), width=18, height=14)

    for (f in seq(length(ranges.group))) { # Multiples of 4
    
        spar(cx=3, mfrow=c(2,2), col='white')
        par(oma=c(3,3,3,3), mar=c(3,3,3,3))
	par(mgp=c(4,2,0))
    
        for (i in seq(length(ranges.group[[f]]))) { # Multiples of 4

	    print(paste0(j, "-", f, "-", i))
            x <- data.matrix(fread(file.group[[f]][i]))

	    Plotting(t, x, jump, ranges.group[[f]][[i]], show.all, id.group[[f]][i], i)

        }
    }
    dev.off()
}

Plotting <- function(t, x, jump, ranges, show.all, id, k) {
    K <- s <- c()

    if (!show.all) {
        # SINGLE jUMP
	eps <- jump[id,4] * 0.1/2
        s <- seq(jump[id,2] - eps, jump[id,2] + jump[id,4] + eps)
        if (s[1] < 0) {
            s <- seq(0, jump[id,2] + jump[id,4] + eps)
        } else if (s[length(s)] > nrow(x)) {
            s <- seq(jump[id,2] - eps, nrow(x))
        } 

	#if (length(which(s > nrow(x))) > 0)
	#    s <- seq(1, nrow(x))
        s <- round(s)
        K <- id
    } else {
        # ALL JUMPS
    	s <- seq(nrow(x))
    	K <- ranges
        if (length(s) > 1e4)
	    s <- round(seq(s[1], s[length(s)], length=1e4))
    }

    #print(paste(s[1], s[length(s)]))
    #s <- seq(697901, 698981)

    matplot(t[s], x[s,], lwd=3.5, lty='solid', type='l', xlab="t", ylab="x(t)")
    ShowJT(t0=jump[K,2], x0=jump[K,3], tau=jump[K,4], jd=jump[K,5], lwd=3, cex=1, jump.line=TRUE)
    
    dy <- abs(max(x[s,]) - min(x[s,]))
    dx <- abs(max(t[s]) - min(t[s]))
    text(x=min(t[s] - 0.0175*dx), y=max(x[s,]) - 0.01*dy, labels=k, cex=3)

}


ShowJT <- function(t0, x0, tau, jd, mpd, jump.line=FALSE, cex=1, lwd=1) {
    for (i in seq(length(t0))) {
	points(t0[i], x0[i], pch=19, col='green', cex=cex)
    	points(t0[i] + tau[i], x0[i] + jd[i], pch=19, col='green', cex=cex)

	# Horizontal 
	ymin <- min(x0[i], x0[i] + jd[i])
	ymax <- max(x0[i], x0[i] + jd[i])
	arrows(x0 = t0[i], x1 = t0[i] + tau[i], y0 = ymin, y1 = ymin, length=0, col='red', lwd=lwd)
	arrows(x0 = t0[i] + tau[i], x1 = t0[i] + tau[i], y0 = ymin, y1 = x0[i] + jd[i], length=0, col='blue', lwd=lwd)
	arrows(x0 = t0[i], x1 = t0[i], y0 = ymin, y1 = x0[i], length=0, col='blue', lwd=lwd)
	if (jump.line) {
	    lines(c(t0[i], t0[i] + tau[i]), c(x0[i], x0[i] + jd[i]), col='magenta', lwd=lwd)
	}
    }
}

#    if (!show.all) { # Draw text
#	c.x <- c.y <- c()
#        if (jump[K,5] < 0) {
#            c.x <- jump[K,2] 
#            c.y <- jump[K,3] + jump[K,5]
#        } else {
#            c.x <- jump[K,2] + jump[K,4] 
#            c.y <- jump[K,3]
#        }
#        c.x <- c.x - jump[K,5]/abs(jump[K,5]) * 0.1 * jump[K,4]
#        c.y <- c.y + 0.15 * abs(jump[K,5])
#        text(x = c.x, y = c.y, labels=paste0("w = ", W[j]), cex=2, col='green')
#    
#        c.x <- c.x - jump[K,5]/abs(jump[K,5]) * 0.1/2 * jump[K,4]
#        c.y <- c.y - 0.15/2 * abs(jump[K,5])
#        text(x = c.x, y = c.y, labels=TeX(paste0("$<D_p> = ", format(jump[K,6], digits=2), " \\pm ", format(jump[K,7],digits=2))), cex=2, col='green')
#    }
