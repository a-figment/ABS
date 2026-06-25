
# Real trajectories
PlotTraj <- function(d, from=0,to=1, s=c(), ...) {
    dx <- dy <- 1/2
    dy <- dx
    PlotChannelBoundary(xVert=c(0, dx, 2*dx, 2*dx, dx, 0), yVert=c(0, dy, 0, d, dy+d,d), from=from,to=to)

    x <- data.matrix(fread(lf("x.dat")))
    y <- data.matrix(fread(lf("y.dat")))
    if (length(s) > 0)
	lines(x[s], y[s], ...)
    else 
	lines(x, y, ...)
    #lines(x[1:5], y[1:5])
    #lines(x[1:2], y[1:2], col='blue')
    points(x[1], y[1], col='green', pch=19, cex=1.1)
    #points(x[nrow(x)], y[nrow(y)], col='red', pch=19, cex=1.1)

}
dev.off()

# Functions for manually computing the crossing heights from a billiard trajectory
int <- function(x) {
    gz <- which(x >= 0)
    id1 <- id2 <- c()
    if (length(gz) == 0) { # All < 0
	val <- ceiling(x)
	#id <- which(val == 0 & x < 0)
    } else {
	val <- c(floor(x[gz]), ceiling(x[-gz]))
	#id <- which(val == 0 & x < 0)
    }

    id1 <- which(val == 0 & x < 0)
    id2 <- which(val == 0 & x > 0)
	
    if (length(id1) > 0)
	val[id1] <- -1
    if (length(id2) > 0)
	val[id2] <- +1

    return(val)
}

int <- function(x) {
    gz <- which(x >= 0)
    if (length(gz) == 0) { # All < 0
	val <- floor(x)
    } else {
	val <- c(ceiling(x[gz]), floor(x[-gz]))
    }
    return(val)
}


#' @x,y 1x2 vectors corresponding to two collision points
getCrossingHeight <- function(x,y) {
    #x <- c(xx[1], xx[2])
    #y <- c(yy[1], yy[2])
    ix <- int(x)
    mid <- yint <- 0
    if (ix[1] != ix[2]) { # A crossing occured

	m <- diff(y)/diff(x) #(y[2]-y[1])/(x[2]-x[1])

	if ((ix[1] == -1 && ix[2] == 1) || (ix[1] == 1 && ix[2] == -1)) {
	    mid <- 0
	} else if (ix[1] > 0 && ix[2] > 0) {
	    mid <- min(ix)
	} else if (ix[1] < 0 && ix[2] < 0) {
	    mid <- max(ix)
	}
	yint <- m * (mid - x[1]) + y[1]
	print(mid)
    }
    return(c(mid, yint))
}

#' @brief Gets the entrance heights from raw trajectories
getHeights <- function() {
    xx <- data.matrix(fread(lf("x.dat")))
    yy <- data.matrix(fread(lf("y.dat")))
    X <- Y <- c()
    for (j in seq(2,nrow(xx)-1)) {
	CH <- getCrossingHeight(xx[j:(j+1)], yy[j:(j+1)])
	X <- c(X, CH[1])
	Y <- c(Y, CH[2])
    }
    return(cbind(X, Y))
}


###
#   COMPLETE SCATTERING MAP: 
###
In <- function(val, bounds) { (val > bounds[1] && val < bounds[2]) }

# Integer functions used
kappa0 <- function(d, h0, t0) { ceiling((d - h0 - tan(t0))/(d*(tan(t0) -1))) }
gamma0 <- function(d, h0, t0) { ceiling((h0 - 1 - d)/(d*(1 - tan(t0)))) }
alpha0 <- function(d, h0, t0) { ceiling((-h0 - tan(t0))/(d*(1 + tan(t0)))) }

#' Updates the position of particles as they cross boundary
UpdatePosition <- function(t1, x0, sym) {
    #print(paste0("t0 = ", t0 * 180/pi, ", t1 = ", t1*180/pi))
    # cos(t0) > 0 always - by the time we enter this function we have already moved into the symmetry branch if the direction was negative

    if (sym && cos(t1) > 0) { #) || (sym && cos(t0) > 0)) { # Moving left
        #print("moving back")
        x1 <- x0 - 1	
    } else if (!sym && cos(t1) > 0) {
        #print("moving forward")
        x1 <- x0 + 1
    } else {
        #print("null move")
        x1 <- x0
    }
}

cot <- function(a) { cos(a)/sin(a) }

{
    #' @brief Scattering map for a first collision on Gamma0
    #' @h0, @t0, @x0 Entrance height, angle and current position
    a0  <- function(d,h0) { -pi/2 + 2*pi }
    a1  <- function(d,h0) { -atan(h0/d) + 2*pi }
    a2  <- function(d,h0) { -atan((h0+d)/2) + 2*pi }
    a3  <- function(d,h0) { -atan(h0) + 2*pi }
    a4  <- function(d,h0) { -atan(2*h0/3) + 2*pi }
    a5  <- function(d,h0) { -atan(h0/2) + 2*pi }   # Note: Crossover between neg/pos angles
    a6  <- function(d,h0) { +atan((d-h0)/2) }
    a7  <- function(d,h0) { +atan(2*(d-h0)/3) }
    a8  <- function(d,h0) { +atan(d-h0) }
    a9  <- function(d,h0) { +atan((2*d-h0)/2) }
    a10 <- function(d,h0) { +atan(2*(2*d-h0)/3) }
    a11 <- function(d,h0) { +atan(2*(3*d-h0)/3) } # Note: angle doesn't coincide with the lower bound for a first bounce on Gamma2

    Chi2 <- function(d, h0, t0) { 
	return(((h0 > 1/3) && In(t0, c(a2(d,h0),  a3(d,h0)))) 
	       || ((h0 < 1/3) && In(t0, c(a2(d,h0),  phi(d,h0)))))
    }
    Chi3 <- function(d, h0, t0) {
	return((h0 < 1/3) && In(t0, c(phi(d,h0),  a3(d,h0))))
    }

    getAnglesInChi.NA <- function(d=1/2, h0, n=1, i) {
	if (i == 0) {
		seq(a0(d,h0), a1(d,h0), length=n)
	} else if (i == 1) {
	    if (a2(d,h0) > a1(d,h0)) {
		seq(a1(d,h0), a2(d,h0), length=n)
	    } else { # Do nothing! No valid angles here
		NA
	    }
	} else if (i == 2) {
	    if (h0 > 1/3) {
		seq(a2(d,h0), a3(d,h0), length=n)
	    } else if (h0 < 1/3) {
		if (a2(d,h0) > a1(d,h0)) {
		    seq(a2(d,h0), phi(d,h0), length=n)
		} else if ((a2(d,h0) < a1(d,h0)) && (phi(d,h0) > a1(d,h0))) {
		    seq(a1(d,h0), phi(d,h0), length=n)
		} else { # Do nothing! No valid angles here
		    NA
		}
	    } else {
		NA
	    }
	} else if (i == 3) {
	    if (h0 < 1/3) {
		if (phi(d,h0) > a1(d,h0)) {
		    seq(phi(d,h0), a3(d,h0), length=n)
		} else {
		    seq(a1(d,h0), a3(d,h0), length=n)
		}
	    } else {
		NA
	    }
	} else if (i == 4) {
	    seq(a3(d,h0), a4(d,h0), length=n)
	} else if (i == 5) {
	    seq(a4(d,h0), a5(d,h0), length=n)
	} else if (i == 6) {
	    c(seq(a5(d,h0), 2*pi-0.00001, length=n/2), seq(0, a6(d,h0), length=n/2))
	} else if (i == 7) {
	    seq(a6(d,h0), a7(d,h0), length=n)
	} else if (i == 8) {
	    seq(a7(d,h0), a8(d,h0), length=n)
	} else if (i == 9) {
	    if (a9(d,h0) < a12(d,h0)) {
		seq(a8(d,h0), a9(d,h0), length=n)
	    } else {
		seq(a8(d,h0), a12(d,h0), length=n)
	    }
	} else if (i == 10) {
	    if (a10(d,h0) < a12(d,h0)) {
		seq(a9(d,h0), a10(d,h0), length=n)
	    } else {
		if (a9(d,h0) < a12(d,h0)) {
		    seq(a9(d,h0), a12(d,h0), length=n)
		} else {
		    NA
		}
	    }
	} else if (i == 11) {
	    if (a10(d,h0) < a12(d,h0)) {
		seq(a10(d,h0), a12(d,h0), length=n)
	    } else {
		NA
	    }
	} else if (i == 12) { # Crossover from G0 to G2
	    NA
	} else if (i == 13) {
		seq(a12(d,h0), a13(d,h0), length=n)
	} else if (i == 14) {
		seq(a13(d,h0), a14(d,h0), length=n)
	} else if (i == 15) {
		seq(a14(d,h0), a15(d,h0), length=n)
	} else if (i == 16) {
		seq(a15(d,h0), a16(d,h0), length=n)
	} else if (i == 17) { #a16 == a17
		seq(a17(d,h0), a18(d,h0), length=n)
	} else if (i == 18) {
		seq(a18(d,h0), a19(d,h0), length=n)
	} else if (i == 19) {
		seq(a19(d,h0), a20(d,h0), length=n)
	} else {
	    print("What else?")
	}
	#} else if (i == 19) {
	#	seq(a18(d,h0), a19(d,h0), length=n)
	#seq(A[[i+1]](d,h0), A[[i+2]](d,h0), length=n)
    }

    getAnglesInChi <- function(d=1/2, h0, n=1, i) {
	if (i == 0) {
		seq(a0(d,h0), a1(d,h0), length=n)
	} else if (i == 1) {
	    if (a2(d,h0) > a1(d,h0)) {
		seq(a1(d,h0), a2(d,h0), length=n)
	    } else { # Do nothing! No valid angles here
		seq(a1(d,h0), a1(d,h0), length=n)
	    }
	} else if (i == 2) {
	    if (h0 > 1/3) {
		seq(a2(d,h0), a3(d,h0), length=n)
	    } else if (h0 < 1/3) {
		if (a2(d,h0) > a1(d,h0)) {
		    seq(a2(d,h0), phi(d,h0), length=n)
		} else if ((a2(d,h0) < a1(d,h0)) && (phi(d,h0) > a1(d,h0))) {
		    seq(a1(d,h0), phi(d,h0), length=n)
		} else { # Do nothing! No valid angles here
		    seq(a1(d,h0), a1(d,h0), length=n)
		}
	    } else {
		seq(a1(d,h0), a1(d,h0), length=n)
	    }
	} else if (i == 3) {
	    if (h0 < 1/3) {
		if (phi(d,h0) > a1(d,h0)) {
		    seq(phi(d,h0), a3(d,h0), length=n)
		} else {
		    seq(a1(d,h0), a3(d,h0), length=n)
		}
	    } else {
		    seq(a3(d,h0), a3(d,h0), length=n)
	    }
	} else if (i == 4) {
	    seq(a3(d,h0), a4(d,h0), length=n)
	} else if (i == 5) {
	    seq(a4(d,h0), a5(d,h0), length=n)
	} else if (i == 6) {
	    seq(a5(d,h0)-2*pi, a6(d,h0), length=n) + 2*pi
	} else if (i == 7) {
	    seq(a6(d,h0), a7(d,h0), length=n)
	} else if (i == 8) {
	    seq(a7(d,h0), a8(d,h0), length=n)
	} else if (i == 9) {
	    if (a9(d,h0) < a12(d,h0)) {
		seq(a8(d,h0), a9(d,h0), length=n)
	    } else {
		seq(a8(d,h0), a12(d,h0), length=n)
	    }
	} else if (i == 10) {
	    if (a10(d,h0) < a12(d,h0)) {
		seq(a9(d,h0), a10(d,h0), length=n)
	    } else {
		if (a9(d,h0) < a12(d,h0)) {
		    seq(a9(d,h0), a12(d,h0), length=n)
		} else {
		    seq(a12(d,h0), a12(d,h0), length=n)
		}
	    }
	} else if (i == 11) {
	    if (a10(d,h0) < a12(d,h0)) {
		seq(a10(d,h0), a12(d,h0), length=n)
	    } else {
		seq(a12(d,h0), a12(d,h0), length=n)
	    }
	} else if (i == 12) { # Crossover from G0 to G2
		seq(a11(d,h0), a11(d,h0), length=n)
	} else if (i == 13) {
		seq(a12(d,h0), a13(d,h0), length=n)
	} else if (i == 14) {
		seq(a13(d,h0), a14(d,h0), length=n)
	} else if (i == 15) {
		seq(a14(d,h0), a15(d,h0), length=n)
	} else if (i == 16) {
		seq(a15(d,h0), a16(d,h0), length=n)
	} else if (i == 17) { #a16 == a17
		seq(a17(d,h0), a18(d,h0), length=n)
	} else if (i == 18) {
		seq(a18(d,h0), a19(d,h0), length=n)
	} else if (i == 19) {
		seq(a19(d,h0), a20(d,h0), length=n)
	} else {
	    print("What else?")
	}
	#} else if (i == 19) {
	#	seq(a18(d,h0), a19(d,h0), length=n)
	#seq(A[[i+1]](d,h0), A[[i+2]](d,h0), length=n)
    }

    # Gamma 2 angles
    a12 <- function(d, h0) { atan(1 - 2*h0) }
    a13 <- function(d, h0) { atan(1 - 2*h0/(1+2*d)) }
    a14 <- function(d, h0) { atan(1 - h0/(d + 1)) }
    a15 <- function(d, h0) { atan(1 + d - h0) }
    a16 <- function(d, h0) { atan(1 + 2*(d - h0)) }

    # Gamma 3 angles
    a17 <- function(d, h0) { atan(1 + 2*(d - h0)) }
    a18 <- function(d, h0) { atan(2*d + 1 - h0) }
    a19 <- function(d, h0) { atan(4*d + 1 - 2*h0) }
    #a20 <- function(d, h0) { atan((2*(d - h0) + 1)/(1-2*d)) }
    a20 <- function(d, h0) { pi/2 }
    phi <- function(d, h0) { -atan(2*(h0 + d)/5) + 2*pi }

}

Sg0 <- function(d,h0,t0,x0, sym) { # BOTTOM LEFT
    X <- c() # = (h0, t0)
    if (In(t0, c(a0(d,h0),  a1(d,h0)))) {			  # Angle checked
	X <- c(-h0*cot(t0), 
	       5*pi/2 - t0)
    } else if ((a2(d,h0) > a1(d,h0)) && In(t0, c(a1(d,h0),  a2(d,h0)))) {		  # Angle checked
	X <- c(h0 + d*(alpha0(d,h0,t0) + 3)*tan(t0) + d*(alpha0(d,h0,t0) + 1), 
	       t0)
    } else if ( ((h0 > 1/3) && In(t0, c(a2(d,h0), a3(d,h0)))) ||
	        ((h0 < 1/3 && a2(d,h0) >= a1(d,h0)) && In(t0, c(a2(d,h0), phi(d,h0)))) ||
		((h0 < 1/3 && a2(d,h0) < a1(d,h0)) && In(t0, c(a1(d,h0), phi(d,h0)))) ) {
	X <- c((h0+d)*cot(t0) + 5/2, 
	       2*pi - t0 + 3*pi/2)
    } else if ( ((h0 < 1/3 && phi(d,h0) >= a1(d,h0)) && In(t0, c(phi(d,h0), a3(d,h0)))) ||
	        ((h0 < 1/3 && phi(d,h0) < a1(d,h0)) && In(t0, c(a1(d,h0), a3(d,h0))))) {
	X <- c(h0 + d + 5*tan(t0)/2, 
	       t0)
    } else if (In(t0, c(a3(d,h0),  a4(d,h0)))) {		  # Angle checked
	X <- c(h0+(1+d)*tan(t0) + d, 
	       t0)
    } else if (In(t0, c(a4(d,h0),  a5(d,h0)))) {		  # Angle checked
	X <- c(h0*cot(t0) + 4*d, 
	       2*pi - t0 + 3*pi/2)	 
    } else if (In(t0, c(a5(d,h0),  2*pi)) || (t0 >= 0 && t0 < a6(d,h0))) {  # Angle checked
	X <- c(h0 + 4*d*tan(t0), 
	       t0)
    } else if (In(t0, c(a6(d,h0),  a7(d,h0)))) {		  # Angle checked
	X <- c(h0 + (2 + d)*tan(t0) - d, 
	       t0)
    } else if (In(t0, c(a7(d,h0),  a8(d,h0)))) {		  # Angle checked
	X <- c(-h0 + 2*(d - tan(t0)), 
	       t0 + pi)
    } else if (In(t0, c(a8(d,h0),  a9(d,h0)))) {		  # Angle checked
	X <- c(-h0 + 1 - 2*tan(t0), 
	       t0 + pi)
    } else if (In(t0, c(a9(d,h0),  a10(d,h0)))) {		  # Angle checked
	X <- c(4*d - (1-h0)/tan(t0), 
	       pi/2 - t0 + pi)
    } else if (In(t0, c(a10(d,h0), a11(d,h0)))) {
	X <- c(-h0 + 3*d - (1 + d)*tan(t0), 
	       t0 + pi)
    } else {
	print("Impossible angle in Sg0.")
    }
    X[3] <- UpdatePosition(X[2], x0, sym)
    return(X)
}

Sg2 <- function(d, h0, t0, x0, sym) { # TOP RIGHT
    X <- c()
    if (In(t0, c(a12(d,h0),  a13(d,h0)))) {
	X <- c(h0 + d*kappa0(d,h0,t0)*(tan(t0)-1) + tan(t0),	# Angle checked
	       t0)
    } else if (In(t0, c(a13(d,h0),  a14(d,h0)))) {		# Angle checked
	X <- c(-h0 + (1 - tan(t0))*(1 + d), 
	       t0 + pi)
    } else if (In(t0, c(a14(d,h0),  a15(d,h0)))) {		# Angle checked
	X <- c(cot(t0)*(h0 + (d + 1)*(tan(t0) - 1)), 
	       -t0 + 3*pi/2)
    } else if (In(t0, c(a15(d,h0),  a16(d,h0)))) {
	X <- c(-h0 + 2*d + 1 - tan(t0),				# Angle checked 
	       t0 + pi)
    } else {
	print("Impossible angle in Sg2.")
    }
    X[3] <- UpdatePosition(X[2], x0, sym)
    return (X)
}

Sg3 <- function(d, h0, t0, x0, sym) { # TOP LEFT
    X <- c()
    if (In(t0, c(a17(d,h0),  a18(d,h0)))) {
	X <- c(2*d + 1 - tan(t0) - h0,				# Angle checked 
	       t0 + pi)
    } else if (In(t0, c(a18(d,h0),  a19(d,h0)))) {		# Angle checked
	X <- c(cot(t0)*(h0 - 1 - 2*d + tan(t0)), 
	       -t0 + 3*pi/2)
    } else if (In(t0, c(a19(d,h0),  a20(d,h0)))) {
	X <- c(cot(t0) * (-h0 + 1 + d*gamma0(d,h0,t0)*(1-tan(t0)) + d*(1 + tan(t0))), 
	       ((-t0 + 5*pi/2) %% (2*pi)))
    } else {
	print("Impossible angle in Sg3.")
    }
    X[3] <- UpdatePosition(X[2], x0, sym)
    return(X)
}

THETA0 <- function(d, h0) { c(3*pi/2, atan(1 - 2*h0)) }
THETA2 <- function(d, h0) { c(atan(1 - 2*h0), atan(1 + 2*d - 2*h0)) }
THETA3 <- function(d, h0) { c(atan(1 + 2*d - 2*h0), pi/2  ) }
#THETA0 <- function(d, h0) { c(-pi/2, atan(1 - 2*h0)) }
#THETA2 <- function(d, h0) { c(atan(1 - 2*h0), atan(1 + 2*d - 2*h0)) }
#THETA3 <- function(d, h0) { c(atan(1 + 2*d - 2*h0), pi/2  ) }

#' @brief Global scattering function
Scatter <- function(d, h0, t0, x0, sym=FALSE) {
    if (In(t0, c(THETA0(d,h0)[1],2*pi)) || (t0 >= 0 && t0 < THETA0(d,h0)[2])) {  # BOTTOM LEFT
	Sg0(d,h0,t0,x0, sym)
    } else if (In(t0, THETA2(d,h0))) { # TOP RIGHT
	Sg2(d,h0,t0,x0, sym)
    } else if (In(t0, THETA3(d,h0))) { # TOP LEFT
	Sg3(d,h0,t0,x0, sym)
    } else {
	#print("Symmetry branch.")
	S <- c(h0, t0)
	if (S[2] > pi && S[2] < 3*pi/2) {
	    S <- Scatter(d,h0, 3*pi-t0, x0, sym=TRUE)
	    S[2] <- 3*pi - S[2] # Revert angle
        #} else if (S[2] > 3*pi/2 && S[2] < 2*pi) {
	#    S <- Scatter(d,h0, 2*pi-t0, x0, sym=TRUE)
	#    S[2] <- 2*pi - S[2]
	} else if (S[2] > pi/2 && S[2] < pi) {
	    S <- Scatter(d,h0, pi-t0, x0, sym=TRUE)
	    S[2] <- pi - S[2]
	} else {
	    print("??????????")
	}
	#print(paste0("Pre-modulo = ", S[2] * 180/pi))
	S[2] <- S[2] %% (2*pi)
	#print(paste0("Post-modulo = ", S[2] * 180/pi))
	S
    }
}


# Positions are correct but something is wrong with symmetry branch
{
    d  <- 1/2
    PlotTraj(d, from=-4,to=3)
    ha <- data.matrix(fread(lf("H0T0")))
    h0 <- ha[1,1] 
    t0 <- ha[2,1]
    #h0 <- 4/10
    #t0 <- pi/6
    S <- Scatter(d, h0, t0, x0=0)
    points(S[3], S[1], col='red', pch=19, cex=0.5)

    n <- function(d,h0,t0) { if (t0 > 0) { floor((d-h0)*cot(t0)/2) } }
    n(d,h0,t0)

    # Computing heights from trajectory
    #H <- getHeights()
    #points(H[,1], H[,2], cex=2.5)
    #H <- as.vector(c(ha[1,1], H[-which(H[,2] == 0),2]))

    # Heights from C++ simulation
    H <- tpose(fread(lf("--H.dat")))
    X <- tpose(fread(lf("--R.dat")))
    points(X[1:30], H[1:30])
 
    mplot(X[,1])

    h <- h0
    X <- 0
    #for (i in seq(length(H))) {
    for (i in seq(501)) {
	h0 <- S[1]
	t0 <- S[2]
	x0 <- S[3]
	S <- Scatter(d, h0, t0, x0)
	points(S[3], S[1], col='red', pch=19, cex=0.7)
	h <- c(h, h0)
	X <- c(X, x0)
    }
}

getTrajectory(d=1/2, h0=ha[1,1], t0=ha[1,2], x0=0, n=5)

m <- min(length(H), length(h))
mplot(H[1:m], h[1:m], type='p')

curve(tan(pi-a4(d,h0))*(x-S[3]) + h0, from = 0, to=10, add=TRUE, lty='dashed', col='red', lwd=1)
curve(tan(pi-a5(d,h0))*(x-S[3]) + h0, from = 0, to=10, add=TRUE, lty='dashed', col='red', lwd=1)


###################################### THE MAP ##########################################
#########################################################################################

#' @brief Gets a trajectory from the scattering map
#' @n The number of iterations of the map
getTrajectory <- function(d, h0, t0, x0, n=1) {
   X <- H <- A <- vector(length=n)
   S <- Scatter(d, h0, t0, x0)
   for (i in seq(n)) {
       H[i] <- S[1]
       A[i] <- S[2]
       X[i] <- S[3]
       S <- Scatter(d, h0=S[1], t0=S[2], x0=S[3])
   }
   return(list(A=A, H=H, X=X))
}

#' @brief Gets a trajectory with randomised initial conditions (inefficient for large studies)
#' @n The number of iterations of the map
getRandomTrajectory <- function(d=1/2, x0=0, n=10, eps=1e-5) {
    t0 <- runif(0, 2*pi-eps, n=1)
    h0 <- runif(eps, d-eps, n=1)
    X <- H <- A <- vector(length=n)
    S <- Scatter(d, h0, t0, x0)
    for (i in seq(n)) {
	H[i] <- S[1]
	A[i] <- S[2]
	X[i] <- S[3]
	S <- Scatter(d, h0=S[1], t0=S[2], x0=S[3])
    }
    return(list(A=A, H=H, X=X))
}



#h0 <- seq(0.01, 0.5-0.01, length=1000)
N <- 10
t0 <- runif(0, 2*pi, n=N)
t0 <- c(runif(0, pi/180, n=N/2), runif(359*pi/180, 2*pi, n=N/2))
h0 <- runif(0.01, 0.5-0.01, n=N)
Traj <- getTrajectory(d, 0.25, t0[1], x0=0, n=1000000)

mplot(Traj$H, type='p', cex=0.1, xlim=c(5e5,5.1e5))
mplot(Traj$X, type='l')
mplot(Traj$X, Traj$H, type='p', cex=0.1, xlim=c(-1000,100))

Heights <- Positions <- Angles <- vector(mode='list', length=length(h0)-1)
Heights[[1]] <- h0
H <- X <- c()
X <- matrix(NA, nrow=length(h0), ncol=3)
H <- matrix(NA, nrow=length(h0), ncol=3)
i <- 2
for (i in seq(2,length(h0))) {
    Traj <- getTrajectory(d, h0[i-1], t0[i-1], x0=0, n=100000)
    Heights[[i]] <- Traj$H
    #H <- c(H, Traj$H)
    #X <- c(X, Traj$X)

    #X <- rbind(X, c(Traj$X[10], Traj$X[100], Traj$X[1000], Traj$X[10000]))
    #X[i-1,] <- c(Traj$X[10], Traj$X[100], Traj$X[1000])
    #Positions[[i-1]] <- Traj$X
    #Angles[[i-1]] <- Traj$A
    
}
Heights <- list.cbind.na(Heights)
Positions <- list.cbind.na(Positions)
Angles <- list.cbind.na(Angles)



spar(mfrow=c(2,2))
Traj <- getRandomTrajecto
mplot(Traj$H[2:length(Traj$H)], Traj$H[1:(length(Traj$H)-1)], type='p', cex=0.05, spa=FALSE, xlab=tex("x_{n}"), ylab=tex("h_n"))#, xlim=c(0,10000))
mplot(Traj$X, Traj$H, type='p', cex=0.05, spa=FALSE, xlab=tex("x_{n}"), ylab=tex("h_n"))#, xlim=c(0,10000))
mplot(Traj$H, type='p', cex=0.1)
mplot(Traj$X, type='l', lwd=0.2)

Traj <- getRandomTrajectory(n=100)
MSD <- Traj$X^2
for (i in seq(1000)) {

    Traj <- getRandomTrajectory(n=100)
    MSD <- MSD + Traj$X^2
    print(i)
}
MSD <- MSD / 1000

MSD <- rowMeans(Positions^2)
mplot(MSD, log='xy')
mu <- getPowerScaling(t=1:length(MSD),	MSD, z=0.9)
mplot(mu, type='o')
mplot(mu, type='o', add=TRUE)


spar(mfrow=c(2,2))
for (i in 10^(seq(3,6))) {
    #s <- seq(1, i, by=10)
    #mplot(s, Traj$X[s], type='l', ylab=tex("x_n"), xlab=tex("n"), lwd=0.75, spa=FALSE)
    mplot(Traj$H[1:i], type='p', cex=0.1, spa=FALSE, log='x')
}
mplot(Traj$X[1:1000000], type='l', ylab=tex("x_n"), xlab=tex("n"), lwd=0.5)
points(rep(0,nrow(Heights)), Heights[,1])

msd <- rowMeans(Positions^2)
mplot(msd, log='xy')

mu <- getPowerScaling(t=1:length(msd),	msd, z=0.8)
mplot(mu)

:
Traj <- getTrajectory(d, h0[1], t0[1], x0, n=1000000)
L <- getTSD_hist(dx = Traj$H, n=100)
mplot(L$x, L$y, type='s', log='', xlim=c(0,1/2))

L <- getTSD_hist(dx = Heights[8,], n=100)
mplot(L$x, L$y, type='s')


H <- tpose(fread(lf("--H.dat")))
A <- tpose(fread(lf("--A.dat")))
X <- tpose(fread(lf("--R.dat")))
T <- tpose(fread(lf("--L.dat")))
which(T < 0)
B <- tpose(fread(lf("--BID.dat")))[-1]

mplot(cumsum(T), X, type='l', lwd=0.4)
L <- getTSD_hist(dx = (H), n=100)
mplot(L$x, L$y, type='s')


### Comparing with C++ implementation


H <- getHeights()
points(H[,1], H[,2])
PlotTraj(1/2, from=-5, to=0, s=1:28, lwd=1)
{
    H <- tpose(fread(lf("--H.dat")))
    A <- tpose(fread(lf("--A.dat")))
    X <- tpose(fread(lf("--R.dat")))

    Traj <- cbind(H, A, X)

    STraj <- getTrajectory(d=1/2, h0=ha[1,1], t0=ha[1,2], x0=0, n=35)
    STraj <- rbind(c(ha[1,1], ha[1,2], 0), cbind(STraj$H, STraj$A, STraj$X))

    Traj
    STraj
    points(Traj[,3], Traj[,1], col='blue')
    points(STraj[,3], STraj[,1])
}

spar(mfrow=c(2,2))
mplot(Traj[,1], STraj[,1], type='p', spa=FALSE)
mplot(Traj[,2], STraj[,2], type='p', spa=FALSE)
mplot(Traj[,3], STraj[,3], type='p', spa=FALSE)

mplot(X, H, pch=19, cex=0.1)

pdf("BranchID.pdf", height=12, width=18)
spar(mfrow=c(2,2))
bid <- tpose(fread(lff("p3", "BID")))[-1,]
hist(bid, breaks=seq(0,20), freq = TRUE, 
     xlab="Branch index", 
     main=tex("(h_0, \\theta_0) = (2/5, \\pi/8)"), 
    xaxt='n')
axis(side=1, labels=seq(0,20), at=seq(0,20)+0.5)

mplot(bid[1:500], type='s', lwd=0.8)

r <- rle(bid)
mplot(r$values, r$lengths, type='s')


seq(-2*pi, 0, length=6) 
seq(-2*pi, 0, length=6) %% (2*pi)



# DEPARTURE FROM A BILLIARD TRAJECTORY

H <- tpose(fread(lf("--H.dat")))
HH <- tpose(fread(lf("RawH.dat")))
X <- tpose(fread(lf("R.dat")))
XX <- tpose(fread(lf("RawX.dat")))
XXX <- getTrajectory(d=1/2, h0=ha[1,1], t0=ha[1,2], x0=0, n=100000)
x <- data.matrix(fread(lf("x.dat")))
y <- data.matrix(fread(lf("y.dat")))

pdf("DepartureCombined.pdf", width=18, height=24)
spar(mfrow=c(2,1))
{
    close <- which(abs(abs(x) - 48.5) < 0.0001)
    id <- (close-10):close
    PlotTraj(d, from=-49.5,to=-47.5, col='blue', lwd=2, s=id)
    st <- 3
    lines(x[(close):(close+st)], y[(close):(close+st)], col='red', lwd=2)
    points(x[close], y[close], col='red', cex=4, lwd=2)
    points(x[id[1]], y[id[1]], col='green', pch=19, cex=2)
    points(x[close+st], y[close+st], col='red', pch=19, cex=2)

    # Extending the trajectory on the other side
    close <- min(which(abs(X[1:36000] - XX[1:36000]) > 0))
    points(X[(close-4):close], H[(close-4):close], cex=2, pch=19)
    points(X[close], H[close], cex=2, pch=19)
    #xx <- data.matrix(fread(lf("x.dat")))
    #yy <- data.matrix(fread(lf("y.dat")))
    lines(xx-49, yy, lty='dashed', lwd=2)
}

    # Long version
    mplot(X, type='l', col='black', lwd=.5, ylab=tex("x_n"), xlab=tex("n"), xlim=c(0,1e5), spa=FALSE)
    mplot(XX, type='l', lwd=.5, add=TRUE, col='blue', spa=FALSE)
    legend('bottomleft', legend=c("Scattering map", "Billiard simulation", "Point of departure"), col=c('black', 'blue', 'red'), lty=c('solid', 'solid', NA), pch=c(NA, NA, 1), cex=2, lwd=2)
    points(34955,-49.5, cex=3, col='red')
dev.off()

lines(XXX$X, lwd=1)
lines(Y)

points(X[34955], H[34955])
points(X[1:34955], H[1:34955], pch=19, cex=0.5, col='red')
points(X[34956:34960], H[34956:34960], pch=19, cex=0.5, col='black')

H <- tpose(fread(lf("--H.dat")))
HH <- tpose(fread(lf("RawH.dat")))
mplot(abs(H[1:36000] - HH[1:36000]), type='l', col='red', lwd=0.8 ,xlim=c(34000, 36000))
)
mplot(cumsum(HH), type='l', col='blue', lwd=1)

X <- tpose(fread(lf("R.dat")))
x <- data.matrix(fread(lf("x.dat")))

# 3D PLOTS
htmlwidgets::saveWidget(rglwidget(width = 800, height = 800), 
                        file = "3DChi.html",
                        libdir = "libs",
                        selfcontained = FALSE
                        )

library(rgl)
plot3d(x=NA, y=NA, z=NA, 
       #xlim=c(0,1/2), ylim=c(0,1/2), zlim=c(0,360), col=cols,
       #xlim=c(0,1/2), ylim=c(0,360), zlim=c(0,360), col=cols,
       #xlim=c(0,1/2), ylim=c(0,360), zlim=c(0,1/2), col=cols,
       xlim=c(0,360), ylim=c(0,1/2), zlim=c(0,360), col=cols,
       xlab=("h0"), ylab=("h1"), zlab=("theta1"))
       #xlab=("theta0"), ylab=("h1"), zlab=("theta1"))
       #xlab=tex("theta_0"), ylab=tex("h_1"), zlab=tex("theta_1"))

eps <- 0.0001
N <- 50
h0 <- seq(eps, 0.5 - eps, length=N)
for (j in seq(0, 19)) {
    for (i in seq(length(h0))) {
	AA <- getAnglesInChi.NA(h0=h0[i], n=length(h0), i=j)
	AA <- AA[-c(1, length(AA))]
	if (length(AA) > 1) {
	    t0 <- AA
	    S1 <- tpose(sapply(t0, function(x) 
			       Scatter(d=1/2, h0=h0[i], t0=x, x0=0, sym = FALSE)))
	    df <- data.frame(x=h0[i], y=S1[,1], z=S1[,2]*180/pi)
	    #df <- data.frame(x=h0[i], y=t0*180/pi, z=S1[,2]*180/pi)
	    #df <- data.frame(x=h0[i], y=t0*180/pi, z=S1[,1])
	    plot3d(x=df$x, y=df$y, z=df$z, col=cols[j+1], add=TRUE)
	}
    }
}

# For each h0, generate angles [0,2pi] and iterate once
for (i in seq(length(h0))) {
    #t0 <- c(seq(3*pi/2 + eps, 2*pi-eps, length=N/2), seq(0, pi/2-eps, length=N/2))
    t0 <- seq(0,2*pi-eps, length=N)
    S1 <- tpose(sapply(t0, function(x) 
		       Scatter(d=1/2, h0=h0[i], t0=x, x0=0, sym = FALSE)))
    df <- data.frame(x=h0[i], y=S1[,1], z=S1[,2]*180/pi)
    plot3d(x=df$x, y=df$y, z=df$z, col='blue', add=TRUE)
}
# For each t0, generate heights [0,h] and iterate once
t0 <- seq(0,2*pi-eps, length=N)
for (i in seq(length(t0))) {
    h0 <- seq(eps,d-eps, length=N)
    S1 <- tpose(sapply(h0, function(x) 
		       Scatter(d=1/2, h0=x, t0=t0[i], x0=0, sym = FALSE)))
    df <- data.frame(x=t0[i]*180/pi, y=S1[,1], z=S1[,2]*180/pi)
    plot3d(x=df$x, y=df$y, z=df$z, radius=200, col='blue', add=TRUE)
}

plot3d(x=NA, y=NA, z=NA, 
       xlim=c(0,1/2), ylim=c(0,1/2), zlim=c(0,360), col=cols,
       xlab="h0", ylab="h1", zlab="theta1")

m <- 1
N <- 100
for (i in seq(length(h0))) {

    t0 <- seq(0,2*pi-eps, length=N)
    S1 <- tpose(sapply(t0, function(x) {
		       Tr <- getTrajectory(d=1/2, h0=h0[i], t0=x, x0=0, n=m)
		       c(Tr$H[m], Tr$A[m], Tr$X[m])
		}))
		       #Scatter(d=1/2, h0=h0[i], t0=x, x0=0, sym = FALSE)))

#    Sprev <- tpose(sapply(t0, function(x) {
#		       Tr <- getTrajectory(d=1/2, h0=h0[i], t0=x, x0=0, n=m-1)
#		       c(Tr$H[m-1], Tr$A[m-1], Tr$X[m-1])
#		}))
    df <- data.frame(x=h0[i], y=S1[,1], z=S1[,2]*180/pi)
    #df <- data.frame(x=Sprev[,1], y=S1[,1], z=S1[,2]*180/pi)
    plot3d(x=df$x, y=df$y, z=df$z, col='blue', add=TRUE)
}

    t1 <- S1[,2]
    h1 <- S1[,1]
    for (j in seq(length(h1))) { 
	S2 <- tpose(sapply(t1, function(x) 
			   Scatter(d=1/2, h0=h1[j], t0=x, x0=S1[j,3], sym = FALSE)))
	#df <- data.frame(x=h1[j], y=S2[,1], z=S2[,2]*180/pi)
	#plot3d(x=df$x, y=df$y, z=df$z, radius=200, col='blue', add=TRUE)

	t2 <- S2[,2]
	h2 <- S2[,1]
	for (k in seq(length(h2))) { 
	    S3 <- tpose(sapply(t2, function(x) 
			       Scatter(d=1/2, h0=h2[k], t0=x, x0=S2[k,3], sym = FALSE)))
	    #df <- data.frame(x=h2[k], y=S3[,1], z=S3[,2]*180/pi)
	    #plot3d(x=df$x, y=df$y, z=df$z, radius=200, col='blue', add=TRUE)

	    t3 <- S3[,2]
	    h3 <- S3[,1]
	    for (s in seq(length(h3))) { 
		S4 <- tpose(sapply(t3, function(x) 
				   Scatter(d=1/2, h0=h3[s], t0=x, x0=S3[k,3], sym = FALSE)))
		df <- data.frame(x=h3[k], y=S4[,1], z=S4[,2]*180/pi)
		plot3d(x=df$x, y=df$y, z=df$z, radius=200, col='blue', add=TRUE)
	    }
	}
    }
}






h0 <- H[1]
t0 <- T0(d,h0)
S1 <- tpose(sapply(t0, function(x) Scatter(d=1/2, h0=h0, t0=x, x0=0, sym = FALSE)))

data <- df <- data.frame(x=h0, y=S1[,1], z=S1[,2]*180/pi)
data <- df <- data.frame(x=h0, y=t0*180/pi, z=S1[,1])
plot3d(x=df$x, y=df$y, z=df$z, radius=200, col='blue',
       xlim=c(0,1/2), ylim=c(0,1/2), zlim=c(0,360),
       xlab="h0", ylab="h1", zlab="theta1")

for (i in seq(2,length(H))) {
    h0 <- H[i]
    t0 <- T0(d,h0)
    S1 <- tpose(sapply(t0, function(x) Scatter(d=1/2, h0=h0, t0=x, x0=0, sym = FALSE)))
    data <- df <- data.frame(x=h0, y=S1[,1], z=S1[,2]*180/pi)
    data <- df <- data.frame(x=h0, y=t0*180/pi, z=S1[,1])
    plot3d(x=df$x, y=df$y, z=df$z, col='blue', add=TRUE)
}


# Run times

Recur <- ReadData(lf("RunTimes"))[,1,]
NRecur <- ReadData(lf("RunTimes"))[,2,]

apply(Recur, 2, mean)
apply(NRecur, 2, mean)

sqrt(apply(Recur, 2, var))
sqrt(apply(NRecur, 2, var))

mplot(Recur, log='xy')
mplot(NRecur, add=TRUE, lty='dashed')
