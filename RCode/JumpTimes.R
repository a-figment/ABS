#' @brief Used for manually inspecting jump times
#' @x Time at which observaitons were made
#' @y Matrix of position or displacement trajectories
#' @w Window half width
#' @i Start index
#' @Data If manual inspection is interrupted, continue where left off
#' 
CommandLoop <- function(x, y, w, sp=0.01, I=1, Data=NA, K = 50) {

    N <- ncol(y)
    data <- vector(mode='list', length=ncol(y)) # Turning points, widths, etc.

    # In case of interruptions, continue from last update
    if (!is.na(Data)) {
	cat("@data list will resume filling from index", I)
        cont <- readline("Continue? [y/n] ")
	if (cont == "y") {
	    data <- Data
	} else {
	    return(Data)
	}
    }

    for (i in seq(I, N)) {
    
        X <- x
        Y <- y[,i]
        start.width <- w

        spar(cx=1.75)
        data[[i]] <- getMM(x=X, y=Y, w=start.width, span=sp)
	JTPlot(x.raw=X, y.raw=Y, y.smooth=data[[i]]$y.hat, ext.idx=data[[i]]$i)

	if (length(data[[i]]$x) > K) {
	    ShowJT(tau=data[[i]]$x[1:K], Y[data[[i]]$i[1:K]], 0)
	} else {
	    ShowJT(tau=data[[i]]$x, Y[data[[i]]$i], 0)
	}
    
        cat("Checking time series", i, "...")
        opt_select <- readline("\n1. Press [Enter] to continue.\n2. Set window width \n3. Discard sample.\nf. Include first.\nl. Include last.\nlf. Include first AND last.\nz. Zoom.\ne. Exit\nInput: ")
	if (opt_select == 'e')
	    return(data)

	while (opt_select != "") {

	    if (opt_select == '2') {		 # Adjust width parameter
	       width_select <- ''
    	       while (!is.na(width_select)) {
    	           width_select <- as.numeric(readline("\nW = "))
    	           if (!is.na(width_select)) {
			data[[i]] <- getMM(x=X, y=Y, w=width_select, span=sp)
			JTPlot(x.raw=X, y.raw=Y, y.smooth=data[[i]]$y.hat, ext.idx=data[[i]]$i)
			if (length(data[[i]]$x) > K) {
			    ShowJT(tau=data[[i]]$x[1:K], Y[data[[i]]$i[1:K]], 0)
			} else {
			    ShowJT(tau=data[[i]]$x, Y[data[[i]]$i], 0)
			}
    	           }
    	       }
	    }


	    if (opt_select == 'f') {
	       data[[i]]$x <- c(x[1], data[[i]]$x)
	       data[[i]]$i <- c(1, data[[i]]$i)
	       ShowJT(tau=data[[i]]$x, Y[data[[i]]$i], 0)
	    } 

	    if (opt_select == 'l') {
	       data[[i]]$x <- c(data[[i]]$x, x[length(x)])
	       data[[i]]$i <- c(data[[i]]$i, length(x))
	       ShowJT(tau=data[[i]]$x, Y[data[[i]]$i], 0)
	    }

	    if (opt_select == 'lf') {
	       data[[i]]$x <- c(x[1], data[[i]]$x, x[length(x)])
	       data[[i]]$i <- c(1, data[[i]]$i, length(x))
	       ShowJT(tau=data[[i]]$x, Y[data[[i]]$i], 0)
	    }

	    if (opt_select == 'z') {
		to <- min(which(data[[i]]$x < X[1000]))
		JTPlot(x.raw=X[1:1000], y.raw=Y[1:1000], y.smooth=data[[i]]$y.hat[1:1000], ext.idx=data[[i]]$i)
		if (length(data[[i]]$x) > 50) {
		    ShowJT(tau=data[[i]]$x[1:50], Y[data[[i]]$i[1:50]], 0)
		} else {
		    ShowJT(tau=data[[i]]$x, Y[data[[i]]$i], 0)
		}
	    }

	    if (opt_select != '3')
		opt_select <- readline("\n1. Press [Enter] to continue.\n2. Set window width \n3. Discard sample.\nf. Include first.\nl. Include last.\nlf. Include first AND last.\nz. Zoom.\ne. Exit\nInput: ")

            if (opt_select == '3') {	         # Discard the sample
	       data[[i]] <- NULL
	       print("Sample discarded!. Press [Enter] to continue.")
	       opt_select = ""
	    }

	    saveRDS(object=data, file="DATA.RData")
	}

    }

    return(data)
}

#' @brief Displays jump times visually. Errors if no jump times found
#' @dxs Pruning parameter - all jump times with displacement < dxs are not considered
#' @tau Jump times
#' @y Min/max
#' 
ShowJT <- function(tau, y, dxs) {

    y.diffs <- diff(y)
    rej <- which(abs(y.diffs) < dxs) 
    if (length(rej) > 0) {
	tau <- insert(vect=tau, pos=rej+1, elems=NA)
    	y <- insert(vect=y, pos=rej+1, elems=NA)
    } 

    if (length(rej) != length(tau)) {

	for (k in seq(length(tau) - 1)) {
    	    yloc <- min(y[k],y[k+1])
    	    #text(labels = TeX(paste0("$\\tau_{",k,"}$ = ", format(tau[k+1]-tau[k], digits=1))), x = (tau[k] + tau[k+1])/2, y=yloc)
    	    #text(labels = TeX(paste0("$\\tau_{",k,"}$")), x = (tau[k] + tau[k+1])/2, y=yloc, offset=0.001)

    	    #print(yloc)
    	    arrows(x0 = tau[k], x1 = tau[k+1], y0 = yloc, y1 = yloc, angle = 90, code=3, col='black', length=0.05, lwd=2)

    	    if (!is.na(tau[k]) && !is.na(tau[k+1])) {
		lines(c(tau[k], tau[k]), c(y[k], yloc), col='black', lty='dashed', lwd=2)
		lines(c(tau[k+1], tau[k+1]), c(y[k+1], yloc), col='black', lty='dashed', lwd=2)
    	        #lines(c(tau[k], tau[k]), c(-1e7, y[k]), lty='dashed', col='black')
    	        #lines(c(tau[k+1], tau[k+1]), c(-1e7, y[k+1]), lty='dashed', col='black')
    	        #points(tau[k], y[k], pch=19, col='black')
    	        #points(tau[k+1], y[k+1], pch=19, col='black')
    	    }

    	}

    	#legend('topleft', 
        #   legend=c(TeX("$\\tau_i$ - ith jump time"), 
    	#	TeX("$x(t)$ - Raw time series"), 
    	#	TeX("Smoothed time series"), 
    	#	TeX("Turning point")), 
        #   col=c('magenta', col='red', col='blue', col='black'), lty=c('solid', 'solid', 'solid', NA), pch=c(NA,NA,NA,19))

    } else {
	cat("All jump times have been pruned!")
	legend('topright', TeX(paste0("No jump times found.", "$\\Delta x^* = $", dxs)))
    }
}


#' @brief Plots raw data (@x, @y), smoothed raw data, and location of turning points
#' @ext.idx Vector of indices which return extreme values of (@x.raw @y.smooth)
#' @show.all Boolean to determine if all max/min found should be shown
JTPlot <- function(x.raw, y.raw, y.smooth, ext.idx, show.all=FALSE) {

    matplot(x.raw, y.raw, lwd=3,
    	type='l', lty='solid', 
    	ylab=TeX("$x(t)$"), 
    	xlab=TeX("$t$"),
    	col='red')
    lines(x.raw, y.smooth, col='blue', lty='solid', lwd=2)

    # This shows all the identified jump times. 
    # Many jump times are often pruned out, so by default we don't show them.
    if (show.all) {
	sapply(ext.idx, function(i) lines(c(x.raw[i],x.raw[i]), c(-1e6, y.smooth[i]), col='black', lty='dashed'))
    	points(x.raw[ext.idx], y.smooth[ext.idx], pch=19, col='black', cex=1)
    }
}

#' @brief 
#' @x.raw Vector of observation times
#' @y.raw Matrix of observtions
#' @save Save each time series to a pdf (one pdf, many pages)
#' 
PlotSample <- function(x.raw, y.raw, id, dxs, sp, w, save=FALSE) {

    if (ncol(data.matrix(x.raw)) == 1) { # Duplicate columns
	x.raw <- matrix(x.raw, nrow=length(x.raw), ncol=ncol(y.raw))	
    }

    y.raw <- data.matrix(y.raw[,id])
    if (save) {
	fname <- ifelse(length(id) == 1, paste0("JTSample_", id, ".pdf"), 
					 paste0("JTSample_", id[1], "-", id[length(id)], ".pdf"))
	cat("Saving", ncol(y.raw), "pages to", fname, "\n")
	pdf(fname, width=10, height=10)
    } else {
	y.raw <- data.matrix(y.raw[,1])
    }

    for (i in seq(ncol(y.raw))) {
	mm <- getMM(x = x.raw[,i], y = y.raw[,i], w = w, sp = sp)
    	spar(cx=1.25)
    	JTPlot(x.raw = x.raw[,i], y.raw = y.raw[,i], y.smooth = mm$y.hat, ext.idx = mm$i, show.all=FALSE)
    	ShowJT(tau = mm$x, y = mm$y.hat[mm$i], dxs=dxs)
    }

    if (save) {
	dev.off()
    }
}

#' @brief Applies a rolling window to a curve (x,y) 
#' @w Window half width
#' @span \in (0,1) is the 'degree of smoothing'. Small @span = tighter fit
#'
argmax <- function(x, y, w, span) {
    #y.smooth <- loess(y ~ x, span=span)$fitted  # Smoothing function 
    y.smooth <- y			        #
    y.max <- zoo::rollapply(zoo::zoo(y.smooth), width=2*w+1, FUN=max, align="center")
    delta <- y.max - y.smooth[-c(1:w, length(y)+1-1:w)]
    i.max <- which(delta <= 0) + w
    return(list(x=x[i.max], i=i.max, y.hat=y.smooth)) # To be wrapped in another list
}

#' @brief Uses @argmax() to find the minima/maxima of a curve
#' @w Window half width
#' @span \in (0,1) is the 'degree of smoothing'. Small @span = tighter fit
#'
getMM <- function(x, y, w, span) {
    Max <- argmax(x=x, y , w=w, span=span)
    Min <- argmax(x=x, -y, w=w, span=span)
   
    MM <- c(Max$x, Min$x)
    SortedMinMax <- sort(MM) # MinMax in ascending order
    IDMM <- order(MM)


    if (length(Max$x) != length(Min$x)) {
	Minmm <- min(length(Max$x), length(Min$x))

	Max$x <- Max$x[1:Minmm]
	Max$i <- Max$i[1:Minmm]

	Min$x <- Min$x[1:Minmm]
	Min$i <- Min$i[1:Minmm]
    }


    # Central object for this study
    ext <- list(x=SortedMinMax, 
		i=c(Max$i, Min$i)[IDMM], 
		y.hat=Max$y.hat, 
		w=w, 
		span=span)
    return(ext)
}

#' @brief Removes trajectories which don't satisfy max(|@X|) < @thresh
#' 
RemoveConfined <- function(X, thresh=1000) {
    IsUnderThresh <- apply(X, 2, function(x) ifelse(max(abs(x)) < thresh, TRUE, FALSE))
    IsConfined <- which(IsUnderThresh == TRUE)    
    nConfined <- length(IsConfined)

    cat("Removing", nConfined, "trajectories which satisfy max(|X|) <", thresh)
    return(X[,-IsConfined])
}


#' @brief Finds min/max for a fixed window width
#' @return A list showing x-coords of min/max 
#' 
AutoLoop <- function(x, y, w, sp) { 

    data <- vector(mode='list', length=ncol(y)) # Turning points, widths, etc.

    if (ncol(data.matrix(x)) == 1) { # Duplicate columns
	x <- matrix(x, nrow=length(x), ncol=ncol(y))	
    }

    for (i in seq(ncol(y))) {
    	data[[i]] <- getMM(x=x[,i], y=y[,i], w=w, span=sp)
    }

    return(data)

}

#' @brief Removes jump times not satisfying a jump DISTANCE of @dxs
#' @tau Vector location of min/max on the horiz. axis
#' @y   Vector location of min/max on the verti. axis
#' @return Vector of jump times
#' 
PruneJT <- function(tau, y, dxs) {
    y.diffs <- diff(y)
    rej <- which(abs(y.diffs) < dxs)
    
    if (length(rej) != 0) { # Some jumps are too small
	new.tau <- insert(vect=tau, pos=rej+1, elems=NA)
	jt <- diff(new.tau)
    	jt <- jt[-which(is.na(jt))]
    } else {
	jt <- diff(tau)
    }

    return(jt)
}

#' @brief Find ballistic speed o
#'
getDiffs <- function(tau, y, dxs) {
    y.diffs <- diff(y)
    rej <- which(abs(y.diffs) < dxs)
    
    if (length(rej) != 0) { # Some jumps are too small
	new.tau <- insert(vect=tau, pos=rej+1, elems=NA)
	new.y <- insert(vect=y, pos=rej+1, elems=NA)

	jt <- diff(new.tau)
    	#jt <- jt[-which(is.na(jt))]

	y <- diff(new.y)
    	#y <- y[-which(is.na(new.y))]
    } else {
	jt <- diff(tau)
	y <- diff(y)
    }

    return(list(x=jt, y=y))
}


#' @brief Finds jump times for severals @dxs
#' @data Central object of this study
#' @dxs Vector of minimum absolute jump displacements
#' @return NA padded matrix of jump times
#' 
getJumpTimes <- function(data, dxs) {
    jt.dxs <- vector(mode='list', length=length(dxs))
    for (j in seq(length(dxs))) {
        jt <- c()

	
        for (i in seq(length(data))) {
	    if (!is.null(data[[i]])) {
		tau <- data[[i]]$x
	    	id <- data[[i]]$i
            	y <- data[[i]]$y.hat[id]
	    }
            
	    njt <- PruneJT(tau=tau, y=y, dxs=dxs[j])
            jt <- c(jt, njt)
	}
        jt.dxs[[j]] <- jt
    }
    return(list.cbind.na(jt.dxs))
}


# Save jump times for later use
setwd("../Studies/S-HPO/Data/")
files <- paste0("p0-a", seq(15,17), "-r1--X.dat")[2]  # IRRATIONALS pi/3, pi/2, 2pi/3
files <- paste0("p0-a", seq(4,6), "-r1--X.dat")[2]    # RATIONALS pi/3, pi/2, 2pi/3
pref <- "JTi-r1-"
for (f in seq(length(files))) {
    t <- seq(1,1e7+1,length=10001)
    X <- tpose(fread(files[f]))
    X <- data.matrix(RemoveConfined(X, 100))
    X <- data.matrix(X)
    data <- AutoLoop(x = t[], y = X[,], w = 10, sp=0.005) # Fixed window width auto-inspect
    dxs <- 100
    jt <- getJumpTimes(data, dxs=dxs)
    fwrite(jt, file=paste0(pref, f, ".dat"), sep=" ", col.names=FALSE) 
    #rm(data)
}

# Load jump times
jt <- list()
for (i in seq(length(files))) 
    jt[[i]] <- data.matrix(fread(paste0(pref, i, ".dat")))
jt <- list.cbind.na(jt)


# ============ STARTING DATA ============
setwd("../Studies/S-HPO/Data/")
{
    t <- seq(1,1e7+1,length=10001)
    X <- tpose(fread("p0-a5-r1--X.dat"))[,1:1000]
    M <- 5000
    #DX <- data.matrix(RemoveConfined(TSD, 5000))
    #t <- apply(tpose(fread("p0-a5-r31--IHFT.dat")), 2, cumsum)
}

W <- 3
S <- seq(2000)
data <- CommandLoop(x = t, y = X, w = W)   # For manual inspection
data <- AutoLoop(x = t[], y = X[,], w = 4, sp=0.005)	     # Fixed window width auto-inspect


DATA <- readRDS("DATA.RData")
data <- CommandLoop(x = t, y = X, w = W, Data=DATA, I=168)   # For manual inspection


#' TESTING REPRODUCABILITY
jt <- list()
k <- c(1, 2000, 4000, 6000)
for (i in seq(2, length(k))) {
    I <- seq(k[i-1], k[i])
    X <- tpose(fread("p0-a5-r1--X.dat"))[,I]

    data <- AutoLoop(x = t[], y = X[,], w = 4, sp=0.005)
    jt[[i-1]] <- getJumpTimes(data, dxs=100)
}


{ #' PLOTTING SAMPLES
    #' Plotting samples of time series with jump times
    r <- rindex(ncol(X), 0)
    r <- 2608
    s <- seq(1,10000)
    PlotSample(x.raw = t[s], y.raw = X[s,], id = r, dxs = 0, w = 2, sp = 0.1, save = FALSE)
    
    # Random sample for a subset of time
    PlotSample(x.raw = t[s], y.raw = DX[s,3366:3377], id = rindex(10, 0), dxs = 1000, w = W, sp = 0.1, save = FALSE)
    
    # Just a random sample 
    PlotSample(x.raw = t, y.raw = X, id = rindex(ncol(X), 0), dxs = 0, w = 3, sp = 0.01, save = FALSE)

    PlotSample(x.raw = t, y.raw = DX, id = rindex(ncol(DX), 9), dxs = 5000, w = W, sp = 0.01, save = TRUE)
    PlotSample(x.raw = t, y.raw = DX, id = 444, dxs = 10000, w = W, sp = 0.01, save = TRUE)
    
    id <- c(994, 126, 6847, 890, 4871, 908, 6028, 3366, 1777, 5561, 5966, 1022, 1792, 4059, 5387, 114, 1011)
    id <- c(5735, 2689, 80, 3417, 3409)
    PlotSample(x.raw = t, y.raw = DX, id = id, dxs = 500, w = W, sp = 0.01, save = TRUE)
}


#' ======================================== JUMP TIME DENSITY START ========================================
pdf("FHJTDist.pdf", width=12, height=12)
{
    dxs <- 100
    jt <- getJumpTimes(data, dxs=dxs)

    spar(cx=1.5)
    L <- getTSD_hist(dx = jt, t = 1, eta = 1, prune = FALSE, n = 200, domain=exp.domain)
    L <- getTSD_hist(dx = jt, t = 1, eta = 1, prune = FALSE, n = 200)
    matplot(L$x, L$y, type='s', log='xy', 
	lty=c('solid'),
    	#lty='solid', lwd=3, 
	lwd=3,
	col=rainbow(10),
	xlim=c(1e4, 1e7),
	#ylim=c(0,max(L$y)),
    	main=("Data for 3 separate ensembles (2000 I.C's each). r = 1/2"),
    	ylab=TeX("$\\Psi(\\tau)$"),
    	xlab=TeX("$\\tau$"))
    LF  <- LinFit(time = L$x[,1], mat = L$y[,1], range=60:120, logs=TRUE,)# show.range=TRUE)
    LF1 <- LinFit(time = L$x[,1], mat = L$y[,1], range=150:193, logs=TRUE)#, show.range=TRUE)
    lines(L$x[,1], exp(LF$b[1]) * L$x[,1]^(LF$ax[1]), col='black', lty='dashed', lwd=2)
    lines(L$x[,1], exp(LF1$b[1]) * L$x[,1]^(LF1$ax[1]), col='blue', lty='dashed', lwd=2)

    legend('topright', 
	   legend=c(TeX(paste0("$\\tau^{", format(LF$ax[1], digits=4), "}$")), 
		    TeX(paste0("$\\tau^{", format(LF1$ax[1], digits=4), "}$"))),
	   lty='dashed', col=c('black', 'blue'), cex=2, lwd=2)
    legend('bottomleft', 
	   legend=c(tex("\\alpha = \\pi / 3"), tex("\\alpha = \\pi / 2"), tex("\\alpha = 2 \\pi / 3")),
	   lty='solid', col=rainbow(10)[1:3], cex=2, lwd=2)

    legend('bottomleft', legend=paste0("Ensemble #", seq(3)), lty='solid', col=rainbow(10), cex=2, lwd=2)

    legend('topright', 
	   legend=c(TeX(paste0("$\\tau^{", format(LF$ax[1], digits=4), "}$")), 
		    TeX(paste0("$\\tau^{", format(LF1$ax[1], digits=4), "}$"))),
	   lty='dashed', col=c('black', 'red'), cex=2, lwd=2)
    legend('bottomleft', 
	   legend=c("Distribution of free paths", "From time series"),
	   lty='solid', col=c('black', 'red'), cex=2, lwd=2)

    # Fitting function and fit
    Psi <- function(tau, tau.0, mu) { ((tau.0^mu) / gamma(mu)) * (exp(-tau.0 / tau) / (tau^(1+mu))) }
    nlf <- nlsLM(y ~ Psi(x, tau.0, mu), data = df, start=list(tau.0 = dxs[1], mu = 1))
    lines(L$x, Psi(L$x, tau.0 = coef(nlf)[1], mu = coef(nlf)[2]), lty='dashed', col='black', lwd=2)

    # From stark paper
    f   <- function(x, mu) { mu / ((1 + x)^(1+mu)) } 
    df <- data.frame(x = L$x[,1], y = L$y[,1])
    nlf <- nlsLM(log(y) ~ log(f(x, mu)), data = log(df), start=list(mu = 1.5), subset=60:95)
    lines(df$x, f(df$x, mu = coef(nlf)[1]), lty='dashed', col='black', lwd=2)

    legend('topleft', legend=c(TeX(paste0("$\\tau_0 = $", format(coef(nlf)[1], digits=3), ", $\\mu = $", format(coef(nlf)[2], digits=3)))), lty='dashed', col='black', cex=2)
    legend('bottomleft', legend=TeX(paste0("$\\Delta x^* =$ ", dxs)), col=rainbow(10), lty='solid', cex=2, lwd=2)
}; 
dev.off()



#' ======================================== JUMP TIME DENSITY FINISH ========================================
#' ======================================== SPEED DISTRIBUTION START ========================================

# Extract jump locations up to different times
S <- c(2500, 5000, 7500, 10000) + 1
mx <- my <- list()
dxs <- 0
for (i in seq(length(S))) {
    s <- seq(S[i])
    data <- AutoLoop(x = t[s], y = X[s,], w = 4, sp=0.01)	     # Fixed window width auto-inspect
    Mx <- My <- c()
    for (j in seq(length(data))) {
        idx <- data[[j]]$i
	tm <- getDiffs(tau=data[[j]]$x, y = data[[j]]$y.hat[idx], dxs=dxs) # Uncomment if dxs > 0
	Mx <- c(Mx, tm$x)
	My <- c(My, tm$y)
    }
    mx[[i]] <- Mx
    my[[i]] <- My
}
mmx  <-  list.cbind.na(mx)[,]
mmy  <-  list.cbind.na(my)[,]
fwrite(mmy, file="Ys_pi2_r05.dat", col.names=FALSE, sep=" ")
fwrite(mmx, file="Xs_pi2_r05.dat", col.names=FALSE, sep=" ")

L <- getTSD_hist(data.matrix(mmy), prune = FALSE, n=200, domain=ts.exp.domain)
L <- getTSD_hist(data.matrix((mmx)), prune = FALSE,  n=150)
spar()
spar(mfrow=c(1,2))
matplot(L$x, L$y, type='s', log='', lwd=2, lty='solid', col=rainbow(10), #xlim=c(1e-6, 1),
	#xlim=c(-750000,750000),
	main=TeX("Distribution of long-jump velocities. $\\alpha = \\pi/2$, $r = 1/2$"),
	ylab=TeX("$\\log (P(|\\Delta x^*|))$"),
        xlab=TeX("$|\\Delta x^*|$"))
abline(v=c(1/2,-1/2), col='black', lty='dashed'); text(x=0.01, y=1e4, labels=TeX(paste0("$|\\Delta x^* / \\tau | = 1/2$", " is the leading front for this configuration -> ")))
abline(v=c(1/2,-1/2), col='black', lty='dashed'); text(x=0.32, y=20, labels=TeX(paste0("$|\\Delta x^* / \\tau | = 1/2$", " is leading front for this config.")))


# Rescaling attempts
T <- (S - 1) * 1000
eta <- 1 - (5/6)
eta <- 1 - 0.5
L <- getTSD_hist(data.matrix(mm), prune = FALSE, n=200, domain=ts.exp.domain, t=T, eta=eta)
L <- getTSD_hist(data.matrix(mmy), prune = FALSE, n=350, t=T, eta=eta)
spar()
matplot(L$x, L$y, type='s', log='y', lwd=2, lty='solid', col=rainbow(10), #xlim=c(1e-6, 10),
	xlim=c(-75000,75000),
	main=TeX("Rescaled distribution of long-jump velocities. $\\alpha = \\pi/2$, $r = 1/2$"),
	ylab=TeX(paste0("$\\log(T^{(1 - 5/6)} \\cdot P(| \\Delta x^* |) )$")),
        xlab=TeX(paste0("$( | \\Delta x^* |) \\cdot T^{-(1 - 5/6)}$")))
	#ylab=TeX(paste0("$\\log(T^{-", eta, "} \\cdot P(| \\Delta x^* / \\tau |) )$")),
        #xlab=TeX(paste0("$( | \\Delta x^* / \\tau |) \\cdot T^{+", eta, "}$")))

legend('topleft', legend=c(TeX("$\\tau$ - Time passed during one jump"), TeX("$\\Delta x^*$ - Displacement after one jump"))) 
legend('topright', legend=paste0("T = ", (S-1)*1000), col=rainbow(10), lty='solid')

#' ======================================== SPEED DISTRIBUTION FINISH ========================================

LF  <- LinFit(time = L$x[,3], mat = L$y[,3], range=160:190, logs=TRUE, show.range=TRUE)
lines(L$x[,3], exp(LF$b[1]) * L$x[,3]^(LF$ax[1]), col='black', lty='dashed', lwd=2)




















#' ======================== LIVE DEMO START ========================
setwd("../Studies/S-HPO/Data/")
id <- c(46, 50, 58, 59, 61, 62, 64, 66, 138, 152)
t <- seq(1,1e7+1,length=10001)
X <- tpose(fread("p0-a5-r1--X.dat"))[,id]

data <- CommandLoop(x = t, y = X, w = 50, K=10000, sp=0.005)   # For manual inspection

#' ======================== LIVE DEMO END ========================

