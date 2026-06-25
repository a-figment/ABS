###
###
#   ROLLING WINDOW ILLUSTRATION
#   The operation of finding a rolling minimum and maximum for
#   some dummy data is illustrated alongside a method for 
#   determining local turning points. The illustration is done
#   by saving a series of images (frames) which show how the 
#   relevant operations are carried out.
#   WARNING: Some hardcoded values and logic sequences are present.
###

setwd("../Studies/S-JumpStatistics/Data/")

y <- c(seq(-1,4), seq(3,-2), seq(-1, 9), seq(8, 3))
x <- seq(length(y))

am <- function(w, y) {
    y.max <- rollmax(zoo(y), k=2*w, align='center')
    delta <- y.max - y[-c(1:(w-1), length(y)+1-1:w)]
    which(delta <= 0) + w - 1
}

w <- 5
s <- -c(1:(w-1), length(y)+1-1:w)
K <- length(y)
y.max.compare <- y.min.compare <- c()
idx.max <- am(w, y)
idx.min <- am(w, -y)
for (i in seq(1, K - 2*w + 1)) {

    png(paste0(i, "-moving_window.png"), width=1480, height=1080)
    spar(cx=2, col='white')
    plot(x, y, pch=19, cex=2, xlab=tex("x_k, k = 1, 2, ..., K"))

    rect(xleft=x[i], xright=x[i + 2*w - 1],
	 ytop=1e6, ybot=-1e6, col=rgb(0,0,1,alpha=0.2))

    if (i >= 2)
	points(x[idx.max[1]], y[idx.max[1]], col='red', pch=19, cex=2)
    if (i >= 19)
	points(x[idx.max[2]], y[idx.max[2]], col='red', pch=19, cex=2)
    if (i >= 8)
	points(x[idx.min], y[idx.min], col='green', pch=19, cex=2)

    j <- seq(i, i+2*w - 1)
    wm <- which.max(-y[j]) + i - 1
    wM <- which.max(y[j]) + i - 1
    y.min.compare <- c(y.min.compare, max(-y[j]))
    y.max.compare <- c(y.max.compare, max(y[j]))

    arrows(x0=x[wm]-1, x1=x[wm], y0=y[wm], y1=y[wm], angle=45, length=0.15, lwd=2, col='red')
    text(x=x[wm]-1.5, y=y[wm]+0.25, 
	 labels="Window min.",
	 cex=2, col='red')

    arrows(x0=x[wM]+1, x1=x[wM], y0=y[wM], y1=y[wM], angle=45, length=0.15, lwd=2, col='red')
    text(x=x[wM]+1.5, y=y[wM]+0.25, 
	 labels="Window max.",
	 cex=2, col='red')

    if (i == 1) {
        text(x=x[wm]+2.6, y=y[wm]-0.28, 
         labels=TeX(paste0("$m_i + y(x_{i+w-1}) = m_{", i, "} + y(x_{", i+w-1, "}) = $", y.min.compare[i] + y[i+w-1])), 
         cex=2, col='blue')
    } else {
        text(x=x[wm]-2.6, y=y[wm]-0.28, 
         labels=TeX(paste0("$m_i + y(x_{i+w-1}) = m_{", i, "} + y(x_{", i+w-1, "}) = $", y.min.compare[i] + y[i+w-1])), 
         cex=2, col='blue')
    }
    text(x=x[wM]-2.6, y=y[wM]-0.28, 
     labels=TeX(paste0("$M_i - y(x_{i+w-1}) = M_{", i, "} - y(x_{", i+w-1, "}) = $",  y.max.compare[i] - y[i+w-1])), 
     cex=2, col='blue')


    legend('topleft', legend=c(TeX("$m_i$ = ith window minimum"), 
			       TeX("$M_i$ = ith window maximum"),
			       "Local minimum",
			       "Local maximum"),
	   col=c(NA, NA, 'green', 'red'), pch=c(NA, NA, 19, 19), cex=2)

   #Sys.sleep(1)
   if(i != K - 2*w + 1)
       dev.off()
}

mm <- getMM(x = x, y=y, w=w+1)
ShowJT(df = mm$df, text=FALSE, eps=100)
mm <- getMM(x = x, y=y, w=w)
ShowJT(df = mm$df, text=FALSE, eps=100)

text(y=-2+0.3, x=9, labels=tex("\\tau_1"), cex=2)
text(y=-2+0.3, x=17, labels=tex("\\tau_2"), cex=2)

text(y=0, x=5+0.5, labels=tex("\\xi_1"), cex=2)
text(y=0, x=24-0.5, labels=tex("\\xi_2"), cex=2)
dev.off()



pdf("panel_illustration.pdf", height=18, width=12)
spar(cx=2.5, col='white', mfrow=c(3,1))
tcx <- 2.2
#par(mfrow = c(2, 2))
##par(cex = 2)
par(mar = c(8, 8, 1.3, 1.3), oma = c(2, 2, 2, 2), mgp=c(1,2,0))
#par(tcl = -0.25)
#par(mgp = c(2, 0.6, 0))

for (i in seq(1, 3)) {

    plot(x, y, pch=19, cex=2, 
	 #xlab=tex("x_k = 1, 2, ..., 30"), ylab=tex("y(x_k)"), 
	 xlab="", ylab="",
	 ylim=c(-3, 9))
    rect(xleft=x[i], xright=x[i + 2*w - 1],
	 ytop=1e6, ybot=-1e6, col=rgb(0,0,1,alpha=0.2))

    text(x = 29.5, y = 8.5, labels=i, cex=4)
    mtext(side=1, text=tex("x_k = 1,2,...,30"), line=5, cex=2)
    mtext(side=2, text=tex("y(x_k)"),		line=5, cex=2)

    if (i >= 2)
	points(x[idx.max[1]], y[idx.max[1]], col='red', pch=19, cex=2)
    if (i >= 19)
	points(x[idx.max[2]], y[idx.max[2]], col='red', pch=19, cex=2)
    if (i >= 8)
	points(x[idx.min], y[idx.min], col='green', pch=19, cex=2)

    j <- seq(i, i+2*w - 1)
    wm <- which.max(-y[j]) + i - 1
    wM <- which.max(y[j]) + i - 1
    y.min.compare <- c(y.min.compare, max(-y[j]))
    y.max.compare <- c(y.max.compare, max(y[j]))


    arrows(x0=x[wM]+2, x1=x[wM], y0=y[wM], y1=y[wM], angle=45, length=0.0, lwd=2, col='black')
    text(x=x[wM]+2.5, y=y[wM]+0.5, 
	 labels="Window max.",
	 cex=tcx, col='black')

    yoff <- -0
    xoff <- +4.5
    if (i == 1) {

        arrows(x0=x[wm]+2, x1=x[wm], y0=y[wm], y1=y[wm], angle=45, length=0.00, lwd=2, col='black')
        text(x=x[wm]+2.5, y=y[wm]+0.5, 
	     labels="Window min.",
    	     cex=tcx, col='black')

        text(x=x[i]+xoff, y=max(y)+yoff, 
	     labels=TeX(paste0("$m_i + y(x_{i+w-1}) = m_{", i, "} + y(x_{", i+w-1, "}) = $", y.min.compare[i] + y[i+w-1])), 
             cex=tcx, col='black')
    } else {

        arrows(x0=x[wm]-2, x1=x[wm], y0=y[wm], y1=y[wm], angle=45, length=0.00, lwd=2, col='black')
        text(x=x[wm]-2.5, y=y[wm]+0.5, 
	     labels="Window min.",
    	     cex=tcx, col='black')

        text(x=x[i]+xoff, y=max(y)+yoff, 
         labels=TeX(paste0("$m_i + y(x_{i+w-1}) = m_{", i, "} + y(x_{", i+w-1, "}) = $", y.min.compare[i] + y[i+w-1])), 
         cex=tcx, col='black')
    }

    text(x=x[i]+xoff, y=max(y)-1+yoff, 
	 labels=TeX(paste0("$M_i - y(x_{i+w-1}) = M_{", i, "} - y(x_{", i+w-1, "}) = $",  y.max.compare[i] - y[i+w-1])), 
     	 cex=tcx, col='black')


    legend('bottomright', legend=c(TeX("$m_i$ = ith window minimum"), 
			       TeX("$M_i$ = ith window maximum"),
			       #"Local minimum",
			       "Local maximum"),
	   col=c(NA, NA, 'red'), pch=c(NA, NA, 19), cex=tcx)
}
dev.off()
