

pdf("THETAnm.pdf", width=14, height=18)
{
    spar(mfrow=c(2,1))
    m.sample <- c(1,2,10,100)
    tmn <- ReadData(lf("THETAn_"))[,1,m.sample]
    tmn[which(tmn == 0, arr.ind=TRUE)] <- NA
    mplot((tmn %% (2*pi)) * 180/pi, yaxt='n', xaxt='n', ylim=c(220, 330), spa=FALSE,
	  ylab=tex("< \\theta_n \\theta_{n+m} > \\ mod \\ 2\\pi"), xlab="n")
    axis(side=2, at=seq(0, 360, by=30))
    axis(side=1, at=seq(0, 500, by=50))
    #legend('bottomright', legend=paste0("m = ", m.sample), lwd=2, cex=2, col=rbow())

    # Converged value as a function o m
    tmn <- ReadData(lf("THETAn_"))[,1,]
    tmn[which(tmn == 0, arr.ind=TRUE)] <- NA
    m <- seq(100)
    limits <- (apply(tmn[-(1:50),], 2, function(x) mean(x, na.rm=TRUE)) %% (2*pi)) * 180/pi
    mplot(m, limits, type='o', xlab="m", ylab=tex("\\bar{< \\theta_n \\theta_{n+m} >} \\ mod \\ 2\\pi"), spa=FALSE)
    abline(h=270, , lty='dashed')
}
