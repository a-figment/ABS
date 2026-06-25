####' 
##   PATH: S-JumpStatistics/Data/ic_dep
##   Considering the jump time as a function of initial velocities 
##   
####' 


setwd("../Studies/S-JumpStatistics/Data/")

{
cnt     <- ReadDataList(lff("", "_jump.dat"), idx=1)
t0      <- ReadDataList(lff("", "_jump.dat"), idx=2)
x0      <- ReadDataList(lff("", "_jump.dat"), idx=3)
jt      <- ReadDataList(lff("", "_jump.dat"), idx=4)
jd      <- ReadDataList(lff("", "_jump.dat"), idx=5)

# Initial conditions
v.x    <- ReadDataList(lff("", "_vel.dat"), idx=1)
v.y    <- ReadDataList(lff("", "_vel.dat"), idx=2)
th     <- atan2(v.y, v.x)
th     <- atan(v.y/v.x)
}

i <- 1
id <- which(cnt[,i] + 1 == 1)
length(id) == nrow(th)

ranges <- sapply(1:(length(id)-1), function(x) seq(id[x], id[x+1]))



ics <- rep(th[1], length(ranges[[1]]))
tau <- jt[ranges[[1]],i]
tau <- jd[ranges[[1]],i] / jt[ranges[[1]],i]
plot(ics, tau, pch=19, cex=0.25,log='y', ylim=c(10, 1e5), xlim=c(-pi/2, pi/2))

pdf("v-vs-th0.pdf", width=14, height=10)
spar(cx=2, col='white')
plot(ics, tau, pch=19, cex=0.25, ylim=c(-0.6, 0.6), xlim=c(-pi/2, pi/2),
     xlab=tex("\\theta_0"), ylab=tex("\\xi / \\tau"))
for (j in seq(2,length(ranges))) {

    s <- seq(1, length(ranges[[j]]), by=200)
    ics <- rep(th[j], length(ranges[[j]]))
    #tau <- jt[ranges[[j]],i]    
    tau <- jd[ranges[[j]],i] / jt[ranges[[j]],i]

    points(ics[s], tau[s], pch=19, cex=0.25)

}
dev.off()

times <- (lapply(ranges, function(x) sum(jt[x,i])))



jump   <- data.matrix(fread(lf("w100_jump.dat")))


plot(jt[1:(nrow(jt)-1),i], jt[2:(nrow(jt)),i], pch=19, log='xy', cex=0.2)
plot(jt[1:(nrow(jt)-5),i], jt[6:(nrow(jt)),i], pch=19, log='xy', cex=0.2)

