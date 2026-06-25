
{
    alpha <- pi/3
    dx <- 1/2
    dy <- (1/2)*(1/(tan(alpha/2)))
    d <- 3.11
    PlotChannelBoundary(xVert=c(-dx,0,dx,dx,0,-dx), yVert=c(0,dy,0,d,d+dy,d), from=-1,to=1
    Traj <- data.matrix(fread(lf("HPO")))
    lines(Traj[,1], Traj[,2])
}

mplot(Traj[,1], Traj[,2])

r <- seq(1/5, 10, length=981)
vb <- data.matrix(fread(lf("VB")))
n <- seq(min(nrow(vb)))

mplot(r[n], vb[n])
mean(vb)
sqrt(var(vb))
