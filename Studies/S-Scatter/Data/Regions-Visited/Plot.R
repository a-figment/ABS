
setwd("Data/Regions-Visited/")
X <- (fread(lf("R.dat")))
H <- tpose(fread(lf("--H.dat")))
A <- tpose(fread(lf("--A.dat")))*180/pi
B <- tpose(fread(lf("--BID.dat")))

BHA <- data.matrix(fread(lf("BHA.dat")))

{
mplot(BHA[,2], BHA[,3], type='p', cex=1, xaxt='n', yaxt='n', #xlim=c(0, 45),
      xlab=tex("\\theta_0"), ylab="Number of regions visisted",)
axis(side=1, at=seq(0,360, by=45), labels=seq(0,360, by=45))
axis(side=2, at=seq(0,10, by=1))#, labels=seq(0,360, by=45))
}
abline(v=180)

BHA[which(BHA[,3] == 7),]
hist(BHA[which(BHA[,3] == 8),2], breaks=200)

mplot(X[,1:10], type='l')
mplot(A,  type='o')
mplot(X[1:500,],A[1:800,], type='o')

mplot(X[-nrow(X),], A[-nrow(X),], type='o')

mplot(H[1:400,],)

df <- data.frame(x=1:(nrow(X)-1), y=X[-nrow(X),], z=A[-nrow(X),])
plot3d(x=df$x, y=df$y, z=df$z, type='p')

x <- data.matrix(fread(lf("--x.dat")))
y <- data.matrix(fread(lf("--y.dat")))
PlotChannelBoundary(xVert, yVert, from=-2, to=2)
lines(x,y, asp=1)
mplot(x)


