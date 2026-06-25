

epx <- tpose(fread("p0-a0-d0--EPx.dat"))
ep <- tpose(fread("p0-a0-d0--EPy.dat"))
ev <- tpose(fread("p0-a0-d0--EscapeVel.dat"))
et <- tpose(fread("p0-a0-d0--EscapeTime.dat"))

x <- tpose(fread("p0-a0-d0--2--xTraj.dat"))
y <- tpose(fread("p0-a0-d0--2--yTraj.dat"))


fz <- apply(ep, 2, function(x) min(which(x == 0))) - 1
fz[which(is.infinite(fz))] <- nrow(ep)

D <- seq(0, 3/4, length=20)
Bpm <- seq(-pi, pi, length=100)
Bp <- seq(-pi/2, pi/2, length=100)

# Associate an integer to the escape locations
i <- 3
K <- seq(fz[i])
ID <- cbind(findInterval(ep[K,i], D), findInterval(ev[K,i], Bp))


# Convert to symbol
s <- apply(ID, 1, function(x) paste0("(",paste(x, collapse=', '), ")"))
barplot(table(s))


for (j in seq(fz[i])) {
    k <- which((s[j] == s[-c(1:j)])) + j
    
    if (length(k) > 0)
	print(k)
}

t <- c(0,cumsum(sqrt(diff(x)^2 + diff(y)^2)))
matplot(t, x, type='l', xlim=c(0,500), ylim=c(-100, 50))
matplot(t, x, type='l', xlim=c(0,400), ylim=c(-10, 0))
matplot(t, x, type='l')


s <- 1:100
lines(x[s],y[s],type='l', asp=1)
abline(v=c(L,-L))
