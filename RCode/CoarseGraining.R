

img.path = "~/Dropbox/PHD/Papers/JumpStatistics/Images/"
setwd("../Studies/S-JumpStatistics/Data/")

alpha <- c(pi/seq(9, 2, by=-1), seq(2, 8)*pi/seq(3,9) )
alpha.lab <- c(tex("\\pi/9"), 
	       tex("\\pi/8"), 
	       tex("\\pi/7"), 
	       tex("\\pi/6"), 
	       tex("\\pi/5"), 
	       tex("\\pi/4"), 
	       tex("\\pi/3"), 
	       tex("\\pi/2"), 
	       tex("2\\pi/3"), 
	       tex("3\\pi/4"), 
	       tex("4\\pi/5"), 
	       tex("5\\pi/6"), 
	       tex("6\\pi/7"), 
	       tex("7\\pi/8"), 
	       tex("8\\pi/9"))
alpha.lab.irr <- c(tex("(\\varphi - 3/5) \\pi/9"), 
	       tex("(\\varphi - 3/5) \\pi/8"), 
	       tex("(\\varphi - 3/5) \\pi/7"), 
	       tex("(\\varphi - 3/5) \\pi/6"), 
	       tex("(\\varphi - 3/5) \\pi/5"), 
	       tex("(\\varphi - 3/5) \\pi/4"), 
	       tex("(\\varphi - 3/5) \\pi/3"), 
	       tex("(\\varphi - 3/5) \\pi/2"), 
	       tex("(\\varphi - 3/5) 2\\pi/3"), 
	       tex("(\\varphi - 3/5) 3\\pi/4"), 
	       tex("(\\varphi - 3/5) 4\\pi/5"), 
	       tex("(\\varphi - 3/5) 5\\pi/6"), 
	       tex("(\\varphi - 3/5) 6\\pi/7"), 
	       tex("(\\varphi - 3/5) 7\\pi/8"), 
	       tex("(\\varphi - 3/5) 8\\pi/9"))


# Step times and step distances
{
    cnt     <- ReadDataList(lff("", "w500_jump.dat"), idx=1) # Index of jumps along each trajectory. 
    t0      <- ReadDataList(lff("", "w500_jump.dat"), idx=2) # Starting time of a jump
    x0      <- ReadDataList(lff("", "w500_jump.dat"), idx=3) # Starting position of a jump
    jt      <- ReadDataList(lff("", "w500_jump.dat"), idx=4) # Jump time
    jd      <- ReadDataList(lff("", "w500_jump.dat"), idx=5) # Jump displacement

    id <- apply(cnt, 2, function(x) which(x == 0))              # diff(id) = number of jumps per trajectory
    num.ics <- apply(cnt, 2, function(x) length(which(x == 0))) # Number of initial conditions 
    ranges <- lapply(list(id), function(x) sapply(1:(length(x)-1), 
			       function(y) seq(x[y], x[y+1]-1)))

    input <- list(NumberOfInitialConditions=num.ics, cbind(id[-1,], diff(id)))
    cat(paste0("\nNumber of initial conditions = ", num.ics, "\n"))
    print(data.frame(StartIndexOfJump=id[-1,], NumJumpsPerIC=diff(id)), max = 20)
}



### 
#   Visualising jump detection
###

# Load trajectory data
i <- 3
x <- data.matrix(fread(lf(paste0("w50_", i, "--x.dat"))))
t <- data.matrix(fread(lf(paste0("w50_", i, "--t.dat"))))

# Plot trajectory
s <- seq(1000000)
mplot(t[s], x[s,], xlab="t", ylab="x(t)")

# Plot jumps
ShowJT(t0=t0[,], x0=x0[,], tau=jt[,], jd=jd[,], cex=0.5)



