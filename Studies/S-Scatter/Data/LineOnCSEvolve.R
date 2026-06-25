




B <- data.matrix(fread(lf("--BID.dat")))
H <- data.matrix(fread(lf("--H.dat")))
A <- data.matrix(fread(lf("--A.dat")))
X <- data.matrix(fread(lf("--R.dat")))
L <- data.matrix(fread(lf("--L.dat")))
dim(H)


mplot(H[,1], rep(1, nrow(H)), type='p',  cex=0.1, col='blue', ylim=c(ncol(H)-100,ncol(H)+1),)

mplot(H[,1], rep(1, nrow(H)), type='p',  cex=0.01, col='blue', ylim=c(0, ncol(H)))
mplot(H[,1], rep(1, nrow(H)), type='p',  cex=0.01, col='blue', ylim=c(100,210))
for (i in seq(2,200, by=1)) {
    id <- which(A[,i] == A[1,1] & A[,i-1] == A[1,1])
    #mplot(H[,i], rep(i, nrow(H)), type='p', cex=0.01, col='blue', add=TRUE)
    mplot(H[id,i], rep(i, length(id)), type='p', cex=0.01, col='blue', add=TRUE)
}




d <- 1
h <- 0.199
theta <- -0.0713
h <- 0.286496180415440004907
theta <- -0.066575290173549184214


#Notice that the following sequence may sometimes not be in our desired interval for ballistic fronts. However, there is a greater probability that it will be in our desired interval since 

t3 <- function(h,a) { ((3/2) + 2*d) * sec(a) }
t4 <- function(h,a) { (d + 1) * sec(a) }
t6 <- function(h,a) { (1 + 2*d) * sec(a) }
t19 <- function(h,a) { sec(a) }

6/(t6(h,theta) + t4(h,theta) + 4*t19(h,theta))


cbind(seq(1,14), 
c(1/t6(h,theta), 
  2/(t6(h,theta) + t3(h,theta)), 
  3/(2*t6(h,theta) + t3(h,theta)),
  4/(3*t6(h,theta) + t3(h,theta)),
  5/(3*t6(h,theta) + t3(h,theta) + t4(h,theta)),
  6/(4*t6(h,theta) + t3(h,theta) + t4(h,theta)),
  7/(5*t6(h,theta) + t3(h,theta) + t4(h,theta)),
  8/(6*t6(h,theta) + t3(h,theta) + t4(h,theta)),
  9/(6*t6(h,theta) + 2*t3(h,theta) + t4(h,theta)),
  10/(7*t6(h,theta) + 2*t3(h,theta) + t4(h,theta)),
  11/(8*t6(h,theta) + 2*t3(h,theta) + t4(h,theta)),
  12/(8*t6(h,theta) + 2*t3(h,theta) + 2*t4(h,theta)),
  13/(9*t6(h,theta) + 2*t3(h,theta) + 2*t4(h,theta)),
  14/(10*t6(h,theta) + 2*t3(h,theta) + 2*t4(h,theta))))



c(2/(t6(h,theta) + t3(h,theta)), 9/(6*t6(h,theta) + 2*t3(h,theta) + t4(h,theta)), 16/(11*t6(h,theta) + 3*t3(h,theta) + 2*t4(h,theta)), 23/(16*t6(h,theta) + 4*t3(h,theta) + 3*t4(h,theta)))
2*7/(2*(5*t6(h,theta) + t3(h,theta) + t4(h,theta)))
1/t3(h,theta)


# Speed of w_2
7/(5*t6(h,theta) + t3(h,theta) + t4(h,theta))
3/(2*t6(h,theta) + t3(h,theta))
10/(7*t6(h,theta) + 2*t3(h,theta) + t4(h,theta))

# Horizontal line bounded by region 3:

d <- 1/2

# Constant angle
a <- -0.07129

xi6 <- -atan((h+d)/5*d)
A <- -d*tan(a) 
xi7 <- -atan(h)
B <- -tan(a)

H <- seq(A, B, length=10)

Sh <- function(const) { (1/2) + const + (5/2)*tan(a) }
Sh <- function(const) { (1/2) + const + (5/2)*tan(a) }

Sh(A)
Sh(B)
Sh(B) - Sh(A)
B - A

4*tan(a) + 1/2
(7/2) * tan(a)+ 1/2


4*tan(a) + 1/2 + 2*tan(a)
(7/2) * tan(a)+ 1/2 + 2*tan(a)


# Analysing periodic sequences
B <- data.matrix(fread(lf("--BID.dat")))
H <- data.matrix(fread(lf("--H.dat")))
A <- data.matrix(fread(lf("--A.dat")))
X <- data.matrix(fread(lf("--R.dat")))
L <- data.matrix(fread(lf("--L.dat")))

mplot(L[1,], X[1,]/L[1,], type='l')

B[1,-1]
tail(L[1,])

is_periodic <- function(vec, sequence, k=1) {
    len_seq <- length(sequence)

    # Iterate through the vector with a step size equal to the length of the sequence
    for (i in seq(1, length(vec) - len_seq + 1, by = len_seq)) {
	v <- as.vector(vec[i:(i + len_seq - 1)])
	if (!identical(v, sequence)) {
	    if (i-15 > 0) {
		print(vec[(i-8):(i+8)])
	    } else  {
		print(vec[(i):(i+15)])
	    }
	    #print(B[k,1:(i+25)])
	    print(v)
	    print(sequence)
	    print(paste0("Diverge time = ", L[k,i]))
	    print(i)
	    return(FALSE)
	}

    }

    return(TRUE)
}

k <- 1
propo <- vector(mode="integer", length=7)
propo[1:7] <- c(6L,3L,6L,6L,4L,6L,6L)

is_periodic(B[k,-1],        propo,k)
is_periodic(B[k,-(1:6388)], propo,k)
is_periodic(B[k,-(1:17325)], propo,k)
is_periodic(B[k,-(1:28262)], propo,k)
is_periodic(B[k,-(1:39199)], propo,k)
L[1,6385] # (6385 - 1) / 6 = 1064 = N
(6385 - 2 + 1) / 7
is_periodic(B[k,-(1:6388)], propo,k)
L[1,17322] - L[1,6389]  
(17322 - 6389 + 1)/7
is_periodic(B[k,-(1:17325)],propo,k)
L[1,28259] - L[1,17326] 
(28259 - 17326+1) / 7
is_periodic(B[k,-(1:28262)],propo,k)
L[1,39196] - L[1,28263] 
(39196 - 28263+1) / 7
is_periodic(B[k,-(1:39199)],propo,k)
(50140-39200+1) / 7
L[1,50140] - L[1,39200]

k <- 2
propo <- vector(mode="integer", length=7)
propo[1:7] <- c(6L,6L,12L,6L,6L,7L,6L)
is_periodic(B[k,-1], propo,k)
L[2,4551]
(4551 -2 + 1)/7
is_periodic(B[k,-(1:4554)], propo,k)
L[2,15488] - L[2,4551] 
(15488 - 4555 + 1)/7
is_periodic(B[k,-(1:15491)], propo,k)
(26425-15492 + 1) / 7
L[2,26425] - L[2,15492] 
is_periodic(B[k,-(1:26428)], propo,k)
(37362-26429 + 1) / 7
is_periodic(B[k,-(1:37365)], propo,k)
(48299-37366 + 1) / 7
is_periodic(B[k,-(1:48302)], propo,k)
(59236-48303 + 1) / 7
is_periodic(B[k,-(1:59239)], propo,k)
(70173 - 59240 + 1) / 7
is_periodic(B[k,-(1:70176)], propo,k)
(70177  + 1) / 7



