
m <- data.matrix(fread(lf("Scatter"),fill = TRUE))
# Re-assign edge ID's
{
    m <- m[,-1]
    m[which(m == 2, arr.ind=TRUE)] <- NA
    m[which(m == 3, arr.ind=TRUE)] <- 2
    m[which(m == 4, arr.ind=TRUE)] <- 3
}

apply(m, 1, function(x) 
      paste0("L", paste(x[-which(is.na(x))], collapse=""), "R"))

X <- apply(m, 1, function(x) paste(x[-which(is.na(x))], collapse=""))

# Remove bouncing ball itineraries?
uX <- unique(X)

# Construct bouncing ball itenaries
BB1 <- sapply(1:10, function(x) paste0(rep("21", x), collapse=""))
BB2 <- sapply(1:10, function(x) paste0(rep("03", x), collapse=""))

which(uX == intersect(uX, BB1))

length(unique(X)) + 2 # Account for removed bouncing balls
