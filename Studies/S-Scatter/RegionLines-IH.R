# ========================= INFINITE HORIZON D = 1 ============================
{
    a0 <- function(d, h0) { rep(rtd * -pi/2, length(h0)) }
    a1 <- function(d, h0) { rtd*(-atan(h0/d)) }
    PROP1 <- function(d,h0,K) { rtd*-atan( (h0 + K*d) / (1 + d*K)) }
    PROP1.id <- function(d,h0,K) { x <- PROP1(d,h0,K); which(x > a1(d,h0) & x < phi0(d,h0)[PHI0.id(d,h0)]); }
    phi0 <- function(d, h0)  { rtd*(-atan((h0+d)/(4*d))) }
    phi1 <- function(d, h0)  { rtd*(-atan((h0+d)/(5*d))) }

    # Gets the range of heights which is valid
    PHI0.id <- function(d,h0) { x <- phi0(d,h0); which(x > a1(d,h0) & x < a2(d,h0)); }
    PHI1.id <- function(d,h0) { x <- phi1(d,h0); which(x > a1(d,h0) & x < a2(d,h0)); }

    #a2 <- function(d, h0) { rtd*(-atan((h0 + d - oh)/(dx2 + d - oh))) }

    a3 <- function(d, h0) { rtd*(-atan(h0/(dx2 + d))) }
    a4 <- function(d, h0) { rtd*(-atan(h0/(dx2 + 2*d))) }
    a5 <- function(d, h0) { rtd*atan((d-h0)/(dx2 + 2*d)) }
    a6 <- function(d, h0) { rtd*atan((dy-h0)/(dx2 + 2*d - oh)) }

    a7 <- function(d, h0) { rtd*atan((dy-h0)/(dx2 + d - oh)) }
    phi2 <- function(d, h0)  { rtd*(atan((dx2 - h0)/(dx2 + 2*d))) }
    phi3 <- function(d, h0)  { rtd*(atan((dx2 - h0)/(dx2 + d))) }
    PHI2.id <- function(d,h0) { x <- phi2(d,h0); which(x > a7(d,h0) & x < a8(d,h0)); }
    PHI3.id <- function(d,h0) { x <- phi3(d,h0); which(x > a7(d,h0) & x < a8(d,h0)); }
    a8 <- function(d, h0) { rtd*atan((d-h0)/dy) }

    phi3 <- function(d, h0) { rtd*atan((dx2 - h0)/(dx2 + d)) }

    # When d = 1/2 this is the angle to the center vertex otherise it is the angle to the upper right vertex
    a12 <- function(d, h0) { rtd*atan((d - h0)/d) } 
    PROP2 <- function(d,h0,K) { rtd*atan( (d*(1 + K) - h0) / (1 + d*K)) }
    PROP2.id <- function(d,h0,K) { x <- PROP2(d,h0,K); which(x > a12(d,h0) & x < a13(d,h0)); }
    a13 <- function(d, h0) { rtd*atan((d + dy - h0)/(dx2 + d - oh)) }
    a14 <- function(d, h0) { rtd*atan((d + dy + dx - h0)/(dx2 + d)) }
    a15 <- function(d, h0) { rtd*atan((d + dx2 - h0)/dx2) }
    a16 <- function(d, h0) { rtd*atan((d + dy - h0) / dx) }

    a17 <- function(d, h0) { rtd*atan((d + dy - h0) / dx) }
    a18 <- function(d, h0) { rtd*atan((2*d + dx2 - h0)/(dx2)) }
    a19 <- function(d, h0) { rtd*atan((2*dy + 2*d - dx - h0)/dx) }
    PROP3    <- function(d,h0,K) { rtd*atan((1 + d*(1 + K) - h0) / (d*K)) }
    PROP3.id <- function(d,h0,K) { x <- PROP3(d,h0,K); which(x > a19(d,h0) & x < a20(d,h0)); }
    a20 <- function(d, h0) { rep(rtd * pi/2, length(h0)) }

    a21 <- function(d,h0) { rtd*atan((dy-h0)/dx) }
    a22 <- function(d,h0) { rtd*atan((d-h0)/dx2) }
}

DrawRegionLines <- function(d, h, K=1, ...) {

    lines(h, a1(d,h), ...)
    lines(h, a2(d,h), ...)
    lines(h, a3(d,h), ...)
    lines(h, a4(d,h), ...)
    lines(h, a5(d,h), ...)
    lines(h, a6(d,h), ...)
    lines(h, a7(d,h), ...)
    lines(h, a8(d,h), ...)
    lines(h, a12(d,h), ...)
}
