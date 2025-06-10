#####################
# Supplementary material for
#   Detecting (non)parallel evolution in multidimensional spaces: 
#   angles, correlations, and eigenanalysis
# by Junya Watanabe
# Version 2021.09.29.
#
# Description:
#   This file includes example R functions concerning
#   the distribution of the angle between a pair of independent vectors
#   that has no preferred directions in the k-dimensional space.
#   Written on R version 3.5.3.
#
# Usage:
#   Similar to the regular R functions for probability distribution,
#   e.g., pnorm(), qnorm(), rnorm(), dnorm().
#####################

#### pang ####
## Distribution function of the random angle.
## This is from the t distribution corresponding to the angle.
pang <- function(q, k = 2, lower.tail = TRUE) {
    T <- sqrt(k - 1) / tan(q)
    pt(T, k - 1, lower.tail = !lower.tail)
}

#### qang ####
## Inverse of distribution function of the random angle.
## Assumed range is [0, pi]
qang <- function(p, k = 2, lower.tail = TRUE) {
    qs <- qt(p, k - 1, lower.tail = !lower.tail)
    ifelse(qs >= 0, atan(sqrt(k - 1) / qs), atan(sqrt(k - 1) / qs) + pi)
}

#### rang ####
## Generate random angles.
## Assumed range is [0, pi].
rang <- function(n, k = 2) {
    rs <- rt(n, k - 1)
    ifelse(rs >= 0, atan(sqrt(k - 1) / rs), atan(sqrt(k - 1) / rs) + pi)
}

#### dang ####
## Density of the random angle.
## Assumed range is [0, pi].
dang <- function(x, k = 2) {
    S <- sin(x) ^ (k - 2)
    b <- beta(1 / 2, (k - 1) / 2)
    return(ifelse((x >= 0) & (x <= pi), S / b, NA))
}
