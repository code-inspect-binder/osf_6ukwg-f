###############################
# Supplementary material for
#    Detecting (non)parallel evolution in multidimensional spaces:
#    angles, correlations, and eigenanalysis
#      by Junya Watanabe
# Version 2021.09.30
# This file involves utility R functions used in the re-analysis in the paper
# Partly depends on the following package:
# - gtools: in kmeans.ex() function
###############################

##### addhull #####
# Function to add convex hull (and points, if specified) on a plot
addhull <- function(x, y, groups = factor(1), hull = TRUE, points = FALSE,
                    off = NULL, lwd.h = 1, lty.h = 1, pch = 20, cex = 1,
                    lwd.p = 1, col = "black", col.h = col, ...) {
    if(missing(y)) {
        if(!is.matrix(x)) stop("Give a pair of vectors or a matrix")
        y <- x[, 2]
        x <- x[, 1]
    }
    if(dim(cbind(x, y))[1] != length(groups))
        groups <- rep(groups, length = dim(cbind(x, y))[1])
    data <- cbind(x, y)
    nc <- length(levels(groups))
    if(length(col.h) < nc) col.h <- rep(col.h, length = nc)
    if(length(lwd.h) < nc) lwd.h <- rep(lwd.h, length = nc)
    if(length(lty.h) < nc) lty.h <- rep(lty.h, length = nc)
    if(length(pch) < nc) pch <- rep(pch, length = nc)
    if(length(cex) < nc) cex <- rep(cex, length = nc)
    if(length(col) < nc) col <- rep(col, length = nc)
    if(length(lwd.p) < nc) lwd.p <- rep(lwd.p, length = nc)
    for(i in 1:nc){
        if(i %in% off) next
        X <- na.omit(subset(data, groups == levels(groups)[i], drop = FALSE))
        if(hull && dim(X)[1] > 1) {
            hp <- chull(X)
            hp <- c(hp, hp[1])
            lines(X[hp, ], lwd = lwd.h[i], lty = lty.h[i], col = col.h[i])
        }
        if(points) points(X, pch = pch[i], cex = cex[i], col = col[i],
                          lwd = lwd.p[i], ...)
    }
}

##### biplot2 #####
# Function for PCA biplot
# Modified from biplot() from the stats package
# The methods are defined below
biplot2 <- function(x, ...) {
    UseMethod("biplot2")
}

##### biplot2.prcomp #####
# prcomp method for biplot2, modified from biplot.prcomp from the stats package
# The modifications are essentially aesthetic, so the usage is the same,
# except the scaling parameter alpha (see text) can be specified.
# The original function takes the argument scale, which is 1 - alpha.
biplot2.prcomp <- function (x, choices = 1L:2L, alpha = 1, scale = 1 - alpha,
                            pc.biplot = FALSE, ...) {
    if (length(choices) != 2L)
        stop("length of choices must be 2")
    if (!length(scores <- x$x))
        stop(gettextf("object '%s' has no scores", deparse(substitute(x))),
            domain = NA)
    if (is.complex(scores))
        stop("biplots are not defined for complex PCA")
    lam <- x$sdev[choices]
    n <- NROW(scores)
    lam <- lam * sqrt(n)
    if (scale < 0 || scale > 1)
        warning("'scale' is outside [0, 1]")
    if (scale != 0)
        lam <- lam^scale
    else lam <- 1
    if (pc.biplot)
        lam <- lam/sqrt(n)
    biplot2.default(t(t(scores[, choices])/lam), t(t(x$rotation[,
        choices]) * lam), ...)
    invisible()
}

##### biplot2.default #####
# default method for biplot2, modified from biplot.prcomp from the stats package
biplot2.default <- function (x, y, var.axes = TRUE, col = c("black", "tomato"),
    cex = rep(par("cex"), 2), show.xlabs = TRUE, show.ylabs = TRUE,
    xlabs = NULL, ylabs = NULL, expand = 1, xlim = NULL, ylim = NULL,
    arrow.len = 0.1, main = NULL, sub = NULL, xlab = NULL, ylab = NULL,
    pch.xpt = 20, cex.xpt = 1, col.xpt = col[1], pos.xtx = 3,
    cex.xtx = 1, col.xtx = col.xpt,
    col.yar = col[2], lwd.yar = 1, cex.ytx = 1, col.ytx = col.yar,
    tcl_in = 0.4, line_inx = -2.5, line_iny = line_inx - 0.5,
    symmetric.limits = TRUE, asp = 1, ...)
{
    n <- nrow(x)
    p <- nrow(y)
    if (missing(xlabs)) {
        xlabs <- dimnames(x)[[1L]]
        if (is.null(xlabs))
            xlabs <- 1L:n
    }
    xlabs <- as.character(xlabs)
    dimnames(x) <- list(xlabs, dimnames(x)[[2L]])
    if (missing(ylabs)) {
        ylabs <- dimnames(y)[[1L]]
        if (is.null(ylabs))
            ylabs <- paste("Var", 1L:p)
    }
    ylabs <- as.character(ylabs)
    dimnames(y) <- list(ylabs, dimnames(y)[[2L]])
    if (length(cex) == 1L)
        cex <- c(cex, cex)
    if (missing(col)) {
        col <- par("col")
        if (!is.numeric(col))
            col <- match(col, palette(), nomatch = 1L)
        col <- c(col, "tomato")
    }
    else if (length(col) == 1L)
        col <- c(col, col)
    unsigned.range <- function(x) c(-abs(min(x, na.rm = TRUE)),
        abs(max(x, na.rm = TRUE)))
    rangx1 <- unsigned.range(x[, 1L])
    rangx2 <- unsigned.range(x[, 2L])
    rangy1 <- unsigned.range(y[, 1L])
    rangy2 <- unsigned.range(y[, 2L])
    if(symmetric.limits) {
        if (missing(xlim) && missing(ylim))
            xlim <- ylim <- rangx1 <- rangx2 <- range(rangx1, rangx2)
        else if (missing(xlim))
            xlim <- rangx1
        else if (missing(ylim))
            ylim <- rangx2
    }
    ratio <- max(rangy1/rangx1, rangy2/rangx2)/expand
    on.exit(par(op))
    op <- par(pty = "s")
    if (!is.null(main))
        op <- c(op, par(mar = par("mar") + c(0, 0, 1, 0)))
    plot(x, type = "n", xlim = xlim, ylim = ylim, col = col[1L],
        xlab = xlab, ylab = ylab, sub = sub, main = main, asp = asp, ...)
    points(x, cex = cex.xpt, col = col.xpt, pch = pch.xpt, ...)
    if(show.xlabs)
        text(x, xlabs, cex = cex.xtx, col = col.xtx, pos = pos.xtx, ...)
    par(new = TRUE)
    dev.hold()
    on.exit(dev.flush(), add = TRUE)
    plot(x / ratio, axes = FALSE, type = "n",
         xlim = if(!is.null(xlim)) xlim * ratio else NULL,
         ylim = if(!is.null(ylim)) ylim * ratio else NULL,
         xlab = "", ylab = "", col = col[1L], asp = asp, ...)
    axis(1, labels = FALSE, col = col[2L], tcl = tcl_in, ...)
    axis(1, tick = FALSE, col = col[2L], line = line_inx, font = 3,
         cex = 0.9, ...)
    axis(2, labels = FALSE, col = col[2L], tcl = tcl_in, ...)
    axis(2, tick = FALSE, col = col[2L], line = line_iny, font = 3,
         cex = 0.9, ...)
    box(col = col[1L])
    if(show.ylabs)
        text(y * 1.2, labels = ylabs, cex = cex.ytx, col = col.ytx, ...)
    if (var.axes)
        arrows(0, 0, y[, 1L], y[, 2L], col = col.yar, lwd = lwd.yar,
            length = arrow.len)
    invisible()
}

##### pcc #####
# Auxiliary function to calculate PC scores
## Arguments (all numeric)
# X: n x p raw data matrix
# R: p x k rotation matrix (eigenvectors)
# M: p vector
## Value
# n x k numeric matrix of PC scores
pcc <- function(X, R, M) {
    result <- sweep(X, 2, M) %*% R
    colnames(result) <- colnames(R)
    return(result)
}

##### rayleigh_test #####
# Rayleigh test for the uniformity of directional data,
# as described by Mardia et al. (1979: chapter 15) and Mardia & Jupp (1999)
# Note that there would probably be better alternatives on CRAN
rayleigh_test <- function(X, correction = TRUE, check = TRUE, convert = TRUE) {
    if(check) {
        if(!isTRUE(all.equal(diag(tcrossprod(X)), 1))) {
            if(convert) {
                X / sqrt(diag(tcrossprod(X)))
                warning("X was converted into directional data")
            } else {
                stop("Directional data (in coordinates) are assumed as X")
            }
        }
    }
    p <- ncol(X)
    n <- nrow(X)
    Mean <- colMeans(X)
    Norm_Mean_2 <- drop(crossprod(Mean))
    S <- p * n * Norm_Mean_2
    if(correction) {
        S <- S * (1 - 1 / 2 / n) + S ^ 2 / (2 * n * (p + 2))
    }
    P_value <- pchisq(S, p, lower.tail = FALSE)
    list(Norm_Mean = sqrt(Norm_Mean_2), n = n, p = p, Statistics = S,
         P_value = P_value)
}


# Clumsy colour palette
cols <- c("tomato", "royalblue4", "orange3", "gray40",
        "lightgreen", "magenta", "darkslategray1",
        "green3", "gold", "mediumorchid","cyan", "pink",
        "brown", "ivory2", "turquoise3", "purple",
        "maroon", "lightblue3", "yellow2", "deeppink")
colt <- apply(col2rgb(cols), 2, function(x)
                                rgb(x[1], x[2], x[3], alpha = 70, max = 255))
ltys <- c("solid", "44", "32", "2353", "1242", "4323")
colg <- cm.colors(200)

#### a2col #####
# Color picker based on angle
a2col <- function(x, s = 1, min.ang = 0, max.ang = pi - min.ang,
                  min.step = 1, max.step = 200) {
    ps <- function(x, p) sign(x) * abs(x) ^ p
    mid <- (min.ang + max.ang) / 2
    gr <- pi / (max.ang - min.ang)
    cd <- (ps((x - mid) * gr, s) + ps(mid, s)) / (2 * ps(mid, s))
    pmin(pmax(ceiling(cd * max.step), min.step), max.step)
}

#### sd.pooled #####
# Ad hoc function to calculate pooled SD from the data
sd.pooled <- function(x, ind) {
    if(length(x) != length(ind)) stop("lengths of x and ind should be the same")
    nl <- length(unique(ind))
    SS <- by(x, ind, function(y) sum((y - mean(y))^2))
    sqrt(sum(SS) / (length(x) - nl))
}

##### sample.str #####
# Resampling function for stratified bootstrap,
# in which within-locality sample size is kept constant
sample.str <- function(data, ind) {
    sss <- function(n) sample(n, n, TRUE)
    if(nrow(data) != length(ind))
        stop("lengths of data and ind should be the same")
    DF <- by(data, ind, function(X) X[sss(nrow(X)), ])
    ans <- DF[[1]]
    for(i in 2:length(DF)) {
        ans <- rbind(ans, DF[[i]])
    }
    ans
}

##### f.mean.diff #####
# Function to calculate phenotypic change vector
# Written after Stuart et al.'s (2017) codes
f.mean.diff <- function(data, id) {
    by.site <- substr(id, 1, 4)
    sites <- unique(by.site)
    i <- 1
    while(i < length(sites)) {
        if(substr(sites[i], 1, 3) == substr(sites[i + 1], 1, 3)) {
            i <- i + 2
        } else {
            sites <- sites[-i]
        }
    }
    by.site <- by.site[by.site %in% sites]
    trait.means <- apply(data, 2, function(x) tapply(x, by.site, mean))
    colnames(trait.means) <- colnames(data)
    by.watershed <- unique(substr(by.site, 1, 3))
    trait.mndiff <- matrix(NA, nrow = length(by.watershed), ncol = ncol(data))
    colnames(trait.mndiff) <- colnames(trait.means)
    indicator <- seq(1, length(unique(by.site)) - 1, 2)
    for (i in 1:length(indicator)) {
        trait.mndiff[i, ] <- trait.means[indicator[i], ] - trait.means[indicator[i] + 1, ]
    }
    rownames(trait.mndiff) <- by.watershed
    return(trait.mndiff)
}

##### kmeans.ex #####
# function to conduct exhaustive search of initial conditions
# with kmeans() function
kmeans.ex <- function(x, k, ...) {
    require(gtools)
    Inds <- combinations(nrow(x), k)
    Ans <- list()
    for(i in seq_len(nrow(Inds))) {
        Ans <- c(Ans, list(kmeans(x, x[Inds[i, ], ], ...)))
    }
    SS <- sapply(Ans, function(obj) obj$tot.withinss)
    ans <- Ans[!duplicated(SS)]
    SS <- SS[!duplicated(SS)]
    invisible(list(best = ans[[which.min(SS)]], SS = SS, all = ans))
}
