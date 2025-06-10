###############################
# Supplementary material for
#    Detecting (non)parallel evolution in multidimensional spaces:
#    angles, correlations, and eigenanalysis
#      by Junya Watanabe
# Version 2021.09.30
# This file involves R codes to reproduce the re-analysis in the paper
# Depends on the following package:
# - car:     For plotting confidence ellipses
# - plotrix: For plotting colour scale rectangles
# - lme4:    For pre-processing (size correction); loaded within "stuart.R"
# - gtools:  Used in kmeans.ex() function; loaded within "utility_functions.R"
###############################

library(car)
library(plotrix)

# Pre-process data from Stuart et al. (2017)
# Ensure the following files from that publication are in the same directory:
#   01.univMorph.csv, 05.geomo_nom.csv, 10.jaw4barCT.csv
source("stuart.R")
# The warning is to be ignored, as per Stuart et al. (2017)

# Load utility functions
source("utility_functions.R")

# Load functions for the probability distribution of random angles
source("functions_random_angle.R")

##############
### Data processing
##############

## Omit incomplete observations
data.comp <- na.omit(morphometrics.nomarine)

## Omit sites without a counterpart
by.site <- substr(data.comp$fishID.univ, 1, 4)
sites <- unique(by.site)
i <- 1
while(i < length(sites)) {
    if(substr(sites[i], 1, 3) == substr(sites[i + 1], 1, 3)) {
        i <- i + 2
    } else {
        sites <- sites[-i]
    }
}
data.comp <- subset(data.comp, by.site %in% sites)
rm(by.site, sites, i)

## Scale data with pooled SD across all localities
## This differs from Stuart et al. (2017), who scaled data with pooled SD for
## each watershed (by t.test())
data.sc <- data.comp
data.sc[, 2:ncol(data.sc)] <- scale(data.sc[, 2:ncol(data.sc)],
                        scale = apply(data.sc[, 2:ncol(data.sc)], 2, sd.pooled,
                                      ind = substr(data.sc[, 1], 1, 4)))

fishID <- data.sc$fishID.univ

# Excluding meristic or density traits
select.col.nm <- select.col[!(select.col %in% c(118, 119, 120, 124))]
data.sc <- data.sc[, trait.select.col][, select.col.nm]

rm(trait.select.col, select.col, select.col.nm)

# Matrix of phenotypic change trajectories X and Z,
# inter-lineage correlation matrix C, and inter-trait cross-product matrix A,
# as well as the number of lineages n (= 13)
# and dimensionality of the vectors p (= 76)
X <- f.mean.diff(data.sc, fishID)
Z <- X / sqrt(diag(tcrossprod(X)))
C <- tcrossprod(Z)
A <- crossprod(Z)
n <- nrow(Z)
p <- ncol(Z) - 4 # -4 accounting for the loss of d.f. from Procrustes alignment


##############
### Analysis of pairwise angles
##############

# Vector of pairwise angles
a <- acos(C[lower.tri(C)])

# Objects for density fills outside 2.5 and 97.5 percentiles (see below)
xl <- seq(to = qang(0.025, p), by = 0.001, length = 200)
xlr <- xl[length(xl):1]
xu <- seq(from = qang(0.975, p), by = 0.001, length = 200)
xur <- xu[length(xu):1]

# Histogram of pairwise angles and scaled density of the null distribution
# in degree (Fig. 2A)
# Note that the density needs to be scaled to accommodate conversion to degree
# ("dang(...) / 180 * pi")
# The density is drawn to have the same area as the histogram,
# hence the additional scaling "* length(a) * 20"
hist(a / pi * 180, ylim = c(0, 93), main = "", xlab = "", ylab = "",
     yaxs = "i", col = "gray40", axes = FALSE)
curve(dang(x / 180 * pi, p) / 180 * pi * length(a) * 20, 20, 160,
      add = TRUE, col = "tomato", lwd = 1.5, lty = 1)
polygon(c(xl, xlr) / pi * 180, c(rep_len(0, length(xl)),
        dang(xlr, p) / 180 * pi * length(a) * 20),
        col = "tomato", border = "tomato")
polygon(c(xu, xur) / pi * 180, c(rep_len(0, length(xu)),
        dang(xur, p) / 180 * pi * length(a) * 20),
        col = "tomato", border = "tomato")
box(bty = "l")
axis(1, at = seq(20, 160, 20), label = FALSE)
axis(1, at = seq(20, 160, 40), label = as.expression(sapply(seq(20, 160, 40),
                                function(x) bquote(.(x)*degree))), tick = FALSE)
axis(2)
mtext("Angle", 1, line = 2.3)


# P-values against the null of no preferred directions,
# and how many of the pairwise angles are outside the 2.5 and 97.5 percentiles
pang(a, p)
sum(pang(a, p) < 0.025)
sum(pang(a, p) > 0.975)


# Monte Carlo null distribution of mean of 78 pairwise angles
# Note that this sampling is valid because the pairwise correlations/angles
# are uncorrelated (strictly) under the null hypothesis
Random_mean <- replicate(100000, mean(rang(length(a), p)))

# The P-value for the Monte Carlo test (P < 1 * 10^-5)
mean(mean(a) > Random_mean)

# Histogram for comparison (Fig. 2B)
hist(Random_mean / pi * 180, xlim = c(85, 94), main = "", xlab = "", ylab = "",
     yaxs = "i", col = "gray40", axes = FALSE)
box(bty = "l")
arrows(mean(a) / pi * 180, 15000, mean(a) / pi * 180, 1000, lwd = 2,
       col = "royalblue4")
text(mean(a) / pi * 180, 15000, "Observed", font = 2, adj = c(0.2, -0.5))
axis(1, at = seq(86, 94, 2), label = FALSE)
axis(1, at = seq(86, 94, 2), label = as.expression(sapply(seq(86, 94, 2),
                                function(x) bquote(.(x)*degree))), tick = FALSE)
axis(2)
mtext("Mean of 78 angles", 1, line = 2.3)



# Eigenanalysis of C and A
eC <- eigen(C)
eA <- eigen(A)
# Positive eigenvalues of A
eAp <- eA$values[eA$values > .Machine$double.eps * 100]
# Trivial sanity check for the equality of non-zero eigenvalues
all.equal(eC$values, eAp)

# Schott's (2005) for independent directions
# The test statistic T, null mean and variance (eq. 16)
(T <- sum(cos(a) ^ 2))
(nullmean <- n * (n - 1) / (2 * p))
(nullvar <- n * (n - 1) * (p - 1) / (p ^ 2 * (p + 2)))
# Z value (effect size against the standard normal distribution) and P-value
(Z_stat <- (T - nullmean) / sqrt(nullvar))
pnorm(Z_stat, lower.tail = FALSE)

# Checking equality of eigenvalue dispersions and sum of squared correlations
# (eq. 15)
all.equal(2 * T, with(eC, sum((values - mean(values))^2)))
all.equal(2 * T, with(eA, sum((values - mean(values))^2)) + n ^ 2 / (p + 4) - n)


# Rayleigh test
# Note that the dimensionality (p) and P-value returned by the function is
# inaccurate, as it is not aware of the true dimensionality
(Rayleigh_result <- rayleigh_test(Z))
# The correct P-value with appropriate dimensionality:
pchisq(Rayleigh_result$Statistics, p, lower.tail = FALSE)



##############
### Ordination
##############
# PCA of the phenotypic change vectors Z, rotation matrix R and mean vector M
# As you would know, signs of eigenvectors are arbitrary and hence may swap,
# depending on specific environments (especially LAPACK libraries used).
# This is non-essential.
PCA.Z <- prcomp(Z, center = FALSE, scale = FALSE)
R <- PCA.Z$rotation
M <- PCA.Z$center

# Labels for graphs
labPC1 <- paste0("PC1 (", sprintf("%.1f",
                                  summary(PCA.Z)$importance[2, 1] * 100), "%)")
labPC2 <- paste0("PC2 (", sprintf("%.1f",
                                  summary(PCA.Z)$importance[2, 2] * 100), "%)")


# Preparing labels for the following plot
# Retaining labels for most highly-loaded traits in the PC1-2
R2 <- R
R2[grepl("gpa", rownames(R2)), ] <- 0
ylabs <- rownames(R2)
ylabs <- gsub("log10.", "", ylabs)
ylabs <- gsub("\\.mm.*", "", ylabs)
ylabs[rank(diag(tcrossprod(R2[, 1:2]))) < 75] <- ""

# PCA biplot (Fig. 2C)
biplot2(PCA.Z, xlab = "", ylab = "", ylabs = ylabs, expand = 0.88,
        cex.xpt = 1.5, cex.xtx = 0.9, cex.ytx = 0.9, arrow.len = 0.05,
        lwd.yar = 0.8,
        col.yar = ifelse(grepl("gpa", rownames(R2)), cols[2], cols[1]),
        col = c("black", cols[1]))
mtext(labPC1, 1, line = 2.3)
mtext(labPC2, 2, line = 2.5)


# PC plot with pairwise angles in colour scale (Fig. 2D)
plot(PCA.Z$x[, 1:2], asp = 1, type = "n", main = "", xlab = "", ylab = "",
     axes = FALSE)
for(i in 1:(n - 1)) {
    for(j in (i + 1):n) {
        lines(PCA.Z$x[c(i, j), 1:2],
              col = colg[a2col(acos(C[i, j]), min.ang = pi / 6)])
    }
}
with(PCA.Z, {
    points(x[, 1:2], pch = 20, cex = 1.5, col = "black")
    text(x[, 1:2], rownames(x), cex = 0.9, pos = 3)
})
box(bty = "o")
axis(1)
axis(2)
mtext(labPC1, 1, line = 2.3)
mtext(labPC2, 2, line = 2.5)
gradient.rect(-0.7, -0.35, -0.0, -0.3, col = colg)
text(-0.7, -0.35, expression(paste("<", 30*degree)), pos = 1)
text(-0.0, -0.35, expression(paste(150*degree, "<")), pos = 1)



##############
### Bootstrap
##############
sites <- substr(fishID, 1, 4)

# Number of individuals involved per site
table(sites)

# Obtaining phenotypic change vectors from bootstrapped observations
# This typically takes several minutes
b <- 5000
set.seed(62208)
Xb <- array(dim = c(dim(X), b))
for(i in seq_len(b))  {
    data.i <- sample.str(data.sc, sites)
    Xb[, , i] <- f.mean.diff(data.i, fishID)
}
Zb <- array(apply(Xb, 3, function(X) X / sqrt(diag(tcrossprod(X)))), dim = c(dim(X), b))


## Ordination + bootstrap clouds (Fig. 2E)
plot(PCA.Z$x[, 1:2], asp = 1, type = "n", main = "", xlab = "", ylab = "",
     axes = FALSE)
# Plotting bootstrap replicates (just 1000 out of 5000 for visual clarity)
for(i in seq_len(nrow(Zb))) {
    S1 <- pcc(t(Zb[i, , 1:1000]), R, M)
    dataEllipse(S1[, 1], S1[, 2], add = TRUE, levels = NA, center.pch = NULL,
          plot.points = TRUE, pch = ".", cex = 0.5, fill = FALSE, col = colt[i])
}
# Drawing confidence ellipses (from the full bootstrap replicates)
for(i in seq_len(nrow(Zb))) {
    S1 <- pcc(t(Zb[i, , ]), R, M)
    dataEllipse(S1[, 1], S1[, 2], add = TRUE, levels = 0.95, center.pch = NULL,
          plot.points = FALSE, fill = FALSE, col = cols[i])
}
with(PCA.Z, {
    points(x[, 1:2], pch = 20, cex = 1.5, col = "black")
    text(x[, 1:2], rownames(x), cex = 0.9, pos = 3)
})
box(bty = "o")
axis(1)
axis(2)
mtext(labPC1, 1, line = 2.3)
mtext(labPC2, 2, line = 2.3)


# Histograms comparing dissimilarity between the original and
# bootstrap replicates between "Moo" and "Pac" (largest and smallest dispersion)
# This shows that the heteroscedasticity is real (not due to visual distortion)
hist(acos(Z[10, ] %*% Zb[10, , ]), breaks = seq(0, pi / 2, 0.1),
     col = "#FF000080", main = "Pac (red) versus Moo (blue)",
     xlab = "Angle between original and bootstrap vectors (radian)")
hist(acos(Z[7, ] %*% Zb[7, , ]), add = TRUE,
     breaks = seq(0, pi / 2, 0.1), col = "#0000FF80")



##############
### k-means clustering
##############

# Repeating exhaustive search for the k-means clustering,
# with a plot of the best clustering for each k
KM2 <- kmeans.ex(PCA.Z$x, 2)
KM2b <- KM2$best
with(PCA.Z, {
    plot(x[, 1:2], pch = 20, asp = 1, col = cols[KM2b$cluster])
    text(x[, 1:2], rownames(x), pos = 3)
    addhull(x[, 1], x[, 2], factor(KM2b$cluster), col.h = cols[1:2])
})
points(KM2b$centers[, 1], KM2b$centers[, 2], pch = 15, col = cols[1:2])

KM3 <- kmeans.ex(PCA.Z$x, 3)
KM3b <- KM3$best
with(PCA.Z, {
    plot(x[, 1:2], pch = 20, asp = 1, col = cols[KM3b$cluster])
    text(x[, 1:2], rownames(x), pos = 3)
    addhull(x[, 1], x[, 2], factor(KM3b$cluster), col.h = cols[1:3])
})
points(KM3b$centers[, 1], KM3b$centers[, 2], pch = 15, col = cols[1:3])

KM4 <- kmeans.ex(PCA.Z$x, 4)
KM4b <- KM4$best
with(PCA.Z, {
    plot(x[, 1:2], pch = 20, asp = 1, col = cols[KM4b$cluster])
    text(x[, 1:2], rownames(x), pos = 3)
    addhull(x[, 1], x[, 2], factor(KM4b$cluster), col.h = cols[1:4])
})
points(KM4b$centers[, 1], KM4b$centers[, 2], pch = 15, col = cols[1:4])

# k = 5 (Fig. 2F)
KM5 <- kmeans.ex(PCA.Z$x, 5)
KM5b <- KM5$best
with(PCA.Z, {
    plot(x[, 1:2], pch = 20, asp = 1, col = cols[KM5b$cluster])
    text(x[, 1:2], rownames(x), pos = 3)
    addhull(x[, 1], x[, 2], factor(KM5b$cluster), col.h = cols[1:5])
})
points(KM5b$centers[, 1], KM5b$centers[, 2], pch = 15, col = cols[1:5])

KM6 <- kmeans.ex(PCA.Z$x, 6)
KM6b <- KM6$best
with(PCA.Z, {
    plot(x[, 1:2], pch = 20, asp = 1, col = cols[KM6b$cluster])
    text(x[, 1:2], rownames(x), pos = 3)
    addhull(x[, 1], x[, 2], factor(KM6b$cluster), col.h = cols[1:6])
})
points(KM6b$centers[, 1], KM6b$centers[, 2], pch = 15, col = cols[1:6])

KM7 <- kmeans.ex(PCA.Z$x, 7)
KM8 <- kmeans.ex(PCA.Z$x, 8)
KM9 <- kmeans.ex(PCA.Z$x, 9)
KM10 <- kmeans.ex(PCA.Z$x, 10)
KM11 <- kmeans.ex(PCA.Z$x, 11)
KM12 <- kmeans.ex(PCA.Z$x, 12)
KMs <- list(NULL, KM2, KM3, KM4, KM5, KM6, KM7, KM8, KM9, KM10, KM11, KM12)

# Plot of the proportion of within-grou sum of squares
plot(cbind(2:12, 1), type = "n", ylim = c(0, 0.45), xlab = "", ylab = "",
     yaxs = "i", axes = FALSE)
par(las = 1)
invisible(sapply(2:12, function(x) points(cbind(x, KMs[[x]]$SS / KM2b$totss),
                                          pch = 20)))
box(bty = "l")
axis(1, at = 2:12, label = FALSE)
axis(1, at = 1:6 * 2, tick = FALSE)
axis(2)
mtext(expression(italic(k)), 1, line = 2)
par(las = 0)
mtext("Unexplained SS / Total SS", 2, line = 2.7)
