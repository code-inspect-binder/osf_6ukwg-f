###############################
# Supplementary material for
#    Detecting (non)parallel evolution in multidimensional spaces:
#    angles, correlations, and eigenanalysis
#      by Junya Watanabe
# Version 2021.09.30
# This file involves R codes adopted (with modifications) from
# Stuart et al. (2017; doi: 10.1038/s41559-017-0158) to pre-process the dataset
# for re-analysis.  Originally retrieved on 2021.08.17 from:
#   http://web.corral.tacc.utexas.edu/Stuart_2017_NatureEE_Data_Code/
# In this file, double hash marks (##) denote Stuart et al's original comments
# Ensure the following files from that publication are in the same directory:
#   01.univMorph.csv, 05.geomo_nom.csv, 10.jaw4barCT.csv
# Partly depends on the following package:
# - lme4: in f.sizecorrect() for size correction
###############################

## === f.size.correct ===
## The output will have all the original columns from dataframe.to.use with the size-corrected columns added on at the end.
f.sizecorrect <- function(standard.length, vector.of.columns, dataframe.to.use, watershed) {
    require(lme4)
  ## Calculate overall mean standard length (Ls)
    Ls <- mean(standard.length, na.rm = TRUE)
  ## Call individual standard length
    L0 <- standard.length
  ## Calculate common, within-group slope, b, for each trait.
  ## Here, I am treating each watershed as a group for random factor, rather than each lake or stream population, because we are comparing parallelism vectors across watersheds.
    b.vector.lmm <- vector()
    for (i in vector.of.columns) {
        abcd <- (dataframe.to.use[i])
        b.model <- lmer((abcd[, ]) ~ (standard.length) + (1 | watershed))
        b <- coef(summary(b.model))[2, 1]
        b.vector.lmm <- c(b.vector.lmm, b)
    }
  ## size correct
    xx <- dataframe.to.use
    columnnames <- colnames(xx)
    j <- 1
    for (i in vector.of.columns) {
        M0 <- xx[, i] ## grab the appropriate column of data
        Ms <- M0 * ((Ls / L0) ^ b.vector.lmm[j]) ## size correction formula
        j <- j + 1
        columnnames <- c(columnnames, paste(colnames(xx[i]), "sc", sep = "."))
        xx <- cbind(xx, Ms)
    }
    colnames(xx) <- columnnames ## Rename the columns in the temporary dataframe xx
    return(xx) ## Output a new dataframe with the name provided in "outputfilename"
}

## 4) Calculating theta, deltaL, and L
## =============================================================================
## 4a) Calculating theta, dL, and centroid distance in MORPHO space for ALL TRAITS, using "01.univMorph.csv" and "05.geomo_nom.csv"
## === 4ai) log 10, size correct univariate traits, then merge with geometric morphometric data. ====
univmorph.t <- read.csv("01.univMorph.csv", header = TRUE)
## Pelvic Spine lengths cause a problem in downstream PCA because of missing data (spines were clipped).
## Create a new column that is the average of left and right pelvic spine lengths for each fish.
meanspine.t <- vector()
for (i in 1:length(univmorph.t$Right.Side.Pelvic.Spine.Length.mm)) {
    temp <- mean(c(univmorph.t$Right.Side.Pelvic.Spine.Length.mm[i], univmorph.t$Left.Side.Pelvic.Spine.Length.mm[i]), na.rm = TRUE)
    meanspine.t <- c(meanspine.t, temp)
}
univmorph.t$mean.pelvic.spine.length <- meanspine.t
## load jaw data
jawctt <- read.csv("10.jaw4barCT.csv", header = TRUE, stringsAsFactors = FALSE)
colnames(jawctt)[6] <- "standard.length.from.calipers.jawct"
jawcttt <- jawctt[, c(1, 6, 7:9, 11, 13:14, 15:19)] ## leaving out sex, gape width and buccal cavity length because they already exist in "morph"
univmorph.t <- merge(univmorph.t, jawcttt, by.x = "fishID.univ", by.y = "fishID.univ", all = TRUE)
univmorph.t <- univmorph.t[, -56] ## remove duplicate standard.length
## log transform and size correct
## log transform
cols.to.log <- c(7:21, 27:36, 41, 42, 43:46, 48, 51:55, 56:66)
columnnames <- colnames(univmorph.t)
for (i in cols.to.log) {
    log.i <- log10(univmorph.t[, i])
    univmorph.t <- cbind(univmorph.t, log.i)
    columnnames <- c(columnnames, paste("log10", colnames(univmorph.t[i]), sep = "."))
}
colnames(univmorph.t) <- columnnames
## size correct
standard.length <- univmorph.t$log10.Standard.Length.mm
vector.of.columns <- c(67:91, 94:114) ## A warning arises for colmn 85: "log10.body.width.midbody.mm". Some sort of failure in convergence for lmer, but still returns a size corrected trait.
dataframe.to.use <- univmorph.t
watershed <- univmorph.t$watershed
univmorph.tsc <- f.sizecorrect(standard.length, vector.of.columns, dataframe.to.use, watershed)

## Merge and define dataframes for downstream use
univmorph <- univmorph.tsc ## this contains columns of untransformed data, log10transformed data, and size-corrected log10data
geomorpho.nomarine <- read.csv("05.geomo_nom.csv", header = TRUE) ## no marine fish
morphometrics.tnm <- merge(univmorph, geomorpho.nomarine, by.x = "fishID.univ", by.y = "fishID.univ", all = FALSE)
morphometrics.nomarine <- morphometrics.tnm[, c(1, 7:21, 27:36, 41, 42, 43:46, 48, 51:55, 56:58, 59:61, 62:66, 67:114, 115:160, 162, 163:200, 201:215)]
## =====

trait.select.col <- c(2:28, 30, 32:49, 50:76, 78, 80:97, 98:122, 124, 126:143, 144:197)
select.col <- c(93:136, 73, 137, 138:175)

# Codes from Stuart et al. (2017) end here
###########################################################

# Just to clean up temporary objects
rm(cols.to.log, columnnames, dataframe.to.use, f.sizecorrect,
   geomorpho.nomarine, i, jawctt, jawcttt, log.i, meanspine.t,
   morphometrics.tnm, standard.length, temp, univmorph, univmorph.t,
   univmorph.tsc, vector.of.columns, watershed)
