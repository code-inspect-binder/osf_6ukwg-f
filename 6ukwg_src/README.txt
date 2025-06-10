###############################
# Supplementary material for
#    Detecting (non)parallel evolution in multidimensional spaces:
#    angles, correlations, and eigenanalysis
#      by Junya Watanabe
# Version 2021.09.30
# This package involves the following files:
# - codes_reanalysis.R:      R codes to reproduce the example re-analysis.
#                            This is the main file; the others are to be loaded
#                            while executing this.
# - function_random_angle.R: R functions for the probability distribution of
#                            random angles (eq. 9).
# - stuart.R:                R codes for data pre-processing from
#                            Stuart et al. (2017).
# - utility_functions.R:     R functions used in the re-analysis.
# The re-analysis of Stuart et al. (2017; doi: 10.1038/s41559-017-0158) dataset
# requires files from that publication available from:
#   http://web.corral.tacc.utexas.edu/Stuart_2017_NatureEE_Data_Code/
#   (accessed 2021.08.17)
# In particular, the following files are to be placed in the same directory:
#   01.univMorph.csv, 05.geomo_nom.csv, 10.jaw4barCT.csv
# These codes were tested in R 3.5.3, with the following packages:
# - car     (version 3.0-2)
# - gtools  (version 3.8.1)
# - lme4    (version 1.1-21)
# - plotrix (version 3.7-4)
###############################
