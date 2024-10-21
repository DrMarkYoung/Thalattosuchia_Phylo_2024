#!/usr/bin/env Rscript

# ======================================================================================= #
#                                                                                         #
#  Master R script - version 2.0 (October, 2024).                                         #
#                                                                                         #
#                                                                                         #
#  Purpose:                                                                               #
#  This is a 'master' script, designed to run all analysis scripts.                       #
#                                                                                         #
#                                                                                         #
#  Version:                                                                               #
#  2.0 supersedes all previous versions. (v1.5 -> v2.0)                                   #
#                                                                                         #
#  It was altered to be GitHub compatible, and the GitHub template created by Daniel      #
#  Morillo. Daniel Morillo's template is licensed under Creative Commons Attribution 4.0  #
#  International license. See: https://github.com/DaniMori/rproj-template  and:           #
#  https://creativecommons.org/licenses/by/4.0/                                           #
#                                                                                         #
#                                                                                         #
#  Related scripts:                                                                       #
#  This script is designed to work in conjunction with the TNT scripts created for        #
#  Young et al. (2024: ZJLS 200:547-617). See: https://doi.org/10.1093/zoolinnean/zlad165 #
#                                                                                         #
#                                                                                         #
#  Origin:                                                                                #
#  Version 1.0 of this script was based on Moon & Stubbs (2020: Communications Biology    #
#  3: 68). See: https://www.nature.com/articles/s42003-020-0779-6                         #
#  It also incorporated some elements from Wright & Lloyd (2020 Palaeontology 63: 997-    #
#  1006). See: https://onlinelibrary.wiley.com/doi/abs/10.1111/pala.12500                 #
#                                                                                         #
#                                                                                         #
#  License and Rights Statement:                                                          #
#  This work is licensed under Creative Commons Attribution 4.0 International license     #
#  (CC BY 4.0). https://creativecommons.org/licenses/by/4.0/                              #
#  Users worldwide are free to copy, redistribute, remix, transform, and build upon the   #
#  material in any medium or format, even commercially. The licensor cannot revoke these  #
#  freedoms as long as you follow the license terms. You must give appropriate credit,    #
#  provide a link to the license, and indicate if changes were made. You may do so in any #
#  reasonable manner, but not in any way that suggests the licensor endorses you or your  #
#  use.                                                                                   #
#                                                                                         #
#                                                                                         #
#  Mark T. Young                                                                          #
#  marktyoung1984@gmail.com                                                               #
# ======================================================================================= #


# ================================================================ #
# 0. Installation of required packages and creation of directories #
# ================================================================ #


# -------------------------------------- #
# 0.1 Required packages and installation #
# -------------------------------------- #

 library(beepr)

### Packages on CRAN:
 list_of_packages <- c("ape", "beepr", "Claddis", "data.table", "devtools", "doParallel", "geoscale", 
                       "geiger", "ggplot2", "graphic", "grid", "foreach", "magrittr", "mailR",
                       "maps", "Morpho", "nlme", "paleotree", "parallel", "phytools", "plot3D", "pscl",
                       "renv","rr2", "strap", "tidyr", "TeachingDemos", "TreeDist", "TreeTools", "vegan")

 new_packages <- list_of_packages[!(list_of_packages %in%
                                   installed.packages()[, "Package"])]


### Change CRAN repository as desired:
 if (length(new_packages)) {
   install.packages(new_packages, repos = "https://www.stats.bris.ac.uk/R/", dependencies = TRUE)
 }

 beep(1)

### IF you have never installed any of these packages, or after an R update, using "dependencies = TRUE" can help if something *odd* occurs ###
### i.e.: install.packages("ape", dependencies = TRUE) ###


 library(devtools)
 devtools::install_github("YuLab-SMU/ggtree")


# ---------------------------------------------- #
# 0.2 Repository file structure and installation #
# ---------------------------------------------- #


### The repository file structure follows the GitHub template created by Daniel Morillo (see above for link to license and template)

### |
### |--- apps         (To store apps, e.g. in Shiny)
### |
### |--- dat          (To store input datasets)
### |
### |--- doc          (To store important documentation of the project)
### |    |
### |    |--- minutes (To store meeting minutes)
### |
### |--- notebooks    (Notebooks to explore data and test processes live here)
### |
### |--- output       (Processing outputs)
### |
### |--- R            (R functions created for this project)
### |
### |--- renv         (System library necessary for `renv` to work)
### |
### |--- src          (Source scripts that implement the main processes)
### |
### |--- www          (Project assets, e.g., images, bibliography files, etc.)


### Check for directories:
 directories <- c("./apps",
                  "./dat", "./dat/Ages", "./dat/Trees", "./dat/Matrices",
                  "./doc", "./doc/minutes",
                  "./notebooks",
                  "./output",
                  "./output/PhyloFigs","./output/PhyloFigs_Geoscaled","./output/PhyloFigs_Geoscaled/Pruned",
                  "./output/StratCongruenceResults_EqW", "./output/StratCongruenceResults_EqW/Boxplots", 
                  "./output/StratCongruenceResults_EIWk15", "./output/StratCongruenceResults_EIWk15/Boxplots", 
                  "./output/TreeDistance",
                  "./R",
                  "./renv",
                  "./src",
                  "./www")


 missing_dir <- directories[!(directories %in% list.dirs())]

### Create new directories:
 for (dir in missing_dir) dir.create(dir)

 beep(1)

# -------------------- #
# 0.3 Function scripts #
# -------------------- #

### loads script:
 source("R/sendEmailNotification.R")



# ===================== #
# 1. Phylogeny figures  #
# ===================== #


### Creates various basic types of phylogeny figures for all phylogenetic analyses. ###

 source("src/1a-Basic_Phylogeny_Images.R")


### Creates various time-scaled strict consensus topologies figures for all phylogenetic analyses. ###

 source("src/1b-Time_Scaled_Phylogeny_Images.R")



# ======================== #
# 2. Tree topology metrics #
# ======================== #

# Runs stratigraphic congruence analyses from the equal weights and the extended           #
# implied weighting (k=15) phylogenetic analyses. Creates boxplots of the various metrics. #

 source("src/2a-StratPhyloCongruence.R")



# Runs various different tree distance metrics for the equal weights and the extended      #
# implied weighting (k=15) phylogenetic analyses. Creates boxplots of the various metrics. #

 source("src/2b-TreeDistance.R")


# ========= #
# 0.3 Done! #
# ========= #

### Noise when script has finished (useful if you are running this in the background):
 library(beepr)

 beep(8)