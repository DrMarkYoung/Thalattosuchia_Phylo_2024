#!/usr/bin/env Rscript

# ======================================================================================= #
#                                                                                         #
#  Tree distance script - version 2.0 (October, 2024).                                    #
#                                                                                         #
#                                                                                         #
#  Purpose:                                                                               #
#  Calculates various tree distance metrics.                                              #
#                                                                                         #
#                                                                                         #
#  Version:                                                                               #
#  2.0 supersedes all previous versions. (v1.5 -> v2.0)                                   #
#                                                                                         #
#                                                                                         #
#  Related scripts:                                                                       #
#  This script is designed to work in conjunction with the TNT scripts created for        #
#  Young et al. (2024: ZJLS 200:547-617). See: https://doi.org/10.1093/zoolinnean/zlad165 #
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



# ========================== #
# 0.1 Parallelisation set-up #
# ========================== #

### Required libraries:
 library(parallel)
 library(foreach)
 library(doParallel)

### Set-up (note, fork is used as I use Mac OS.)
 core_number <- detectCores() - 1
 cluster <- makeForkCluster(core_number)
 registerDoParallel(cluster)


# ===================== #
# 1. Format input files #
# ===================== #

# ------------------------------------------------------------------------------- #
#   This script uses the same input files as 2a-StratPhyloCongruence.R            #
#   If you want to run only the tree distance script (and not the stratigraphic   #
#   congruence script), please copy and paste section 1 from the 2a script here   #
# ------------------------------------------------------------------------------- #


# ============================================================== #
# 2. Tree distance/similarity metric results for 100 random MPCs #
# ============================================================== #

# --------------------------------------------------------------- #
# 2.1 The distance/similarity metric analyses for 100 random MPCs #
# --------------------------------------------------------------- #

### Required libraries:
 library(beepr)
 library(TreeDist)

clusterCall(cluster, function () {
library(TreeDist)
})

### Smith (2020) Shared phylogenetic information tree distance metric:

 SPI_EqW_C0vC1 <- SharedPhylogeneticInfo (EqW_C0_Trees_subsample, EqW_C1_Trees_subsample, normalize = TRUE)
 SPI_EqW_C0vC2 <- SharedPhylogeneticInfo (EqW_C0_Trees_subsample, EqW_C2_Trees_subsample, normalize = TRUE)
 SPI_EqW_C0vC3 <- SharedPhylogeneticInfo (EqW_C0_Trees_subsample, EqW_C3_Trees_subsample, normalize = TRUE)
 SPI_EqW_C0vC4 <- SharedPhylogeneticInfo (EqW_C0_Trees_subsample, EqW_C4_Trees_subsample, normalize = TRUE)
 SPI_EIWk15_C0vC1 <- SharedPhylogeneticInfo (EIWk15_C0_Trees_subsample, EIWk15_C1_Trees_subsample, normalize = TRUE)
 SPI_EIWk15_C0vC2 <- SharedPhylogeneticInfo (EIWk15_C0_Trees_subsample, EIWk15_C2_Trees_subsample, normalize = TRUE)
 SPI_EIWk15_C0vC3 <- SharedPhylogeneticInfo (EIWk15_C0_Trees_subsample, EIWk15_C3_Trees_subsample, normalize = TRUE)
 SPI_EIWk15_C0vC4 <- SharedPhylogeneticInfo (EIWk15_C0_Trees_subsample, EIWk15_C4_Trees_subsample, normalize = TRUE)


### Smith (2020) Mutual clustering information tree distance metric:

 MCI_EqW_C0vC1 <- MutualClusteringInfo (EqW_C0_Trees_subsample, EqW_C1_Trees_subsample, normalize = TRUE)
 MCI_EqW_C0vC2 <- MutualClusteringInfo (EqW_C0_Trees_subsample, EqW_C2_Trees_subsample, normalize = TRUE)
 MCI_EqW_C0vC3 <- MutualClusteringInfo (EqW_C0_Trees_subsample, EqW_C3_Trees_subsample, normalize = TRUE)
 MCI_EqW_C0vC4 <- MutualClusteringInfo (EqW_C0_Trees_subsample, EqW_C4_Trees_subsample, normalize = TRUE)
 MCI_EIWk15_C0vC1 <- MutualClusteringInfo (EIWk15_C0_Trees_subsample, EIWk15_C1_Trees_subsample, normalize = TRUE)
 MCI_EIWk15_C0vC2 <- MutualClusteringInfo (EIWk15_C0_Trees_subsample, EIWk15_C2_Trees_subsample, normalize = TRUE)
 MCI_EIWk15_C0vC3 <- MutualClusteringInfo (EIWk15_C0_Trees_subsample, EIWk15_C3_Trees_subsample, normalize = TRUE)
 MCI_EIWk15_C0vC4 <- MutualClusteringInfo (EIWk15_C0_Trees_subsample, EIWk15_C4_Trees_subsample, normalize = TRUE)


### Nye et al. (2006) similarity metric, normalised against the average number of splits in the two topologies (default):

 Nye1_EqW_C0vC1 <- NyeSimilarity (EqW_C0_Trees_subsample, EqW_C1_Trees_subsample, normalize = TRUE)
 Nye1_EqW_C0vC2 <- NyeSimilarity (EqW_C0_Trees_subsample, EqW_C2_Trees_subsample, normalize = TRUE)
 Nye1_EqW_C0vC3 <- NyeSimilarity (EqW_C0_Trees_subsample, EqW_C3_Trees_subsample, normalize = TRUE)
 Nye1_EqW_C0vC4 <- NyeSimilarity (EqW_C0_Trees_subsample, EqW_C4_Trees_subsample, normalize = TRUE)
 Nye1_EIWk15_C0vC1 <- NyeSimilarity (EIWk15_C0_Trees_subsample, EIWk15_C1_Trees_subsample, normalize = TRUE)
 Nye1_EIWk15_C0vC2 <- NyeSimilarity (EIWk15_C0_Trees_subsample, EIWk15_C2_Trees_subsample, normalize = TRUE)
 Nye1_EIWk15_C0vC3 <- NyeSimilarity (EIWk15_C0_Trees_subsample, EIWk15_C3_Trees_subsample, normalize = TRUE)
 Nye1_EIWk15_C0vC4 <- NyeSimilarity (EIWk15_C0_Trees_subsample, EIWk15_C4_Trees_subsample, normalize = TRUE)


### Nye et al. (2006) similarity metric, normalised against the maximum similarity possible for the unconstrained topology, using NSplits:
### (i.e. the best possible match to the unconstrained topology) ###

 Nye2_EqW_C0vC1 <- NyeSimilarity (EqW_C0_Trees_subsample, EqW_C1_Trees_subsample, normalize = TreeTools::NSplits(EqW_C0_Trees_subsample))
 Nye2_EqW_C0vC2 <- NyeSimilarity (EqW_C0_Trees_subsample, EqW_C2_Trees_subsample, normalize = TreeTools::NSplits(EqW_C0_Trees_subsample))
 Nye2_EqW_C0vC3 <- NyeSimilarity (EqW_C0_Trees_subsample, EqW_C3_Trees_subsample, normalize = TreeTools::NSplits(EqW_C0_Trees_subsample))
 Nye2_EqW_C0vC4 <- NyeSimilarity (EqW_C0_Trees_subsample, EqW_C4_Trees_subsample, normalize = TreeTools::NSplits(EqW_C0_Trees_subsample))
 Nye2_EIWk15_C0vC1 <- NyeSimilarity (EIWk15_C0_Trees_subsample, EIWk15_C1_Trees_subsample, normalize = TreeTools::NSplits(EIWk15_C0_Trees_subsample))
 Nye2_EIWk15_C0vC2 <- NyeSimilarity (EIWk15_C0_Trees_subsample, EIWk15_C2_Trees_subsample, normalize = TreeTools::NSplits(EIWk15_C0_Trees_subsample))
 Nye2_EIWk15_C0vC3 <- NyeSimilarity (EIWk15_C0_Trees_subsample, EIWk15_C3_Trees_subsample, normalize = TreeTools::NSplits(EIWk15_C0_Trees_subsample))
 Nye2_EIWk15_C0vC4 <- NyeSimilarity (EIWk15_C0_Trees_subsample, EIWk15_C4_Trees_subsample, normalize = TreeTools::NSplits(EIWk15_C0_Trees_subsample))

 beep(1)


# -------------------------------------------------------------------- #
# 2.2 Unpack the data from the tree distance/similarity metric results #
# -------------------------------------------------------------------- #

### This unpacks the data - depreciated but this works at the moment ###
### Required libraries:
 library(tidyr)

 SPI_EqW_C0vC1_unpacked <- extract_numeric(SPI_EqW_C0vC1)
 SPI_EqW_C0vC2_unpacked <- extract_numeric(SPI_EqW_C0vC2)
 SPI_EqW_C0vC3_unpacked <- extract_numeric(SPI_EqW_C0vC3)
 SPI_EqW_C0vC4_unpacked <- extract_numeric(SPI_EqW_C0vC4)
 SPI_EIWk15_C0vC1_unpacked <- extract_numeric(SPI_EIWk15_C0vC1)
 SPI_EIWk15_C0vC2_unpacked <- extract_numeric(SPI_EIWk15_C0vC2)
 SPI_EIWk15_C0vC3_unpacked <- extract_numeric(SPI_EIWk15_C0vC3)
 SPI_EIWk15_C0vC4_unpacked <- extract_numeric(SPI_EIWk15_C0vC4)

 MCI_EqW_C0vC1_unpacked <- extract_numeric(MCI_EqW_C0vC1)
 MCI_EqW_C0vC2_unpacked <- extract_numeric(MCI_EqW_C0vC2)
 MCI_EqW_C0vC3_unpacked <- extract_numeric(MCI_EqW_C0vC3)
 MCI_EqW_C0vC4_unpacked <- extract_numeric(MCI_EqW_C0vC4)
 MCI_EIWk15_C0vC1_unpacked <- extract_numeric(MCI_EIWk15_C0vC1)
 MCI_EIWk15_C0vC2_unpacked <- extract_numeric(MCI_EIWk15_C0vC2)
 MCI_EIWk15_C0vC3_unpacked <- extract_numeric(MCI_EIWk15_C0vC3)
 MCI_EIWk15_C0vC4_unpacked <- extract_numeric(MCI_EIWk15_C0vC4)

 Nye1_EqW_C0vC1_unpacked <- extract_numeric(Nye1_EqW_C0vC1)
 Nye1_EqW_C0vC2_unpacked <- extract_numeric(Nye1_EqW_C0vC2)
 Nye1_EqW_C0vC3_unpacked <- extract_numeric(Nye1_EqW_C0vC3)
 Nye1_EqW_C0vC4_unpacked <- extract_numeric(Nye1_EqW_C0vC4)
 Nye1_EIWk15_C0vC1_unpacked <- extract_numeric(Nye1_EIWk15_C0vC1)
 Nye1_EIWk15_C0vC2_unpacked <- extract_numeric(Nye1_EIWk15_C0vC2)
 Nye1_EIWk15_C0vC3_unpacked <- extract_numeric(Nye1_EIWk15_C0vC3)
 Nye1_EIWk15_C0vC4_unpacked <- extract_numeric(Nye1_EIWk15_C0vC4)

 Nye2_EqW_C0vC1_unpacked <- extract_numeric(Nye2_EqW_C0vC1)
 Nye2_EqW_C0vC2_unpacked <- extract_numeric(Nye2_EqW_C0vC2)
 Nye2_EqW_C0vC3_unpacked <- extract_numeric(Nye2_EqW_C0vC3)
 Nye2_EqW_C0vC4_unpacked <- extract_numeric(Nye2_EqW_C0vC4)
 Nye2_EIWk15_C0vC1_unpacked <- extract_numeric(Nye2_EIWk15_C0vC1)
 Nye2_EIWk15_C0vC2_unpacked <- extract_numeric(Nye2_EIWk15_C0vC2)
 Nye2_EIWk15_C0vC3_unpacked <- extract_numeric(Nye2_EIWk15_C0vC3)
 Nye2_EIWk15_C0vC4_unpacked <- extract_numeric(Nye2_EIWk15_C0vC4)

 beep(1)


# ------------------------------------------------------------------------------- #
# 2.3 PDFs of the equal weights analyses tree distance/similarity metric boxplots #
# ------------------------------------------------------------------------------- #

### Creates boxplots comparing the tree distance/similarity metrics of the different constraint analyses against the unconstrained analysis, creates PDFs (colours chosen so that people with colourblindness can differentiate the boxes easily):

## Shared phylogenetic Information tree distance metric boxplots:
 boxplot(SPI_EqW_C0vC1_unpacked, SPI_EqW_C0vC2_unpacked, SPI_EqW_C0vC3_unpacked, SPI_EqW_C0vC4_unpacked, xlab= "Phylogenetic hypotheses", ylab = "Shared phylogenetic Information  tree distance metric", col=c("#E69F00", "#56B4E9", "#009E73", "#F0E442"), names= c("C0vsC1", "C0vC2", "C0vC3", "C0vC4"))

pdf("output/TreeDistance/Shared_Phylogenetic_Information_Metric_EqW_Results_C0vC1vC3vC4.pdf", width=17, height=9) 
 boxplot(SPI_EqW_C0vC1_unpacked, SPI_EqW_C0vC2_unpacked, SPI_EqW_C0vC3_unpacked, SPI_EqW_C0vC4_unpacked, xlab= "Phylogenetic hypotheses", ylab = "Shared phylogenetic Information tree distance metric", col=c("#E69F00", "#56B4E9", "#009E73", "#F0E442"), names= c("C0vsC1", "C0vC2", "C0vC3", "C0vC4"))

 dev.off()


## Mutual clustering information tree distance metric boxplots:
 boxplot(MCI_EqW_C0vC1_unpacked, MCI_EqW_C0vC2_unpacked, MCI_EqW_C0vC3_unpacked, MCI_EqW_C0vC4_unpacked, xlab= "Phylogenetic hypotheses", ylab = "Mutual clustering information tree distance metric", col=c("#E69F00", "#56B4E9", "#009E73", "#F0E442"), names= c("C0vsC1", "C0vC2", "C0vC3", "C0vC4"))

 pdf("output/TreeDistance/Mutal_Clustering_Information_Metric_EqW_Results_C0vC1vC3vC4.pdf", width=17, height=9) 
 boxplot(MCI_EqW_C0vC1_unpacked, MCI_EqW_C0vC2_unpacked, MCI_EqW_C0vC3_unpacked, MCI_EqW_C0vC4_unpacked, xlab= "Phylogenetic hypotheses", ylab = "Mutual clustering information tree distance metric", col=c("#E69F00", "#56B4E9", "#009E73", "#F0E442"), names= c("C0vsC1", "C0vC2", "C0vC3", "C0vC4"))

 dev.off()


## Shared phylogenetic Information tree distance metric boxplots:
 boxplot(Nye1_EqW_C0vC1_unpacked, Nye1_EqW_C0vC2_unpacked, Nye1_EqW_C0vC3_unpacked, Nye1_EqW_C0vC4_unpacked, xlab= "Phylogenetic hypotheses", ylab = "Nye et al. (2006) similarity metric, normalised against the average number of splits in the two topologies", col=c("#E69F00", "#56B4E9", "#009E73", "#F0E442"), names= c("C0vsC1", "C0vC2", "C0vC3", "C0vC4"))

 pdf("output/TreeDistance/Nye_et_al_1_Metric_Results_EqW_C0vC1vC3vC4.pdf", width=17, height=9) 
 boxplot(Nye1_EqW_C0vC1_unpacked, Nye1_EqW_C0vC2_unpacked, Nye1_EqW_C0vC3_unpacked, Nye1_EqW_C0vC4_unpacked, xlab= "Phylogenetic hypotheses", ylab = "Nye et al. (2006) similarity metric, normalised against the average number of splits in the two topologies", col=c("#E69F00", "#56B4E9", "#009E73", "#F0E442"), names= c("C0vsC1", "C0vC2", "C0vC3", "C0vC4"))

 dev.off()


## Shared phylogenetic Information tree distance metric boxplots:
 boxplot(Nye2_EqW_C0vC1_unpacked, Nye2_EqW_C0vC2_unpacked, Nye2_EqW_C0vC3_unpacked, Nye2_EqW_C0vC4_unpacked, xlab= "Phylogenetic hypotheses", ylab = "Nye et al. (2006) similarity metric, normalised against the maximum similarity possible for the unconstrained topology", col=c("#E69F00", "#56B4E9", "#009E73", "#F0E442"), names= c("C0vsC1", "C0vC2", "C0vC3", "C0vC4"))

 pdf("output/TreeDistance/Nye_et_al_2_Metric_Results_EqW_C0vC1vC3vC4.pdf", width=17, height=9) 
 boxplot(Nye2_EqW_C0vC1_unpacked, Nye2_EqW_C0vC2_unpacked, Nye2_EqW_C0vC3_unpacked, Nye2_EqW_C0vC4_unpacked, xlab= "Phylogenetic hypotheses", ylab = "Nye et al. (2006) similarity metric, normalised against the maximum similarity possible for the unconstrained topology", col=c("#E69F00", "#56B4E9", "#009E73", "#F0E442"), names= c("C0vsC1", "C0vC2", "C0vC3", "C0vC4"))

 dev.off()

 beep(1)


# ------------------------------------------------------------------------------------------------- #
# 2.4 PDFs of the extended implied weights (k=15) analyses tree distance/similarity metric boxplots #
# ------------------------------------------------------------------------------------------------- #

### Creates boxplots comparing the tree distance/similarity metrics of the different constraint analyses against the unconstrained analysis, creates PDFs (colours chosen so that people with colourblindness can differentiate the boxes easily):

## Shared phylogenetic Information tree distance metric boxplots:
 boxplot(SPI_EIWk15_C0vC1_unpacked, SPI_EIWk15_C0vC2_unpacked, SPI_EIWk15_C0vC3_unpacked, SPI_EIWk15_C0vC4_unpacked, xlab= "Phylogenetic hypotheses", ylab = "Shared phylogenetic Information  tree distance metric", col=c("#E69F00", "#56B4E9", "#009E73", "#F0E442"), names= c("C0vsC1", "C0vC2", "C0vC3", "C0vC4"))

pdf("output/TreeDistance/Shared_Phylogenetic_Information_Metric_EIWk15_Results_C0vC1vC3vC4.pdf", width=17, height=9) 
 boxplot(SPI_EIWk15_C0vC1_unpacked, SPI_EIWk15_C0vC2_unpacked, SPI_EIWk15_C0vC3_unpacked, SPI_EIWk15_C0vC4_unpacked, xlab= "Phylogenetic hypotheses", ylab = "Shared phylogenetic Information  tree distance metric", col=c("#E69F00", "#56B4E9", "#009E73", "#F0E442"), names= c("C0vsC1", "C0vC2", "C0vC3", "C0vC4"))

 dev.off()


## Mutual clustering information tree distance metric boxplots:
 boxplot(MCI_EIWk15_C0vC1_unpacked, MCI_EIWk15_C0vC2_unpacked, MCI_EIWk15_C0vC3_unpacked, MCI_EIWk15_C0vC4_unpacked, xlab= "Phylogenetic hypotheses", ylab = "Mutual clustering information tree distance metric", col=c("#E69F00", "#56B4E9", "#009E73", "#F0E442"), names= c("C0vsC1", "C0vC2", "C0vC3", "C0vC4"))

 pdf("output/TreeDistance/Mutal_Clustering_Information_Metric_EIWk15_Results_C0vC1vC3vC4.pdf", width=17, height=9) 
 boxplot(MCI_EIWk15_C0vC1_unpacked, MCI_EIWk15_C0vC2_unpacked, MCI_EIWk15_C0vC3_unpacked, MCI_EIWk15_C0vC4_unpacked, xlab= "Phylogenetic hypotheses", ylab = "Mutual clustering information tree distance metric", col=c("#E69F00", "#56B4E9", "#009E73", "#F0E442"), names= c("C0vsC1", "C0vC2", "C0vC3", "C0vC4"))

 dev.off()


## Shared phylogenetic Information tree distance metric boxplots:
 boxplot(Nye1_EIWk15_C0vC1_unpacked, Nye1_EIWk15_C0vC2_unpacked, Nye1_EIWk15_C0vC3_unpacked, Nye1_EIWk15_C0vC4_unpacked, xlab= "Phylogenetic hypotheses", ylab = "Nye et al. (2006) similarity metric, normalised against the average number of splits in the two topologies", col=c("#E69F00", "#56B4E9", "#009E73", "#F0E442"), names= c("C0vsC1", "C0vC2", "C0vC3", "C0vC4"))

 pdf("output/TreeDistance/Nye_et_al_1_Metric_Results_EIWk15_C0vC1vC3vC4.pdf", width=17, height=9) 
 boxplot(Nye1_EIWk15_C0vC1_unpacked, Nye1_EIWk15_C0vC2_unpacked, Nye1_EIWk15_C0vC3_unpacked, Nye1_EIWk15_C0vC4_unpacked, xlab= "Phylogenetic hypotheses", ylab = "Nye et al. (2006) similarity metric, normalised against the average number of splits in the two topologies", col=c("#E69F00", "#56B4E9", "#009E73", "#F0E442"), names= c("C0vsC1", "C0vC2", "C0vC3", "C0vC4"))

 dev.off()


## Shared phylogenetic Information tree distance metric boxplots:
 boxplot(Nye2_EIWk15_C0vC1_unpacked, Nye2_EIWk15_C0vC2_unpacked, Nye2_EIWk15_C0vC3_unpacked, Nye2_EIWk15_C0vC4_unpacked, xlab= "Phylogenetic hypotheses", ylab = "Nye et al. (2006) similarity metric, normalised against the maximum similarity possible for the unconstrained topology", col=c("#E69F00", "#56B4E9", "#009E73", "#F0E442"), names= c("C0vsC1", "C0vC2", "C0vC3", "C0vC4"))

 pdf("output/TreeDistance/Nye_et_al_2_Metric_Results_EIWk15_C0vC1vC3vC4.pdf", width=17, height=9) 
 boxplot(Nye2_EIWk15_C0vC1_unpacked, Nye2_EIWk15_C0vC2_unpacked, Nye2_EIWk15_C0vC3_unpacked, Nye2_EIWk15_C0vC4_unpacked, xlab= "Phylogenetic hypotheses", ylab = "Nye et al. (2006) similarity metric, normalised against the maximum similarity possible for the unconstrained topology", col=c("#E69F00", "#56B4E9", "#009E73", "#F0E442"), names= c("C0vsC1", "C0vC2", "C0vC3", "C0vC4"))

 dev.off()

 beep(1)


# ======================== #
# 0.1 Parallelisation ends #
# ======================== #

# Stop cluster
 stopImplicitCluster()


# ========= #
# 0.2 Done! #
# ========= #

### Noise when script has finished (useful if you are running this in the background):
 library(beepr)

 beep(8)