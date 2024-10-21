#!/usr/bin/env Rscript

# ======================================================================================= #
#                                                                                         #
#  Stratigraphic phylogenetic congruence script - version 2.0 (October, 2024).            #
#                                                                                         #
#                                                                                         #
#  Purpose:                                                                               #
#  Calculates the stratigraphic congruence of competing topologies.                       #
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
#  Origin:                                                                                #
#  Version 1.0 was a modification of the StatCongruenceCode.R script created by Wright &  #
#  Lloyd (2020) Palaeontology 63: 997-1006.                                               #
#  See: https://onlinelibrary.wiley.com/doi/abs/10.1111/pala.12500                        #
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
#  If you plan to use or amend this script, please ensure that Wright & Lloyd (2020) are  #
#  appropriately credited as the originators.                                             #
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

### Required libraries:
 library(ape)
 library(beepr)
 library(geoscale)
 library(TreeTools)

clusterCall(cluster, function () {
library(ape)
})

### Loads the resultant cladograms from the phylogenetic analyses, and creates a Nexus format using Ape:
 EqW_C0_Trees <- ReadTntTree("dat/Trees/xiw_equal_C0_RC_AllTrees.tre")
 ape::write.tree(EqW_C0_Trees, file='dat/Trees/xiw_equal_C0_RC_AllTrees.txt') 

 EqW_C1_Trees <- ReadTntTree("dat/Trees/xiw_equal_C1_RC_AllTrees.tre")
 ape::write.tree(EqW_C1_Trees, file='dat/Trees/xiw_equal_C1_RC_AllTrees.txt') 

 EqW_C2_Trees <- ReadTntTree("dat/Trees/xiw_equal_C2_RC_AllTrees.tre")
 ape::write.tree(EqW_C2_Trees, file='dat/Trees/xiw_equal_C2_RC_AllTrees.txt')

 EqW_C3_Trees <- ReadTntTree("dat/Trees/xiw_equal_C3_RC_AllTrees.tre")
 ape::write.tree(EqW_C3_Trees, file='dat/Trees/xiw_equal_C3_RC_AllTrees.txt') 

 EqW_C4_Trees <- ReadTntTree("dat/Trees/xiw_equal_C4_RC_AllTrees.tre")
 ape::write.tree(EqW_C4_Trees, file='dat/Trees/xiw_equal_C4_RC_AllTrees.txt') 

 EIWk15_C0_Trees <- ReadTntTree("dat/Trees/xiw_k15_C0_RC_AllTrees.tre")
 ape::write.tree(EIWk15_C0_Trees, file='dat/Trees/xiw_k15_C0_RC_AllTrees.txt') 

 EIWk15_C1_Trees <- ReadTntTree("dat/Trees/xiw_k15_C1_RC_AllTrees.tre")
 ape::write.tree(EIWk15_C1_Trees, file='dat/Trees/xiw_k15_C1_RC_AllTrees.txt') 

 EIWk15_C2_Trees <- ReadTntTree("dat/Trees/xiw_k15_C2_RC_AllTrees.tre")
 ape::write.tree(EIWk15_C2_Trees, file='dat/Trees/xiw_k15_C2_RC_AllTrees.txt') 

 EIWk15_C3_Trees <- ReadTntTree("dat/Trees/xiw_k15_C3_RC_AllTrees.tre")
 ape::write.tree(EIWk15_C3_Trees, file='dat/Trees/xiw_k15_C3_RC_AllTrees.txt') 

 EIWk15_C4_Trees <- ReadTntTree("dat/Trees/xiw_k15_C4_RC_AllTrees.tre")
 ape::write.tree(EIWk15_C4_Trees, file='dat/Trees/xiw_k15_C4_RC_AllTrees.txt') 


### Randomly sample 100 cladograms from the input sample, and creates NEXUS format using Ape:
 EqW_C0_Trees_subsample<-sample(EqW_C0_Trees,size=100)
 ape::write.tree(EqW_C0_Trees_subsample, file='dat/Trees/xiw_equal_C0_RC_Subsample_Trees.txt') 

 EqW_C1_Trees_subsample<-sample(EqW_C1_Trees,size=100)
 ape::write.tree(EqW_C1_Trees_subsample, file='dat/Trees/xiw_equal_C1_RC_Subsample_Trees.txt') 

 EqW_C2_Trees_subsample<-sample(EqW_C2_Trees,size=100)
 ape::write.tree(EqW_C2_Trees_subsample, file='dat/Trees/xiw_equal_C2_RC_Subsample_Trees.txt') 

 EqW_C3_Trees_subsample<-sample(EqW_C3_Trees,size=100)
 ape::write.tree(EqW_C3_Trees_subsample, file='dat/Trees/xiw_equal_C3_RC_Subsample_Trees.txt') 

 EqW_C4_Trees_subsample<-sample(EqW_C4_Trees,size=100)
 ape::write.tree(EqW_C4_Trees_subsample, file='dat/Trees/xiw_equal_C4_RC_Subsample_Trees.txt')

 EIWk15_C0_Trees_subsample<-sample(EIWk15_C0_Trees,size=100)
 ape::write.tree(EIWk15_C0_Trees_subsample, file='dat/Trees/xiw_k15_C0_RC_Subsample_Trees.txt') 

 EIWk15_C1_Trees_subsample<-sample(EIWk15_C1_Trees,size=100)
 ape::write.tree(EIWk15_C1_Trees_subsample, file='dat/Trees/xiw_k15_C1_RC_Subsample_Trees.txt') 

 EIWk15_C2_Trees_subsample<-sample(EIWk15_C2_Trees,size=100)
 ape::write.tree(EIWk15_C2_Trees_subsample, file='dat/Trees/xiw_k15_C2_RC_Subsample_Trees.txt')

 EIWk15_C3_Trees_subsample<-sample(EIWk15_C3_Trees,size=100)
 ape::write.tree(EIWk15_C3_Trees_subsample, file='dat/Trees/xiw_k15_C3_RC_Subsample_Trees.txt') 

 EIWk15_C4_Trees_subsample<-sample(EIWk15_C4_Trees,size=100)
 ape::write.tree(EIWk15_C4_Trees_subsample, file='dat/Trees/xiw_k15_C4_RC_Subsample_Trees.txt') 


### Ladderise the random 100 cladogram sample using Ape, and convert them into a multiphylo object:
 EqW_C0_Trees_subsample <- lapply(EqW_C0_Trees_subsample, function(x) {ladderize(x)})
 class(EqW_C0_Trees_subsample) <- "multiPhylo"

 EqW_C1_Trees_subsample <- lapply(EqW_C1_Trees_subsample, function(x) {ladderize(x)})
 class(EqW_C1_Trees_subsample) <- "multiPhylo"

 EqW_C2_Trees_subsample <- lapply(EqW_C2_Trees_subsample, function(x) {ladderize(x)})
 class(EqW_C2_Trees_subsample) <- "multiPhylo"

 EqW_C3_Trees_subsample <- lapply(EqW_C3_Trees_subsample, function(x) {ladderize(x)})
 class(EqW_C3_Trees_subsample) <- "multiPhylo"

 EqW_C4_Trees_subsample <- lapply(EqW_C4_Trees_subsample, function(x) {ladderize(x)})
 class(EqW_C4_Trees_subsample) <- "multiPhylo"

 EIWk15_C0_Trees_subsample <- lapply(EIWk15_C0_Trees_subsample, function(x) {ladderize(x)})
 class(EIWk15_C0_Trees_subsample) <- "multiPhylo"

 EIWk15_C1_Trees_subsample <- lapply(EIWk15_C1_Trees_subsample, function(x) {ladderize(x)})
 class(EIWk15_C1_Trees_subsample) <- "multiPhylo"

 EIWk15_C2_Trees_subsample <- lapply(EIWk15_C2_Trees_subsample, function(x) {ladderize(x)})
 class(EIWk15_C2_Trees_subsample) <- "multiPhylo"

 EIWk15_C3_Trees_subsample <- lapply(EIWk15_C3_Trees_subsample, function(x) {ladderize(x)})
 class(EIWk15_C3_Trees_subsample) <- "multiPhylo"

 EIWk15_C4_Trees_subsample <- lapply(EIWk15_C4_Trees_subsample, function(x) {ladderize(x)})
 class(EIWk15_C4_Trees_subsample) <- "multiPhylo"

 beep(1)


# ==================================== #
# 2. Stratigraphic congruence analyses #
# ==================================== #

### Required libraries:
 library(beepr)
 library(geiger)
 library(TeachingDemos)
 library(graphics)
 library(strap)

clusterCall(cluster, function () {
library(strap)
})

### Loads in Age file:
 CrocAges <- read.table("dat/Ages/CrocAges.txt", header=T) 


# ---------------------------------------------------------------------- #
# 2.1 Comparison of the equal weights analyses, using the 'equal' method #
# ---------------------------------------------------------------------- #

### Stratigraphic congruence metrics on the 100 random no-constraints analysis input cladograms plus 1000 randomly generated trees:
 StratPhyloCongruenceResults_EqW_C0 <- StratPhyloCongruence(trees = c(EqW_C0_Trees_subsample), ages = CrocAges, rlen = 1, method = "equal", samp.perm = 1000, hard = TRUE, randomly.sample.ages = FALSE, fix.topology = TRUE, fix.outgroup = TRUE, outgroup.taxon = "Gracilisuchus_stipanicicorum")


### writes C0 results to file:
 txtStart("output/StratCongruenceResults_EqW/EqW_C0_StratCongResults_Subsample_Trees.txt")
 StratPhyloCongruenceResults_EqW_C0
 txtStop()


### Stratigraphic congruence metrics on the 100 random constraint one analysis input cladograms plus 1000 randomly generated trees:
 StratPhyloCongruenceResults_EqW_C1 <- StratPhyloCongruence(trees = c(EqW_C1_Trees_subsample), ages = CrocAges, rlen = 1, method = "equal", samp.perm = 1000, rand.perm = 1000, hard = TRUE, randomly.sample.ages = FALSE, fix.topology = TRUE, fix.outgroup = TRUE, outgroup.taxon = "Gracilisuchus_stipanicicorum")

### writes C1 results to file:
 txtStart("output/StratCongruenceResults_EqW/EqW_C1_StratCongResults_Subsample_Trees.txt")
 StratPhyloCongruenceResults_EqW_C1
 txtStop()


### Stratigraphic congruence metrics on the 100 random constraint two analysis input cladograms plus 1000 randomly generated trees:
 StratPhyloCongruenceResults_EqW_C3 <- StratPhyloCongruence(trees = c(EqW_C3_Trees_subsample), ages = CrocAges, rlen = 1, method = "equal", samp.perm = 1000, rand.perm = 1000, hard = TRUE, randomly.sample.ages = FALSE, fix.topology = TRUE, fix.outgroup = TRUE, outgroup.taxon = "Gracilisuchus_stipanicicorum")

### writes C3 results to file:
 txtStart("output/StratCongruenceResults_EqW/EqW_C3_StratCongResults_Subsample_Trees.txt")
 StratPhyloCongruenceResults_EqW_C3
 txtStop()


### Stratigraphic congruence metrics on the 100 random constraint two analysis input cladograms plus 1000 randomly generated trees:
 StratPhyloCongruenceResults_EqW_C4 <- StratPhyloCongruence(trees = c(EqW_C4_Trees_subsample), ages = CrocAges, rlen = 1, method = "equal", samp.perm = 1000, rand.perm = 1000, hard = TRUE, randomly.sample.ages = FALSE, fix.topology = TRUE, fix.outgroup = TRUE, outgroup.taxon = "Gracilisuchus_stipanicicorum")

### writes C4 results to file:
 txtStart("output/StratCongruenceResults_EqW/EqW_C4_StratCongResults_Subsample_Trees.txt")
 StratPhyloCongruenceResults_EqW_C4
 txtStop()



### Creates boxplots comparing the different constraint analyses and creates PDFs (colours chosen so that people with colourblindness can differentiate the boxes easily):

## SCI boxplots:
boxplot(StratPhyloCongruenceResults_EqW_C0$input.tree.results[, "SCI"], StratPhyloCongruenceResults_EqW_C1$input.tree.results[, "SCI"], StratPhyloCongruenceResults_EqW_C3$input.tree.results[, "SCI"], StratPhyloCongruenceResults_EqW_C4$input.tree.results[, "SCI"], xlab= "Phylogenetic hypotheses", ylab = "Stratigraphic Consistency Index", col=c("#E69F00", "#56B4E9", "#009E73", "#F0E442"), names= c("No constraints", "Crocodyliform constraint", "Neosuchian constraint", "Tethysuchian constraint"))

pdf("output/StratCongruenceResults_EqW/Boxplots/EqW_SCI_Results_C0vC1vC3vC4.pdf", width=17, height=9) 
boxplot(StratPhyloCongruenceResults_EqW_C0$input.tree.results[, "SCI"], StratPhyloCongruenceResults_EqW_C1$input.tree.results[, "SCI"], StratPhyloCongruenceResults_EqW_C3$input.tree.results[, "SCI"], StratPhyloCongruenceResults_EqW_C4$input.tree.results[, "SCI"], xlab= "Phylogenetic hypotheses", ylab = "Stratigraphic Consistency Index", col=c("#E69F00", "#56B4E9", "#009E73", "#F0E442"), names= c("No constraints", "Crocodyliform constraint", "Neosuchian constraint", "Tethysuchian constraint"))
dev.off()

## RCI boxplots:
boxplot(StratPhyloCongruenceResults_EqW_C0$input.tree.results[, "RCI"], StratPhyloCongruenceResults_EqW_C1$input.tree.results[, "RCI"], StratPhyloCongruenceResults_EqW_C3$input.tree.results[, "RCI"], StratPhyloCongruenceResults_EqW_C4$input.tree.results[, "RCI"], xlab= "Phylogenetic hypotheses", ylab = "Relative Completeness Index", col=c("#E69F00", "#56B4E9", "#009E73", "#F0E442"), names= c("No constraints", "Crocodyliform constraint", "Neosuchian constraint", "Tethysuchian constraint"))

pdf("output/StratCongruenceResults_EqW/Boxplots/EqW_RCI_Results_C0vC1vC3vC4.pdf", width=17, height=9) 
boxplot(StratPhyloCongruenceResults_EqW_C0$input.tree.results[, "RCI"], StratPhyloCongruenceResults_EqW_C1$input.tree.results[, "RCI"], StratPhyloCongruenceResults_EqW_C3$input.tree.results[, "RCI"], StratPhyloCongruenceResults_EqW_C4$input.tree.results[, "RCI"], xlab= "Phylogenetic hypotheses", ylab = "Relative Completeness Index", col=c("#E69F00", "#56B4E9", "#009E73", "#F0E442"), names= c("No constraints", "Crocodyliform constraint", "Neosuchian constraint", "Tethysuchian constraint"))
dev.off()

## GER boxplots:
boxplot(StratPhyloCongruenceResults_EqW_C0$input.tree.results[, "GER"], StratPhyloCongruenceResults_EqW_C1$input.tree.results[, "GER"], StratPhyloCongruenceResults_EqW_C3$input.tree.results[, "GER"], StratPhyloCongruenceResults_EqW_C4$input.tree.results[, "GER"], xlab= "Phylogenetic hypotheses", ylab = "Gap Excess Ratio (GER)", col=c("#E69F00", "#56B4E9", "#009E73", "#F0E442"), names= c("No constraints", "Crocodyliform constraint", "Neosuchian constraint", "Tethysuchian constraint"))

pdf("output/StratCongruenceResults_EqW/Boxplots/EqWGER_Results_C0vC1vC3vC4.pdf", width=17, height=9) 
boxplot(StratPhyloCongruenceResults_EqW_C0$input.tree.results[, "GER"], StratPhyloCongruenceResults_EqW_C1$input.tree.results[, "GER"], StratPhyloCongruenceResults_EqW_C3$input.tree.results[, "GER"], StratPhyloCongruenceResults_EqW_C4$input.tree.results[, "GER"], xlab= "Phylogenetic hypotheses", ylab = "Gap Excess Ratio (GER)", col=c("#E69F00", "#56B4E9", "#009E73", "#F0E442"), names= c("No constraints", "Crocodyliform constraint", "Neosuchian constraint", "Tethysuchian constraint"))
dev.off()

## MSM* boxplots:
boxplot(StratPhyloCongruenceResults_EqW_C0$input.tree.results[, "MSM*"], StratPhyloCongruenceResults_EqW_C1$input.tree.results[, "MSM*"], StratPhyloCongruenceResults_EqW_C3$input.tree.results[, "MSM*"], StratPhyloCongruenceResults_EqW_C4$input.tree.results[, "MSM*"], xlab= "Phylogenetic hypotheses", ylab = "Manhattan Stratigraphic Measure*", col=c("#E69F00", "#56B4E9", "#009E73", "#F0E442"), names= c("No constraints", "Crocodyliform constraint", "Neosuchian constraint", "Tethysuchian constraint"))

pdf("output/StratCongruenceResults_EqW/Boxplots/EqW_MSM*_Results_C0vC1vC3vC4.pdf", width=17, height=9) 
boxplot(StratPhyloCongruenceResults_EqW_C0$input.tree.results[, "MSM*"], StratPhyloCongruenceResults_EqW_C1$input.tree.results[, "MSM*"], StratPhyloCongruenceResults_EqW_C3$input.tree.results[, "MSM*"], StratPhyloCongruenceResults_EqW_C4$input.tree.results[, "MSM*"], xlab= "Phylogenetic hypotheses", ylab = "Manhattan Stratigraphic Measure*", col=c("#E69F00", "#56B4E9", "#009E73", "#F0E442"), names= c("No constraints", "Crocodyliform constraint", "Neosuchian constraint", "Tethysuchian constraint"))
dev.off()

## GER* boxplots:
boxplot(StratPhyloCongruenceResults_EqW_C0$input.tree.results[, "GER*"], StratPhyloCongruenceResults_EqW_C1$input.tree.results[, "GER*"], StratPhyloCongruenceResults_EqW_C3$input.tree.results[, "GER*"], StratPhyloCongruenceResults_EqW_C4$input.tree.results[, "GER*"], xlab= "Phylogenetic hypotheses", ylab = "Gap Excess Ratio (GER*)", col=c("#E69F00", "#56B4E9", "#009E73", "#F0E442"), names= c("No constraints", "Crocodyliform constraint", "Neosuchian constraint", "Tethysuchian constraint"))

pdf("output/StratCongruenceResults_EqW/Boxplots/EqW_GER*_Results_C0vC1vC3vC4.pdf", width=17, height=9) 
boxplot(StratPhyloCongruenceResults_EqW_C0$input.tree.results[, "GER*"], StratPhyloCongruenceResults_EqW_C1$input.tree.results[, "GER*"], StratPhyloCongruenceResults_EqW_C3$input.tree.results[, "GER*"], StratPhyloCongruenceResults_EqW_C4$input.tree.results[, "GER*"], xlab= "Phylogenetic hypotheses", ylab = "Gap Excess Ratio (GER*)", col=c("#E69F00", "#56B4E9", "#009E73", "#F0E442"), names= c("No constraints", "Crocodyliform constraint", "Neosuchian constraint", "Tethysuchian constraint"))
dev.off()

## GERt boxplots:
boxplot(StratPhyloCongruenceResults_EqW_C0$input.tree.results[, "GERt"], StratPhyloCongruenceResults_EqW_C1$input.tree.results[, "GERt"], StratPhyloCongruenceResults_EqW_C3$input.tree.results[, "GERt"], StratPhyloCongruenceResults_EqW_C4$input.tree.results[, "GERt"], xlab= "Phylogenetic hypotheses", ylab = "Gap Excess Ratio (GERt)", col=c("#E69F00", "#56B4E9", "#009E73", "#F0E442"), names= c("No constraints", "Crocodyliform constraint", "Neosuchian constraint", "Tethysuchian constraint"))

pdf("output/StratCongruenceResults_EqW/Boxplots/EqW_GERt_Results_C0vC1vC3vC4.pdf", width=17, height=9) 
boxplot(StratPhyloCongruenceResults_EqW_C0$input.tree.results[, "GERt"], StratPhyloCongruenceResults_EqW_C1$input.tree.results[, "GERt"], StratPhyloCongruenceResults_EqW_C3$input.tree.results[, "GERt"], StratPhyloCongruenceResults_EqW_C4$input.tree.results[, "GERt"], xlab= "Phylogenetic hypotheses", ylab = "Gap Excess Ratio (GERt)", col=c("#E69F00", "#56B4E9", "#009E73", "#F0E442"), names= c("No constraints", "Crocodyliform constraint", "Neosuchian constraint", "Tethysuchian constraint"))
dev.off()

## MIG boxplots:
boxplot(StratPhyloCongruenceResults_EqW_C0$input.tree.results[, "MIG"], StratPhyloCongruenceResults_EqW_C1$input.tree.results[, "MIG"], StratPhyloCongruenceResults_EqW_C3$input.tree.results[, "MIG"], StratPhyloCongruenceResults_EqW_C4$input.tree.results[, "MIG"], xlab= "Phylogenetic hypotheses", ylab = "Minimum Implied Gap", col=c("#E69F00", "#56B4E9", "#009E73", "#F0E442"), names= c("No constraints", "Crocodyliform constraint", "Neosuchian constraint", "Tethysuchian constraint"))

pdf("output/StratCongruenceResults_EqW/Boxplots/EqW_MIG_Results_C0vC1vC3vC4.pdf", width=17, height=9) 
boxplot(StratPhyloCongruenceResults_EqW_C0$input.tree.results[, "MIG"], StratPhyloCongruenceResults_EqW_C1$input.tree.results[, "MIG"], StratPhyloCongruenceResults_EqW_C3$input.tree.results[, "MIG"], StratPhyloCongruenceResults_EqW_C4$input.tree.results[, "MIG"], xlab= "Phylogenetic hypotheses", ylab = "Minimum Implied Gap", col=c("#E69F00", "#56B4E9", "#009E73", "#F0E442"), names= c("No constraints", "Crocodyliform constraint", "Neosuchian constraint", "Tethysuchian constraint"))
dev.off()

 beep(1)


# ---------------------------------------------------------------------------------------- #
# 2.2 Comparison of the extended implied weights (k=15) analyses, using the 'equal' method #
# ---------------------------------------------------------------------------------------- #


### Stratigraphic congruence metrics on the 100 random no-constraints analysis input cladograms plus 1000 randomly generated trees:
 StratPhyloCongruenceResults_EIWk15_C0 <- StratPhyloCongruence(trees = c(EIWk15_C0_Trees_subsample), ages = CrocAges, rlen = 1, method = "equal", samp.perm = 1000, rand.perm = 1000, hard = TRUE, randomly.sample.ages = FALSE, fix.topology = TRUE, fix.outgroup = TRUE, outgroup.taxon = "Gracilisuchus_stipanicicorum")

### writes results to file:
 txtStart("output/StratCongruenceResults_EIWk15/EIWk15_C0_StratCongResults_Subsample_Trees.txt")
 StratPhyloCongruenceResults_EIWk15_C0
 txtStop()


### Stratigraphic congruence metrics on the 100 random constraint one analysis input cladograms plus 1000 randomly generated trees:
 StratPhyloCongruenceResults_EIWk15_C1 <- StratPhyloCongruence(trees = c(EIWk15_C1_Trees_subsample), ages = CrocAges, rlen = 1, method = "equal", samp.perm = 1000, rand.perm = 1000, hard = TRUE, randomly.sample.ages = FALSE, fix.topology = TRUE, fix.outgroup = TRUE, outgroup.taxon = "Gracilisuchus_stipanicicorum")

### writes results to file:
 txtStart("output/StratCongruenceResults_EIWk15/EIWk15_C1_StratCongResults_Subsample_Trees.txt")
 StratPhyloCongruenceResults_EIWk15_C1
 txtStop()


### Stratigraphic congruence metrics on the 100 random constraint two analysis input cladograms plus 1000 randomly generated trees:
 StratPhyloCongruenceResults_EIWk15_C3 <- StratPhyloCongruence(trees = c(EIWk15_C3_Trees_subsample), ages = CrocAges, rlen = 1, method = "equal", samp.perm = 1000, rand.perm = 1000, hard = TRUE, randomly.sample.ages = FALSE, fix.topology = TRUE, fix.outgroup = TRUE, outgroup.taxon = "Gracilisuchus_stipanicicorum")

### writes results to file:
 txtStart("output/StratCongruenceResults_EIWk15/EIWk15_C3_StratCongResults_Subsample_Trees.txt")
 StratPhyloCongruenceResults_EIWk15_C3
 txtStop()


### Stratigraphic congruence metrics on the 100 random constraint two analysis input cladograms plus 1000 randomly generated trees:
 StratPhyloCongruenceResults_EIWk15_C4 <- StratPhyloCongruence(trees = c(EIWk15_C4_Trees_subsample), ages = CrocAges, rlen = 1, method = "equal", samp.perm = 1000, rand.perm = 1000, hard = TRUE, randomly.sample.ages = FALSE, fix.topology = TRUE, fix.outgroup = TRUE, outgroup.taxon = "Gracilisuchus_stipanicicorum")

### writes results to file:
 txtStart("output/StratCongruenceResults_EIWk15/EIWk15_C4_StratCongResults_Subsample_Trees.txt")
 StratPhyloCongruenceResults_EIWk15_C4
 txtStop()


### Creates boxplots comparing the different constraint analyses and creates PDFs (colours chosen so that people with colourblindness can differentiate the boxes easily):

## SCI boxplots:
boxplot(StratPhyloCongruenceResults_EIWk15_C0$input.tree.results[, "SCI"], StratPhyloCongruenceResults_EIWk15_C1$input.tree.results[, "SCI"], StratPhyloCongruenceResults_EIWk15_C3$input.tree.results[, "SCI"], StratPhyloCongruenceResults_EIWk15_C4$input.tree.results[, "SCI"], xlab= "Phylogenetic hypotheses", ylab = "Stratigraphic Consistency Index", col=c("#E69F00", "#56B4E9", "#009E73", "#F0E442"), names= c("No constraints", "Crocodyliform constraint", "Neosuchian constraint", "Tethysuchian constraint"))

pdf("output/StratCongruenceResults_EIWk15/Boxplots/EIWk15_SCI_Results_C0vC1vC3vC4.pdf", width=17, height=9) 
boxplot(StratPhyloCongruenceResults_EIWk15_C0$input.tree.results[, "SCI"], StratPhyloCongruenceResults_EIWk15_C1$input.tree.results[, "SCI"], StratPhyloCongruenceResults_EIWk15_C3$input.tree.results[, "SCI"], StratPhyloCongruenceResults_EIWk15_C4$input.tree.results[, "SCI"], xlab= "Phylogenetic hypotheses", ylab = "Stratigraphic Consistency Index", col=c("#E69F00", "#56B4E9", "#009E73", "#F0E442"), names= c("No constraints", "Crocodyliform constraint", "Neosuchian constraint", "Tethysuchian constraint"))
dev.off()

## RCI boxplots:
boxplot(StratPhyloCongruenceResults_EIWk15_C0$input.tree.results[, "RCI"], StratPhyloCongruenceResults_EIWk15_C1$input.tree.results[, "RCI"], StratPhyloCongruenceResults_EIWk15_C3$input.tree.results[, "RCI"], StratPhyloCongruenceResults_EIWk15_C4$input.tree.results[, "RCI"], xlab= "Phylogenetic hypotheses", ylab = "Relative Completeness Index", col=c("#E69F00", "#56B4E9", "#009E73", "#F0E442"), names= c("No constraints", "Crocodyliform constraint", "Neosuchian constraint", "Tethysuchian constraint"))

pdf("output/StratCongruenceResults_EIWk15/Boxplots/EIWk15_RCI_Results_C0vC1vC3vC4.pdf", width=17, height=9) 
boxplot(StratPhyloCongruenceResults_EIWk15_C0$input.tree.results[, "RCI"], StratPhyloCongruenceResults_EIWk15_C1$input.tree.results[, "RCI"], StratPhyloCongruenceResults_EIWk15_C3$input.tree.results[, "RCI"], StratPhyloCongruenceResults_EIWk15_C4$input.tree.results[, "RCI"], xlab= "Phylogenetic hypotheses", ylab = "Relative Completeness Index", col=c("#E69F00", "#56B4E9", "#009E73", "#F0E442"), names= c("No constraints", "Crocodyliform constraint", "Neosuchian constraint", "Tethysuchian constraint"))
dev.off()

## GER boxplots:
boxplot(StratPhyloCongruenceResults_EIWk15_C0$input.tree.results[, "GER"], StratPhyloCongruenceResults_EIWk15_C1$input.tree.results[, "GER"], StratPhyloCongruenceResults_EIWk15_C3$input.tree.results[, "GER"], StratPhyloCongruenceResults_EIWk15_C4$input.tree.results[, "GER"], xlab= "Phylogenetic hypotheses", ylab = "Gap Excess Ratio (GER)", col=c("#E69F00", "#56B4E9", "#009E73", "#F0E442"), names= c("No constraints", "Crocodyliform constraint", "Neosuchian constraint", "Tethysuchian constraint"))

pdf("output/StratCongruenceResults_EIWk15/Boxplots/EIWk15_GER_Results_C0vC1vC3vC4.pdf", width=17, height=9) 
boxplot(StratPhyloCongruenceResults_EIWk15_C0$input.tree.results[, "GER"], StratPhyloCongruenceResults_EIWk15_C1$input.tree.results[, "GER"], StratPhyloCongruenceResults_EIWk15_C3$input.tree.results[, "GER"], StratPhyloCongruenceResults_EIWk15_C4$input.tree.results[, "GER"], xlab= "Phylogenetic hypotheses", ylab = "Gap Excess Ratio (GER)", col=c("#E69F00", "#56B4E9", "#009E73", "#F0E442"), names= c("No constraints", "Crocodyliform constraint", "Neosuchian constraint", "Tethysuchian constraint"))
dev.off()

## MSM* boxplots:
boxplot(StratPhyloCongruenceResults_EIWk15_C0$input.tree.results[, "MSM*"], StratPhyloCongruenceResults_EIWk15_C1$input.tree.results[, "MSM*"], StratPhyloCongruenceResults_EIWk15_C3$input.tree.results[, "MSM*"], StratPhyloCongruenceResults_EIWk15_C4$input.tree.results[, "MSM*"], xlab= "Phylogenetic hypotheses", ylab = "Manhattan Stratigraphic Measure*", col=c("#E69F00", "#56B4E9", "#009E73", "#F0E442"), names= c("No constraints", "Crocodyliform constraint", "Neosuchian constraint", "Tethysuchian constraint"))

pdf("output/StratCongruenceResults_EIWk15/Boxplots/EIWk15_MSM*_Results_C0vC1vC3vC4.pdf", width=17, height=9) 
boxplot(StratPhyloCongruenceResults_EIWk15_C0$input.tree.results[, "MSM*"], StratPhyloCongruenceResults_EIWk15_C1$input.tree.results[, "MSM*"], StratPhyloCongruenceResults_EIWk15_C3$input.tree.results[, "MSM*"], StratPhyloCongruenceResults_EIWk15_C4$input.tree.results[, "MSM*"], xlab= "Phylogenetic hypotheses", ylab = "Manhattan Stratigraphic Measure*", col=c("#E69F00", "#56B4E9", "#009E73", "#F0E442"), names= c("No constraints", "Crocodyliform constraint", "Neosuchian constraint", "Tethysuchian constraint"))
dev.off()

## GER* boxplots:
boxplot(StratPhyloCongruenceResults_EIWk15_C0$input.tree.results[, "GER*"], StratPhyloCongruenceResults_EIWk15_C1$input.tree.results[, "GER*"], StratPhyloCongruenceResults_EIWk15_C3$input.tree.results[, "GER*"], StratPhyloCongruenceResults_EIWk15_C4$input.tree.results[, "GER*"], xlab= "Phylogenetic hypotheses", ylab = "Gap Excess Ratio (GER*)", col=c("#E69F00", "#56B4E9", "#009E73", "#F0E442"), names= c("No constraints", "Crocodyliform constraint", "Neosuchian constraint", "Tethysuchian constraint"))

pdf("output/StratCongruenceResults_EIWk15/Boxplots/EIWk15_GER*_Results_C0vC1vC3vC4.pdf", width=17, height=9) 
boxplot(StratPhyloCongruenceResults_EIWk15_C0$input.tree.results[, "GER*"], StratPhyloCongruenceResults_EIWk15_C1$input.tree.results[, "GER*"], StratPhyloCongruenceResults_EIWk15_C3$input.tree.results[, "GER*"], StratPhyloCongruenceResults_EIWk15_C4$input.tree.results[, "GER*"], xlab= "Phylogenetic hypotheses", ylab = "Gap Excess Ratio (GER*)", col=c("#E69F00", "#56B4E9", "#009E73", "#F0E442"), names= c("No constraints", "Crocodyliform constraint", "Neosuchian constraint", "Tethysuchian constraint"))
dev.off()

## GERt boxplots:
boxplot(StratPhyloCongruenceResults_EIWk15_C0$input.tree.results[, "GERt"], StratPhyloCongruenceResults_EIWk15_C1$input.tree.results[, "GERt"], StratPhyloCongruenceResults_EIWk15_C3$input.tree.results[, "GERt"], StratPhyloCongruenceResults_EIWk15_C4$input.tree.results[, "GERt"], xlab= "Phylogenetic hypotheses", ylab = "Gap Excess Ratio (GERt)", col=c("#E69F00", "#56B4E9", "#009E73", "#F0E442"), names= c("No constraints", "Crocodyliform constraint", "Neosuchian constraint", "Tethysuchian constraint"))

pdf("output/StratCongruenceResults_EIWk15/Boxplots/EIWk15_GERt_Results_C0vC1vC3vC4.pdf", width=17, height=9) 
boxplot(StratPhyloCongruenceResults_EIWk15_C0$input.tree.results[, "GERt"], StratPhyloCongruenceResults_EIWk15_C1$input.tree.results[, "GERt"], StratPhyloCongruenceResults_EIWk15_C3$input.tree.results[, "GERt"], StratPhyloCongruenceResults_EIWk15_C4$input.tree.results[, "GERt"], xlab= "Phylogenetic hypotheses", ylab = "Gap Excess Ratio (GERt)", col=c("#E69F00", "#56B4E9", "#009E73", "#F0E442"), names= c("No constraints", "Crocodyliform constraint", "Neosuchian constraint", "Tethysuchian constraint"))
dev.off()

## MIG boxplots:
boxplot(StratPhyloCongruenceResults_EIWk15_C0$input.tree.results[, "MIG"], StratPhyloCongruenceResults_EIWk15_C1$input.tree.results[, "MIG"], StratPhyloCongruenceResults_EIWk15_C3$input.tree.results[, "MIG"], StratPhyloCongruenceResults_EIWk15_C4$input.tree.results[, "MIG"], xlab= "Phylogenetic hypotheses", ylab = "Minimum Implied Gap", col=c("#E69F00", "#56B4E9", "#009E73", "#F0E442"), names= c("No constraints", "Crocodyliform constraint", "Neosuchian constraint", "Tethysuchian constraint"))

pdf("output/StratCongruenceResults_EIWk15/Boxplots/EIWk15_MIG_Results_C0vC1vC3vC4.pdf", width=17, height=9) 
boxplot(StratPhyloCongruenceResults_EIWk15_C0$input.tree.results[, "MIG"], StratPhyloCongruenceResults_EIWk15_C1$input.tree.results[, "MIG"], StratPhyloCongruenceResults_EIWk15_C3$input.tree.results[, "MIG"], StratPhyloCongruenceResults_EIWk15_C4$input.tree.results[, "MIG"], xlab= "Phylogenetic hypotheses", ylab = "Minimum Implied Gap", col=c("#E69F00", "#56B4E9", "#009E73", "#F0E442"), names= c("No constraints", "Crocodyliform constraint", "Neosuchian constraint", "Tethysuchian constraint"))
dev.off()

 beep(1)


# ======================== #
# 0.1 Parallelisation ends #
# ======================== #

# stop cluster
 stopImplicitCluster()


# ========= #
# 0.2 Done! #
# ========= #

### Noise when script has finished (useful if you are running this in the background):
 library(beepr)

 beep(8)