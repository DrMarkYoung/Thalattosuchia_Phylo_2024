#!/usr/bin/env Rscript

# ======================================================================================= #
#                                                                                         #
#  Basic tree figure script - version 2.0 (October, 2024).                                #
#                                                                                         #
#                                                                                         #
#  Purpose:                                                                               #
#  This script creates 'basic' phylogeny images for the strict consensus topology.        #
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


# ============================== #
# 1. Input and format tree files #
# ============================== #

### Required libraries:
 library(ape)
 library(beepr)
 library(geoscale)
 library(TreeTools)
 library(magrittr)
 library(phytools)


### Loads in the strict consensus files from the TNT phylogenetic analyses, ladderises, and then converts them into nexus format using Ape. Plots them using Phytools. The %>% requires magrittr:

 EqW_C0_Strict <- ReadTntTree("dat/Trees/xiw_equal_C0_RC_nelsen.tre") %>% ladderize
 ape::write.tree(EqW_C0_Strict, file='dat/Trees/xiw_equal_C0_RC_Nelsen.txt') 
 plotTree(EqW_C0_Strict,ftype="i",fsize=0.3,lwd=1)

 EqW_C1_Strict <- ReadTntTree("dat/Trees/xiw_equal_C1_RC_nelsen.tre") %>% ladderize
 ape::write.tree(EqW_C1_Strict, file='dat/Trees/xiw_equal_C1_RC_Nelsen.txt') 
 plotTree(EqW_C1_Strict,ftype="i",fsize=0.3,lwd=1)

 EqW_C2_Strict <- ReadTntTree("dat/Trees/xiw_equal_C2_RC_nelsen.tre") %>% ladderize
 ape::write.tree(EqW_C2_Strict, file='dat/Trees/xiw_equal_C2_RC_Nelsen.txt')
 plotTree(EqW_C2_Strict,ftype="i",fsize=0.3,lwd=1)

 EqW_C3_Strict <- ReadTntTree("dat/Trees/xiw_equal_C3_RC_nelsen.tre") %>% ladderize
 ape::write.tree(EqW_C3_Strict, file='dat/Trees/xiw_equal_C3_RC_Nelsen.txt')
 plotTree(EqW_C3_Strict,ftype="i",fsize=0.3,lwd=1)

 EqW_C4_Strict <- ReadTntTree("dat/Trees/xiw_equal_C4_RC_nelsen.tre") %>% ladderize
 ape::write.tree(EqW_C4_Strict, file='dat/Trees/xiw_equal_C4_RC_Nelsen.txt')
 plotTree(EqW_C4_Strict,ftype="i",fsize=0.3,lwd=1)

 EIW15_C0_Strict <- ReadTntTree("dat/Trees/xiw_k15_C0_RC_nelsen.tre") %>% ladderize
 ape::write.tree(EIW15_C0_Strict, file='dat/Trees/xiw_k15_C0_RC_Nelsen.txt') 
 plotTree(EIW15_C0_Strict,ftype="i",fsize=0.3,lwd=1)

 EIW15_C1_Strict <- ReadTntTree("dat/Trees/xiw_k15_C1_RC_nelsen.tre") %>% ladderize
 ape::write.tree(EIW15_C1_Strict, file='dat/Trees/xiw_k15_C1_RC_Nelsen.txt') 
 plotTree(EIW15_C1_Strict,ftype="i",fsize=0.3,lwd=1)

 EIW15_C2_Strict <- ReadTntTree("dat/Trees/xiw_k15_C2_RC_nelsen.tre") %>% ladderize
 ape::write.tree(EIW15_C2_Strict, file='dat/Trees/xiw_k15_C2_RC_Nelsen.txt')
 plotTree(EIW15_C2_Strict,ftype="i",fsize=0.3,lwd=1)

 EIW15_C3_Strict <- ReadTntTree("dat/Trees/xiw_k15_C3_RC_nelsen.tre") %>% ladderize
 ape::write.tree(EIW15_C3_Strict, file='dat/Trees/xiw_k15_C3_RC_Nelsen.txt')
 plotTree(EIW15_C3_Strict,ftype="i",fsize=0.3,lwd=1)

 EIW15_C4_Strict <- ReadTntTree("dat/Trees/xiw_k15_C4_RC_nelsen.tre") %>% ladderize
 ape::write.tree(EIW15_C4_Strict, file='dat/Trees/xiw_k15_C4_RC_Nelsen.txt')
 plotTree(EIW15_C4_Strict,ftype="i",fsize=0.3,lwd=1)


# -------------------------------------------------------------------- #
#   NOTE: herein the script uses the common TNT notation I created.    #  
#   C0= no constraints. C1= constraint 1, etc.                         #
#   EqW=equal weights analysis. EIWk15=extended implied weights k=15   #
# -------------------------------------------------------------------- #

 beep(1)


# ========================================================================= #
# 2. Creates standard and fan topologies using Phytools, and generates PDFs #
# ========================================================================= #

### Required libraries:
 library(beepr)
 library(maps)
 library(phytools)


### Create a vector with all the different topologies:
all_trees <- vector("list", 10)
 names(all_trees) <- c("EqW_C0_Strict", "EqW_C1_Strict", "EqW_C2_Strict", "EqW_C3_Strict", "EqW_C4_Strict", "EIW15_C0_Strict","EIW15_C1_Strict","EIW15_C2_Strict","EIW15_C3_Strict", "EIW15_C4_Strict")

# Assign the phylogenies to the vector:
 all_trees[[1]] <- EqW_C0_Strict
 all_trees[[2]] <- EqW_C1_Strict
 all_trees[[3]] <- EqW_C2_Strict
 all_trees[[4]] <- EqW_C3_Strict
 all_trees[[5]] <- EqW_C4_Strict
 all_trees[[6]] <- EIW15_C0_Strict
 all_trees[[7]] <- EIW15_C1_Strict
 all_trees[[8]] <- EIW15_C2_Strict
 all_trees[[9]] <- EIW15_C3_Strict
 all_trees[[10]] <- EIW15_C4_Strict


### For loop to create PDFs of strict consensus topologies (standard, fan, and half fan):
 for (i in 1:length(all_trees)) {
   pdf(paste("output/PhyloFigs/Strict_Consensus - " ,names(all_trees[i]),".pdf", sep=""), width=13, height=9)
   plotTree(all_trees[[i]], ftype="i", fsize=0.3, lwd=1)
   dev.off()
   
   pdf(paste("output/PhyloFigs/Strict_Fan_Consensus - " ,names(all_trees[i]),".pdf", sep=""), width=20, height=20)
   plotTree(all_trees[[i]], ftype="i",type="fan",fsize=0.3,lwd=1.5)
   dev.off()
   
   pdf(paste("output/PhyloFigs/Strict_Half_Fan_Consensus - " ,names(all_trees[i]),".pdf", sep=""), width=20, height=20)
   plotTree(all_trees[[i]],ftype="i",type="fan",part=0.5,fsize=0.3,lwd=1.5)
   dev.off()
 }

 beep(1)


# ================================================================================ #
# 3. Creates phylogenies using Phytools with 'branch length=1', and generates PDFs #
# ================================================================================ #

### required libraries:
 library(beepr)
 library(maps)
 library(phytools)


# Computes phylogeny with branch length = 1:
 EqW_C0_Strict_BL1 <- compute.brlen(EqW_C0_Strict,1)
 EqW_C1_Strict_BL1 <- compute.brlen(EqW_C1_Strict,1)
 EqW_C2_Strict_BL1 <- compute.brlen(EqW_C2_Strict,1)
 EqW_C3_Strict_BL1 <- compute.brlen(EqW_C3_Strict,1)
 EqW_C4_Strict_BL1 <- compute.brlen(EqW_C4_Strict,1)
 EIW15_C0_Strict_BL1 <- compute.brlen(EIW15_C0_Strict,1)
 EIW15_C1_Strict_BL1 <- compute.brlen(EIW15_C1_Strict,1)
 EIW15_C2_Strict_BL1 <- compute.brlen(EIW15_C2_Strict,1)
 EIW15_C3_Strict_BL1 <- compute.brlen(EIW15_C3_Strict,1)
 EIW15_C4_Strict_BL1 <- compute.brlen(EIW15_C4_Strict,1)


# Create a vector with all the different topologies:
 all_trees_BL1 <- vector("list", 10)
 names(all_trees_BL1) <- c("EqW_C0_Strict_BL1", "EqW_C1_Strict_BL1", "EqW_C2_Strict_BL1", "EqW_C3_Strict_BL1", "EqW_C4_Strict_BL1", "EIW15_C0_Strict_BL1","EIW15_C1_Strict_BL1","EIW15_C2_Strict_BL1","EIW15_C3_Strict_BL1", "EIW15_C4_Strict_BL1")


# Assign the phylogenies to the vector:
 all_trees_BL1[[1]] <- EqW_C0_Strict_BL1
 all_trees_BL1[[2]] <- EqW_C1_Strict_BL1
 all_trees_BL1[[3]] <- EqW_C2_Strict_BL1
 all_trees_BL1[[4]] <- EqW_C3_Strict_BL1
 all_trees_BL1[[5]] <- EqW_C4_Strict_BL1
 all_trees_BL1[[6]] <- EIW15_C0_Strict_BL1
 all_trees_BL1[[7]] <- EIW15_C1_Strict_BL1
 all_trees_BL1[[8]] <- EIW15_C2_Strict_BL1
 all_trees_BL1[[9]] <- EIW15_C3_Strict_BL1
 all_trees_BL1[[10]] <- EIW15_C4_Strict_BL1
 
 
# For loop to create PDFs of strict consensus topologies (standard) for branch length of 1:
 for (i in 1:length(all_trees_BL1)) {
   pdf(paste("output/PhyloFigs/Strict_Consensus_BL1 - " ,names(all_trees_BL1[i]),".pdf", sep=""), width=15, height=15)
   plotTree(all_trees_BL1[[i]], ftype="i", fsize=0.3, lwd=1)
   dev.off()
 }

 beep(1)


# ======== #
# 0. Done! #
# ======== #

### Noise when script has finished (useful if you are running this in the background):
library(beepr)

 beep(8)