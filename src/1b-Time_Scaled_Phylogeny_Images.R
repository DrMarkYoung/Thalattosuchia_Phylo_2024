#!/usr/bin/env Rscript

# ======================================================================================= #
#                                                                                         #
#  Time-scaled phylogeny script - version 2.0 (October, 2024).                            #
#                                                                                         #
#                                                                                         #
#  Purpose:                                                                               #
#  This script creates time-calibrated strict consensus topologies using the Strap and    #
#  Paleotree R packages.                                                                  #
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
 EqW_C0_Strict <- ReadTntTree("dat/Trees/xiw_equal_C0_RC_Nelsen.tre") %>% ladderize
 ape::write.tree(EqW_C0_Strict, file='dat/Trees/xiw_equal_C0_RC_Nelsen.txt') 
 plotTree(EqW_C0_Strict,ftype="i",fsize=0.3,lwd=1)

 EqW_C1_Strict <- ReadTntTree("dat/Trees/xiw_equal_C1_RC_Nelsen.tre") %>% ladderize
 ape::write.tree(EqW_C1_Strict, file='dat/Trees/xiw_equal_C1_RC_Nelsen.txt') 
 plotTree(EqW_C1_Strict,ftype="i",fsize=0.3,lwd=1)

 EqW_C2_Strict <- ReadTntTree("dat/Trees/xiw_equal_C2_RC_Nelsen.tre") %>% ladderize
 ape::write.tree(EqW_C2_Strict, file='dat/Trees/xiw_equal_C2_RC_Nelsen.txt')
 plotTree(EqW_C2_Strict,ftype="i",fsize=0.3,lwd=1)

 EqW_C3_Strict <- ReadTntTree("dat/Trees/xiw_equal_C3_RC_Nelsen.tre") %>% ladderize
 ape::write.tree(EqW_C3_Strict, file='dat/Trees/xiw_equal_C3_RC_Nelsen.txt')
 plotTree(EqW_C3_Strict,ftype="i",fsize=0.3,lwd=1)

 EqW_C4_Strict <- ReadTntTree("dat/Trees/xiw_equal_C4_RC_Nelsen.tre") %>% ladderize
 ape::write.tree(EqW_C4_Strict, file='dat/Trees/xiw_equal_C4_RC_Nelsen.txt')
 plotTree(EqW_C4_Strict,ftype="i",fsize=0.3,lwd=1)

 EIWk15_C0_Strict <- ReadTntTree("dat/Trees/xiw_k15_C0_RC_Nelsen.tre") %>% ladderize
 ape::write.tree(EIWk15_C0_Strict, file='dat/Trees/xiw_k15_C0_RC_Nelsen.txt') 
 plotTree(EIWk15_C0_Strict,ftype="i",fsize=0.3,lwd=1)

 EIWk15_C1_Strict <- ReadTntTree("dat/Trees/xiw_k15_C1_RC_Nelsen.tre") %>% ladderize
 ape::write.tree(EIWk15_C1_Strict, file='dat/Trees/xiw_k15_C1_RC_Nelsen.txt') 
 plotTree(EIWk15_C1_Strict,ftype="i",fsize=0.3,lwd=1)

 EIWk15_C2_Strict <- ReadTntTree("dat/Trees/xiw_k15_C2_RC_Nelsen.tre") %>% ladderize
 ape::write.tree(EIWk15_C2_Strict, file='dat/Trees/xiw_k15_C2_RC_Nelsen.txt')
 plotTree(EIWk15_C2_Strict,ftype="i",fsize=0.3,lwd=1)

 EIWk15_C3_Strict <- ReadTntTree("dat/Trees/xiw_k15_C3_RC_Nelsen.tre") %>% ladderize
 ape::write.tree(EIWk15_C3_Strict, file='dat/Trees/xiw_k15_C3_RC_Nelsen.txt')
 plotTree(EIWk15_C3_Strict,ftype="i",fsize=0.3,lwd=1)

 EIWk15_C4_Strict <- ReadTntTree("dat/Trees/xiw_k15_C4_RC_Nelsen.tre") %>% ladderize
 ape::write.tree(EIWk15_C4_Strict, file='dat/Trees/xiw_k15_C4_RC_Nelsen.txt')
 plotTree(EIWk15_C4_Strict,ftype="i",fsize=0.3,lwd=1)


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


### For loop to check if the topologies are rooted and binary:
 for (i in 1:length(all_trees)) {
   tree <- all_trees[[i]]
   names <- names(all_trees[i])
   print(paste(names, "is rooted:", is.rooted(tree)))
   print(paste(names, "is binary:", is.binary(tree)))
}


# -------------------------------------------------------------------- #
#   NOTE: herein the script uses the common TNT notation I created.    #  
#   C0= no constraints. C1= constraint 1, etc.                         #
#   EqW=equal weights analysis. EIWk15=extended implied weights k=15   #
# -------------------------------------------------------------------- #

 beep(1)


# =================== #
# 2. Input data files #
# =================== #

### required libraries:
 library(beepr)
 library(geiger)


### Loads in Age file:
 Ages <- read.table("dat/Ages/CrocAges.txt", header=T)


### For loop to name check (using geiger) that the pruned phylogeny and age files match:
 for (i in 1:length(all_trees)) {
   tree <- all_trees[[i]]
   names <- names(all_trees[i])
   print(paste(names, "is", name.check(tree, Ages)))
}

 beep(1)


# ================================================================================= #
# 3. Creates time-dated phylogenies using Strap, 'equal method', and generates PDFs #
# ================================================================================= #

### Required libraries:
 library(beepr)
 library(strap)


### Creates time-dated phylogeny using the equal method:
 EqW_C0_Strict_EqualMethod <- DatePhylo(EqW_C0_Strict, Ages, method="equal",rlen=1)
 EqW_C1_Strict_EqualMethod <- DatePhylo(EqW_C1_Strict, Ages, method="equal",rlen=1)
 EqW_C2_Strict_EqualMethod <- DatePhylo(EqW_C2_Strict, Ages, method="equal",rlen=1)
 EqW_C3_Strict_EqualMethod <- DatePhylo(EqW_C3_Strict, Ages, method="equal",rlen=1)
 EqW_C4_Strict_EqualMethod <- DatePhylo(EqW_C4_Strict, Ages, method="equal",rlen=1)
 EIWk15_C0_Strict_EqualMethod <- DatePhylo(EIWk15_C0_Strict, Ages, method="equal",rlen=1)
 EIWk15_C1_Strict_EqualMethod <- DatePhylo(EIWk15_C1_Strict, Ages, method="equal",rlen=1)
 EIWk15_C2_Strict_EqualMethod <- DatePhylo(EIWk15_C2_Strict, Ages, method="equal",rlen=1)
 EIWk15_C3_Strict_EqualMethod <- DatePhylo(EIWk15_C3_Strict, Ages, method="equal",rlen=1)
 EIWk15_C4_Strict_EqualMethod <- DatePhylo(EIWk15_C4_Strict, Ages, method="equal",rlen=1)


### Create a vector with all the different time-scaled topologies:
all_timescaled_EM_trees <- vector("list", 10)
 names(all_timescaled_EM_trees) <- c("EqW_C0_Strict_EqualMethod", "EqW_C1_Strict_EqualMethod", "EqW_C2_Strict_EqualMethod", "EqW_C3_Strict_EqualMethod", "EqW_C4_Strict_EqualMethod", "EIWk15_C0_Strict_EqualMethod","EIWk15_C1_Strict_EqualMethod","EIWk15_C2_Strict_EqualMethod","EIWk15_C3_Strict_EqualMethod", "EIWk15_C4_Strict_EqualMethod")

# Assign the phylogenies to the vector:
 all_timescaled_EM_trees[[1]] <- EqW_C0_Strict_EqualMethod
 all_timescaled_EM_trees[[2]] <- EqW_C1_Strict_EqualMethod
 all_timescaled_EM_trees[[3]] <- EqW_C2_Strict_EqualMethod
 all_timescaled_EM_trees[[4]] <- EqW_C3_Strict_EqualMethod
 all_timescaled_EM_trees[[5]] <- EqW_C4_Strict_EqualMethod
 all_timescaled_EM_trees[[6]] <- EIWk15_C0_Strict_EqualMethod
 all_timescaled_EM_trees[[7]] <- EIWk15_C1_Strict_EqualMethod
 all_timescaled_EM_trees[[8]] <- EIWk15_C2_Strict_EqualMethod
 all_timescaled_EM_trees[[9]] <- EIWk15_C3_Strict_EqualMethod
 all_timescaled_EM_trees[[10]] <- EIWk15_C4_Strict_EqualMethod


### For loop to create equal method-dated strict consensus topology PDFs, and plots against the ICS stratigraphic chart (in different orientations and time bins):

 for (i in 1:length(all_timescaled_EM_trees)) {
   pdf(paste("output/PhyloFigs_Geoscaled/Time_Tree_EM - " ,names(all_timescaled_EM_trees[i]),".pdf", sep=""), width=20, height=20)
   plotTree(all_timescaled_EM_trees[[i]], ftype="i", fsize=0.3, lwd=1)
   dev.off()


# Geoscaled, Period units and oriented right:
   pdf(paste("output/PhyloFigs_Geoscaled/Geoscaled_Tree_EM_ICS_Period_Right - " ,names(all_timescaled_EM_trees[i]),".pdf", sep=""), width=17, height=9)
	geoscalePhylo(all_timescaled_EM_trees[[i]], Ages, cex.age=0.8, cex.ts=0.7, cex.tip=0.3, quat.rm=TRUE, units = c("Period"), boxes="Period", direction="rightwards") 
   dev.off()


# Geoscaled, Era, Period, Epoch and Age units and oriented right:
 pdf(paste("output/PhyloFigs_Geoscaled/Geoscaled_Tree_EM_ICS_EraAge_Right - " ,names(all_timescaled_EM_trees[i]),".pdf", sep=""), width=17, height=9)
 	geoscalePhylo(all_timescaled_EM_trees[[i]], Ages, cex.age=0.8, cex.ts=0.7, cex.tip=0.2, quat.rm=TRUE, units = c("Era","Period","Epoch","Age"), boxes="Age", direction="rightwards") 
 dev.off()


# Geoscaled, Period units and oriented upwards:
 pdf(paste("output/PhyloFigs_Geoscaled/Geoscaled_Tree_EM_ICS_Period_Up - " ,names(all_timescaled_EM_trees[i]),".pdf", sep=""), width=17, height=9)
	geoscalePhylo(all_timescaled_EM_trees[[i]], Ages, cex.age=0.8, cex.ts=0.7, cex.tip=0.2, quat.rm=TRUE, units = c("Period"), boxes="Age", direction="upwards") 
 dev.off()


# Geoscaled, Era, Period, Epoch and Age units and oriented upward:
 pdf(paste("output/PhyloFigs_Geoscaled/Geoscaled_Tree_EM_ICS_EraAge_Up - " ,names(all_timescaled_EM_trees[i]),".pdf", sep=""), width=17, height=9)
	geoscalePhylo(all_timescaled_EM_trees[[i]], Ages, cex.age=0.8, cex.ts=0.7, cex.tip=0.2, quat.rm=TRUE, units = c("Era","Period","Epoch","Age"), boxes="Age", direction="upwards") 
 dev.off()

 }

 beep(1)


# =================================================================================================== #
# 4. Creates time-dated phylogenies using Paleotree, '1Ma min branch length method', & generates PDFs #
# =================================================================================================== #

### Required libraries:
 library(beepr)
 library(paleotree)

### Creates time-dated phylogeny using the minimum branch length (= 1 Ma) method:
 EqW_C0_Strict_mbl1 <- timePaleoPhy(EqW_C0_Strict, Ages, "mbl", 1)
 EqW_C1_Strict_mbl1 <- timePaleoPhy(EqW_C1_Strict, Ages, "mbl", 1)
 EqW_C2_Strict_mbl1 <- timePaleoPhy(EqW_C2_Strict, Ages, "mbl", 1)
 EqW_C3_Strict_mbl1 <- timePaleoPhy(EqW_C3_Strict, Ages, "mbl", 1)
 EqW_C4_Strict_mbl1 <- timePaleoPhy(EqW_C4_Strict, Ages, "mbl", 1)
 EIWk15_C0_Strict_mbl1 <- timePaleoPhy(EIWk15_C0_Strict, Ages, "mbl", 1)
 EIWk15_C1_Strict_mbl1 <- timePaleoPhy(EIWk15_C1_Strict, Ages, "mbl", 1)
 EIWk15_C2_Strict_mbl1 <- timePaleoPhy(EIWk15_C2_Strict, Ages, "mbl", 1)
 EIWk15_C3_Strict_mbl1 <- timePaleoPhy(EIWk15_C3_Strict, Ages, "mbl", 1)
 EIWk15_C4_Strict_mbl1 <- timePaleoPhy(EIWk15_C4_Strict, Ages, "mbl", 1)


### Create a vector with all the different time-scaled topologies:
all_timescaled_mbl1_trees <- vector("list", 10)
 names(all_timescaled_mbl1_trees) <- c("EqW_C0_Strict_mbl1", "EqW_C1_Strict_mbl1", "EqW_C2_Strict_mbl1", "EqW_C3_Strict_mbl1", "EqW_C4_Strict_mbl1", "EIWk15_C0_Strict_mbl1","EIWk15_C1_Strict_mbl1","EIWk15_C2_Strict_mbl1","EIWk15_C3_Strict_mbl1", "EIWk15_C4_Strict_mbl1")

# Assign the phylogenies to the vector:
 all_timescaled_mbl1_trees[[1]] <- EqW_C0_Strict_mbl1
 all_timescaled_mbl1_trees[[2]] <- EqW_C1_Strict_mbl1
 all_timescaled_mbl1_trees[[3]] <- EqW_C2_Strict_mbl1
 all_timescaled_mbl1_trees[[4]] <- EqW_C3_Strict_mbl1
 all_timescaled_mbl1_trees[[5]] <- EqW_C4_Strict_mbl1
 all_timescaled_mbl1_trees[[6]] <- EIWk15_C0_Strict_mbl1
 all_timescaled_mbl1_trees[[7]] <- EIWk15_C1_Strict_mbl1
 all_timescaled_mbl1_trees[[8]] <- EIWk15_C2_Strict_mbl1
 all_timescaled_mbl1_trees[[9]] <- EIWk15_C3_Strict_mbl1
 all_timescaled_mbl1_trees[[10]] <- EIWk15_C4_Strict_mbl1


### For loop to create mbl1-dated strict consensus topology PDFs, and plots against the ICS stratigraphic chart (in different orientations and time bins):

 for (i in 1:length(all_timescaled_mbl1_trees)) {
   pdf(paste("output/PhyloFigs_Geoscaled/Time_Tree_mbl1 - " ,names(all_timescaled_mbl1_trees[i]),".pdf", sep=""), width=20, height=20)
   plotTree(all_timescaled_mbl1_trees[[i]], ftype="i", fsize=0.3, lwd=1)
   dev.off()


# Geoscaled, Period units and oriented right:
   pdf(paste("output/PhyloFigs_Geoscaled/Geoscaled_Tree_mbl1_ICS_Period_Right - " ,names(all_timescaled_mbl1_trees[i]),".pdf", sep=""), width=17, height=9)
	geoscalePhylo(all_timescaled_mbl1_trees[[i]], Ages, cex.age=0.8, cex.ts=0.7, cex.tip=0.3, quat.rm=TRUE, units = c("Period"), boxes="Period", direction="rightwards") 
   dev.off()


# Geoscaled, Era, Period, Epoch and Age units and oriented right:
 pdf(paste("output/PhyloFigs_Geoscaled/Geoscaled_Tree_mbl1_ICS_EraAge_Right - " ,names(all_timescaled_mbl1_trees[i]),".pdf", sep=""), width=17, height=9)
 	geoscalePhylo(all_timescaled_mbl1_trees[[i]], Ages, cex.age=0.8, cex.ts=0.7, cex.tip=0.3, quat.rm=TRUE, units = c("Era","Period","Epoch","Age"), boxes="Age", direction="rightwards") 
 dev.off()


# Geoscaled, Period units and oriented upwards:
 pdf(paste("output/PhyloFigs_Geoscaled/Geoscaled_Tree_mbl1_ICS_Period_Up - " ,names(all_timescaled_mbl1_trees[i]),".pdf", sep=""), width=17, height=9)
	geoscalePhylo(all_timescaled_mbl1_trees[[i]], Ages, cex.age=0.8, cex.ts=0.7, cex.tip=0.2, quat.rm=TRUE, units = c("Period"), boxes="Age", direction="upwards") 
 dev.off()


# Geoscaled, Era, Period, Epoch and Age units and oriented upward:
 pdf(paste("output/PhyloFigs_Geoscaled/Geoscaled_Tree_mbl1_ICS_EraAge_Up - " ,names(all_timescaled_mbl1_trees[i]),".pdf", sep=""), width=17, height=9)
	geoscalePhylo(all_timescaled_mbl1_trees[[i]], Ages, cex.age=0.8, cex.ts=0.7, cex.tip=0.2, quat.rm=TRUE, units = c("Era","Period","Epoch","Age"), boxes="Age", direction="upwards") 
 dev.off()

 }

 beep(1)


# =========================================================================== #
# 5. Creates time-dated phylogenies for Thalattosuchia-only, & generates PDFs #
# =========================================================================== #


# ------------------------------------- #
# 5.1 Pruning dataset and name checking #
# ------------------------------------- #

### Required libraries:
 library(beepr)
 library(strap)
 library(geiger)

### set species to be included in phylogeny image:
thalattosuchian.species<-c("Moroccan_teleosauroid","Turnersuchus_hingleyae","Plagiophthalmosuchus_gracilirostris","Mystriosaurus_laurillardi","Chinese_teleosauroid","Indosinosuchus_potamosiamensis","Platysuchus_multiscrobiculatus","Teleosaurus_cadomensis","Seldsienean_megistorhynchus","Mycterosuchus_nasutus","Bathysuchus_megarhinus","Aeolodon_priscus","Sericodon_jugleri","Macrospondylus_bollensis","Charitomenosuchus_leedsi","Deslongchampsina_larteti","Andrianavoay_baroni","Proexochokefalos_heberti","Proexochokefalos_cf_bouchardi","Neosteneosaurus_edwardsi","Yvridiosuchus_boutilieri","Lemmysuchus_obtusidens","Machimosaurus_buffetauti","Machimosaurus_mosae","Machimosaurus_hugii","Machimosaurus_rex","Pelagosaurus_typus","Pelagosaurus_tomarensis","Opisuchus_meieri","Teleidosaurus_calvadosii","Peipehsuchus_teleorhinus","Metriorhynchoid_indet","Magyarosuchus_fitosi","Eoneustes_bathonicus","Eoneustes_gaudryi","Chile_metriorhynchoid","Zoneait_nargorum","Cricosaurus_albersdoerferi","Cricosaurus_bambergensis","Cricosaurus_elegans","Cricosaurus_suevicus","Cricosaurus_schroederi","Cricosaurus_araucanensis","Cricosaurus_vignaudi","Cricosaurus_lithographicus","Cricosaurus_rauhuti","Cricosaurus_puelchorum","Dakosaurus_lissocephalus","Rhacheosaurus_gracilis","Cretaceous_rhacheosaurin","Swiss_rhacheosaurin","Cricosaurus_saltillensis","Metriorhynchus_palpebrosus","USNM_419640_Cuba","MNHNH_P3009_Cuba","Maledictosuchus_nuyivijanan","Maledictosuchus_riclaensis","Gracilineustes_acutus","Gracilineustes_leedsi","Metriorhynchus_brevirostris","Metriorhynchus_cf_brevirostris","Thalattosuchus_superciliosus","Metriorhynchus_casamiquelai","Metriorhynchus_westermanni","Neptunidraco_ammoniticus","Metriorhynchus_brachyrhynchus","Tyrannoneustes_lythrodectikos","English_rostrum","Swiss_rostrum","Chouquet_cf_hastifer","Druegendorf_merged","Mr_Passmores_specimen","cf_Torvoneustes","Torvoneustes_coryphaeus","Torvoneustes_carpenteri","Torvoneustes_mexicanus","Torvoneustes_sp_England","Torvoneustes_sp_Czech","Purranisaurus_potens","Ieldraan_melkshamensis","Geosaurus_giganteus","Geosaurus_grandis","Geosaurus_lapparenti","Suchodus_durobrivensis","Plesiosuchus_manselii","Plesiosuchina_indet_France","Plesiosuchina_indet_Sicily","Dakosaurus_andiniensis","Dakosaurus_maximus","Mr_Leeds_dakosaur","Vaches_Noires_dakosaur", "Enaliosuchus_macrospondylus","Neustosaurus_gigondarum","Sinemurian_taxon")


### Loads in Age file:
 ThalattoAges <- read.table("dat/Ages/ThalattosuchianAges.txt", header=T)


### Prunes the dataset, NOTE: use drop.tip to only keep selected OTUs
 EqW_C0_Strict_Pruned <- drop.tip(EqW_C0_Strict,setdiff(EqW_C0_Strict $tip.label,thalattosuchian.species))
 EqW_C1_Strict_Pruned <- drop.tip(EqW_C1_Strict,setdiff(EqW_C1_Strict $tip.label,thalattosuchian.species))
 EqW_C2_Strict_Pruned <- drop.tip(EqW_C2_Strict,setdiff(EqW_C2_Strict $tip.label,thalattosuchian.species))
 EqW_C3_Strict_Pruned <- drop.tip(EqW_C3_Strict,setdiff(EqW_C3_Strict $tip.label,thalattosuchian.species))
 EqW_C4_Strict_Pruned <- drop.tip(EqW_C4_Strict,setdiff(EqW_C4_Strict $tip.label,thalattosuchian.species))
 EIWk15_C0_Strict_Pruned <- drop.tip(EIWk15_C0_Strict,setdiff(EIWk15_C0_Strict $tip.label,thalattosuchian.species))
 EIWk15_C1_Strict_Pruned <- drop.tip(EIWk15_C1_Strict,setdiff(EIWk15_C1_Strict $tip.label,thalattosuchian.species))
 EIWk15_C2_Strict_Pruned <- drop.tip(EIWk15_C2_Strict,setdiff(EIWk15_C2_Strict $tip.label,thalattosuchian.species))
 EIWk15_C3_Strict_Pruned <- drop.tip(EIWk15_C3_Strict,setdiff(EIWk15_C3_Strict $tip.label,thalattosuchian.species))
 EIWk15_C4_Strict_Pruned <- drop.tip(EIWk15_C4_Strict,setdiff(EIWk15_C4_Strict $tip.label,thalattosuchian.species))


### Create a vector with all the different topologies:
all_thalattosuchian_trees <- vector("list", 10)
 names(all_thalattosuchian_trees) <- c("EqW_C0_Strict_Pruned", "EqW_C1_Strict_Pruned", "EqW_C2_Strict_Pruned", "EqW_C3_Strict_Pruned", "EqW_C4_Strict_Pruned", "EIWk15_C0_Strict_Pruned","EIWk15_C1_Strict_Pruned","EIWk15_C2_Strict_Pruned","EIWk15_C3_Strict_Pruned", "EIWk15_C4_Strict_Pruned")

# Assign the phylogenies to the vector:
 all_thalattosuchian_trees[[1]] <- EqW_C0_Strict_Pruned
 all_thalattosuchian_trees[[2]] <- EqW_C1_Strict_Pruned
 all_thalattosuchian_trees[[3]] <- EqW_C2_Strict_Pruned
 all_thalattosuchian_trees[[4]] <- EqW_C3_Strict_Pruned
 all_thalattosuchian_trees[[5]] <- EqW_C4_Strict_Pruned
 all_thalattosuchian_trees[[6]] <- EIWk15_C0_Strict_Pruned
 all_thalattosuchian_trees[[7]] <- EIWk15_C1_Strict_Pruned
 all_thalattosuchian_trees[[8]] <- EIWk15_C2_Strict_Pruned
 all_thalattosuchian_trees[[9]] <- EIWk15_C3_Strict_Pruned
 all_thalattosuchian_trees[[10]] <- EIWk15_C4_Strict_Pruned


### For loop to name check (using geiger) that the pruned phylogeny and age files match:
 for (i in 1:length(all_thalattosuchian_trees)) {
   thalatto_trees <- all_thalattosuchian_trees[[i]]
   thalatto_names <- names(all_thalattosuchian_trees[i])
   print(paste(thalatto_names, "is", name.check(thalatto_trees, ThalattoAges)))
}


# ---------------- #
# 5.2 Time-scaling #
# ---------------- #

### Creates time-dated phylogeny using the equal method for only thalattosuchians:
 EqW_C0_Strict_EqualMethod_Pruned <- DatePhylo(EqW_C0_Strict_Pruned, ThalattoAges, method="equal",rlen=1)
 EqW_C1_Strict_EqualMethod_Pruned <- DatePhylo(EqW_C1_Strict_Pruned, ThalattoAges, method="equal",rlen=1)
 EqW_C2_Strict_EqualMethod_Pruned <- DatePhylo(EqW_C2_Strict_Pruned, ThalattoAges, method="equal",rlen=1)
 EqW_C3_Strict_EqualMethod_Pruned <- DatePhylo(EqW_C3_Strict_Pruned, ThalattoAges, method="equal",rlen=1)
 EqW_C4_Strict_EqualMethod_Pruned <- DatePhylo(EqW_C4_Strict_Pruned, ThalattoAges, method="equal",rlen=1)
 EIWk15_C0_Strict_EqualMethod_Pruned <- DatePhylo(EIWk15_C0_Strict_Pruned, ThalattoAges, method="equal",rlen=1)
 EIWk15_C1_Strict_EqualMethod_Pruned <- DatePhylo(EIWk15_C1_Strict_Pruned, ThalattoAges, method="equal",rlen=1)
 EIWk15_C2_Strict_EqualMethod_Pruned <- DatePhylo(EIWk15_C2_Strict_Pruned, ThalattoAges, method="equal",rlen=1)
 EIWk15_C3_Strict_EqualMethod_Pruned <- DatePhylo(EIWk15_C3_Strict_Pruned, ThalattoAges, method="equal",rlen=1)
 EIWk15_C4_Strict_EqualMethod_Pruned <- DatePhylo(EIWk15_C4_Strict_Pruned, ThalattoAges, method="equal",rlen=1)


### Create a vector with all the different time-scaled topologies:
all_thalattosuchian_timescaled_EM_trees <- vector("list", 10)
 names(all_thalattosuchian_timescaled_EM_trees) <- c("EqW_C0_Strict_EqualMethod_Pruned", "EqW_C1_Strict_EqualMethod_Pruned", "EqW_C2_Strict_EqualMethod_Pruned", "EqW_C3_Strict_EqualMethod_Pruned", "EqW_C4_Strict_EqualMethod_Pruned", "EIWk15_C0_Strict_EqualMethod_Pruned","EIWk15_C1_Strict_EqualMethod_Pruned","EIWk15_C2_Strict_EqualMethod_Pruned","EIWk15_C3_Strict_EqualMethod_Pruned", "EIWk15_C4_Strict_EqualMethod_Pruned")

# Assign the phylogenies to the vector:
 all_thalattosuchian_timescaled_EM_trees[[1]] <- EqW_C0_Strict_EqualMethod_Pruned
 all_thalattosuchian_timescaled_EM_trees[[2]] <- EqW_C1_Strict_EqualMethod_Pruned
 all_thalattosuchian_timescaled_EM_trees[[3]] <- EqW_C2_Strict_EqualMethod_Pruned
 all_thalattosuchian_timescaled_EM_trees[[4]] <- EqW_C3_Strict_EqualMethod_Pruned
 all_thalattosuchian_timescaled_EM_trees[[5]] <- EqW_C4_Strict_EqualMethod_Pruned
 all_thalattosuchian_timescaled_EM_trees[[6]] <- EIWk15_C0_Strict_EqualMethod_Pruned
 all_thalattosuchian_timescaled_EM_trees[[7]] <- EIWk15_C1_Strict_EqualMethod_Pruned
 all_thalattosuchian_timescaled_EM_trees[[8]] <- EIWk15_C2_Strict_EqualMethod_Pruned
 all_thalattosuchian_timescaled_EM_trees[[9]] <- EIWk15_C3_Strict_EqualMethod_Pruned
 all_thalattosuchian_timescaled_EM_trees[[10]] <- EIWk15_C4_Strict_EqualMethod_Pruned



### For loop to create equal method-dated strict consensus topology PDFs, and plots against the ICS stratigraphic chart (in different orientations and time bins):

 for (i in 1:length(all_thalattosuchian_timescaled_EM_trees)) {
   pdf(paste("output/PhyloFigs_Geoscaled/Pruned/Time_Tree_EM - " ,names(all_thalattosuchian_timescaled_EM_trees[i]),".pdf", sep=""), width=20, height=20)
   plotTree(all_thalattosuchian_timescaled_EM_trees[[i]], ftype="i", fsize=0.3, lwd=1)
   dev.off()


# Geoscaled, Period units and oriented right:
   pdf(paste("output/PhyloFigs_Geoscaled/Pruned/Geoscaled_Tree_EM_ICS_Period_Right - " ,names(all_thalattosuchian_timescaled_EM_trees[i]),".pdf", sep=""), width=17, height=9)
	geoscalePhylo(all_thalattosuchian_timescaled_EM_trees[[i]], Ages, cex.age=0.8, cex.ts=0.7, cex.tip=0.5, quat.rm=TRUE, units = c("Period"), boxes="Period", direction="rightwards") 
   dev.off()


# Geoscaled, Era, Period, Epoch and Age units and oriented right:
 pdf(paste("output/PhyloFigs_Geoscaled/Pruned/Geoscaled_Tree_EM_ICS_EraAge_Right - " ,names(all_thalattosuchian_timescaled_EM_trees[i]),".pdf", sep=""), width=17, height=9)
 	geoscalePhylo(all_thalattosuchian_timescaled_EM_trees[[i]], Ages, cex.age=0.8, cex.ts=0.7, cex.tip=0.5, quat.rm=TRUE, units = c("Era","Period","Epoch","Age"), boxes="Age", direction="rightwards") 
 dev.off()


# Geoscaled, Period units and oriented upwards:
 pdf(paste("output/PhyloFigs_Geoscaled/Pruned/Geoscaled_Tree_EM_ICS_Period_Up - " ,names(all_thalattosuchian_timescaled_EM_trees[i]),".pdf", sep=""), width=17, height=9)
	geoscalePhylo(all_thalattosuchian_timescaled_EM_trees[[i]], Ages, cex.age=0.8, cex.ts=0.7, cex.tip=0.5, quat.rm=TRUE, units = c("Period"), boxes="Age", direction="upwards") 
 dev.off()


# Geoscaled, Era, Period, Epoch and Age units and oriented upward:
 pdf(paste("output/PhyloFigs_Geoscaled/Pruned/Geoscaled_Tree_EM_ICS_EraAge_Up - " ,names(all_thalattosuchian_timescaled_EM_trees[i]),".pdf", sep=""), width=17, height=9)
	geoscalePhylo(all_thalattosuchian_timescaled_EM_trees[[i]], Ages, cex.age=0.8, cex.ts=0.7, cex.tip=0.5, quat.rm=TRUE, units = c("Era","Period","Epoch","Age"), boxes="Age", direction="upwards") 
 dev.off()

 }

 beep(1)


# ======== #
# 0. Done! #
# ======== #

### Noise when script has finished (useful if you are running this in the background):
 library(beepr)

 beep(8)