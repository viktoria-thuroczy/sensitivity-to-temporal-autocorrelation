rm(list = ls())
#load packages
library(Rcompadre)
library(here)
library(readxl)
library(dplyr)

library(Rage)
library(tidyverse)
library(popbio)

library(rotl)
library(ape)
library(phylotools)


setwd("C:/Users/vikto/OneDrive/Documents/stochastic simulation paper/code/data")

load("COMADRE_v.4.23.3.1.RData")
load("COMPADRE_v.6.23.5.0.RData")

compadre <- as_cdb(compadre)
comadre <- as_cdb(comadre)

###Data Filtering ----

## flag matrices with issues (so able to calculate growth rates)
# filter: irreducibility ergodicity and primitivity
comadre <- cdb_flag(comadre)
compadre <- cdb_flag(compadre)
comadre<- subset(comadre, check_NA_A==FALSE & check_NA_U==FALSE & check_ergodic==TRUE
                 & check_primitive==TRUE & check_irreducible==TRUE)
compadre<- subset(compadre, check_NA_A==FALSE & check_NA_U==FALSE & check_ergodic==TRUE
                  & check_primitive==TRUE & check_irreducible==TRUE)


## subset to unmanipulated, wild populations and have 3 or more life stages
comadre<- subset(comadre, MatrixTreatment=="Unmanipulated" & 
                   MatrixCaptivity == "W" &
                   MatrixDimension >= 3)
compadre<- subset(compadre, MatrixTreatment=="Unmanipulated" &
                    MatrixCaptivity == "W" &
                    MatrixDimension >= 3)


## only individual matrixes
comadre<- subset(comadre, MatrixComposite=="Individual")
compadre<- subset(compadre, MatrixComposite=="Individual")


## remove zero fecundity and clonal species:
I<- which(comadre$MatrixFec=="Yes")
comadre<- comadre[I,]
I<- which(compadre$MatrixFec=="Yes")
compadre<- compadre[I,]


## remove any matrices that didn't OBSERVE any fecundity:
sum_fec<- unlist(lapply(matF(comadre), function(x){sum(sum(x))}))
I<- which(sum_fec>0)
length(I)
# we want to drop just these matrices, but other matrices from the same species/author might be fine
comadre<- comadre[I,]

sum_fec<- unlist(lapply(matF(compadre), function(x){sum(sum(x))}))
I<- which(sum_fec>0)
length(I)
# we want to drop just these matrices, but other matrices from the same species/author might be fine
compadre<- compadre[I,]


## remove species with clonal reproduction
sum_clonal<- unlist(lapply(matC(compadre), function(x){sum(sum(x))}))
I<- which(sum_clonal>0)
length(I)
clonal<- unique(cdb_flatten(compadre)[I, "SpeciesAuthor"])
length(clonal$SpeciesAuthor)
# we want to drop all matrices from those SpeciesAuthor's, not just the ones where clonality was observed:
I<- which(compadre$SpeciesAuthor %in% clonal$SpeciesAuthor)
compadre<- compadre[-I,] # drop these
# There is clonality observed in 213 of the COMPADRE matrices, from 47 species/studies
sum_clonal<- unlist(lapply(matC(comadre), function(x){sum(sum(x))}))
I<- which(sum_clonal>0)
length(I)
clonal<- unique(cdb_flatten(comadre)[I, "SpeciesAuthor"])
# we want to drop all matrices from those SpeciesAuthor's, not just the ones where clonality was observed:
I<- which(comadre$SpeciesAuthor %in% clonal$SpeciesAuthor)
comadre<- comadre[-I,] # drop these


## restrict to only divisible matrices (need to be able to separate A into F and U)
I<- which(comadre$MatrixSplit=="Divided")
comadre<- comadre[I,]
I<- which(compadre$MatrixSplit=="Divided")
compadre<- compadre[I,]


## filter one population per species - population with the most matrixes (longest study period)
population_counts <- comadre %>%
  group_by(SpeciesAccepted, MatrixPopulation) %>%
  summarize(Count = n()) %>%
  ungroup()
# find the matrix population that appears most for each species
max_populations <- population_counts %>%
  group_by(SpeciesAccepted) %>%
  filter(Count == max(Count)) %>%
  select(SpeciesAccepted, MatrixPopulation)
# join with the original dataframe to keep only the desired rows
comadre <- comadre %>%
  inner_join(max_populations, by = c("SpeciesAccepted", "MatrixPopulation"))

population_counts <- compadre %>%
  group_by(SpeciesAccepted, MatrixPopulation) %>%
  summarize(Count = n()) %>%
  ungroup()
# find the matrix population that appears most for each species
max_populations <- population_counts %>%
  group_by(SpeciesAccepted) %>%
  filter(Count == max(Count)) %>%
  select(SpeciesAccepted, MatrixPopulation)
# join with the original dataframe to keep only the desired rows
compadre <- compadre %>%
  inner_join(max_populations, by = c("SpeciesAccepted", "MatrixPopulation"))

## only keep data with same matrix dimensions - population with the most available matrices
population_counts <- comadre %>%
  group_by(SpeciesAccepted, MatrixDimension) %>%
  summarize(Count = n()) %>%
  ungroup()
# find the matrix population that appears most for each species
max_populations <- population_counts %>%
  group_by(SpeciesAccepted) %>%
  filter(Count == max(Count)) %>%
  select(SpeciesAccepted, MatrixDimension)
# join with the original dataframe to keep only the desired rows
comadre <- comadre %>%
  inner_join(max_populations, by = c("SpeciesAccepted", "MatrixDimension"))

#compadre
population_counts <- compadre %>%
  group_by(SpeciesAccepted, MatrixDimension) %>%
  summarize(Count = n()) %>%
  ungroup()
# find the matrix population that appears most for each species
max_populations <- population_counts %>%
  group_by(SpeciesAccepted) %>%
  filter(Count == max(Count)) %>%
  select(SpeciesAccepted, MatrixDimension)
# join with the original dataframe to keep only the desired rows
compadre <- compadre %>%
  inner_join(max_populations, by = c("SpeciesAccepted", "MatrixDimension"))

## if species have same number of matrices in multiple dimensions, keep the ones with higher dimension
comadre <- comadre %>%
  group_by(SpeciesAccepted) %>%
  filter(MatrixDimension == max(MatrixDimension)) %>%
  ungroup()

compadre <- compadre %>%
  group_by(SpeciesAccepted) %>%
  filter(MatrixDimension == max(MatrixDimension)) %>%
  ungroup()

## only keep species which have 5 or more matrices
species_counts <- table(comadre$SpeciesAccepted)
# initialize an empty vector to store species names that appear more than 4 times
species_more_than_4 <- character()
# loop through the names and check if they appear more than 4 times
for (species_name in names(species_counts)) {
  if (species_counts[species_name] > 4) {
    species_more_than_4 <- c(species_more_than_4, species_name)
  }
}
comadre <- comadre[comadre$SpeciesAccepted %in% species_more_than_4, ]

species_counts <- table(compadre$SpeciesAccepted)
# initialize an empty vector to store species names that appear more than 4 times
species_more_than_4 <- character()
# loop through the names and check if they appear more than 4 times
for (species_name in names(species_counts)) {
  if (species_counts[species_name] > 4) {
    species_more_than_4 <- c(species_more_than_4, species_name)
  }
}
compadre <- compadre[compadre$SpeciesAccepted %in% species_more_than_4, ]

## exclude species without good data
comadre <- comadre %>%
  filter(!comadre$SpeciesAccepted == "Scolytus ventralis")

compadre <- compadre %>%
  filter(!SpeciesAuthor == "Primula_veris_4")


save(comadre, file = "comadre.RData")
save(compadre, file = "compadre.RData")

### Calculating LHT ----

# collapsing MPMs into one mean MPM
comadre$id_stage <- cdb_id_stages(comadre, "MatrixClassOrganized")
comadre_collapsed <- cdb_collapse(comadre, columns = c("id_stage"))

compadre$id_stage <- ifelse(
  compadre$SpeciesAccepted == "Trillium grandiflorum",
  cdb_id_stages(compadre, "MatrixClassAuthor"),
  cdb_id_stages(compadre, "MatrixClassOrganized"))

compadre_collapsed <- cdb_collapse(compadre, columns = c("id_stage"))

table(comadre_collapsed$SpeciesAccepted, comadre_collapsed$MatrixComposite)
table(compadre_collapsed$SpeciesAccepted, compadre_collapsed$MatrixComposite)

#save(comadre_collapsed, file = "comadre_collapsed.RData")
#save(compadre_collapsed, file = "compadre_collapsed.RData")


## comadre
# make dataframe to load LHTs into
comadre_LHT <- data.frame(
  SpeciesAccepted = character(),
  OrganismType = character(),
  Gen_T = numeric(),
  Net_Rep_Rate = numeric(),
  Do_Iteroparity = numeric(),
  Ro_Senescence = numeric(),
  Age_Maturity = numeric(),
  stringsAsFactors = FALSE
)

# calculating LHTs loop through each species
for (species in comadre_collapsed$SpeciesAccepted) {
  ind <- cdb_check_species(comadre_collapsed, species, return_db = TRUE) 
  
  organism_type <- ind$OrganismType
  
  Gen_T <- gen_time(matU(ind[[1]])[[1]], matF(ind[[1]])[[1]])
  Net_Rep_Rate <- net_repro_rate(matU(ind[[1]])[[1]], matF(ind[[1]])[[1]])
  
  lx <- mpm_to_lx(matU(ind[[1]])[[1]])
  mx <- mpm_to_mx(matU(ind[[1]])[[1]], matF(ind[[1]])[[1]])
  
  Do_Iteroparity <- entropy_d(lx, mx)
  
  Ro_Senescence <- entropy_k_stage(matU(ind[[1]])[[1]])
  Age_Maturity <- mature_age(matU(ind[[1]])[[1]], matF(ind[[1]])[[1]])
  
  # Append the results to the data frame
  comadre_LHT <- rbind(comadre_LHT, data.frame(
    SpeciesAccepted = species,
    OrganismType = organism_type,
    Gen_T = Gen_T,
    Net_Rep_Rate = Net_Rep_Rate,
    Do_Iteroparity = Do_Iteroparity,
    Ro_Senescence = Ro_Senescence,
    Age_Maturity = Age_Maturity,
    stringsAsFactors = FALSE
  ))
}

## compadre 
# make dataframe to load LHTs into
compadre_LHT <- data.frame(
  SpeciesAccepted = character(),
  ORganismType = character(),
  Gen_T = numeric(),
  Net_Rep_Rate = numeric(),
  Do_Iteroparity = numeric(),
  Ro_Senescence = numeric(),
  Age_Maturity = numeric(),
  stringsAsFactors = FALSE
)

# calculating LHTs loop through each species
for (species in compadre_collapsed$SpeciesAccepted) {
  ind <- cdb_check_species(compadre_collapsed, species, return_db = TRUE) 
  
  organism_type <- ind$OrganismType
  
  Gen_T <- gen_time(matU(ind[[1]])[[1]], matF(ind[[1]])[[1]])
  Net_Rep_Rate <- net_repro_rate(matU(ind[[1]])[[1]], matF(ind[[1]])[[1]])
  
  lx <- mpm_to_lx(matU(ind[[1]])[[1]])
  mx <- mpm_to_mx(matU(ind[[1]])[[1]], matF(ind[[1]])[[1]])
  
  Do_Iteroparity <- entropy_d(lx, mx)
  Ro_Senescence <- entropy_k(lx)
  
  Age_Maturity <- mature_age(matU(ind[[1]])[[1]], matF(ind[[1]])[[1]])
  
  # Append the results to the data frame
  compadre_LHT <- rbind(compadre_LHT, data.frame(
    SpeciesAccepted = species,
    OrganismType = organism_type, 
    Gen_T = Gen_T,
    Net_Rep_Rate = Net_Rep_Rate,
    Do_Iteroparity = Do_Iteroparity,
    Ro_Senescence = Ro_Senescence,
    Age_Maturity = Age_Maturity,
    stringsAsFactors = FALSE
  ))
}

LHT <- rbind(compadre_LHT, comadre_LHT)

#make adjustments so it fits phylogenetic analysis:
LHT$SpeciesAccepted[LHT$SpeciesAccepted=="Chamaecrista lineata var. keyensis"]="Chamaecrista lineata"
LHT$SpeciesAccepted[LHT$SpeciesAccepted=="Echinospartum ibericum subsp. algibicum"]="Echinospartum ibericum"
LHT$SpeciesAccepted[LHT$SpeciesAccepted=="Eriogonum longifolium var. gnaphalifolium"]="Eriogonum longifolium"
LHT$SpeciesAccepted[LHT$SpeciesAccepted=="Oenothera coloradensis subsp. coloradensis"]="Oenothera coloradensis"
LHT$SpeciesAccepted[LHT$SpeciesAccepted=="Petrocoptis pyrenaica subsp. pseudoviscosa"]="Petrocoptis pyrenaica"

write.csv(LHT, "LifeHistoryTraits.csv", row.names = FALSE)

## Load Phylogeny ----
#assembly tree based on a larger database of species, to provide more robust phylogeny
phyl <- read.csv("Phylogeny_species.csv")

#make adjustments so it fits phylogenetic analysis:
phyl$SpeciesAccepted[phyl$SpeciesAccepted=="Chamaecrista lineata var. keyensis"]="Chamaecrista lineata"
phyl$SpeciesAccepted[phyl$SpeciesAccepted=="Echinospartum ibericum subsp. algibicum"]="Echinospartum ibericum"
phyl$SpeciesAccepted[phyl$SpeciesAccepted=="Eriogonum longifolium var. gnaphalifolium"]="Eriogonum longifolium"
phyl$SpeciesAccepted[phyl$SpeciesAccepted=="Oenothera coloradensis subsp. coloradensis"]="Oenothera coloradensis"
phyl$SpeciesAccepted[phyl$SpeciesAccepted=="Petrocoptis pyrenaica subsp. pseudoviscosa"]="Petrocoptis pyrenaica"

resolved_names <- tnrs_match_names(phyl$SpeciesAccepted)

#change phylogeny ID to close relative which is resolved (so not to loose data)
resolved_names$ott_id[resolved_names$unique_name == "Polemonium van-bruntiae" & resolved_names$ott_id == 5738006] <- 3902247
resolved_names$ott_id[resolved_names$unique_name == "Asplenium x adulterinum" & resolved_names$ott_id == 16952] <- 235392

rn<- na.omit(resolved_names$ott_id)
tree <- tol_induced_subtree(ott_ids = rn)

tree$tip.label <- gsub("_ott.*$", "", tree$tip.label)
tree$tip.label <- gsub("_", " ", tree$tip.label)

# resolve multichomoty
tree  <- multi2di(tree)
# calculate branch length
tree <- compute.brlen(tree)

# resolve differences in species names (names synonymous)
tree$tip.label[tree$tip.label == "Triadica sebifera"] <- "Sapium sebiferum"
tree$tip.label[tree$tip.label == "Phyllanthus emblica (species in kingdom Archaeplastida)"] <- "Phyllanthus emblica"
tree$tip.label[tree$tip.label == "Erythranthe cardinalis"] <- "Mimulus cardinalis"
tree$tip.label[tree$tip.label == "Erythranthe lewisii"] <- "Mimulus lewisii"
tree$tip.label[tree$tip.label == "Polemonium vanbruntiae"] <- "Polemonium van-bruntiae"
tree$tip.label[tree$tip.label == "Polygonum basiramium"] <- "Polygonella basiramia"
tree$tip.label[tree$tip.label == "Asplenium viride"] <- "Asplenium adulterinum"
tree$tip.label[tree$tip.label == "Anser caerulescens"]  <- "Chen caerulescens"

write.tree(tree, "Tree.tre")
