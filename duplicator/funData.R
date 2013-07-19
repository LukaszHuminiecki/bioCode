###############################################################################
## this file is called: funData.R
###############################################################################

# print this file name (useful for slurm debugging)
print("funData.R")

# initialize environment if used to store results generated here (but loaded env may be reused)
env_duplicator <- new.env()

# load all required packages: library, install.packages(), R CMD INSTALL
library(hgu95a.db)

# source all R functions
source("funBasic.R")
source("pacGraph.R")

# define all constants, filenames, parameters, read datasets, etc
pre_output_file <- "env_duplicator"

# process args: args <- commandArgs(); print (args); if (!is.na(args[3])) data_file  <-  args[3];

# fallback arguments if no args: if (is.na(args[3])) data_file  <- "env_data_hs_litendata_v1";

# define directory prefix, use as: data_file <- paste0(directory_path, pre_data_file);
directory_path <- "~/nb/DATA"
#data_file <- paste0(directory_path, pre_data_file);
output_file <- paste0(directory_path, pre_output_file);
allaffypem <- read.table("data/expression/allaffypem")
FULL <- TRUE

# load all dependencies (other R environments)
load("data/cancerMutationLists/env_cancerMutationLists")
attach(env_cancerMutationLists)
load("data/networks/hcsm")
attach(hcsm)
load("data/duplicator/env_duplicator_vectors")
attach(env_duplicator_vectors)
load("data/duplicator/env_duplicator_base")
attach(env_duplicator_base)
print("Loaded envs")
load("data/cancerMutationLists/env_cancerMutationLists")
attach(env_cancerMutationLists)
load("data/networks/hcsm")
attach(hcsm)
load("data/duplicator/env_duplicator_vectors")
attach(env_duplicator_vectors)
load("data/duplicator/env_duplicator_base")
attach(env_duplicator_base)

# other: perhaps something which depends on one or more of the above

###############################################################################
## Routine code starts here:
###############################################################################

# check timing
t1 <- Sys.time()
print(t1)

# prepare arguments
duplicationWave <- duplicationVector2duplicationWave(dupEvent12_genes_LL_hs)
duplicationWaveHUMAN <- unique(duplicationVector2duplicationWave(dupEvent12_genes_LL_hs[dupEvent12_genes_LL_hs$taxon=="HUMAN",]))
duplicationWaveHomoPanGorilla  <- unique(duplicationVector2duplicationWave(dupEvent12_genes_LL_hs[dupEvent12_genes_LL_hs$taxon=="Homo/Pan/Gorilla",]))
duplicationWaveCatarrhini <- unique(duplicationVector2duplicationWave(dupEvent12_genes_LL_hs[dupEvent12_genes_LL_hs$taxon=="Catarrhini",]))
duplicationWaveEutheria <- unique(duplicationVector2duplicationWave(dupEvent12_genes_LL_hs[dupEvent12_genes_LL_hs$taxon=="Eutheria",]))
duplicationWaveTheria <- unique(duplicationVector2duplicationWave(dupEvent12_genes_LL_hs[dupEvent12_genes_LL_hs$taxon=="Theria",]))
duplicationWaveAmniota <- unique(duplicationVector2duplicationWave(dupEvent12_genes_LL_hs[dupEvent12_genes_LL_hs$taxon=="Amniota",]))
duplicationWaveTetrapoda <- unique(duplicationVector2duplicationWave(dupEvent12_genes_LL_hs[dupEvent12_genes_LL_hs$taxon=="Tetrapoda",]))
duplicationWaveEuteleostomi <- unique(duplicationVector2duplicationWave(dupEvent12_genes_LL_hs[dupEvent12_genes_LL_hs$taxon=="Euteleostomi",]))
duplicationWaveChordata <- unique(duplicationVector2duplicationWave(dupEvent12_genes_LL_hs[dupEvent12_genes_LL_hs$taxon=="Chordata",]))
duplicationWaveDeuterostomia <- unique(duplicationVector2duplicationWave(dupEvent12_genes_LL_hs[dupEvent12_genes_LL_hs$taxon=="Deuterostomia",]))
duplicationWaveBilateria <- unique(duplicationVector2duplicationWave(dupEvent12_genes_LL_hs[dupEvent12_genes_LL_hs$taxon=="Bilateria",]))
duplicationWaveEumetazoa <- unique(duplicationVector2duplicationWave(dupEvent12_genes_LL_hs[dupEvent12_genes_LL_hs$taxon=="Eumetazoa",]))
duplicationWaveMetazoa <- unique(duplicationVector2duplicationWave(dupEvent12_genes_LL_hs[dupEvent12_genes_LL_hs$taxon=="Metazoa",]))
duplicationWaveFungiMetazoa <- unique(duplicationVector2duplicationWave(dupEvent12_genes_LL_hs[dupEvent12_genes_LL_hs$taxon=='Fungi\\Metazoa',]))
duplicationWaveEukaryota <- unique(duplicationVector2duplicationWave(dupEvent12_genes_LL_hs[dupEvent12_genes_LL_hs$taxon=="Eukaryota",]))

duplicationWave2 <- (dupEvent12_genes_LL_hs2)
duplicationWaveCatarrhini <- unique(duplicationVector2duplicationWave(dupEvent12_genes_LL_hs[dupEvent12_genes_LL_hs$taxon=="Catarrhini",]))
duplicationWaveEutheria <- unique(duplicationVector2duplicationWave(dupEvent12_genes_LL_hs[dupEvent12_genes_LL_hs$taxon=="Eutheria",]))
duplicationWaveTheria <- unique(duplicationVector2duplicationWave(dupEvent12_genes_LL_hs[dupEvent12_genes_LL_hs$taxon=="Theria",]))
duplicationWaveAmniota <- unique(duplicationVector2duplicationWave(dupEvent12_genes_LL_hs[dupEvent12_genes_LL_hs$taxon=="Amniota",]))
duplicationWaveTetrapoda <- unique(duplicationVector2duplicationWave(dupEvent12_genes_LL_hs[dupEvent12_genes_LL_hs$taxon=="Tetrapoda",]))
duplicationWaveEuteleostomi <- unique(duplicationVector2duplicationWave(dupEvent12_genes_LL_hs[dupEvent12_genes_LL_hs$taxon=="Euteleostomi",]))
duplicationWaveChordata <- unique(duplicationVector2duplicationWave(dupEvent12_genes_LL_hs[dupEvent12_genes_LL_hs$taxon=="Chordata",]))
duplicationWaveDeuterostomia <- unique(duplicationVector2duplicationWave(dupEvent12_genes_LL_hs[dupEvent12_genes_LL_hs$taxon=="Deuterostomia",]))
duplicationWaveBilateria <- unique(duplicationVector2duplicationWave(dupEvent12_genes_LL_hs[dupEvent12_genes_LL_hs$taxon=="Bilateria",]))
duplicationWaveEumetazoa <- unique(duplicationVector2duplicationWave(dupEvent12_genes_LL_hs[dupEvent12_genes_LL_hs$taxon=="Eumetazoa",]))
duplicationWaveMetazoa <- unique(duplicationVector2duplicationWave(dupEvent12_genes_LL_hs[dupEvent12_genes_LL_hs$taxon=="Metazoa",]))
duplicationWaveFungiMetazoa <- unique(duplicationVector2duplicationWave(dupEvent12_genes_LL_hs[dupEvent12_genes_LL_hs$taxon=='Fungi\\Metazoa',]))
duplicationWaveEukaryota <- unique(duplicationVector2duplicationWave(dupEvent12_genes_LL_hs[dupEvent12_genes_LL_hs$taxon=="Eukaryota",]))

duplicationWave2 <- (dupEvent12_genes_LL_hs2)
duplicationWaveHUMAN2 <- unique((dupEvent12_genes_LL_hs2[dupEvent12_genes_LL_hs2$taxon=="HUMAN",]))
duplicationWaveHomoPanGorilla2  <- unique((dupEvent12_genes_LL_hs2[dupEvent12_genes_LL_hs2$taxon=="Homo/Pan/Gorilla",]))
duplicationWaveCatarrhini2 <- unique((dupEvent12_genes_LL_hs2[dupEvent12_genes_LL_hs2$taxon=="Catarrhini",]))
duplicationWaveEutheria2 <- unique((dupEvent12_genes_LL_hs2[dupEvent12_genes_LL_hs2$taxon=="Eutheria",]))
duplicationWaveTheria2 <- unique((dupEvent12_genes_LL_hs2[dupEvent12_genes_LL_hs2$taxon=="Theria",]))
duplicationWaveAmniota2 <- unique((dupEvent12_genes_LL_hs2[dupEvent12_genes_LL_hs2$taxon=="Amniota",]))
duplicationWaveTetrapoda2 <- unique((dupEvent12_genes_LL_hs2[dupEvent12_genes_LL_hs2$taxon=="Tetrapoda",]))
duplicationWaveEuteleostomi2 <- unique((dupEvent12_genes_LL_hs2[dupEvent12_genes_LL_hs2$taxon=="Euteleostomi",]))
duplicationWaveChordata2 <- unique((dupEvent12_genes_LL_hs2[dupEvent12_genes_LL_hs2$taxon=="Chordata",]))
duplicationWaveDeuterostomia2 <- unique((dupEvent12_genes_LL_hs2[dupEvent12_genes_LL_hs2$taxon=="Deuterostomia",]))
duplicationWaveBilateria2 <- unique((dupEvent12_genes_LL_hs2[dupEvent12_genes_LL_hs2$taxon=="Bilateria",]))
duplicationWaveEumetazoa2 <- unique((dupEvent12_genes_LL_hs2[dupEvent12_genes_LL_hs2$taxon=="Eumetazoa",]))
duplicationWaveMetazoa2 <- unique((dupEvent12_genes_LL_hs2[dupEvent12_genes_LL_hs2$taxon=="Metazoa",]))
duplicationWaveFungiMetazoa2 <- unique((dupEvent12_genes_LL_hs2[dupEvent12_genes_LL_hs2$taxon=='Fungi\\Metazoa',]))
t2 <- Sys.time()
print(t2)
print("Prepared waves")

#
# orderedTaxa, orderedTaxaColors, orderedTaxaShape, orderedTaxaList
#
 
env_duplicator$orderedTaxa <- orderedTaxa <- c("HUMAN", "Homo/Pan/Gorilla", "Catarrhini", "Eutheria" , "Theria" ,"Amniota", "Tetrapoda", "Euteleostomi", "Chordata", "Deuterostomia", "Bilateria","Eumetazoa", "Metazoa", "Fungi/Metazoa", "Eukaryota")
env_duplicator$orderedTaxaColors <- orderedTaxaColors <- c("red", "red", "orange", "orange", "yellow", "yellow", "white", "green", "green", "blue", "blue","brown", "brown", "grey", "grey")
env_duplicator$orderedTaxaShape <- orderedTaxaShape <- c("circle", "box", "circle", "box" , "circle", "box", "circle", "circle", "box", "circle", "box", "circle", "box", "circle", "box")
orderedTaxaList <- list()
orderedTaxaList[[1]] <- orderedTaxa[1:7]
orderedTaxaList[[2]] <- orderedTaxa[8]
orderedTaxaList[[3]] <- orderedTaxa[9:10]
orderedTaxaList[[4]] <- orderedTaxa[11]
env_duplicator$orderedTaxaList <- orderedTaxaList
print("Prepared orderedTaxa")
