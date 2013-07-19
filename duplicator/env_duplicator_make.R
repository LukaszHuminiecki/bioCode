###############################################################################
## this file is called: env_duplicator_make.R
###############################################################################

# print this file name (useful for slurm debugging)
print("env_duplicator_make.R")

# initialize environment if used to store results generated here (but loaded env may be reused)
env_duplicator <- new.env()

# load all required packages: library, install.packages(), R CMD INSTALL
library(hgu95a.db)

# source all R functions
source("~/fantom5/R/duplicator/funBasic.R")
source("~/fantom5/R/duplicator/pacGraph.R")

# define all constants, filenames, parameters, read datasets, etc
pre_output_file <- "env_duplicator"

# process args: args <- commandArgs(); print (args); if (!is.na(args[3])) data_file  <-  args[3];

# fallback arguments if no args: if (is.na(args[3])) data_file  <- "env_data_hs_litendata_v1";

# define directory prefix, use as: data_file <- paste0(directory_path, pre_data_file);
directory_path <- "~/fantom5/R/nb/DATA"
#data_file <- paste0(directory_path, pre_data_file);
output_file <- paste0(directory_path, pre_output_file);
allaffypem <- read.table("~/fantom5/assets/allaffypem")
FULL <- TRUE

# load all dependencies (other R environments)
load("~/fantom5/BILSCancerProject/Duplicator/data/cancerMutationLists/env_cancerMutationLists")
attach(env_cancerMutationLists)
load("~/fantom5/BILSCancerProject/Duplicator/data/networks/hcsm")
attach(hcsm)
load("~/nb/DATA/env_fantom5_vectors")
attach(env_fantom5_vectors)
load("~/nb/DATA/env_fantom5_base")
attach(env_fantom5_base)
print("Loaded envs")
load("~/fantom5/BILSCancerProject/Duplicator/data/cancerMutationLists/env_cancerMutationLists")
attach(env_cancerMutationLists)
load("~/fantom5/BILSCancerProject/Duplicator/data/networks/hcsm")
attach(hcsm)
load("~/nb/DATA/env_fantom5_vectors")
attach(env_fantom5_vectors)
load("~/nb/DATA/env_fantom5_base")
attach(env_fantom5_base)

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

# prepare igraphs from graphNELs

igraph <- igraph.from.graphNEL(hcsm1)
graph <- hcsm1

vn <- V(igraph)$name
x <- duplicationWaveEuteleostomi2$primary_acc.x
y <- duplicationWaveEuteleostomi2$primary_acc.y
xx <- x[x %in% vn & y %in% vn]
yy <- y[x %in% vn & y %in% vn]
t3 <- Sys.time()
print("Loaded graph objects")

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

# must convert functions to feed on parallel vectors rather than data.frame because this is trouble

if(FULL) {
CancerMutationGeneList_pre <- Broad2000
CancerMutationGeneListHUGOids <- convertHUGO2EG(CancerMutationGeneList_pre[,2])
CancerMutationGeneList <- as.data.frame(cbind(CancerMutationGeneList_pre[,1], CancerMutationGeneListHUGOids))
colnames(CancerMutationGeneList) <- c("score", "id")
print("Loaded cancer mutations")
t4 <- Sys.time()
print(t4)
}

#
# networkAnnotatedWithDuplication
#

print("networkAnnotatedWithDuplication")
networkAnnotatedWithDuplication <- annotateAttributeDuplication(graph, duplicationWave, orderedTaxa)
# networkAnnotatedWithDuplicationMutation <- annotateAttributeMutation(networkAnnotatedWithDuplication, CancerMutationGeneList)
# networkAnnotatedWithDuplicationMutationTargets <- annotateAttributeDrugTarget(networkAnnotatedWithDuplicationMutation, CancerDrugTargets)
#inetworkAnnotatedWithDuplicationMutation <- igraph.from.graphNEL(networkAnnotatedWithDuplicationMutation)
t5 <- Sys.time()
print(t5)

#
# wholeNetworkAnalysisWrapper
#

print("wholeNetworkAnalysis calls")
igraph_wholeNetworkAnalysisWrapper <- wholeNetworkAnalysisWrapper_4igraph(igraph)
igraph_allNeighbourhoodsTriadCensus <- allNeighbourhoodsTriadCensus(igraph)
igraph_allNodeTriadCensus <- allNodeTriadCensus(igraph)
igraph_allTriadCensus <- allTriadCensus(igraph)
t6 <- Sys.time()
print(t6)

#
# rAll and rAllList
#

# SLOW

if(FALSE) {
print("rAll")
rAll_igraph <- list()
rAll_graph <- list()
for(txTarget in orderedTaxa) {
	print(txTarget)
	duplicationWaveTxTarget <- unique(duplicationVector2duplicationWave(dupEvent12_genes_LL_hs[dupEvent12_genes_LL_hs$taxon %in% txTarget,]))
	rAll_igraph[[txTarget[1]]] <- duplicationWavesAnalysisWrapper_4igraph(duplicationWaveTxTarget, igraph)
	rAll_graph[[txTarget[1]]] <- duplicationWavesAnalysisWrapper_4graph(duplicationWaveTxTarget, graph)
}
t7 <- Sys.time()
print(t7)

print("rAllList")
rAllList_igraph <- list()
rAllList_graph <- list()
for(txTargetList in orderedTaxaList) {
	print(txTarget)
	duplicationWaveTxTarget <- unique(duplicationVector2duplicationWave(dupEvent12_genes_LL_hs[dupEvent12_genes_LL_hs$taxon %in% txTargetList,]))
	rAllList_igraph[[txTargetList[1]]] <- duplicationWavesAnalysisWrapper_4igraph(duplicationWaveTxTarget, igraph)
	rAllList_graph[[txTargetList[1]]] <- duplicationWavesAnalysisWrapper_4graph(duplicationWaveTxTarget, graph)
	}
}
t8 <- Sys.time()
print(t8)

#
# get layout with colors/shapes for taxa and then draw
#

if(FULL) {
print("initializeLayoutWithTaxaShapeColor")
#lg=initializeLayoutWithTaxaShapeColor(networkAnnotatedWithDuplication, "neato", orderedTaxaColors, orderedTaxaShape, orderedTaxa)
#renderGraph(lg)
}
t9 <- Sys.time()
print(t9)

#
# Make dupEvent12_genes_LL_hs2_Annotated by calling
# duplicationPairsAnalysisWrapper
#

# SLOW

if(FALSE) {

ancestralGraph <- list()
print("dupEvent12_genes_LL_hs2_Annotated")
dupEvent12_genes_LL_hs2_Annotated <- NULL
x <- intersect(V(igraph), dupEvent12_genes_LL_hs2$primary_acc.x)
y <- intersect(V(igraph), dupEvent12_genes_LL_hs2$primary_acc.y)
for (i in 1:length(dupEvent12_genes_LL_hs2[,1])) {
	print(c(i, "out of", length(dupEvent12_genes_LL_hs2[,1])))
	xx <- dupEvent12_genes_LL_hs2$primary_acc.x[i]
	yy <- dupEvent12_genes_LL_hs2$primary_acc.y[i]
	dupEvent12_genes_LL_hs2_Annotated_DuplicationPairsAnalysisWrapper <- list()
	if ((xx %in% x) && (yy %in% y)) {
		dupEvent12_genes_LL_hs2_Annotated_DuplicationPairsAnalysisWrapper <- duplicationPairsAnalysisWrapper(xx, yy, igraph)
		}	
# now call wrapper to contract.vertices and get ancestral graph for each taxon 
ancestralGraph <- list()
for(txTarget in orderedTaxa) {
	duplicationDoubleWaveTxTarget <- unique(dupEvent12_genes_LL_hs2[dupEvent12_genes_LL_hs$taxon %in% txTarget,])
	ancestralNetwork[[txTarget]] <- getAncestralNetwork(duplicationDoubleWaveTxTargetOverlap$primary_acc.x, duplicationDoubleWaveTxTargetOverlap$primary_acc.y, igraph)
	}

t11 <- Sys.time()
print(t11)

#
# annotateEdgesWithPC
#

if(FULL) {
print("annotateEdgesWithPC")
affyidExpressionMatrix <- allaffypem[1:48]
igraphannotateEdgesWithPC <- annotateEdgesWithPC(igraph, affyidExpressionMatrix)
}
t12 <- Sys.time()
print(t12)

#
#
#

if(FULL) {
env_duplicator$dupEvent12_genes_LL_hs2_Annotated <- dupEvent12_genes_LL_hs2_Annotated 
save(env_duplicator, file = env_duplicator_file_name)
}
t13 <- Sys.time()
print(t13)

#
#
#

print("Ended env_duplicator_make.R")
