#   
#  ##   #  # ###  #   #  ###  ###  ###  ##   ###
#  # #  #  # #  # #   # #    #   #  #  #  #  #  #
#  # #  #  # ###  #   # #    #####  #  #  #  ###
#  # #  #  # #    #   # #    #   #  #  #  #  # #
#  ##    ##  #    ### #  ### #   #  #   ##   #  #
#

# print this file name (useful for slurm debugging)
print("funBasic.R")

#
# grep "function" to see all functions
#

#
# New functions for dealing with Tfbs conservation, Tfbs distribution by taxa and loci from ENCODE
#
# can be used with Tfbs_x_sum_taxon; Tfbs_y_sum_taxon; Tfbs_union_sum_taxon; Tfbs_intersection_sum_taxon
# usage example, return percentage of genes from given taxon associated with given Tfbs
# top_taxa_named_vector(Tfbs_x_sum_taxon, "CTCF")
#

top_taxa_named_vector <- function(named_vector, target) {
x <- sort(table(names(named_vector[named_vector == target])))    # show me top taxa for given Tfbs
y <- sort(table(names(named_vector)))                            # background
result <- sort(x[intersect(names(x), names(y))] / y[intersect(names(x), names(y))] *100) 
result
}

top_tfbs_per_taxon <- function(named_vector, target) {
x <- sort(table(named_vector[names(named_vector) == target]))    # show me top taxa for given Tfbs
y <- sort(table(named_vector))                                   # background
result <- sort(x[intersect(names(x), names(y))])
result
}

top_tfbs_per_taxon_normalized <- function(named_vector, target) {
x <- sort(table(named_vector[names(named_vector) == target]))    # show me top taxa for given Tfbs
y <- sort(table(named_vector))  	               	         # background
result <- sort(x[intersect(names(x), names(y))] / y[intersect(names(x), names(y))] *100)
result
}

tfbs_per_taxon <- function(named_vector, target) {
named_vector[names(named_vector) == target]
}

# can be used with Tfbs_x_sum; Tfbs_y_sum; Tfbs_union_sum; Tfbs_intersection_sum
# usage example, return percentage of genes from given taxon associated with given Tfbs
# top_loci_named_vector(Tfbs_x_sum, "Pol2")
top_loci_named_vector <- function(named_vector, target) {
x <- sort(table(names(named_vector[named_vector == target])))    # show me top taxa for given Tfbs
y <- sort(table(names(named_vector)))                            # background
result <- sort(x[intersect(names(x), names(y))])
result
}

top_loci_normalized_named_vector <- function(named_vector, target) {
x <- sort(table(names(named_vector[named_vector == target])))    # show me top taxa for given Tfbs
y <- sort(table(names(named_vector)))                            # background
result <- sort(x[intersect(names(x), names(y))] / y[intersect(names(x), names(y))] *100) # variant return
result
}

#
# Tfbs_usage_per_name_in_named_vector return a heatmap ready data frame where
# Tfbs usage has been cross-tabulated against names (taxa or loci)
# This could however may be needed to be normalized by number of duplications in each taxon
# or perhaps even number of Tfbs in each locus
#

Tfbs_usage_per_name_in_named_vector <- function(named_vector) {
named_vector_copy <- named_vector
names(named_vector) <- NULL
named_vector_as_data_frame <- table(as.data.frame(cbind(named_vector, names(named_vector_copy))))
}

#
# annotateEdgesWithPC annotates interaction network with PC values inferred from Affy
#

annotateEdgesWithPC <- function(g, affyidExpressionMatrix) {
map <- as.data.frame(hgu95aENTREZID)
PCAll <- NULL
edgeList<-get.edgelist(g)
for(i in 1:length(edgeList[,1]))
{
n1 <- edgeList[i,1]
n2 <- edgeList[i,2]
m1 <- map[map$gene_id==n1,]
m2 <- map[map$gene_id==n2,]
#only take unique mappers at this stage!
if(length(m1[,1])==1 && length(m2[,1])==1)
 {
 # not get PC from precalc expression and put that as the edge!
 x <- affyidExpressionMatrix[affyidExpressionMatrix[,1]==m1$probe_id,][2:length(affyidExpressionMatrix)]
 y <- affyidExpressionMatrix[affyidExpressionMatrix[,1]==m2$probe_id,][2:length(affyidExpressionMatrix)]
 xn <- (as.numeric(x))
 yn <- (as.numeric(y))
 PCAll[i] <- PC <- cor(xn, yn)
 if(is.na(PC)) {next}
 g <- set.edge.attribute(g, "weight", c(n1, n2), PC)
 if(PC > 0.8){g <- set.edge.attribute(g, "color", c(n1, n2), "red")}
 else if(PC > 0.6){g <- set.edge.attribute(g, "color", c(n1, n2), "orange")}
 else if(PC > 0.4){g <- set.edge.attribute(g, "color", c(n1, n2), "yellow")}
 else if(PC < 0.0){g <- set.edge.attribute(g, "color", c(n1, n2), "green")}
 }}
return(g)
}

#
# getAncestralNetwork()
# x <- c("5319", "4744")
# y <- c("9525", "6908")
# getAncestralNetwork(xx, yy, iNetwork)
#
# still problem with multiple mappers after grep 
# the condition has length > 1 and only the first element will be used
#
# in general, this caused a lot of problems because of mapping vector
# mapping vector is index from 1:max. always max nodes will be created even of no mapping provided. 
# this was very flaky. decided to first contract and then select for all but empty nodes
# difficult function to get to work.
#
#

getAncestralNetwork <- function(x, y, iNetwork) {
x <- as.character(x)
y <- as.character(y)
verts <- V(iNetwork)
vn <- verts$name  
mappingVector <- 1:length(vn)
for(i in 1:length(x)) {
 xPos <- yPos <- NULL
 xPos <- grep(x[i], vn)
 yPos <- grep(y[i], vn)
 if(length(xPos) == 0 | length(yPos) == 0) {next;}
 if(xPos > yPos) {mappingVector[xPos] <- mappingVector[yPos]}
   else {mappingVector[yPos] <- mappingVector[xPos]}}
 iNetwork <- contract.vertices(iNetwork, mappingVector)
 v <- 1:length(index)
 index <- mappingVector == 1:length(mappingVector)
vec <- 1:length(mappingVector)
selector <- vec[index]
 return(subgraph(iNetwork, selector))
}

#
#
#

duplicationPairsAnalysisWrapper <- function(x, y, iNetwork) {
res <- list()
	res[["shortestPathsAll"]] <- shortest.paths(iNetwork, v=x, to=y, mode="all") # inf?
	res[["shortestPathsOut"]] <- shortest.paths(iNetwork, v=x, to=y, mode="out")
	res[["shortestPathsIn"]] <- shortest.paths(iNetwork, v=x, to=y, mode="in")
	res[["xBetweenness"]] <- betweenness(iNetwork, v=x, directed = TRUE, weights = NULL)
	res[["yBetweenness"]] <- betweenness(iNetwork, v=y, directed = TRUE, weights = NULL)
	res[["xDegreeAll"]] <- degree(iNetwork, v=x, loops = TRUE, mode = "all")
	res[["xDegreeOut"]] <- degree(iNetwork, v=x, loops = TRUE, mode = "out")
	res[["xDegreeIn"]]  <- degree(iNetwork, v=x, loops = TRUE, mode = "in")
	res[["yDegreeAll"]] <- degree(iNetwork, v=y, loops = TRUE, mode = "all")
	res[["yDegreeOut"]] <- degree(iNetwork, v=y, loops = TRUE, mode = "out")
	res[["yDegreeIn"]]  <- degree(iNetwork, v=y, loops = TRUE, mode = "in")
	res[["xEccentricityAll"]] <- eccentricity(iNetwork, vids=x, mode="all")
	res[["xEccentricityOut"]] <- eccentricity(iNetwork, vids=x, mode="out")
	res[["xEccentricityIn"]]  <- eccentricity(iNetwork, vids=x, mode="in")
	res[["yEccentricityAll"]] <- eccentricity(iNetwork, vids=y, mode="all")
	res[["yEccentricityOut"]] <- eccentricity(iNetwork, vids=y, mode="out")
	res[["yEccentricityIn"]]  <- eccentricity(iNetwork, vids=y, mode="in")
	res[["similarityJaccardAll"]] <- similarity.jaccard(iNetwork, vids = c(x, y), mode = "all", loops = FALSE)[2]
	res[["similarityJaccardOut"]] <- similarity.jaccard(iNetwork, vids = c(x, y), mode = "out", loops = FALSE)[2]
	res[["similarityJaccardIn"]] <- similarity.jaccard(iNetwork, vids = c(x, y), mode = "in", loops = FALSE)[2]
	res[["similarityDiceAll"]] <- similarity.dice(iNetwork, vids = c(x, y), mode = "all", loops = FALSE)[2]
	res[["similarityDiceOut"]] <- similarity.dice(iNetwork, vids = c(x, y), mode = "out", loops = FALSE)[2]
	res[["similarityDiceIn"]] <- similarity.dice(iNetwork, vids = c(x, y), mode = "in", loops = FALSE)[2]
	res[["similarityInvlogweightedAll"]] <- similarity.invlogweighted(iNetwork, vids = c(x, y), mode = "all")[2]
	res[["similarityInvlogweightedOut"]] <- similarity.invlogweighted(iNetwork, vids = c(x, y), mode = "out")[2]
	res[["similarityInvlogweightedIn"]] <- similarity.invlogweighted(iNetwork, vids = c(x, y), mode = "in")[2]
return(res)
}

#
# nodesWithAttribute returns subset of graphs nodes with desired attribute value
# graphNEL network, character attributeType, character attributeName
# for example, attributeType="duplicationAge", attributeName="HUMAN"
# returns: character vector
#

getNodesWithAttribute <- function(network, attributeType, attributeName) {
nodeDataNetwork <- as.vector(unlist(nodeData(network, nodes(network), "duplicationAge")))
index <- nodeDataNetwork == attributeName # why all this NA
x <- nodes(network)[index]
x[!is.na(x)] # omit N values: nodes without attribute of attributeType
}

#
# initializeLayoutWithTaxaShapeColor API
# initialize graphLayout with layoutType, orderedTaxaColors, orderedTaxaShape, orderedTaxa
# layout, layoutType=""neato"" or "dot" or "twopi"
#

initializeLayoutWithTaxaShapeColor <- function(network, layoutType="neato", orderedTaxaColors, orderedTaxaShape, orderedTaxa) {
lg=layoutGraph(x=network, layoutType=layoutType)
lgNodeData <- nodeData(lg, nodes(lg), "duplicationAge")
lgNodeDataVector <- unlist(lgNodeData)
for (i in 1:length(lgNodeDataVector)) {
        nodeRenderInfo(lg)$fill[i] <- orderedTaxaColors[orderedTaxa==lgNodeDataVector[i]]
        nodeRenderInfo(lg)$shape[i] <- orderedTaxaShape[orderedTaxa==lgNodeDataVector[i]]
        }
return(lg)
}

#
# subigraph()
#

subIgraph <- function(g, duplicationWaveTxTarget) {
sg <- induced.subgraph(g, V(g)[V(g) %in% duplicationWaveTxTarget[,2]])
return(sg)
}

#
# Wrapper for duplicationWave analysis
#

duplicationWavesAnalysisWrapper_4graph <- function(duplicationWaveTxTarget, graph) {
res <- list()
res[["subsetDegree"]] <- subsetDegree(duplicationWaveTxTarget, graph)
# res[["subsetGraphBetween"]] <- subsetGraphBetween(duplicationWaveTxTarget, graph) # Error in rep(1:nN, elem) : invalid 'times' argument
# res[["subsetGraphDegree"]] <- subsetGraphDegree(duplicationWaveTxTarget, graph) # object 'mean_Degree' not found
# res[["subsetBetween001"]] <- subsetBetween001(duplicationWaveTxTarget, graph) # object 'rand_between_top' not found
# res[["subNetwork"]] <- subNetwork(duplicationWaveTxTarget, graph) # could not find function "subNetwork"
res[["subsetBetween"]] <- subsetBetween(duplicationWaveTxTarget, graph)
res[["subsetBetweenTop"]] <- subsetBetweenTop(duplicationWaveTxTarget, graph)
res[["subsetGraph"]] <- subsetGraph(duplicationWaveTxTarget, graph)
# res[["subsetGraphCC"]] <- subsetGraphCC(duplicationWaveTxTarget, graph) # Error in rep(1:nN, elem) : invalid 'times' argument
# res[["subsetTransl"]] <- subsetTransl(duplicationWaveTxTarget, graph) # Error in intersect(transl[, 1], testset)
#res[["dyadicity"]] <- dyadicity(duplicationWaveTxTarget, graph) # object 'm11' not found
#res[["heterophilicity"]] <- heterophilicity(duplicationWaveTxTarget, graph) # object 'm11' not found
#res[["nodeDistribution"]] <- nodeDistribution(duplicationWaveTxTarget, graph) # object 'm11' not found
}

duplicationWavesAnalysisWrapper_4igraph <- function(duplicationWaveTxTarget, igraph) {
sg <- subIgraph(igraph, duplicationWaveTxTarget)
res <- wholeNetworkAnalysisWrapper(sg)
return(res)
}

#
# wholeNetworkAnalysisWrapper
#

wholeNetworkAnalysisWrapper_4igraph <- function(igraph) {

if(!is.igraph(igraph)) {print("Usage: wholeNetworkAnalysisWrapper works with igraph!")}

res <- list()
print("DegreeAll")
res[["DegreeAll"]] <- degree(igraph, loops = TRUE, mode = "all")
res[["DegreeOut"]] <- degree(igraph, loops = TRUE, mode = "out")
res[["DegreeIn"]]  <- degree(igraph, loops = TRUE, mode = "in")
print("averagePathLengthFF")
res[["averagePathLengthFF"]] <- average.path.length(igraph, directed=FALSE, unconnected=FALSE)
res[["averagePathLengthFT"]] <- average.path.length(igraph, directed=FALSE, unconnected=TRUE)
res[["averagePathLengthTF"]] <- average.path.length(igraph, directed=TRUE, unconnected=FALSE)
res[["averagePathLengthTT"]] <- average.path.length(igraph, directed=TRUE, unconnected=TRUE)
print("edgeBetweennessCommunity")
res[["edgeBetweennessCommunity"]] <- edge.betweenness.community(igraph)
res[["triadCensus"]] <- triad.census(igraph)
res[["graphMotifs3"]] <- graph.motifs(igraph, size = 3)
res[["graphMotifs4"]] <- graph.motifs(igraph, size = 4)
res[["graphMotifsNo3"]] <- graph.motifs.no(igraph, size = 3)
res[["graphMotifsNo4"]] <- graph.motifs.no(igraph, size = 4)
res[["graphMotifsEst3"]] <- graph.motifs.est(igraph,size = 3)
res[["graphMotifsEst4"]] <- graph.motifs.est(igraph,size = 4)
return(res)
}

#
# allNeighbourhoodsTriadCensus
# for each vertex, return triads for its neighbours in relation to that vertex

allNeighbourhoodsTriadCensus <- function(graph) {
res <- list()
for(i in 1:length(V(graph))) {
	subGraph = graph.neighborhood(graph, order = 1, V(graph)[i], mode = 'all')[[1]]
	allMotifs = triad.census(subGraph)
	removeNode = delete.vertices(subGraph, V(subGraph)[1])
	node1Motifs = allMotifs - triad.census(removeNode)
	res[[i]] <- node1Motifs
	}
return(res)
}

#
# allNodeTriadCensus
# for each vertex, return its triads
#

allNodeTriadCensus <- function(graph) {
res <- list()
for(i in 1:length(V(graph))) {
        subGraph = graph.neighborhood(graph, order = 1, V(graph)[i], mode = 'all')[[1]]
        allMotifs = triad.census(subGraph)
        removeNode = delete.vertices(subGraph, V(subGraph)[1])
        res[[i]] <- triad.census(removeNode)
        }
return(res)
}

#
# allTriadCensus
# for each vertex, return its allMotifs triads
#

allTriadCensus <- function(graph) {
res <- list()
for(i in 1:length(V(graph))) {
        subGraph = graph.neighborhood(graph, order = 1, V(graph)[i], mode = 'all')[[1]]
        allMotifs = triad.census(subGraph)
        removeNode = delete.vertices(subGraph, V(subGraph)[1])
        res[[i]] <- allMotifs
        } 
return(res)
}

# \name{nodeDistribution}
# \description{calculates dyadicity, heterophilicity for a subset of nodes in a larger network}
# \usage{nodeDistribution(graphNEL hcsm1, character testset)}
# \examples{nodeDistribution(hcsm1, nodes(hcsm1))}

nodeDistribution <- function (testset, hcsm) {
uhcsm <- ugraph(hcsm)
testset_isnot <- setdiff(nodes(uhcsm), testset)
testset_is <- intersect(nodes(uhcsm), testset)
test_subGraph <- subGraph(as.character(testset_is), uhcsm)
negation_subGraph <- subGraph(as.character(testset_isnot), uhcsm)
M <- length(unlist(edges(uhcsm)))/2
N <- length(uhcsm@nodes)
m11observed <- length(unlist(edges(test_subGraph)))/2
n1          <- length(test_subGraph@nodes)
m00observed <- length(unlist(edges(negation_subGraph)))/2
n0          <- length(negation_subGraph@nodes)
m10observed <- M - m11 - m00
p <- 2*M/(N*(N-1))
m11expected <- n1*(n1-1)*p/2
expected <- n1*(N-n1)*p
dyadicity       <- m11observed / m11expected
heterophilicity <- m10observed / m10expected
c(dyadicity, heterophilicity)}

# \name{dyadicity}
# \description{calculates dyadicity for a subset of nodes in a larger network}
# \usage{dyadicity(graphNEL hcsm1, character testset)}
# \examples{dyadicity(hcsm1, nodes(hcsm1))}

dyadicity <- function (testset, hcsm) {
uhcsm <- ugraph(hcsm)
testset_isnot <- setdiff(nodes(uhcsm), testset)
testset_is <- intersect(nodes(uhcsm), testset)
test_subGraph <- subGraph(as.character(testset_is), uhcsm)
negation_subGraph <- subGraph(as.character(testset_isnot), uhcsm)
M <- length(unlist(edges(uhcsm)))/2
N <- length(uhcsm@nodes)
m11observed <- length(unlist(edges(test_subGraph)))/2
n1          <- length(test_subGraph@nodes)
m00observed <- length(unlist(edges(negation_subGraph)))/2
n0          <- length(negation_subGraph@nodes)
m10observed <- M - m11 - m00
p <- 2*M/(N*(N-1))
m11expected <- n1*(n1-1)*p/2
m10expected <- n1*(N-n1)*p
dyadicity       <- m11observed / m11expected
heterophilicity <- m10observed / m10expected
c(dyadicity)}

# \name{heterophilicity}
# \description{calculates heterophilicity for a subset of nodes in a larger network}
# \usage{heterophilicity(graphNEL hcsm1, character testset)}
# \examples{heterophilicity(hcsm1, nodes(hcsm1))}

heterophilicity <- function (testset, hcsm) {
uhcsm <- ugraph(hcsm)
testset_isnot <- setdiff(nodes(uhcsm), testset)
testset_is <- intersect(nodes(uhcsm), testset)
test_subGraph <- subGraph(as.character(testset_is), uhcsm)
negation_subGraph <- subGraph(as.character(testset_isnot), uhcsm)
M <- length(unlist(edges(uhcsm)))/2
N <- length(uhcsm@nodes)
m11observed <- length(unlist(edges(test_subGraph)))/2
n1          <- length(test_subGraph@nodes)
m00observed <- length(unlist(edges(negation_subGraph)))/2
n0          <- length(negation_subGraph@nodes)
m10observed <- M - m11 - m00
p <- 2*M/(N*(N-1))
m11expected <- n1*(n1-1)*p/2
m10expected <- n1*(N-n1)*p
dyadicity       <- m11observed / m11expected
heterophilicity <- m10observed / m10expected
c(heterophilicity)}

# \name{dyadicity_p}
# \description{calculates dyadicity with a p-value for a subset of nodes in a larger network}
# \usage{dyadicity_p(graphNEL hcsm1, character testset)}
# \examples{heterophilicity(hcsm1, nodes(hcsm1))}

dyadicity_p <- function (testset, hcsm, R = 100) {
size <- length(nodes(hcsm)); 
testset_is <- intersect(nodes(hcsm), testset)
xx <- length(testset_is)
testset_dyadicity <- dyadicity(hcsm, testset) 
mean_dyadicity  <- numeric(R);
for (i in 1:R) {
mean_dyadicity[i] <- dyadicity(hcsm, nodes(hcsm)[sample(1:size, xx)])}
p  <- length(mean_dyadicity[mean_dyadicity < testset_dyadicity])/R
c(testset_dyadicity, mean(mean_dyadicity), p)}

# \name{youngestFromduplicationWave}
# \description{given duplicationWave data.frame and orderedTaxa, sorts according to taxa age and returns 1st row}
# \usage{youngestFromduplicationWave(data.frame duplicationWave, characterVector orderedTaxa)}
# \examples{youngestFromduplicationWave(duplicationWave, orderedTaxa)}

youngestFromduplicationWave <- function(duplicationWave, orderedTaxa) {
	errMsg= "(data.frame) duplicationWave, (character) orderedTaxa"
        if(!(class(duplicationWave)=="data.frame")) {
        return(errMsg)
        }
        if(!(class(orderedTaxa)=="character")) {
        return(errMsg)
        }
duplicationWave$taxon <- factor(duplicationWave$taxon, levels=orderedTaxa)
duplicationToken <- duplicationWave[order(duplicationWave$taxon), ][1,]
return(duplicationToken)
}

# \name{annotateAttributeDuplication}
# \description{given network and duplicationWave, annotates network nodes with taxon of youngest duplication}
# \usage{Usage: (graphNEL) network, (data.frame) duplicationWave}

# Example annotateAttributeDuplication(hcsm, duplicationWave, orderedTaxa)

annotateAttributeDuplication <- function(network, duplicationWave, orderedTaxa) {
	errMsg="Usage: (graphNEL) network, (data.frame) duplicationWave, (character) orderedTaxa"
	if(!(class(network)[1]=="graphNEL" && class(duplicationWave)[1]=="data.frame")) {
	return(errMsg)}

	networkWithDuplicationAge <- initNodeAttribute (graph=network, attribute.name='duplicationAge',attribute.type='char',default.value='undefined')
	ids <- as.character(duplicationWave$"id")
	idsNode <- nodes(network)
	for (node in nodes(network)) {
		index <- which(as.character(duplicationWave$"id") %in% node)
		duplicationWave2 <- duplicationWave[index,]
		duplicationToken <- youngestFromduplicationWave(duplicationWave2, orderedTaxa) # only return youngest duplication
		nodeData(networkWithDuplicationAge, node, "duplicationAge") <- duplicationToken$taxon
		}
	return(networkWithDuplicationAge)
}

#
# 
#

annotateAttributeMutation <- function(network, CancerMutationGeenList) {
        errMsg="Usage: (graphNEL) network, (data.frame) CancerMutationGeenList"
        if(!(class(network)[1]=="graphNEL" && class(CancerMutationGeenList)[1]=="data.frame")) {
        return(errMsg)}

        networkWithMutations <- initNodeAttribute (graph=network, attribute.name="CancerMutationGeenList",attribute.type="integer",default.value='undefined')
        ids <- as.character(CancerMutationGeenList$"id")
        idsNode <- nodes(network)
        for (node in nodes(network)) {
                index <- which(as.character(CancerMutationGeenList$"id") %in% node)
                CancerMutationGeenList2 <- CancerMutationGeenList[index,]
		CancerMutationGeenList3 <- CancerMutationGeenList2[1,] # take only first row
                nodeData(networkWithMutations, node, "CancerMutationGeenList") <- CancerMutationGeenList$score # score not taxon!
                }
        return(networkWithMutations)
}

#
# CancerDrugTargets lists which drugs target which proteins (network nodes)
#

annotateAttributeDrugTarget <- function(network, CancerDrugTargets) {
        errMsg="Usage: (graphNEL) network, (data.frame) CancerDrugTargets"
        if(!(class(network)[1]=="graphNEL" && class(CancerDrugTargets)[1]=="data.frame")) {
        return(errMsg)}

        networkWithTargets <- initNodeAttribute (graph=network, attribute.name="CancerDrugTargets",attribute.type='char',default.value='undefined')
        ids <- as.character(CancerDrugTargets$"id")
        idsNode <- nodes(network)
        for (node in nodes(network)) {
                index <- which(as.character(CancerDrugTargets$"id") %in% node)
                CancerDrugTargets2 <- CancerDrugTargets[index,] # TODO: not sure how it copes with multiple drugTargeting
                nodeData(networkWithTargets, node, "CancerDrugTargets") <- duplicationToken$taxon
                }
        return(networkWithTargets)
}

#
# duplicationVector is the 6-column data.frame of TreeFam gene duplication data 
# [1] "family"      "familySide"  "gene"        "node"        "primary_acc"       "taxon"
# duplicationWave needs only two columns: "taxon", "primary_acc"
#

# \name{duplicationVector2duplicationWave}
# \description{converts 6-column duplicationVector to 2-column duplicationWave}
# \usage{Usage: (graphNEL) network, (data.frame) duplicationWave}

duplicationVector2duplicationWave <- function(duplicationVector) {
	# API: (data.frame) duplicationVector; check arg type
        if(!(class(duplicationVector)=="data.frame")) {
        return("duplicationVector needs to be a data.frame of 6 columns")
        }
        if(!(length(duplicationVector)==6)) {
        return("duplicationVector needs to be a data.frame of 6 columns")
        }
duplicationWave <-  unique(duplicationVector[,c(2,6)])
colnames(duplicationWave) <- c("taxon", "id")
return(duplicationWave)
}

#
# converts vector of HUGOid to a vector of EGids preserving order
#

convertHUGO2EG <- function (HUGOid) {
EGid <- NULL
for (i in 1:length(HUGOid))
	{
	aHUGOid <- HUGOid[i]
	aEGid <- as.data.frame(org.Hs.egSYMBOL2EG[aHUGOid])$gene_id # will ignore multiple mapping?
	EGid[i] <- aEGid
	}
return(EGid)
}

#
#
#

mergeIGraphList <- function(IGraphList) {  # merge all graphs in a list into 1
	mergeIGraphListAll <- IGraphList[[1]]  # not sure if this works yet
	for (i in length(IGraphList)) {
	mergeIGraphListAll = mergeIGraphListAll + IGraphList[[1]]
	}
	mergeIGraphListAll
}

#
# metacompare is a wrapper function for a number of vector comparison methods
# metacompare(c("af", "sf", "sf"), c("sf", "ff"))
# a and b should be vectors, although it seems to "work" for data.frames and matrices as well

metacompare <- function (a,b) {
print(union(a,b))
print(intersect(a,b))
print(length(union(a,b)))
print(length(intersect(a,b)))
print(setdiff(a,b))
print(setdiff(b,a))
print(setequal(a,b))
}

#
# Uses eapply (environment apply: apply to each obj in env)
#

des <- function (env) {
c(eapply(env, dim), eapply(env, ls))}
dem <- function (df) {
df[1:10,]
}

#
# Example: subsetBetween(nodes(hcsm)[1:10], hcsm, R=100)
# Returns: character vector with following slots: testSetbetweennessCentrality, meanRandomizedBetweennessCentrality, p-value 
#

subsetBetween <- function (testset, hcsm, R=100) {
size <- length(nodes(hcsm))
between <- brandes.betweenness.centrality(hcsm)
testset_is <- intersect(nodes(hcsm), testset)
testset_all <- is.element(nodes(hcsm), testset_is)
testset_between <- sum(between$relative.betweenness.centrality.vertices[testset_all])
rand_between  <- numeric(R)
for (i in 1:R) {rand_between[i] <-
sum(between$relative.betweenness.centrality.vertices[sample(1:size,length(testset_is))])}
c(testset_between, mean(rand_between), length(rand_between[rand_between <
testset_between])/R)}

#
# Example: subsetBetweenTop(nodes(hcsm)[1:10], hcsm, R=100)
#

subsetBetweenTop <- function (testset, hcsm, R=100) {
size <- length(nodes(hcsm))
between <- brandes.betweenness.centrality(hcsm)
testset_is <- intersect(nodes(hcsm), testset)
testset_all <- is.element(nodes(hcsm), testset_is)
testset_between <- between$relative.betweenness.centrality.vertices[testset_all]
testset_between_top <-
sum(sort(testset_between)[length(testset_is):(length(testset_is)-length(testset_is)/5)])
rand_between_top  <- numeric(R);
for (i in 1:R) {rand_between <-
between$relative.betweenness.centrality.vertices[sample(1:size,length(testset_is))]
rand_between_top[i] <-
sum(sort(rand_between)[length(testset_is):(length(testset_is)-length(testset_is)/5)])}
c(length(rand_between_top[rand_between_top <
testset_between_top])/R, testset_between_top, mean(rand_between_top))}

#
# Example: subsetBetween001(nodes(hcsm)[1:10], hcsm, R=100)
#

subsetBetween001 <- function (testset, hcsm, R=100) {
size <- length(nodes(hcsm))
between <- brandes.betweenness.centrality(hcsm)
testset_is <- intersect(nodes(hcsm), testset)
testset_all <- is.element(nodes(hcsm), testset_is)
testset_between <- between$relative.betweenness.centrality.vertices[testset_all]
testset_between_top <- length(testset_between[testset_between>0.01])
rand_between  <- numeric(R)
for (i in 1:R) {rand_between <-
between$relative.betweenness.centrality.vertices[sample(1:size,length(testset_is))]
rand_between_top[i] <- length(rand_between[rand_between>0.01])}
c(length(rand_between_top[rand_between_top <
testset_between_top])/R, testset_between_top, mean(rand_between_top))}

#
# Example: subsetGraphBetween(nodes(hcsm)[1:10], hcsm, R=100)
#

subsetGraphBetween <- function (testset, hcsm, R=100) {
size <- length(nodes(hcsm))
uhcsm   <- ugraph(hcsm); testset_is <- intersect(nodes(uhcsm), testset)
test_subGraph <- subGraph(as.character(testset_is), uhcsm)
between <- brandes.betweenness.centrality(test_subGraph)
testset_between <- sum(between$betweenness.centrality.vertices)
rand_between  <- numeric(R)
for (i in 1:R) {rand_subGraph <- subGraph(as.character(sample(1:size,
length(testset_is))), uhcsm)
rand_between[i] <-
sum(brandes.betweenness.centrality(rand_subGraph)$betweenness.centrality.vertices)}
c(length(rand_between[rand_between <
testset_between])/R, testset_between, mean(rand_between, na.rm = T))}

#
# Example: bringEdges(hcsm, nodes(hcsm)[4])   # bring edges for 1 node
#          bringEdges(hcsm, nodes(hcsm)[4:6]) # bring edges for 3 consecutive nodes

bringEdges <- function (hcsm, node) {
index <- is.element(nodes(hcsm), node)
edges <- edges(hcsm)[index]
edges}

#
#
#

bringUndirectedEdges <- function (hcsm, node) {
uhcsm <- ugraph(hcsm)
index <- is.element(nodes(uhcsm), node)
edges <- as.vector(unlist(edges(uhcsm)[index]))
edges}

#
#
#

bringOutEdges <- function (hcsm, node) {
index <- is.element(nodes(hcsm), node)
edges <- as.vector(unlist(edges(hcsm)[index]))
edges}

#
#
#

bringInEdges <- function (hcsm, node) {
uhcsm <- ugraph(hcsm)
index <- is.element(nodes(hcsm), node)
edges1 <- as.vector(unlist(edges(uhcsm)[index]))
edges2 <- as.vector(unlist(edges(hcsm)[index]))
edges  <- setdiff(edges1, edges2)
edges}

#
#
#

bringCommonUndirectedEdges <- function (hcsm, node1, node2) {
edges1 <- bringUndirectedEdges(hcsm, node1)
edges2 <- bringUndirectedEdges(hcsm, node2)
edges <- intersect(edges1, edges2)
edges}

#
#
#

bringCommonOutEdges <- function (hcsm, node1, node2) {
edges1 <- bringOutEdges(hcsm, node1)
edges2 <- bringOutEdges(hcsm, node2)
edges <- intersect(edges1, edges2)
edges}

#
#
#

bringCommonInEdges <- function (hcsm, node1, node2) {
edges1 <- bringInEdges(hcsm, node1)
edges2 <- bringInEdges(hcsm, node2)
edges <- intersect(edges1, edges2)
edges}

#
#
#

bringCommonScrambledEdges <- function (hcsm, node1, node2) {
edges1 <- bringInEdges(hcsm, node1)
edges2 <- bringOutEdges(hcsm, node1)
edges3 <- bringInEdges(hcsm, node2)
edges4 <- bringOutEdges(hcsm, node2)
edges <- union(intersect(edges1, edges4), intersect(edges2, edges3))
edges}

#
#
#

bringAllUndirectedEdges <- function (hcsm, node1, node2) {
edges1 <- bringUndirectedEdges(hcsm, node1)
edges2 <- bringUndirectedEdges(hcsm, node2)
edges <- union(edges1, edges2)
edges}

#
#
#

subsetGraphCC <- function (testset, hcsm, R=100) {
size <- length(nodes(hcsm))
testset_is <- intersect(nodes(hcsm), testset)
test_subGraph <- subGraph(as.character(testset_is), ugraph(hcsm))
cc <- connectedComp(test_subGraph)
test_subGraph_max_cc <- length(cc$'1')
max_cc <- numeric(R)
t1 <- new("GOHyperGParams", geneIds = cc$'1', annotation =
"hgu95av2.db", ontology="MF", pvalueCutoff = 0.0001, conditional = TRUE,
testDirection = "over")
over <- summary(hyperGTest(t1))
for (i in 1:R) {
rand_subGraph <- subGraph(nodes(hcsm)[sample(1:size,length(testset_is))],ugraph(hcsm));
rand_cc <- connectedComp(rand_subGraph);
max_cc[i] <- length(rand_cc$'1'); }
p <- length(max_cc[max_cc < test_subGraph_max_cc])/R
c(p, test_subGraph_max_cc, mean(max_cc))}

#
#
#

subsetGraphDegree <- function (testset, hcsm, R=100) {
size <- length(nodes(hcsm))
testset_is <- intersect(nodes(hcsm), testset)
test_subGraph <- subGraph(as.character(testset_is), ugraph(hcsm))
test_subGraph_degree_mean <- mean(degree(test_subGraph))
mean_degree <- numeric(R)
for (i in 1:R) {
rand_subGraph <- subGraph(nodes(hcsm)[sample(1:size, length(testset_is))],
ugraph(hcsm)); mean_Degree[i] <- mean(degree(rand_subGraph))}
p <- length(mean_Degree[mean_Degree < test_subGraph_degree_mean])/R
c(p, test_subGraph_degree_mean, mean(mean_Degree))}

#
#
#

subsetDegree <- function (testset, hcsm, R=100) {

if(class(hcsm)[1] != "graphNEL") {return("Usage subsetDegree(testset, graphNEL, repeat)")}

size <- length(nodes(hcsm)); testset_is <- intersect(nodes(hcsm), testset)
uhcsm   <- ugraph(hcsm); xx <- length(testset_is)
testset_out <- is.element(names(degree((hcsm))$outDegree), testset_is)
testset_in  <- is.element(names(degree((hcsm))$inDegree), testset_is)
testset_all <- is.element(names(degree(uhcsm)), testset_is)
testset_out_degree <- degree((hcsm))$outDegree[testset_out]
testset_in_degree  <- degree((hcsm))$inDegree[testset_in]
testset_all_degree <- degree(uhcsm)[testset_all]
testset_out_degree_mean <- mean(testset_out_degree)
testset_in_degree_mean  <- mean(testset_in_degree)
testset_degree_mean <- mean(testset_all_degree)
mean_outDegree  <- numeric(R); mean_inDegree   <- numeric(R)
mean_Degree     <- numeric(R); for (i in 1:R) {
mean_outDegree[i] <- mean(degree((hcsm))$outDegree[sample(1:size, xx)])
mean_inDegree[i] <- mean(degree((hcsm))$inDegree[sample(1:size, xx)])
mean_Degree[i] <- mean(degree(uhcsm)[sample(1:size, xx)])
p_out  <- length(mean_outDegree[mean_outDegree < testset_out_degree_mean])/R
p_in   <- length(mean_inDegree[mean_inDegree < testset_in_degree_mean])/R
p      <- length(mean_Degree[mean_Degree < testset_degree_mean])/R}
c(testset_out_degree_mean, mean(mean_outDegree), p_out, testset_in_degree_mean,
mean(mean_inDegree), p_in, testset_degree_mean, mean(mean_Degree), p)}

#
#
#

subsetGraph <- function (testset, hcsm) {
size <- length(nodes(hcsm))
testset_is <- intersect(nodes(hcsm), testset)
test_subGraph <- subGraph(as.character(testset_is), ugraph(hcsm))
test_subGraph}

#
#
#

bringDataSubTaxon <- function(geneset, data) {
set <- t(geneset)
nn <- bringRefSeq(set)
mm <- rownames(data)
index <- is.element(mm, nn)
index2 = 1:length(index)
positions <- index2[index]
datasubset <- data[positions,]
datasubset
}

#
#
#

bringRefSeq <- function(geneset) {
  nn <- NULL
  for (i in 1:length(geneset))  {
  n <- as.character(geneset[i])
  nn <- c(nn, org.Hs.egREFSEQ[[n]])
  }
  nn
}

#
#
#

bringExpressionVector <- function(refseq, data)  {
expressionVector <- NA
if (refseq %in% (rownames(data)))
  expressionVector <- as.vector(unlist(data[refseq,]))
  expressionVector
} 

#
#
#

subsetTransl <- function (testset, transl, R=10000) {
testset_is <- intersect(transl[,1], testset)
testset_all <- is.element(transl[,1], testset_is)
testset_tAI <- mean(transl[,3][testset_is], na.rm = TRUE)
testset_FE  <- mean(transl[,4][testset_is], na.rm = TRUE)
testset_ATG <- mean(transl[,5][testset_is], na.rm = TRUE)
size <- length(testset_is)
rand_tAI  <- numeric(R)
rand_FE   <- numeric(R)
rand_ATG  <- numeric(R)
for (i in 1:R) {rand_is <- sample(transl[,1], size)
rand_tAI[i] <- mean(transl[,3][rand_is], na.rm = TRUE)
rand_FE[i]  <- mean(transl[,4][rand_is], na.rm = TRUE)
rand_ATG[i] <- mean(transl[,5][rand_is], na.rm = TRUE)}
c(testset_tAI, mean(rand_tAI), length(rand_tAI[rand_tAI < testset_tAI])/R,
   testset_FE, mean(rand_FE), length(rand_FE[rand_FE < testset_FE])/R,
testset_ATG, mean(rand_ATG), length(rand_ATG[rand_ATG < testset_ATG])/R)
}

#
#
#

pairsTransl <- function (pairs, transl) {
limit <- length(pairs[,1])
tAI_A <- NULL; FE_A <- NULL; ATG_A <- NULL; tAI_B <- NULL; FE_B <- NULL; ATG_B <- NULL;
for (i in 1:limit) {
  if (length(transl[,1][is.element(transl[,1], pairs[i,1])]) == 0) next
  if (length(transl[,1][is.element(transl[,1], pairs[i,2])]) == 0) next
  tAI_A[i] <- transl[,3][is.element(transl[,1], pairs[i,1])]
  FE_A[i]  <- transl[,4][is.element(transl[,1], pairs[i,1])]
  ATG_A[i] <- transl[,5][is.element(transl[,1], pairs[i,1])]
  tAI_B[i] <- transl[,3][is.element(transl[,1], pairs[i,2])]
  FE_B[i]  <- transl[,4][is.element(transl[,1], pairs[i,2])]
  ATG_B[i] <- transl[,5][is.element(transl[,1], pairs[i,2])]
  }
  cbind(tAI_A, FE_A, ATG_A, tAI_B, FE_B, ATG_B)
}

#
# 
#

bind <- function(x) {
load(x)
attach(env_data)
ls(env_data)
}

#
#
#

check_GOBP_over <- function (question, U, p) {
answer <- new("GOHyperGParams", geneIds = question, annotation = "hgu95av2.db", ontology="BP", 
pvalueCutoff = p, conditional = TRUE, testDirection = "over", universeGeneIds = U)
a2  <- hyperGTest(answer)}
check_GOBP_under <- function (question, U, p) {
answer <- new("GOHyperGParams", geneIds = question, annotation = "hgu95av2.db", ontology="BP", 
pvalueCutoff = p, conditional = TRUE, testDirection = "under", universeGeneIds = U)
a2  <- hyperGTest(answer)}

#
#
#

check_GOMF_over <- function (question, U, p) {
answer <- new("GOHyperGParams", geneIds = question, annotation = "hgu95av2.db", ontology="MF", 
pvalueCutoff = p, conditional = TRUE, testDirection = "over", universeGeneIds = U)
a2  <- hyperGTest(answer)}
check_GOMF_under <- function (question, U, p) {
answer <- new("GOHyperGParams", geneIds = question, annotation = "hgu95av2.db", ontology="MF", 
pvalueCutoff = p, conditional = TRUE, testDirection = "under", universeGeneIds = U)
a2  <- hyperGTest(answer)}

#
#
#

check_GOCC_over <- function (question, U, p) {
answer <- new("GOHyperGParams", geneIds = question, annotation = "hgu95av2.db", ontology="CC", 
pvalueCutoff = p, conditional = TRUE, testDirection = "over") 
a2  <- hyperGTest(answer)}
check_GOCC_under <- function (question, U, p) {
answer <- new("GOHyperGParams", geneIds = question, annotation = "hgu95av2.db", ontology="CC", 
pvalueCutoff = p, conditional = TRUE, testDirection = "under")
a2  <- hyperGTest(answer)}

#
#
#

check_KEGG_over <- function (question, U, p) {
answer <- new("KEGGHyperGParams", geneIds = question, annotation = "hgu95av2.db", pvalueCutoff = p, 
 testDirection = "over")
a2  <- hyperGTest(answer)}
check_KEGG_under <- function (question, U, p) {
answer <- new("GOHyperGParams", geneIds = question, annotation = "hgu95av2.db", pvalueCutoff = p, 
 testDirection = "under")
a2  <- hyperGTest(answer)}

#
#
#

check_PFAM_over <- function (question, U, p) {
answer <- new("PFAMHyperGParams", geneIds = question, annotation = "hgu95av2.db", pvalueCutoff = p, 
 testDirection = "over")
a2  <- hyperGTest(answer)}
check_PFAM_under <- function (question, U, p) {
answer <- new("PFAMHyperGParams", geneIds = question, annotation = "hgu95av2.db", pvalueCutoff = p, 
 testDirection = "under")
a2  <- hyperGTest(answer)}

#
# bringLocationEG
#

bringLocationEG <- function(id) {
gene <- transcripts(txdb, list(gene_id = id))
gene <- reduce(gene)
chr <- as.character(seqnames(gene))
pos <- start(gene)
str <- as.character(strand(gene))
return(c(chr, pos, str))
}

#
# bringLocationRF
#

bringLocationRF <- function(id) {
refseq <- refseq_TTSes[refseq_TTSes$refseq==id,]
chr <- refseq[2]
str <- refseq[3]
pos <- refseq[4]
return(c(chr, pos, str))
}

#
# jaccardIndexEG2("891", "9133")
#

jaccardIndexEG2 <- function(id1, id2, range=500, ignoreStrand=FALSE) {
res1 <- bringLocationEG(id1)
chr1 <- res1[1]
pos1 <- as.numeric(res1[2]) 
str1 <- res1[3]
res2 <- bringLocationEG(id2)
chr2 <- res2[1]
pos2 <- as.numeric(res2[2])
str2 <- res2[3]
res <- list()
res <- jaccardIndex2(chr1, pos1, str1, chr2, pos2, str2, range, ignoreStrand)
return(res)
}

#
# jaccardIndexRF2("NM_032291", "NM_032291")
#

jaccardIndexRF2 <- function(id1, id2, range=500, ignoreStrand=FALSE) {
res1 <- bringLocationRF(id1)
chr1 <- res1[0]
pos1 <- as.numeric(res1[1]) 
str1 <- res1[2]
res2 <- bringLocationRF(id2)
chr2 <- res2[0]
pos2 <- as.numeric(res2[1])
str2 <- res2[2]
res <- list()
res <- jaccardIndex2(chr1, pos1, str1, chr2, pos2, str2, range, ignoreStrand)
return(res)
}

#
# jaccardIndex2
#

jaccardIndex2 <- function(chr1, pos1, str1, chr2, pos2, str2, range, ignoreStrand) {
pos1FetchTfbs <- fetchTfbs(chr1, pos1, str1, range, ignoreStrand)
pos2FetchTfbs <- fetchTfbs(chr2, pos2, str2, range, ignoreStrand)
u <- union(pos1FetchTfbs, pos2FetchTfbs)
i <- intersect(pos1FetchTfbs, pos2FetchTfbs)
if(length(u)==0) {Jaccard <- NaN}
else {Jaccard <- length(i)/length(u)}
res <- list()
res[["Jaccard"]] <- Jaccard
res[["union"]] <- u
res[["intersection"]] <- i
res[["id1FetchTfbs"]] <- pos1FetchTfbs
res[["id2FetchTfbs"]] <- pos2FetchTfbs
return(res)
}

#
# fetchTfbs("chr1", 23232, "+")
#

fetchTfbs <- function(chr, pos, str, range=500, ignoreStrand=FALSE) {
query <- GRanges(seqnames=chr,ranges=IRanges(start=as.numeric(pos)-range,end=as.numeric(pos)+range),strand=str)
overlaps <- findOverlaps(query, Tfbs, ignore.strand = ignoreStrand)
sh <- subjectHits(overlaps)
x <- Tfbs[sh,]
yes <- as.data.frame(x[,1])
detectedTfbs <- as.vector(yes$name)
return(detectedTfbs)
}

#
# fetch_ENCODE_Tfbs
#

fetch_ENCODE_Tfbs <- function(query, Jaccard_window = 500) {
print(query)
gene <- transcripts(txdb, list(gene_id = query))
gene <- reduce(gene)  # Emergency hack! This function can deal with only one promoter so alt. promoters merged into earliest TSS!
gene_promoter <- flank(gene, Jaccard_window, start = TRUE, both = TRUE, use.names = TRUE, ignore.strand=FALSE)
genes_Tfbs <- findOverlaps(gene_promoter, Tfbs)
sh <- subjectHits(genes_Tfbs)
x <- Tfbs[sh,]
yes <- as.data.frame(x[,1])
final <- as.vector(yes$name)
chr <- as.vector(seqnames(gene_promoter))
location <- paste0("Promoter=", chr, sep=":", start(gene_promoter), sep="-", end(gene_promoter))
final_elements <- paste0(final, sep=", ",collapse="")
return <- paste0(location, sep=", includes: ", final_elements)
return
}

#
# fetch_ENCODE_Tfbs_750
#

fetch_ENCODE_Tfbs_750 <- function(query, Jaccard_window = 500) {
print(query)
gene <- transcripts(txdb, list(gene_id = query))
gene <- reduce(gene)  # Emergency hack! This function can deal with only one promoter so alt. promoters merged into earliest TSS!
gene_promoter <- flank(gene, Jaccard_window, start = TRUE, both = TRUE, use.names = TRUE, ignore.strand=FALSE)
genes_Tfbs <- findOverlaps(gene_promoter, Tfbs)
sh <- subjectHits(genes_Tfbs)
x <- Tfbs[sh,]
yes <- as.data.frame(x[,1])
final <- as.vector(yes$name)
width_score <- as.data.frame(x[,2])[,4:5]
width_score_750 <- width_score[width_score$score > 750,]
more_than_750 <- rownames(width_score_750)
more_than_750_index <- 1:length(final) %in% more_than_750
final_more_than_750 <- final[more_than_750_index]
final <- final_more_than_750
chr <- as.vector(seqnames(gene_promoter))
location <- paste0("Promoter=", chr, sep=":", start(gene_promoter), sep="-", end(gene_promoter))
final_elements <- paste0(final, sep=", ",collapse="")
return <- paste0(location, sep=", includes (at 0.75 quantile cut-off): ", final_elements)
return
}

#
# Jaccard_ENCODE_Tfbs("891", "9133")
#

Jaccard_ENCODE_Tfbs <- function(query1, query2, Jaccard_window = 1000) {
if(!is.character(query1)) {print("Usage: (character) query1, (character) query2, (int) Jaccard_window")}
if(!is.character(query1)) {print("Usage: (character) query1, (character) query2, (int) Jaccard_window")}
if(!is.numeric(Jaccard_window)) {print("Usage: (character) query1, (character) query2, (int) Jaccard_window")}

#print(paste0("Calculating Jaccard index for pair: ", query1, sep = " ", query2))
gene <- transcripts(txdb, list(gene_id = query1))
print("done gene <- transcripts(txdb, list(gene_id = query1))")
gene <- reduce(gene)  # Emergency hack! This function can deal with only one promoter so alt. promoters merged into earliest TSS!
gene_promoter <- NULL
try(gene_promoter <- flank(gene, Jaccard_window, start = TRUE, both = TRUE, use.names = TRUE, ignore.strand=FALSE))
if(is.null(gene_promoter)) { return(NaN) }
print("done: if(is.null(gene_promoter)) { return(NaN) }")
genes_Tfbs <- findOverlaps(gene_promoter, Tfbs)
sh <- subjectHits(genes_Tfbs)
x <- Tfbs[sh,]
yes <- as.data.frame(x[,1])
final <- as.vector(yes$name)
print("done: final <- as.vector(yes$name)")
chr <- as.vector(seqnames(gene_promoter))
#location <- paste0("Promoter=", chr, sep=":", start(gene_promoter), sep="-", end(gene_promoter))
final_elements <- paste0(final, sep=", ",collapse="")
Tfbs_x <- final
print(paste0("Done: ", query1))

gene <- transcripts(txdb, list(gene_id = query2))
gene <- reduce(gene)  # Emergency hack! This function can deal with only one promoter so alt. promoters merged into earliest TSS!
gene_promoter <- NULL
try(gene_promoter <- flank(gene, Jaccard_window, start = TRUE, both = TRUE, use.names = TRUE, ignore.strand=FALSE))
if(is.null(gene_promoter)) { return(NaN) }
genes_Tfbs <- findOverlaps(gene_promoter, Tfbs)
sh <- subjectHits(genes_Tfbs)
x <- Tfbs[sh,]
yes <- as.data.frame(x[,1])
final <- as.vector(yes$name)
chr <- as.vector(seqnames(gene_promoter))
#location <- paste0("Promoter=", chr, sep=":", start(gene_promoter), sep="-", end(gene_promoter))
final_elements <- paste0(final, sep=", ",collapse="")
Tfbs_y <- final
print(paste0("Done: ", query2))
u <- union(Tfbs_x, Tfbs_y)
i <- intersect(Tfbs_x, Tfbs_y)
Jaccard <- length(i)/length(u)
res <- list()
res[["Jaccard"]] <- Jaccard
res[["union"]] <- u
res[["intersection"]] <- i
res[["Tfbs_x"]] <- Tfbs_x
res[["Tfbs_y"]] <- Tfbs_y
res[["returnsFromat"]] <- "Jaccard, union, intersection, Tfbs_x, Tfbs_y"
return(res)
}

#
# Jaccard_ENCODE_Tfbs("891", "9133") should give the same result as Jaccard_ENCODE_Tfbs_bed("NM_031966", "NM_004701", 250)
#

load("~/nb/DATA/env_promoter_bed_500_res_env")
resJaccard_ENCODE_Tfbs_substrate_500 <- env_bed_jaccard$resJaccard_ENCODE_Tfbs_substrate
hg19_refseq <- read.table("data/expression/hg19_refseq.bed")

#
# function(query1, query2, Jaccard_window = 500)
# Jaccard_ENCODE_Tfbs_bed()
#
# Jaccard_ENCODE_Tfbs_bed("NM_198471", "NM_015158", 500)
#

Jaccard_ENCODE_Tfbs_bed <- function(query1, query2, Jaccard_window = 500, taxon) {
if(!is.character(query1)) 
{print("Usage: (character) query1, (character) query2, (int) Jaccard_window"); return()}
if(!is.character(query2)) 
{print("Usage: (character) query1, (character) query2, (int) Jaccard_window"); return()}
if(!is.numeric(Jaccard_window)) 
{print("Usage: (character) query1, (character) query2, (int) Jaccard_window"); return()}
if(is.na(query1)) {print("Usage: (character) query1, (character) query2, (int) Jaccard_window"); return()}
if(is.na(query2)) {print("Usage: (character) query1, (character) query2, (int) Jaccard_window"); return()}

#print(paste0("Calculating Jaccard index for pair: ", query1, sep = " ", query2))

if(Jaccard_window == 250) {
Tfbs_x <- resJaccard_ENCODE_Tfbs_substrate_250[[query1]]
Tfbs_y <- resJaccard_ENCODE_Tfbs_substrate_250[[query2]]
} else if(Jaccard_window == 500) {
Tfbs_x <- resJaccard_ENCODE_Tfbs_substrate_500[[query1]]
Tfbs_y <- resJaccard_ENCODE_Tfbs_substrate_500[[query2]]
} else if(Jaccard_window == 1000) {
Tfbs_x <- resJaccard_ENCODE_Tfbs_substrate_1000[[query1]]
Tfbs_y <- resJaccard_ENCODE_Tfbs_substrate_1000[[query2]]
} else if(Jaccard_window == 1500) {
Tfbs_x <- resJaccard_ENCODE_Tfbs_substrate_1500[[query1]]
Tfbs_y <- resJaccard_ENCODE_Tfbs_substrate_1500[[query2]]
} else if(Jaccard_window == 2000) {
Tfbs_x <- resJaccard_ENCODE_Tfbs_substrate_2000[[query1]]
Tfbs_y <- resJaccard_ENCODE_Tfbs_substrate_2000[[query2]]
}

u <- union(Tfbs_x, Tfbs_y)
i <- intersect(Tfbs_x, Tfbs_y) 
Jaccard <- length(i)/length(u)
res <- list()

query1_hg19_refseq <- hg19_refseq[hg19_refseq$V4 == query1,]
query2_hg19_refseq <- hg19_refseq[hg19_refseq$V4 == query2,]

if(dim(query1_hg19_refseq)[1]==0) {return()}
if(dim(query2_hg19_refseq)[1]==0) {return()}

if(query1_hg19_refseq$V6 == "+") {
query1_loc <- paste0(query1_hg19_refseq$V1, ":", query1_hg19_refseq$V2-Jaccard_window, "-",query1_hg19_refseq$V2+Jaccard_window)
}
if(query1_hg19_refseq$V6 == "-") {
query1_loc <- paste0(query1_hg19_refseq$V1, ":", query1_hg19_refseq$V3-Jaccard_window, "-",query1_hg19_refseq$V3+Jaccard_window)
}
if(query2_hg19_refseq$V6 == "+") {
query2_loc <- paste0(query2_hg19_refseq$V1, ":", query2_hg19_refseq$V2-Jaccard_window, "-",query2_hg19_refseq$V2+Jaccard_window)
}
if(query2_hg19_refseq$V6 == "-") {
query2_loc <- paste0(query2_hg19_refseq$V1, ":", query2_hg19_refseq$V3-Jaccard_window, "-",query2_hg19_refseq$V3+Jaccard_window)
}

res[["taxon"]] <- taxon
res[["query1_TSS_window"]] <- query1_loc
res[["query2_TSS_window"]] <- query2_loc
res[["query1"]] <- query1
res[["query2"]] <- query2
res[["Jaccard"]] <- Jaccard
res[["union"]] <- u
res[["intersection"]] <- i
res[["Tfbs_x"]] <- Tfbs_x
res[["Tfbs_y"]] <- Tfbs_y
res[["returnsFromat"]] <- "Jaccard, union, intersection, Tfbs_x, Tfbs_y"
return(res)
}

#
#
#
