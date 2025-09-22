library(QFeatures)

data(subR25pept)
protID <- parentProtId(subR25pept[[2]])
X <- QFeatures::adjacencyMatrix(subR25pept[[2]])

X.split <- DaparToolshed::splitAdjacencyMat(X)
X.shared <- X.split$Xshared
X.unique <- X.split$Xspec


#adjacencyMatrix(subR25pept[[2]]) <- X.unique
#rowdata.pep <- rowData(subR25pept[[2]])


# subR25pept <- aggregateFeatures4Prostar(
#   object = subR25pept,
#   i = length(subR25pept),
#   name = 'aggregated',
#   fcol = 'adjacencyMatrix',
#   fun = 'colSumsMat')
# 
# 
# .names <- "Sequence"
# 
# proteinNames <- rownames(subR25pept[[2]])
# data <- rowData(subR25pept[[1]])
# 
# new.col <- BuildColumnToProteinDataset(
#   peptideData = rowData(subR25pept[[1]]), 
#   matAdj = adjacencyMatrix(subR25pept[[2]]), 
#   columnName = "Sequence",
#   proteinNames = rownames(subR25pept[[2]]))