data(Exp1_R25_pept, package="DaparToolshedData")
obj <- Exp1_R25_pept[1:10]
protID <- parentProtId(obj[[2]])
X <- adjacencyMatrix(obj[[2]])

X.split <- DaparToolshed::splitAdjacencyMat(X)
X.shared <- X.split$Xshared
X.unique <- X.split$Xspec


#adjacencyMatrix(obj[[2]]) <- X.unique
#rowdata.pep <- rowData(obj[[2]])


# obj <- aggregateFeatures4Prostar(
#   object = obj,
#   i = length(obj),
#   name = 'aggregated',
#   fcol = 'adjacencyMatrix',
#   fun = 'colSumsMat')
# 
# 
# .names <- "Sequence"
# 
# proteinNames <- rownames(obj[[length(obj)]])
# data <- rowData(obj[[length(obj)-1]])
# 
# new.col <- BuildColumnToProteinDataset(
#   peptideData = rowData(obj[[length(obj)-1]]), 
#   matAdj = adjacencyMatrix(obj[[2]]), 
#   columnName = "Sequence",
#   proteinNames = rownames(obj[[length(obj)]]))