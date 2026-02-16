library(QFeatures)

data(subR25pept)
protID <- parentProtId(subR25pept[[1]])
X <- QFeatures::adjacencyMatrix(subR25pept[[1]])

X.split <- DaparToolshed::splitAdjacencyMat(X)
X.shared <- X.split$Xshared
X.unique <- X.split$Xspec
