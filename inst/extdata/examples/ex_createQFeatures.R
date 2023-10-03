data.file <- system.file("extdata", "Exp1_R25_pept.txt", package = "DaparToolshedData")

data <- read.table(data.file, header = TRUE, sep = "\t", as.is = TRUE, stringsAsFactors = FALSE)

sample.file <- system.file("extdata", "samples_Exp1_R25.txt", package = "DaparToolshedData")

sample <- read.table(sample.file, header = TRUE, sep = "\t", as.is = TRUE, stringsAsFactors = FALSE)


metacell.indices <- 43:48
metacell.names <- colnames(data)[metacell.indices]

qdata.indices <- 56:61
qdata.names <- colnames(data)[qdata.indices]


ft <- createQFeatures(data = data, 
                      sample = sample,
                      indQData = 56:61,
                      keyId = "Sequence",
                      analysis = "test",
                      indexForMetacell = 43:48,
                      typeDataset = "peptide",
                      parentProtId = "Protein_group_IDs",
                      force.na = TRUE,
                      software = "maxquant"
                      )

ft <- createQFeatures(data = data, 
                      sample = sample,
                      indQData = qdata.names,
                      keyId = "Sequence",
                      analysis = "test",
                      indexForMetacell = metacell.names,
                      typeDataset = "peptide",
                      parentProtId = "Protein_group_IDs",
                      force.na = TRUE,
                      software = "maxquant")
