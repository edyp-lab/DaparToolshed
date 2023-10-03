file <- system.file("extdata", "Exp1_R25_pept.txt", package="DaparToolshedData")
data <- read.table(file, header=TRUE, sep="\t",stringsAsFactors = FALSE)
metadataFile <- system.file("extdata", "samples_Exp1_R25.txt",
package="DaparToolshedData")
metadata <- read.table(metadataFile, header=TRUE, sep="\t", as.is=TRUE,
stringsAsFactors = FALSE)
conds <- metadata$Condition
qdata <- data[,56:61]
df <- data[ , 43:48]
df <- BuildMetacell('maxquant', 'peptide', qdata, conds, df)
df <- BuildMetacell('proline', 'peptide', qdata, conds, df)
