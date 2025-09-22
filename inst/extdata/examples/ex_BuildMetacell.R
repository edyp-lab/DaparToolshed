data(subR25pept)
conds <- design.qf(subR25pept)$Condition
qdata <- SummarizedExperiment::assay(subR25pept[[2]])
metacell1 <- BuildMetacell('maxquant', 'peptide', qdata, conds)
metacell2 <- BuildMetacell('proline', 'peptide', qdata, conds)
