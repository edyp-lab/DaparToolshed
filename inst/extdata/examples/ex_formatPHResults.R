

library(SummarizedExperiment)
data(subR25prot)
obj <- subR25prot

obj <- NAIsZero(obj, 1)
obj <- NAIsZero(obj, 2)
qdata <- SummarizedExperiment::assay(obj[[2]])
conds <- design.qf(obj)$Condition
anova_tests <- apply(qdata, 1, classic1wayAnova, conditions = as.factor(conds))
anova_tests <- t(anova_tests)

names(anova_tests) <- rownames(qdata)
tms <- lapply(
  anova_tests,
  function(x) {
    summary(multcomp::glht(x,
      linfct = multcomp::mcp(conditions = "Tukey")
    ),
      test = multcomp::adjusted("none")
    )
  }
)
res <- formatPHResults(tms)