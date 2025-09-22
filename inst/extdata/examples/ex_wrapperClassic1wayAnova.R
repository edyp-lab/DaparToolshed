library(QFeatures)
data(subR25prot)
obj <- subR25prot
filter <- FunctionFilter('qMetacellOnConditions',
cmd = 'delete',
mode = 'AtLeastOneCond',
pattern = c("Missing POV", "Missing MEC"),
conds = design.qf(obj)$Condition,
percent = TRUE,
th = 0.8,
operator = '<')

obj <- filterFeaturesOneSE(obj, filter)

anovatest <- wrapperClassic1wayAnova(obj, 2)




anova_tests <- t(apply(assay(obj[[2]]), 1, classic1wayAnova,
  conditions = as.factor(design.qf(obj)$Condition)
))