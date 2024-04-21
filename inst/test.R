library(DaparToolshed)
library(omXplore)
library(MagellanNTK)
data(ft_na)
obj <- ft_na[[1]]
conds <- get_group(ft_na)
list_tags <- c("None" = "None",metacell.def(typeDataset(ft_na[[1]]))$node)
operator = setNames(nm = SymFilteringOperators())
keep_vs_remove <- setNames(nm = c("delete", "keep"))
val_vs_percent <- setNames(nm = c("Count", "Percentage"))
value = NULL

shiny::runApp(
mod_FunctionFilter_Generator(
obj, conds, list_tags, keep_vs_remove, val_vs_percent, operator))


shiny::runApp(mod_VariableFilter_Generator(obj, keep_vs_remove, value, operator))

shiny::runApp(mod_Metacell_Filtering(ft_na, 1))

shiny::runApp(mod_Variable_Filtering(ft_na, 1))




data(ft_na, package='DaparToolshed')
obj <- ft_na[[1]]
conds <- get_group(ft_na)
list_tags <- c("None" = "None", metacell.def(typeDataset(ft_na[[1]]))$node)
operator = setNames(nm = SymFilteringOperators())
keep_vs_remove <- setNames(nm = c("delete", "keep"))
val_vs_percent <- setNames(nm = c("Count", "Percentage"))

shiny::runApp(mod_build_qMetacell_FunctionFilter(obj, 
  conds, list_tags, keep_vs_remove, val_vs_percent, operator))







data(ft_na)
obj <- ft_na[[1]]
shiny::runApp(
  mod_VariableFilter_Generator(
    obj, conds, list_tags, keep_vs_remove, val_vs_percent, operator))
)

library(omXplore)
library(DaparToolshed)
library(MagellanNTK)
data(ft_na)
path <- system.file('workflow/PipelineProtein', package = 'DaparToolshed')
shiny::runApp(workflowApp("PipelineProtein_Filtering", path, dataIn = ft_na))

