# Display a CC

Display a CC

## Usage

``` r
display.CC.visNet(
  g,
  layout = NULL,
  obj = NULL,
  prot.tooltip = NULL,
  pept.tooltip = NULL
)
```

## Arguments

- g:

  A cc (a list)

- layout:

  xxxxx

- obj:

  xxx

- prot.tooltip:

  xxx

- pept.tooltip:

  xxx

## Value

A plot

## Author

Thomas Burger, Samuel Wieczorek

## Examples

``` r
data(subR25pept)
X <- QFeatures::adjacencyMatrix(subR25pept[[1]])
ll <- get.pep.prot.cc(X)
g <- buildGraph(ll[[1]], X)
display.CC.visNet(g)

{"x":{"nodes":{"id":[1,2,3],"group":["shared.peptide","protein","protein"],"label":["AADAIIIK","116","115"],"title":["<p>1<br>Tooltip !<\/p>","<p>2<br>Tooltip !<\/p>","<p>3<br>Tooltip !<\/p>"],"size":[10,20,20]},"edges":{"from":[1,1],"to":[2,3]},"nodesToDataframe":true,"edgesToDataframe":true,"options":{"width":"100%","height":"100%","nodes":{"shape":"dot"},"manipulation":{"enabled":false},"groups":{"spec.peptide":{"color":"#5CA3F7"},"shared.peptide":{"color":"#0EA513"},"useDefaultGroups":true,"protein":{"color":"#ECB57C","shape":"dot"}},"edges":{"width":2,"color":"#A9A9A9"}},"groups":["shared.peptide","protein"],"width":"100%","height":"100%","idselection":{"enabled":false,"style":"width: 150px; height: 26px","useLabels":true,"main":"Select by id"},"byselection":{"enabled":false,"style":"width: 150px; height: 26px","multiple":false,"hideColor":"rgba(200,200,200,0.5)","highlight":false},"main":null,"submain":null,"footer":null,"background":"rgba(0, 0, 0, 0)","highlight":{"enabled":false,"hoverNearest":false,"degree":1,"algorithm":"all","hideColor":"rgba(200,200,200,0.5)","labelOnly":true},"collapse":{"enabled":false,"fit":false,"resetHighlight":true,"clusterOptions":null,"keepCoord":true,"labelSuffix":"(cluster)"}},"evals":[],"jsHooks":[]}
```
