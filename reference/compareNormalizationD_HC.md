# Builds a plot from a dataframe. Same as compareNormalizationD but uses the library `highcharter`

Plot to compare the quantitative proteomics data before and after
normalization using the package `highcharter`

## Usage

``` r
compareNormalizationD_HC(
  qDataBefore,
  qDataAfter,
  keyId = NULL,
  conds = NULL,
  pal = NULL,
  subset.view = NULL,
  n = 100,
  type = "scatter"
)
```

## Arguments

- qDataBefore:

  A dataframe that contains quantitative data before normalization.

- qDataAfter:

  A dataframe that contains quantitative data after normalization.

- keyId:

  xxx

- conds:

  A vector of the conditions (one condition per sample).

- pal:

  xxx

- subset.view:

  xxx

- n:

  An integer that is equal to the maximum number of displayed points.
  This number must be less or equal to the size of the dataset. If it is
  less than it, it is a random selection

- type:

  scatter or line

## Value

A plot

## Author

Samuel Wieczorek

## Examples

``` r
library(SummarizedExperiment)
data(subR25prot)
qDataBefore <- assay(subR25prot[[2]])
conds <- design.qf(subR25prot)$Condition
id <- rowData(subR25prot[[2]])[, idcol(subR25prot[[2]])]
# pal <- ExtendPalette(2)
qDataAfter <- LOESS(qDataBefore, conds, type = "overall")

n <- 1
compareNormalizationD_HC(
qDataBefore = qDataBefore,
qDataAfter = qDataAfter, 
keyId = id, 
pal = NULL, 
n = n,
subset.view = seq_len(n),
conds = conds)
#> Warning: Color palette set to default.

{"x":{"hc_opts":{"chart":{"reflow":true,"type":"scatter","zoomType":"None","showAxes":true,"resetZoomButton":{"position":{"align":"left","verticalAlign":"top"}}},"title":{"text":null},"yAxis":{"title":{"text":null}},"credits":{"enabled":false},"exporting":{"enabled":true,"filename":"compareNormalization","buttons":{"contextButton":{"menuItems":["downloadPNG","downloadSVG","downloadPDF"]}}},"boost":{"enabled":false},"plotOptions":{"series":{"label":{"enabled":false},"turboThreshold":0},"treemap":{"layoutAlgorithm":"squarified"}},"series":[{"name":"Intensity_C_R1","data":[{"x":21.32980947458864,"y":0.993965784906908,"name":"CON__A2I7N1;CON__A2I7N0"}]},{"name":"Intensity_C_R2","data":[{"x":21.29457139231285,"y":0.9925397606508786,"name":"CON__A2I7N1;CON__A2I7N0"}]},{"name":"Intensity_C_R3","data":[{"x":21.36591657001674,"y":0.9902907203565321,"name":"CON__A2I7N1;CON__A2I7N0"}]},{"name":"Intensity_D_R1","data":[{"x":21.31561837611934,"y":0.9983696923938775,"name":"CON__A2I7N1;CON__A2I7N0"}]},{"name":"Intensity_D_R2","data":[{"x":20.89947493961432,"y":1.007609625100951,"name":"CON__A2I7N1;CON__A2I7N0"}]},{"name":"Intensity_D_R3","data":[{"x":21.37283955976277,"y":0.9801693975900653,"name":"CON__A2I7N1;CON__A2I7N0"}]}],"colors":["#E41A1C","#E41A1C","#E41A1C","#377EB8","#377EB8","#377EB8"],"tooltip":{"headerFormat":"","pointFormat":"Id: {point.name}"}},"theme":{"chart":{"backgroundColor":"transparent"},"colors":["#7cb5ec","#434348","#90ed7d","#f7a35c","#8085e9","#f15c80","#e4d354","#2b908f","#f45b5b","#91e8e1"]},"conf_opts":{"global":{"Date":null,"VMLRadialGradientURL":"http =//code.highcharts.com/list(version)/gfx/vml-radial-gradient.png","canvasToolsURL":"http =//code.highcharts.com/list(version)/modules/canvas-tools.js","getTimezoneOffset":null,"timezoneOffset":0,"useUTC":true},"lang":{"contextButtonTitle":"Chart context menu","decimalPoint":".","downloadCSV":"Download CSV","downloadJPEG":"Download JPEG image","downloadPDF":"Download PDF document","downloadPNG":"Download PNG image","downloadSVG":"Download SVG vector image","downloadXLS":"Download XLS","drillUpText":"◁ Back to {series.name}","exitFullscreen":"Exit from full screen","exportData":{"annotationHeader":"Annotations","categoryDatetimeHeader":"DateTime","categoryHeader":"Category"},"hideData":"Hide data table","invalidDate":null,"loading":"Loading...","months":["January","February","March","April","May","June","July","August","September","October","November","December"],"noData":"No data to display","numericSymbolMagnitude":1000,"numericSymbols":["k","M","G","T","P","E"],"printChart":"Print chart","resetZoom":"Reset zoom","resetZoomTitle":"Reset zoom level 1:1","shortMonths":["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"],"shortWeekdays":["Sat","Sun","Mon","Tue","Wed","Thu","Fri"],"thousandsSep":" ","viewData":"View data table","viewFullscreen":"View in full screen","weekdays":["Sunday","Monday","Tuesday","Wednesday","Thursday","Friday","Saturday"]}},"type":"chart","fonts":[],"debug":false},"evals":[],"jsHooks":[]}
```
