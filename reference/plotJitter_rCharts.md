# Display a a jitter plot for CC

Display a a jitter plot for CC

## Usage

``` r
plotJitter_rCharts(df, clickFunction = NULL)
```

## Arguments

- df:

  xxxx

- clickFunction:

  xxxx

## Value

A plot

## Author

Thomas Burger, Samuel Wieczorek

## Examples

``` r
data(subR25pept)
X <- BuildAdjacencyMatrix(subR25pept[[1]])
ll <- get.pep.prot.cc(X)[1:4]
n.prot <- unlist(lapply(ll, function(x) {length(x$proteins)}))
n.pept <- unlist(lapply(ll, function(x) {length(x$peptides)}))
df <- tibble::tibble(
x = jitter(n.pept),
y = jitter(n.prot),
index = seq_len(length(ll))
)
plotJitter_rCharts(df)
#> Warning: There is no tooltip in the object.

{"x":{"hc_opts":{"chart":{"reflow":true,"type":"scatter","zoomType":"xy","showAxes":true,"resetZoomButton":{"position":{"align":"left","verticalAlign":"top"}}},"title":{"text":null},"yAxis":{"title":{"text":"Nb of proteins ic CC"}},"credits":{"enabled":false},"exporting":{"enabled":true,"filename":"plotCC","buttons":{"contextButton":{"menuItems":["downloadPNG","downloadSVG","downloadPDF"]}}},"boost":{"enabled":false},"plotOptions":{"series":{"label":{"enabled":false},"turboThreshold":0,"animation":{"duration":100},"cursor":"pointer","point":{"events":{"click":"function(event) \n                {\n                Shiny.onInputChange('eventPointClicked', \n                [this.index]+'_'+ [this.series.name]);\n                }"}}},"treemap":{"layoutAlgorithm":"squarified"}},"series":[{"data":[{"x":0.998103682314977,"y":2.027599583845586,"index":1},{"x":0.9915935055259615,"y":1.872761585377157,"index":2},{"x":1.000320117501542,"y":1.844665393792093,"index":3},{"x":1.00815483382903,"y":0.8807650652714074,"index":4}],"type":"scatter"}],"legend":{"enabled":false},"xAxis":{"title":{"text":"Nb of peptides ic CC"}},"tooltip":{"headerFormat":"","pointFormat":null}},"theme":{"chart":{"backgroundColor":"transparent"},"colors":["#7cb5ec","#434348","#90ed7d","#f7a35c","#8085e9","#f15c80","#e4d354","#2b908f","#f45b5b","#91e8e1"]},"conf_opts":{"global":{"Date":null,"VMLRadialGradientURL":"http =//code.highcharts.com/list(version)/gfx/vml-radial-gradient.png","canvasToolsURL":"http =//code.highcharts.com/list(version)/modules/canvas-tools.js","getTimezoneOffset":null,"timezoneOffset":0,"useUTC":true},"lang":{"contextButtonTitle":"Chart context menu","decimalPoint":".","downloadCSV":"Download CSV","downloadJPEG":"Download JPEG image","downloadPDF":"Download PDF document","downloadPNG":"Download PNG image","downloadSVG":"Download SVG vector image","downloadXLS":"Download XLS","drillUpText":"◁ Back to {series.name}","exitFullscreen":"Exit from full screen","exportData":{"annotationHeader":"Annotations","categoryDatetimeHeader":"DateTime","categoryHeader":"Category"},"hideData":"Hide data table","invalidDate":null,"loading":"Loading...","months":["January","February","March","April","May","June","July","August","September","October","November","December"],"noData":"No data to display","numericSymbolMagnitude":1000,"numericSymbols":["k","M","G","T","P","E"],"printChart":"Print chart","resetZoom":"Reset zoom","resetZoomTitle":"Reset zoom level 1:1","shortMonths":["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"],"shortWeekdays":["Sat","Sun","Mon","Tue","Wed","Thu","Fri"],"thousandsSep":" ","viewData":"View data table","viewFullscreen":"View in full screen","weekdays":["Sunday","Monday","Tuesday","Wednesday","Thursday","Friday","Saturday"]}},"type":"chart","fonts":[],"debug":false},"evals":["hc_opts.plotOptions.series.point.events.click"],"jsHooks":[]}
```
