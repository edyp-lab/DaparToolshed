# Compare two assays

This plot compares the quantitative proteomics data between two assays.
It can be used for example to compare the effect of the normalization
process.

The comparison is made with the division operator.

## Usage

``` r
plotCompareAssays(
  obj,
  i,
  j,
  info = NULL,
  pal.name = "Set1",
  subset.view = NULL,
  n = 100,
  type = "scatter",
  FUN = NULL
)
```

## Arguments

- obj:

  An instance of the class

- i:

  A numeric matrix containing quantitative data after normalization.

- j:

  A numeric matrix containing quantitative data after normalization

- info:

  xxx

- pal.name:

  xxx

- subset.view:

  xxx

- n:

  xxx

- type:

  The type of plot. Available values are 'scatter' (default) or 'line'

- FUN:

  xxx

## Value

A plot

## Author

Samuel Wieczorek, Enora Fremy

## Examples

``` r
data(subR25prot)
plotCompareAssays(subR25prot, 1, 2, n = 5)

{"x":{"hc_opts":{"chart":{"reflow":true,"type":"scatter","zoomType":"None","showAxes":true,"width":0,"height":0,"resetZoomButton":{"position":{"align":"left","verticalAlign":"top"}}},"title":{"text":null},"yAxis":{"title":{"text":null}},"credits":{"enabled":false},"exporting":{"enabled":false},"boost":{"enabled":false},"plotOptions":{"series":{"label":{"enabled":false},"turboThreshold":0},"treemap":{"layoutAlgorithm":"squarified"}},"series":[{"name":"Intensity_C_R1","data":[{"x":290790000,"y":9.668627675042858e-08,"name":null},{"x":13332000,"y":1.775306772185144e-06,"name":null},{"x":390980000,"y":7.300250542125265e-08,"name":null},{"x":null,"y":null,"name":null},{"x":163940000,"y":1.66454755777164e-07,"name":null}]},{"name":"Intensity_C_R2","data":[{"x":284320000,"y":9.877230003368785e-08,"name":null},{"x":null,"y":null,"name":null},{"x":387480000,"y":7.362843657965466e-08,"name":null},{"x":null,"y":null,"name":null},{"x":125990000,"y":2.135783315237225e-07,"name":null}]},{"name":"Intensity_C_R3","data":[{"x":288030000,"y":9.756498921526693e-08,"name":null},{"x":15418000,"y":1.54871658433968e-06,"name":null},{"x":460460000,"y":6.249945835564968e-08,"name":null},{"x":null,"y":null,"name":null},{"x":137990000,"y":1.959561452220847e-07,"name":null}]},{"name":"Intensity_D_R1","data":[{"x":103900000,"y":2.56310109849288e-07,"name":null},{"x":null,"y":null,"name":null},{"x":163210000,"y":1.671598199505637e-07,"name":null},{"x":2721500,"y":7.8544812071972e-06,"name":null},{"x":44454000,"y":5.715078486469354e-07,"name":null}]},{"name":"Intensity_D_R2","data":[{"x":95622000,"y":2.772462325110676e-07,"name":null},{"x":null,"y":null,"name":null},{"x":193480000,"y":1.422762518192853e-07,"name":null},{"x":null,"y":null,"name":null},{"x":42653000,"y":5.942405899074297e-07,"name":null}]},{"name":"Intensity_D_R3","data":[{"x":97103000,"y":2.732460637946269e-07,"name":null},{"x":null,"y":null,"name":null},{"x":125390000,"y":2.145453948097931e-07,"name":null},{"x":null,"y":null,"name":null},{"x":36597000,"y":6.8653775047633e-07,"name":null}]}],"colors":["#E41A1C","#377EB8"],"tooltip":{"enabled":false}},"theme":{"chart":{"backgroundColor":"transparent"},"colors":["#7cb5ec","#434348","#90ed7d","#f7a35c","#8085e9","#f15c80","#e4d354","#2b908f","#f45b5b","#91e8e1"]},"conf_opts":{"global":{"Date":null,"VMLRadialGradientURL":"http =//code.highcharts.com/list(version)/gfx/vml-radial-gradient.png","canvasToolsURL":"http =//code.highcharts.com/list(version)/modules/canvas-tools.js","getTimezoneOffset":null,"timezoneOffset":0,"useUTC":true},"lang":{"contextButtonTitle":"Chart context menu","decimalPoint":".","downloadCSV":"Download CSV","downloadJPEG":"Download JPEG image","downloadPDF":"Download PDF document","downloadPNG":"Download PNG image","downloadSVG":"Download SVG vector image","downloadXLS":"Download XLS","drillUpText":"◁ Back to {series.name}","exitFullscreen":"Exit from full screen","exportData":{"annotationHeader":"Annotations","categoryDatetimeHeader":"DateTime","categoryHeader":"Category"},"hideData":"Hide data table","invalidDate":null,"loading":"Loading...","months":["January","February","March","April","May","June","July","August","September","October","November","December"],"noData":"No data to display","numericSymbolMagnitude":1000,"numericSymbols":["k","M","G","T","P","E"],"printChart":"Print chart","resetZoom":"Reset zoom","resetZoomTitle":"Reset zoom level 1:1","shortMonths":["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"],"shortWeekdays":["Sat","Sun","Mon","Tue","Wed","Thu","Fri"],"thousandsSep":" ","viewData":"View data table","viewFullscreen":"View in full screen","weekdays":["Sunday","Monday","Tuesday","Wednesday","Thursday","Friday","Saturday"]}},"type":"chart","fonts":[],"debug":false},"evals":[],"jsHooks":[]}
```
