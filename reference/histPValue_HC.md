# Plots a histogram of p-values

Plots a histogram of p-values

## Usage

``` r
histPValue_HC(pval_ll, bins = 80, pi0 = 1)
```

## Arguments

- pval_ll:

  A vector of the p-values.

- bins:

  A integer indicating the number of cells for the histogram

- pi0:

  A `float` between 0 and 1 corresponding to the proportion of true null
  hypotheses.

## Value

A histogram of the p-values with pi0

## Author

Samuel Wieczorek

## Examples

``` r
data(subR25prot)
obj <- subR25prot
# Simulate imputation
obj <- NAIsZero(obj, 1)
obj <- NAIsZero(obj, 2)
allComp <- limmaCompleteTest(SummarizedExperiment::assay(obj[[length(obj)]]), design.qf(obj), comp.type="OnevsOne")
histPValue_HC(allComp$P_Value[1])

{"x":{"hc_opts":{"chart":{"reflow":true,"type":"column"},"title":{"text":"P-value histogram"},"yAxis":{"title":{"text":"Density"},"plotLines":[{"color":"blue","width":2,"value":1,"zIndex":5}]},"credits":{"enabled":false},"exporting":{"enabled":true,"filename":"histPVal","buttons":{"contextButton":{"menuItems":["downloadPNG","downloadSVG","downloadPDF"]}}},"boost":{"enabled":false},"plotOptions":{"series":{"label":{"enabled":false},"turboThreshold":0,"stacking":"normal","animation":{"duration":100},"connectNulls":true,"marker":{"enabled":false}},"treemap":{"layoutAlgorithm":"squarified"},"column":{"groupPadding":0,"pointPadding":0,"borderWidth":0}},"series":[{"data":[49,0,2.220446049250313e-16,0.9999999999999996,0,0,0,4.440892098500626e-16,1.000000000000001,0,0,0,0,0.9999999999999982,0,0,0,0,0,0,1.998401444325282e-15,0,2.999999999999996,1.998401444325282e-15,0,0,1.999999999999997,0,0,0,0,0,0,0,0,0,0,0,0,0.9999999999999982,0,0,0,0,0,0,1.999999999999997,1.000000000000009,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.9999999999999982,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.9999999999999982],"name":"p-value density"},{"data":[1,1,1,1,0,0,0.9999999999999991,1,1,0,0,0,0,1,0,0,0,0,0,0,1,0,1,1,0,0,1,0.9999999999999991,0,0,0,0.9999999999999991,0,0.9999999999999991,0,0,0,0,0.9999999999999991,1,0,0,0,0,0,0,1,1,0,0.9999999999999991,0.9999999999999991,0,0.9999999999999991,0,0,0,0,0,0.9999999999999991,0,0,0,0,0,0,0,0,0,0,0.9999999999999991,0,0,0,0,0,0,0.9999999999999991,1,0.9999999999999991,0,0.9999999999999991,0,0,0,0.9999999999999991,0,0,0.9999999999999991,0,0,0.9999999999999991,0.9999999999999991,0.9999999999999991,0,0.9999999999999991,0,0.9999999999999991,0.9999999999999991,1],"name":"p-value density"}],"legend":{"enabled":false},"colors":["green","red"],"xAxis":{"title":{"text":"P-value"},"categories":[0,0.01,0.02,0.03,0.04,0.05,0.06,0.07000000000000001,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.3,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.4,0.41,0.42,0.43,0.44,0.45,0.46,0.47,0.48,0.49,0.5,0.51,0.52,0.53,0.54,0.55,0.5600000000000001,0.5700000000000001,0.58,0.59,0.6,0.61,0.62,0.63,0.64,0.65,0.66,0.67,0.68,0.6900000000000001,0.7000000000000001,0.71,0.72,0.73,0.74,0.75,0.76,0.77,0.78,0.79,0.8,0.8100000000000001,0.8200000000000001,0.8300000000000001,0.84,0.85,0.86,0.87,0.88,0.89,0.9,0.91,0.92,0.93,0.9400000000000001,0.9500000000000001,0.96,0.97,0.98,0.99]},"tooltip":{"headerFormat":"","pointFormat":"<b> {series.name} <\/b>: {point.y} ","valueDecimals":2},"annotations":[{"labelOptions":{"backgroundColor":"transparent","verticalAlign":"top","y":-30,"borderWidth":0,"x":20,"style":{"fontSize":"1.5em","color":"blue"}},"labels":[{"point":{"xAxis":0,"yAxis":0,"x":80,"y":1},"text":"pi0=1"}]}]},"theme":{"chart":{"backgroundColor":"transparent"},"colors":["#7cb5ec","#434348","#90ed7d","#f7a35c","#8085e9","#f15c80","#e4d354","#2b908f","#f45b5b","#91e8e1"]},"conf_opts":{"global":{"Date":null,"VMLRadialGradientURL":"http =//code.highcharts.com/list(version)/gfx/vml-radial-gradient.png","canvasToolsURL":"http =//code.highcharts.com/list(version)/modules/canvas-tools.js","getTimezoneOffset":null,"timezoneOffset":0,"useUTC":true},"lang":{"contextButtonTitle":"Chart context menu","decimalPoint":".","downloadCSV":"Download CSV","downloadJPEG":"Download JPEG image","downloadPDF":"Download PDF document","downloadPNG":"Download PNG image","downloadSVG":"Download SVG vector image","downloadXLS":"Download XLS","drillUpText":"◁ Back to {series.name}","exitFullscreen":"Exit from full screen","exportData":{"annotationHeader":"Annotations","categoryDatetimeHeader":"DateTime","categoryHeader":"Category"},"hideData":"Hide data table","invalidDate":null,"loading":"Loading...","months":["January","February","March","April","May","June","July","August","September","October","November","December"],"noData":"No data to display","numericSymbolMagnitude":1000,"numericSymbols":["k","M","G","T","P","E"],"printChart":"Print chart","resetZoom":"Reset zoom","resetZoomTitle":"Reset zoom level 1:1","shortMonths":["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"],"shortWeekdays":["Sat","Sun","Mon","Tue","Wed","Thu","Fri"],"thousandsSep":" ","viewData":"View data table","viewFullscreen":"View in full screen","weekdays":["Sunday","Monday","Tuesday","Wednesday","Thursday","Friday","Saturday"]}},"type":"chart","fonts":[],"debug":false},"evals":[],"jsHooks":[]}
```
