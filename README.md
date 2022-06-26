# PCGPier
Calculate the gene-pair based prognostic signature PCGPier for prostate cancer samples. This prognosis prediction study established a gene-pairs based RFS prediction classifier for prostate cancer patients, which is universal applicable in multiple cohorts with different gene expression detection methods. 

Installation
=
You may install this package with:
```r
if (!require("devtools")) 
    install.packages("devtools")
devtools::install_github("AHMUJia/PCGPier")
```

Example
=
```r
# load R package and internal data set
library(PCGPier)
load(system.file("extdata", "demo.RData", package = "PCGPier", mustWork = TRUE)) 

# calculate PCGPier
riskscore  <- PCGPier(expr = demo)

# print
head(PCGPier)
#     SampleID  Riskscore
# 1 GSM1133136 -1.0774821
# 2 GSM1133137 -1.6531444
# 3 GSM1133138 -3.1462112
# 4 GSM1133139  0.6813608
# 5 GSM1133140 -0.9929728
# 6 GSM1133141 -1.3352130
