### setwd("your work dictionary")
setwd("/home/zcfang/")
### Calculate GSVA scores of six gene signatures， and the asscociation with os and pri
source("./TLS_project/code/signatures_TCGA/compute_gs_signatures.R")
### Calculate immune infiltration degree， and the asscociation with os and pfi
source("./TLS_project/code/signatures_TCGA/compute_infiltration.R")
### Calculate Calculate four scores related to immunity ， and the asscociation with os and pfi
source("./TLS_project/code/signatures_TCGA/compute_infiltration.R")
### Calculate TMB ， and the asscociation with os and pfi
source("./TLS_project/code/signatures_TCGA/compute_TMB.R")
### get_meta-analysis
source("./TLS_project/code/signatures_TCGA/get_meta-analysis.R")
