
suppressMessages(library(survcomp))
suppressMessages(library(tidyverse))

create_directory <- function(signature){
    dir <- paste0("./TLS_project/results/siganatures_TCGA/", signature)
    if (file.exists(dir)){unlink(dir, recursive = TRUE)}
    dir.create(dir, recursive = TRUE)
    dir.create(paste0(dir, "/KMPlot"))
    dir.create(paste0(dir, "/KMPlot/OS"))
    dir.create(paste0(dir, "/KMPlot/PFI"))
    dir.create(paste0(dir, "/Overall"))
    dir.create(paste0(dir, "/Overall/OS"))
    dir.create(paste0(dir, "/Overall/PFI"))
}

source("./TLS_project/code/summary_figures/Get_KMplot.R")
source("./TLS_project/code/meta_analysis/Get_Association.R")

load("./TLS_project/data/TCGA_pheno.Rdata")
load("./TLS_project/data/TCGA_infiltration.Rdata")

create_directory("infiltration")

cox_os <- NULL
cox_pfi <- NULL

for (i in 1:length(names(TCGA_pheno))){
    pheno <- TCGA_pheno[[names(TCGA_pheno)[i]]]
    infiltration_data <- TCGA_infiltration[[names(TCGA_pheno)[i]]]
    
        
    pheno$score <- as.numeric(infiltration_data[pheno$sample_id, "infiltration"])
    pheno <- pheno[which(!is.na(pheno$score)),]
        
    OS_KMplot_dir <- paste0("./TLS_project/results/siganatures_TCGA/", "infiltration", "/KMPlot/OS/", names(TCGA_pheno)[i], ".pdf")
    PFI_KMplot_dir <- paste0("./TLS_project/results/siganatures_TCGA/", "infiltration", "/KMPlot/PFI/", names(TCGA_pheno)[i], ".pdf")
        
    Get_KMplot(cancer_type = names(TCGA_pheno)[i], status = pheno$os.status, time = pheno$os, score = pheno$score, data = pheno, dir = OS_KMplot_dir)
    Get_KMplot(cancer_type = names(TCGA_pheno)[i], status = pheno$pfi.status, time = pheno$pfi, score = pheno$score, data = pheno, dir = PFI_KMplot_dir)
        
    cox_os <- rbind(cox_os, Get_HR_continous(status = pheno$os.status, time = pheno$os, score = pheno$score, data = pheno))
    cox_pfi <- rbind(cox_pfi, Get_HR_continous(status = pheno$pfi.status, time = pheno$pfi, score = pheno$score, data = pheno))

    }

cox_os <- cbind(names(TCGA_pheno), cox_os)
cox_pfi <- cbind(names(TCGA_pheno), cox_pfi)
colnames(cox_os) <- c("study", "HR", "SE", "95di_low", "95di_high", "Pval")
colnames(cox_pfi) <- c("study", "HR", "SE", "95di_low", "95di_high", "Pval")
save(cox_os, file = paste0("./TLS_project/results/siganatures_TCGA/", "infiltration", "/COX_OS.Rdata"))
save(cox_pfi, file = paste0("./TLS_project/results/siganatures_TCGA/", "infiltration", "/COX_PFI.Rdata"))
