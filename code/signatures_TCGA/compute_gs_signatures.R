
suppressMessages(library(GSVA))
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

load("./TLS_project/data/TCGA_expr.Rdata")
load("./TLS_project/data/TCGA_pheno.Rdata")
load("./TLS_project/data/gene_signatures.Rdata")

for (i in 1:length(names(gene_signatures))){
    create_directory(names(gene_signatures)[i])
    geneset <- gene_signatures[[names(gene_signatures)[i]]] %>% data.frame()
    cox_os <- NULL
    cox_pfi <- NULL
    for (j in 1:length(names(TCGA_expr))){
        expr <- TCGA_expr[[names(TCGA_expr)[j]]]
        pheno <- TCGA_pheno[[names(TCGA_expr)[j]]]
        
        ssgsea <- GSVA::gsva(unique(as.matrix(expr)),
                         geneset,
                         method='ssgsea',
                         mx.diff=TRUE,
                         kcdf='Gaussian',
                         parallel.sz=1)
        df_ssgsea=as.data.frame(t(ssgsea))
        colnames(df_ssgsea)="ssgsea"
        
        pheno$score <- df_ssgsea[match(pheno$sample_id,rownames(df_ssgsea)),]
        pheno <- pheno[which(!is.na(pheno$score)),]
        
        OS_KMplot_dir <- paste0("./TLS_project/results/siganatures_TCGA/", names(gene_signatures)[i], "/KMPlot/OS/", names(TCGA_expr)[j], ".pdf")
        PFI_KMplot_dir <- paste0("./TLS_project/results/siganatures_TCGA/", names(gene_signatures)[i], "/KMPlot/PFI/", names(TCGA_expr)[j], ".pdf")
        
        Get_KMplot(cancer_type = names(TCGA_expr)[j], status = pheno$os.status, time = pheno$os, score = pheno$score, data = pheno, dir = OS_KMplot_dir)
        Get_KMplot(cancer_type = names(TCGA_expr)[j], status = pheno$pfi.status, time = pheno$pfi, score = pheno$score, data = pheno, dir = PFI_KMplot_dir)
        
        cox_os <- rbind(cox_os, Get_HR_continous(status = pheno$os.status, time = pheno$os, score = pheno$score, data = pheno))
        cox_pfi <- rbind(cox_pfi, Get_HR_continous(status = pheno$pfi.status, time = pheno$pfi, score = pheno$score, data = pheno))

    }
    cox_os <- cbind(names(TCGA_expr), cox_os)
    cox_pfi <- cbind(names(TCGA_expr), cox_pfi)
    colnames(cox_os) <- c("study", "HR", "SE", "95di_low", "95di_high", "Pval")
    colnames(cox_pfi) <- c("study", "HR", "SE", "95di_low", "95di_high", "Pval")
    save(cox_os, file = paste0("./TLS_project/results/siganatures_TCGA/", names(gene_signatures)[i], "/COX_OS.Rdata"))
    save(cox_pfi, file = paste0("./TLS_project/results/siganatures_TCGA/", names(gene_signatures)[i], "/COX_PFI.Rdata"))
}
