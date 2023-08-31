
suppressMessages(library(survcomp))
suppressMessages(library(tidyverse))
suppressMessages(library(plotROC))

create_directory <- function(signature){
    dir <- paste0("./TLS_project/results/siganatures_ICB/", signature)
    if (file.exists(dir)){unlink(dir, recursive = TRUE)}
    dir.create(dir, recursive = TRUE)
    dir.create(paste0(dir, "/KMPlot"))
    dir.create(paste0(dir, "/KMPlot/OS"))
    dir.create(paste0(dir, "/Overall"))
    dir.create(paste0(dir, "/Overall/OS"))
    dir.create(paste0(dir, "/Overall/Response"))
}

source("./TLS_project/code/summary_figures/Get_KMplot.R")
source("./TLS_project/code/meta_analysis/Get_Association.R")

load("./TLS_project/data/ICB_expr.Rdata")
load("./TLS_project/data/ICB_pheno.Rdata")
load("./TLS_project/data/ICB_infiltration.Rdata")

create_directory("infiltration")

cox_os <- NULL
log_response <- NULL
auc <- NULL

for (i in 1:length(names(ICB_expr))){
    pheno <- ICB_pheno[[names(ICB_expr)[i]]]
    infiltration_data <- ICB_infiltration[[names(ICB_expr)[i]]]
    
    pheno$score <- as.numeric(infiltration_data[pheno$sample_id, "infiltration"])
    pheno <- pheno[which(!is.na(pheno$score)),]
    
    d = as.data.frame( cbind( pheno$response , pheno$score ) )
    colnames(d) = c( "response" , "sig" )
    d = d[ !is.na( d$response ) , ]
    d$response = as.numeric( as.character( d$response ) )
    d$sig = as.numeric( as.character( d$sig ) )
    basicplot <- ggplot(d, aes(d = response, m = sig )) + geom_roc(n.cuts = 100, labels = FALSE)
    auc <- c(auc, as.numeric(calc_auc(basicplot)$AUC))
        
    if (!names(ICB_expr)[i] %in% c("Auslander_Melanoma_GSE115821", "Kim_GC_PRJEB25780")){
        OS_KMplot_dir <- paste0("./TLS_project/results/siganatures_ICB/", "infiltration", "/KMPlot/OS/", names(ICB_expr)[i], ".pdf")
        Get_KMplot(cancer_type = names(ICB_expr)[i], 
             status = as.numeric(pheno$status), 
             time = as.numeric(pheno$os), 
             score = as.numeric(pheno$score), 
             data = pheno, 
             dir = OS_KMplot_dir)
        cox_os <- rbind(cox_os, Get_HR_continous(status = as.numeric(pheno$status), 
                                           as.numeric(pheno$os), 
                                           score = as.numeric(pheno$score), 
                                           data = pheno))
    } else{
        cox_os <- rbind(cox_os, rep(NA, 5))
    }

    #cox_os <- rbind(cox_os, Get_HR_continous(status = pheno$os.status, time = pheno$os, score = pheno$score, data = pheno))
    log_response <- rbind(log_response, Get_coef_continous(score = pheno$score, response = pheno$response))

}

cox_os <- cbind(names(ICB_expr), cox_os) %>% as.data.frame
log_response <- cbind(names(ICB_expr), log_response) %>% as.data.frame
auc <- cbind(names(ICB_expr), auc) %>% as.data.frame
colnames(cox_os) <- c("study", "HR", "SE", "95di_low", "95di_high", "Pval")
colnames(log_response) <- c("study", "coef", "SE", "95di_low", "95di_high", "Pval")
colnames(auc) <- c("study", "AUC")
save(cox_os, file = paste0("./TLS_project/results/siganatures_ICB/", "infiltration", "/COX_OS.Rdata"))
save(log_response, file = paste0("./TLS_project/results/siganatures_ICB/", "infiltration", "/Log_Response.Rdata"))
save(auc, file = paste0("./TLS_project/results/siganatures_ICB/", "infiltration", "/AUC.Rdata"))
