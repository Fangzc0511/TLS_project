
suppressMessages(library(tidyverse))

source("./TLS_project/code/meta_analysis/meta-analysis.R")

signature_list <- list.files("./TLS_project/results/siganatures_ICB/")

for (signatureID in signature_list){
    prefix <- paste0("./TLS_project/results/siganatures_ICB/", signatureID, "/")
    cox_os_Rdata <- paste0(prefix, "COX_OS.Rdata")
    log_RData <- paste0(prefix, "Log_Response.Rdata")
    load(cox_os_Rdata)
    cox_os <- as.data.frame(cox_os)
    cox_os = cox_os[ as.numeric(cox_os$SE) <= 10 & !is.na( as.numeric(cox_os$Pval) ) , ]
    if( nrow( cox_os[ !is.na( cox_os[ , 1 ] ) , ] ) ){
      Get_Cox_Forestplot( data = cox_os , prefix = prefix , signatureID = signatureID ) 
    }
    
    load( log_RData )
    log_response <- as.data.frame(log_response)
    ## Meta-analysis of the Log Regression models (Response) Continous
    log_response = log_response[ as.numeric(log_response$SE) <= 10 & !is.na( as.numeric(log_response$Pval) ) , ]
    if( nrow( log_response[ !is.na( log_response[ , 1 ] ) , ] ) ){
      Get_LogReg_Forestplot( data = log_response , prefix= prefix , signatureID = signatureID ) 
    }
    
    
}

######per_cancer

create_directory <- function(signature, cancer_type){
    dir <- paste0("./TLS_project/results/siganatures_ICB_per_cancer/", signature, "/", cancer_type)
    if (file.exists(dir)){unlink(dir, recursive = TRUE)}
    dir.create(dir, recursive = TRUE)
    dir.create(paste0(dir, "/Overall"))
    dir.create(paste0(dir, "/Overall/OS"))
    dir.create(paste0(dir, "/Overall/Response"))
}

load("./TLS_project/data/ICB_pheno.Rdata")

all_data <- names(ICB_pheno)

tumorID_data <- list(RCC = all_data[grep("RCC", all_data)],
                     Melanoma = all_data[grep("Melanoma", all_data)],
                     NSCLC = all_data[grep("NSCLC", all_data)],
                     UC = all_data[grep("UC", all_data)],
                     GBM = all_data[grep("GBM", all_data)],
                     GC = all_data[grep("GC", all_data)]
)

for (signatureID in signature_list){
    for (tumorID in names(tumorID_data)){
        create_directory(signature = signatureID, cancer_type = tumorID)
        prefix <- paste0("./TLS_project/results/siganatures_ICB_per_cancer/",signatureID , "/", tumorID, "/")
        
        cox_os_Rdata <- paste0("./TLS_project/results/siganatures_ICB/", signatureID, "/COX_OS.Rdata")
        log_RData <- paste0("./TLS_project/results/siganatures_ICB/", signatureID, "/Log_Response.Rdata")
        
        load(cox_os_Rdata)
        cox_os <- as.data.frame(cox_os)
        cox_os <- cox_os[which(cox_os$study %in% intersect(cox_os$study, tumorID_data[[tumorID]])),]
        cox_os = cox_os[ as.numeric(cox_os$SE) <= 10 & !is.na( as.numeric(cox_os$Pval) ) , ]
        if( nrow( cox_os[ !is.na( cox_os[ , 1 ] ) , ] ) ){
          Get_Cox_Forestplot( data = cox_os , prefix= prefix  , signatureID = signatureID ) 
        }
        
        load( log_RData )
        log_response <- as.data.frame(log_response)
        log_response <- log_response[which(log_response$study %in% intersect(log_response$study, tumorID_data[[tumorID]])),]
        ## Meta-analysis of the Log Regression models (Response) Continous
        log_response = log_response[ as.numeric(log_response$SE) <= 10 & !is.na( as.numeric(log_response$Pval) ) , ]
        if( nrow( log_response[ !is.na( log_response[ , 1 ] ) , ] ) ){
          Get_LogReg_Forestplot( data = log_response , prefix= prefix , signatureID = signatureID ) 
        }
    }
}


