
source("./TLS_project/code/meta_analysis/meta-analysis.R")

signature_list <- list.files("./TLS_project/results/siganatures_TCGA/")

for (signatureID in signature_list){
    prefix <- paste0("./TLS_project/results/siganatures_TCGA/", signatureID, "/")
    cox_os_Rdata <- paste0(prefix, "COX_OS.Rdata")
    cox_pfi_Rdata <- paste0(prefix, "COX_PFI.Rdata")
    load(cox_os_Rdata)
    cox_os <- as.data.frame(cox_os)
    cox_os = cox_os[ cox_os$SE <= 10 & !is.na( cox_os$Pval ) , ]
    if( nrow( cox_os[ !is.na( cox_os[ , 1 ] ) , ] ) ){
      # cancer <- Get_Cancer( cancer=cox_os )
      # seq <- Get_Seq( seq=cox_os )
      Get_Cox_Forestplot( data = cox_os , prefix = prefix , signatureID = signatureID) 
    }
    load(cox_pfi_Rdata)
    cox_pfi <- as.data.frame(cox_pfi)
    if( nrow( cox_pfi[ !is.na( cox_os[ , 1 ] ) , ] ) ){
      # cancer <- Get_Cancer( cancer=cox_os )
      # seq <- Get_Seq( seq=cox_os )
      Get_Cox_Forestplot_Pfi( data = cox_pfi , prefix = prefix , signatureID = signatureID) 
    }
    
    
}
