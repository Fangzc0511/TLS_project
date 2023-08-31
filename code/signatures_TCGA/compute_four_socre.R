
suppressMessages(library(GSVA))
suppressMessages(library(GSEABase))
suppressMessages(library(ggplot2))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))
suppressMessages(library(tidyverse))

load("./TLS_project/data/TCGA_expr.Rdata")
load("./TLS_project/data/TCGA_pheno.Rdata")

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

ipsmap <- function (x) {
    if (x<=0) {
      ips<-0
    } else {
      if (x>=3) {
        ips<-10
      } else {
        ips<-round(x*10/3, digits=0)
      }
    }
    return(ips)
  }

TCGA_for_score <- list()

for (dataset in 1:length(names(TCGA_expr))){
    cancertype <- names(TCGA_expr)[dataset]
    mydata <- TCGA_expr[[cancertype]]
    
    mydata_matrix <- as.matrix(mydata)
    
    hallmark <- GSEABase::getGmt("./TLS_project/data/IPRES/msdb+other_1.gmt")
    es <- GSVA::gsva(mydata_matrix,
                   hallmark,
                   method='gsva',
                   mx.diff=TRUE,
                   kcdf='Gaussian',
                   parallel.sz=1)
    
    es_scale=t(apply(t(es),2,scale))
    colnames(es_scale)=colnames(es)
    df_IPRES_score=as.data.frame(apply(t(es_scale),1,mean))
    colnames(df_IPRES_score)="IPRES_score"
    
    #####IFN_score
    IFN_expanded = c('CD3D','IDO1', 'CIITA','CD3E', 'CCL5','GZMK', 'CD2','HLA-DRA', 'CXCL13', 'IL2RG','NKG7','HLA-E','CXCR6','LAG3', 'TAGAP','CXCL10','STAT1','GZMB')
    mydata_matrix_IFN=mydata_matrix[intersect(rownames(mydata_matrix),IFN_expanded),]
    df_IFN_score=as.data.frame(apply(mydata_matrix_IFN,2,mean))
    colnames(df_IFN_score)="IFN_score"
    
    ########COX_IS
    postive_gene=c("VEGFA","CCL2","CXCL8","CXCL1","CXCL2","CSF3","IL6","IL1B","IL1A")
    #postive_gene=c("VEGFA","CCL2","CXCL8","CXCL1","CXCL2","CSF3","IL6","IL1B","IL1A")
    negtive_gene=c("CCL5","CXCL9","CXCL10","CXCL11","IL12A","IL12B","IFNG","CD8A","CD8B","GZMA","GZMB","EOMES","PRF1","STAT1","TBX21")
    mydata_matrix_pos=mydata_matrix[intersect(rownames(mydata_matrix),postive_gene),]
    df_COX_pos=as.data.frame(apply(mydata_matrix_pos,2,mean))
    colnames(df_COX_pos)="COX_pos_score"
    mydata_matrix_neg=mydata_matrix[intersect(rownames(mydata_matrix),negtive_gene),]
    df_COX_neg=as.data.frame(apply(mydata_matrix_neg,2,mean))
    colnames(df_COX_neg)="COX_neg_score"
    df_COX_IS=df_COX_pos/df_COX_neg
    colnames(df_COX_IS)="COX_IS"
    
    gene_expression <- TCGA_expr[[cancertype]]
    sample_names <- names(gene_expression)
    
    IPSG<-read.table("./TLS_project/data/Immunophenogram/IPS_genes.txt",header=TRUE, sep="\t", dec = ".",check.names=FALSE)
    unique_ips_genes<-as.vector(unique(IPSG$NAME))
    
    IPS<-NULL
    MHC<-NULL
    CP<-NULL
    EC<-NULL
    SC<-NULL
    AZ<-NULL
    
    GVEC<-row.names(gene_expression)
    VEC<-as.vector(IPSG$GENE)
    ind<-which(is.na(match(VEC,GVEC)))
    MISSING_GENES<-VEC[ind]
    dat<-IPSG[ind,]
    if (length(MISSING_GENES)>0){
      cat("differently named or missing genes: ",MISSING_GENES,"\n")
    }
    for (x in 1:length(ind)) {
      print(IPSG[ind,])
    }
    
    for (i in 1:length(sample_names)){
        GE<-gene_expression[[i]]
        mGE<-mean(GE)
        sGE<-sd(GE)
        Z1<-(gene_expression[as.vector(IPSG$GENE),i]-mGE)/sGE
        W1<-IPSG$WEIGHT
        WEIGHT<-NULL
        MIG<-NULL
        k<-1
        for (gen in unique_ips_genes){
          MIG[k]<- mean(Z1[which (as.vector(IPSG$NAME)==gen)],na.rm=TRUE)
          WEIGHT[k]<- mean(W1[which (as.vector(IPSG$NAME)==gen)])
          k<-k+1
        }
        WG<-MIG*WEIGHT
        MHC[i]<-mean(WG[1:10])
        CP[i]<-mean(WG[11:20])
        EC[i]<-mean(WG[21:24])
        SC[i]<-mean(WG[25:26])
        AZ[i]<-sum(MHC[i],CP[i],EC[i],SC[i])
        IPS[i]<-ipsmap(AZ[i])
  }
    DF<-data.frame(SAMPLE=sample_names,MHC=MHC,EC=EC,SC=SC,CP=CP,AZ=AZ,IPS=IPS)
    data_4_score=cbind(DF,df_IPRES_score,df_IFN_score,df_COX_IS)
    TCGA_for_score[[cancertype]] <- data_4_score
    ######################################################################################################################
    
    
}

for (signatureID in c("IPS", "IPRES_score", "IFN_score", "COX_IS")){
    cox_os <- NULL
    cox_pfi <- NULL
    create_directory(signatureID)
    for (i in 1:length(names(TCGA_pheno))){
        pheno <- TCGA_pheno[[names(TCGA_pheno)[i]]]
        data_4_score <- TCGA_for_score[[names(TCGA_pheno)[i]]]
        
        pheno$score <- NA
        for (sample in intersect(pheno$sample_id, data_4_score$SAMPLE)){
            pheno[which(pheno$sample_id==sample), "score"] <- as.numeric(data_4_score[which(data_4_score$SAMPLE==sample), signatureID])
        }
        pheno <- pheno[which(!is.na(pheno$score)),]
        
        OS_KMplot_dir <- paste0("./TLS_project/results/siganatures_TCGA/", signatureID, "/KMPlot/OS/", names(TCGA_pheno)[i], ".pdf")
        PFI_KMplot_dir <- paste0("./TLS_project/results/siganatures_TCGA/", signatureID, "/KMPlot/PFI/", names(TCGA_pheno)[i], ".pdf")
        
        Get_KMplot(cancer_type = names(TCGA_pheno)[i], status = pheno$os.status, time = pheno$os, score = pheno$score, data = pheno, dir = OS_KMplot_dir)
        Get_KMplot(cancer_type = names(TCGA_pheno)[i], status = pheno$pfi.status, time = pheno$pfi, score = pheno$score, data = pheno, dir = PFI_KMplot_dir)
        
        cox_os <- rbind(cox_os, Get_HR_continous(status = pheno$os.status, time = pheno$os, score = pheno$score, data = pheno))
        cox_pfi <- rbind(cox_pfi, Get_HR_continous(status = pheno$pfi.status, time = pheno$pfi, score = pheno$score, data = pheno))
        
    }
    cox_os <- cbind(names(TCGA_pheno), cox_os)
    cox_pfi <- cbind(names(TCGA_pheno), cox_pfi)
    colnames(cox_os) <- c("study", "HR", "SE", "95di_low", "95di_high", "Pval")
    colnames(cox_pfi) <- c("study", "HR", "SE", "95di_low", "95di_high", "Pval")
    save(cox_os, file = paste0("./TLS_project/results/siganatures_TCGA/", signatureID, "/COX_OS.Rdata"))
    save(cox_pfi, file = paste0("./TLS_project/results/siganatures_TCGA/", signatureID, "/COX_PFI.Rdata"))
    

}
