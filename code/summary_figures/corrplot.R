
suppressMessages(library(GSVA))
suppressMessages(library(tidyverse))
suppressMessages(library(survival))
suppressMessages(library(reshape2))
suppressMessages(library(corrplot))

load("./TLS_project/data/ICB_expr.Rdata")
load("./TLS_project/data/ICB_pheno.Rdata")
load("./TLS_project/data/gene_signatures.Rdata")
load("./TLS_project/data/ICB_four_score.Rdata")
load("./TLS_project/data/ICB_infiltration.Rdata")

file = "./TLS_project/results/summary_figures/corrplot/"
if( dir.exists( file ) ) {
   unlink( file , recursive = TRUE )
}
dir.create( file )

cancer_type <- names(ICB_expr)
score_list <- list()
all_signatures <- c(names(gene_signatures),"infiltration", "IPS", "IPRES_score", "IFN_score", "COX_IS")
for (i in 1:length(cancer_type)){
    df_expression <- ICB_expr[[cancer_type[i]]]
    signature_score <- matrix(data = NA, nrow = ncol(df_expression), ncol = length(all_signatures))
    rownames(signature_score) <- colnames(df_expression)
    colnames(signature_score) <- all_signatures
    
    #########gene signatures
    for (j in 1:length(names(gene_signatures))){
        gene_list <- gene_signatures[[j]]

        geneset=data.frame(gene_list)
        ssgsea <- GSVA::gsva(unique(as.matrix(df_expression)),
                             geneset,
                             method='ssgsea',
                             mx.diff=TRUE,
                             kcdf='Gaussian',
                             parallel.sz=1)
        df_ssgsea=as.data.frame(t(ssgsea))
        colnames(df_ssgsea)="ssgsea"

        signature_score[, names(gene_signatures)[j]] <- df_ssgsea[match(rownames(df_ssgsea), rownames(signature_score)), "ssgsea"]
    }
    
    ########infiltration
    infiltration_score <- ICB_infiltration[[cancer_type[i]]]
    signature_score[, "infiltration"] <- infiltration_score[match(infiltration_score$sample_id, rownames(signature_score)), "infiltration"]
    
    ########four score
    four_score <- subset(ICB_four_score[[cancer_type[i]]], select = c("SAMPLE", "IPS", "IPRES_score", "IFN_score", "COX_IS"))
    for (k in c("IPS", "IPRES_score", "IFN_score", "COX_IS")){
        signature_score[, k] <- four_score[match(four_score$SAMPLE, rownames(signature_score)), k]
    }
    score_list[[i]] <- signature_score

}
names(score_list) <- cancer_type

cor <- NULL
for (i in 1:length(cancer_type)) {
  cor_percancer = cor( score_list[[i]] , method="s" )
  cor_percancer = cbind( cancer_type[i], nrow( score_list[[i]] ) , melt( cor_percancer ) )
  cor_percancer = cor_percancer[ cor_percancer$Var1 != cor_percancer$Var2 , ]
  cor <- rbind(cor, cor_percancer)
}
colnames( cor ) = c( "dataset", "N" , "var1" , "var2" , "value" )
cor = cbind( cor , apply( cor , 1 , function(x){ paste( x[ "var1" ] , x[ "var2" ] , sep="_" ) } ) )
colnames( cor ) = c( "dataset" , "N" , "var1" , "var2" , "value" , "id" )
cor = cor[ , c( "dataset" , "N" , "var1" , "var2" , "id" , "value" ) ]

mat = matrix( NA , nrow= length(all_signatures) , ncol= length(all_signatures) )
colnames(mat) = rownames(mat) = all_signatures

for( i in 1:length( all_signatures ) ){
  for( j in 1:length( all_signatures ) ){
    if( i != j ){
      c = cor[ cor$var1 %in% all_signatures[i] & cor$var2 %in% all_signatures[j] , ]
      
      mat[ all_signatures[i] , all_signatures[j] ] = round( median( c$value[ !is.na( c$value ) ] , na.rm=TRUE ) , 2 )
      
    } else{
      mat[ all_signatures[i] , all_signatures[j] ] = 1
    }
  }
}

pdf( "./TLS_project/results/summary_figures/corrplot/CorrPlot_Overall.pdf" , height=10 , width=10 , bg="transparent" , family="Times" )
corrplot(mat , tl.cex = 1.3,
         method = "color" , type = "lower" , lower.col = "black", number.cex = .4 ,  tl.col = "black" , tl.srt = 45 , diag = FALSE )
text(5,13.6,"Overall",cex=2,family="serif")
dev.off()

tumorID_data <- list(RCC = cancer_type[grep("RCC", cancer_type)],
                     Melanoma = cancer_type[grep("Melanoma", cancer_type)],
                     NSCLC = cancer_type[grep("NSCLC", cancer_type)],
                     UC = cancer_type[grep("UC", cancer_type)],
                     GBM = cancer_type[grep("GBM", cancer_type)],
                     GC = cancer_type[grep("GC", cancer_type)]
)

for( k in 1:length(tumorID_data) ){
  
  mat= matrix( NA , nrow= length(all_signatures) , ncol= length(all_signatures) )
  colnames(mat) = rownames(mat) = all_signatures
  
  for( i in 1:length( all_signatures ) ){
      for( j in 1:length( all_signatures ) ){
        if( i != j ){
          c = cor[ cor$var1 %in% all_signatures[i] & cor$var2 %in% all_signatures[j] , ]

          mat[ all_signatures[i] , all_signatures[j] ] = round( median( c$value[ !is.na( c$value ) ] , na.rm=TRUE ) , 2 )

        } else{
          mat[ all_signatures[i] , all_signatures[j] ] = 1
        }
      }
    }
  
  pdf( paste( "./TLS_project/results/summary_figures/corrplot//CorrPlot_", names(tumorID_data)[k] , ".pdf" , sep="" ) , height=10 , width=10 , bg="transparent", family="Times" )
  corrplot(mat,tl.cex = 1.3, 
           method = "color" , type = "lower" , lower.col = "black", number.cex = .4 , tl.col = "black", tl.srt = 45 , diag = FALSE )
  text(5,13.6,names(tumorID_data)[k],cex=2.)
  dev.off()
}
