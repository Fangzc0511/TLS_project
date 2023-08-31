
all_sample <- c(3,5,12,13,15)
cutoff <- 0.145


suppressMessages(library(tidyverse))
suppressMessages(library(GSVA))
suppressMessages(library(meta))
suppressMessages(library(metafor))
suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(beeswarm))
suppressMessages(library(plotROC))
suppressMessages(library(survcomp))

load("./TLS_project/data/ICB_expr.Rdata")
load("./TLS_project/data/ICB_pheno.Rdata")
load("./TLS_project/data/gene_signatures.Rdata")

#file = "./TLS_project/results/denovo_Single_Gene"
#if( dir.exists( file ) ) {
#   unlink( file , recursive = TRUE )
#}
#dir.create( file )

total_gene <- c()
for (i in 1:length(names(gene_signatures))){
    total_gene <- union(total_gene, gene_signatures[[names(gene_signatures)[i]]][,1])
}
candidate_gene <- total_gene
for (i in 1:length(names(ICB_expr))){
    candidate_gene <- intersect(candidate_gene, rownames(ICB_expr[[names(ICB_expr)[[i]]]]))
}

ntest <- length(all_sample)

result_per_test <- NULL

cancer_type_test <- names(ICB_expr)[all_sample]
cancer_type_train <- names(ICB_expr)[setdiff(1:length(names(ICB_expr)), all_sample)]

pval = matrix( NA , nrow=length(candidate_gene) , ncol=length(cancer_type_train) )
colnames(pval) = cancer_type_train
rownames(pval) = candidate_gene
coef = matrix( NA , nrow=length(candidate_gene) , ncol=length(cancer_type_train) )
colnames(coef) = cancer_type_train
rownames(coef) = candidate_gene
se = matrix( NA , nrow=length(candidate_gene) , ncol=length(cancer_type_train) )
colnames(se) = cancer_type_train
rownames(se) = candidate_gene


for (per_train in 1:length(cancer_type_train)){
    exp <- ICB_expr[[cancer_type_train[per_train]]][rownames(ICB_expr[[cancer_type_train[per_train]]]) %in% candidate_gene, ]
    pheno <- ICB_pheno[[cancer_type_train[per_train]]]
    for (k in 1:length(rownames(exp))){
        g <- exp[k, ] %>% as.numeric() %>% scale()
        g[is.na(g)] <- 0
        rownames(g) <- colnames(exp[k, ])
        pheno$score <- g[match(pheno$sample_id,rownames(g)),]
        fit <- glm( pheno$response ~ pheno$score , family=binomial( link="logit" ) )
        if (nrow(summary(fit)$coefficients)==1){break}
        pval[rownames(exp)[k], cancer_type_train[per_train]] = summary(fit)$coefficients[ 2 , 4 ]
        coef[rownames(exp)[k], cancer_type_train[per_train]] = round( summary(fit)$coefficients[ 2 , 1  ] , 2 )
        se[rownames(exp)[k], cancer_type_train[per_train]] = round( summary(fit)$coefficients[ 2 , 2 ] , 2 ) 
    }
}
fdr <- matrix( p.adjust( pval , method= "fdr" ) , ncol= ncol( pval ) , nrow= nrow( pval ) , dimnames= dimnames( pval ) ) 


meta_res = NULL
for(i in 1:nrow(coef)){
        data = as.data.frame( cbind( colnames(coef) , coef[ i , ] , se[ i , ] , fdr[ i , ] ) )
        colnames(data) = c( "study" , "coef" , "se" , "fdr" )
        data$coef = as.numeric(as.character( data$coef ))
        data$se = as.numeric(as.character( data$se ))
        data$fdr = as.numeric(as.character( data$fdr )) 
        data = data[ !is.na(data$coef) , ]

        meta <- metagen( TE = coef,
                       seTE = se,
                       data = data,
                       studlab = study ,
                       fixed = FALSE ,
                       random = TRUE ,
                       control = list( maxiter = 10000 , stepadj=0.5 ) )

        meta_res <- rbind( meta_res , c( rownames(coef)[i] , 
                                       meta$TE.random ,  
                                       meta$seTE.random ,  
                                       meta$pval.random , 
                                       meta$I2 , 
                                       meta$pval.Q ) )
        }

meta_res = as.data.frame(meta_res)
colnames(meta_res) = c( "gene" , "coef" , "se" , "pval" , "I2" , "I2_pval" )
meta_res$gene = as.character( meta_res$gene )
meta_res$coef = as.numeric(as.character( meta_res$coef ))
meta_res$se = as.numeric(as.character( meta_res$se ))
meta_res$pval = as.numeric(as.character( meta_res$pval )) 
meta_res$I2 = as.numeric(as.character( meta_res$I2 ))
meta_res$I2_pval = as.numeric(as.character( meta_res$I2_pval )) 
meta_res$fdr <- p.adjust( meta_res$pval , method= "fdr" ) 

meta_res_filter <- meta_res[which(abs(meta_res$coef) > cutoff),]
response_gene <- meta_res_filter[,1][which(meta_res_filter[,2] > 0)] %>% as.data.frame()
write.table(response_gene, file = "./TLS_project/results/denovo_Single_Gene/response_gene.txt", row.names = F, quote = F, col.names = F)


signature_list <- list.files("./TLS_project/results/siganatures_ICB/")
signature_list <- signature_list[-which(signature_list=="TMB")]
data_auc <- matrix(data = NA, nrow = length(cancer_type_test), ncol = length(signature_list) + 1 )
rownames(data_auc) <- cancer_type_test
colnames(data_auc) <- c(signature_list, "TLSscore")

for (signatureID in signature_list){
    prefix <- paste0("./TLS_project/results/siganatures_ICB/", signatureID)
    auc_RData <- paste0(prefix, "/AUC.Rdata")
    load( auc_RData )
    auc <- as.data.frame(auc)
    data_auc[, signatureID] <- auc[match(rownames(data_auc), auc$study) ,]$AUC

}

data_auc[which(is.na(data_auc))] <- 0

for (per_test in 1:length(cancer_type_test)){
    df_expression <- ICB_expr[[cancer_type_test[per_test]]]
    response_geneset=data.frame(response_gene)
    response_ssgsea <- GSVA::gsva(unique(as.matrix(df_expression)),
                                  response_geneset,
                                  method='ssgsea',
                                  mx.diff=TRUE,
                                  kcdf='Gaussian',
                                  parallel.sz=1)
    df_response_ssgsea=as.data.frame(t(response_ssgsea))
    colnames(df_response_ssgsea)="response_ssgsea"

    pheno_test <- ICB_pheno[[cancer_type_test[per_test]]]
    pheno_test$score <- df_response_ssgsea[match(pheno_test$sample_id, rownames(df_response_ssgsea)),]
    
    d = as.data.frame( cbind( pheno_test$response , pheno_test$score ) )
    colnames(d) = c( "response" , "sig" )
    d = d[ !is.na( d$response ) , ]
    d$response = as.numeric( as.character( d$response ) )
    d$sig = as.numeric( as.character( d$sig ) )
    basicplot <- ggplot(d, aes(d = response, m = sig )) + geom_roc(n.cuts = 100, labels = FALSE)
    data_auc[cancer_type_test[per_test], "TLSscore"] <- calc_auc(basicplot)$AUC

}

data_auc <- data_auc[,names(sort(colMeans(apply(data_auc, 2, as.numeric))))]



##########################################################################
##########################################################################
df1 = melt(t(data_auc))

df2 <- data.frame(cohort=df1$Var2,sig=df1$Var1,auc=as.numeric(df1$value))

df2 <- df2[which(!is.na(df2$auc)),]

study_color = brewer.pal( n = length(cancer_type_test) ,name = "Dark2")
names( study_color ) = cancer_type_test



pdf( "./TLS_project/results/summary_figures/boxplot_auc/auc.pdf" , height=9,width=15,bg="transparent", family="Times")
zones=matrix(c(1,2), ncol=2, byrow=TRUE)
layout(zones, widths=c(3.8,1))


################################################
################################################
## Plot
par(mar=c(12,3,3,0))
xLabels <- c(signature_list, "TLSscore")
yLabels <- seq( round( min( df2$auc , na.rm=TRUE ) , 1 ) ,  round( max( df2$auc , na.rm=TRUE ) , 1 ) + 0.2  , by=0.2 )
boxplot( auc ~ sig , data= df2 , ylab="auc value" , xlab="" , main="" , 
         col= "white" , 
         boxlty = 1 , outline= FALSE , axes= FALSE , ylim= c( min( df2$auc ) ,  max(df2$auc)))

beeswarm( auc ~ sig , data= df2 , 
          pwpch =  rep(19,times=length(df2$cohort))  , corral= "wrap" , cex= .8 , 
          pwcol = adjustcolor( study_color[ df2$cohort ] , alpha.f = .9 ) , 
          add = TRUE)

axis(side = 2, at=yLabels, labels= as.character( yLabels ), las= 2,
     cex.axis=1,tick=1 , col="black")

mtext( "AUC Distribution across the validation cohorts" ,cex = 2)
for (i in 1:length(signature_list)){
  t <- t.test(df2[which(df2$sig==colnames(data_auc)[i]),"auc"],as.numeric(data_coef[,"TLSscore"]))
    #t <- t.test(df2[which(df2$sig==colnames(data_auc)[i]),"auc"],as.numeric(data_coef[,"TLSscore"]),alternative = "less")

  label <- ifelse(t$p.value < 0.001, "***", 
                  ifelse(t$p.value < 0.01, "**", 
                         ifelse(t$p.value < 0.05, "*", paste0("P=", as.character(round(t$p.value, 3))))))

  text(i,0.83,label,cex=1.1)
  
}


Map(axis, side = 1, at=seq( 1 , length(xLabels) , 1 ) , 
    col.axis = ifelse( xLabels %in% "TLSscore" , "#f9a825" ,  "#1976d2" ), 
    padj=.5, labels=xLabels, las= 2,
    cex.axis=1.2 , tick=1 , col= "black" )

axis( 1 , at= seq( 1 , length(xLabels) , 1 ) , labels= FALSE )

#abline( h= 0 , lty= 2 , lwd= 1 )

################################################
################################################
## Legend

par(mar=c(0,0,1,0))
plot.new()
legend( x=-0.08 , y=.95 ,title="Study",
        legend = cancer_type_test,
        pt.bg= adjustcolor( study_color[ sort( unique( df2$cohort ) ) ] , alpha.f = .9 ) ,
        col= study_color[ sort( unique( df2$cohort ) ) ]  ,
        pt.cex= 1 ,
        y.intersp=1.2, cex=1.25, pch= c( 21 , 21 , 21 , 21 , 21 , 24 ) , lwd=1, lty=0, bty='n', ncol=1)
dev.off()
