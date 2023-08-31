
suppressMessages(library(RColorBrewer))
suppressMessages(library(beeswarm))
suppressMessages(library(tidyverse))

all_sample <- c(3,5,12,13,15)

file = "./TLS_project/results/summary_figures/boxplot/"
if( dir.exists( file ) ) {
   unlink( file , recursive = TRUE )
}
dir.create( file )

load("./TLS_project/data/ICB_expr.Rdata")
load("./TLS_project/data/ICB_pheno.Rdata")
response_geneset=read.table("./TLS_project/results/denovo_Single_Gene/response_gene.txt",sep="\t",header = F)

cancer_type <- names(ICB_expr)[all_sample]
for (per_test in 1:length(cancer_type)){
    df_expression <- ICB_expr[[cancer_type[per_test]]]
    response_geneset=data.frame(response_geneset)
    response_ssgsea <- GSVA::gsva(unique(as.matrix(df_expression)),
                                  response_geneset,
                                  method='ssgsea',
                                  mx.diff=TRUE,
                                  kcdf='Gaussian',
                                  parallel.sz=1)
    df_response_ssgsea=as.data.frame(t(response_ssgsea))
    colnames(df_response_ssgsea)="response_ssgsea"

    pheno_test <- ICB_pheno[[cancer_type[per_test]]]
    pheno_test$score <- df_response_ssgsea[match(pheno_test$sample_id, rownames(df_response_ssgsea)),]
    
    split_data <- split(pheno_test$score, pheno_test$response)
    
    re_score <- split_data[[which(names(split_data)==1)]]
    nonre_score <- split_data[[which(names(split_data)==0)]]
    df1 <- data.frame(response=rep("R", length(re_score)), TLSscore=re_score)
    df2 <- data.frame(response=rep("NR", length(nonre_score)), TLSscore=nonre_score)
    data <- rbind(df1,df2)

    pdf( paste( "./TLS_project/results/summary_figures/boxplot/", cancer_type[per_test] , ".pdf" , sep = "" ) , height=3.5 ,width=3.5 , bg="transparent", family="Times")
    xLabels <- paste( names(table(data$response)) , "\n(N=" , table(data$response) , ")" , sep = "" )
    yLabels <- seq( round( min( data$TLSscore , na.rm=TRUE ) , 1 ) ,  round( max( data$TLSscore , na.rm=TRUE ) , 1 ) , by=.25 ) 
    boxplot( TLSscore ~ response , data= data , ylab= "TLSscore" , xlab="Response" , main="" , 
           col=adjustcolor( c( "#f44336" , "#00bcd4" ), alpha.f = .4), 	    	 
           boxlty = 1 ,outline=FALSE, axes=FALSE, ylim = c( min( data$TLSscore , na.rm = TRUE ) ,  max( data$TLSscore , na.rm = TRUE ) ) )


    beeswarm( TLSscore ~ response , data= data , 
            pch = 19, corral="wrap", cex=.2, 
            col=adjustcolor( c( "#f44336" , "#00bcd4" ), alpha.f = .5), 
            add = TRUE)

    axis(side = 2, at=yLabels, labels=yLabels, las= 2,
       cex.axis=1,tick=1,col="black")

    axis(side = 1, at=seq(1,length(xLabels),1) , padj=.5, labels=xLabels, las= 1,
       cex.axis=1,tick=1,col="black")

    t = t.test( TLSscore ~ response , data= data ,alternative = "less")
    mtext( paste( as.character(cancer_type[per_test]) , "\n", "p.value",
                ifelse( t$p.value <= .001 , "≤ 0.001" , paste( "=" , format.pval( t$p.value , 1 ) ) ) , sep="" ) , 
         col="#6D6D6D" )

    dev.off()


}
