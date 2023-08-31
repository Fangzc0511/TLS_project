
suppressMessages(library(pheatmap))
suppressMessages(library(RColorBrewer))

all_signatures <- c("gs_12_Chemokine", "gs_Clubb_James_H_A_etal", "gs_Meylan_Maxime_etal", "gs_PagliaruloFabio_etal", "gs_Wang_Bin_etal", "gs_Wu_Zhenghao_etal", "TMB", "infiltration", "IPS", "IPRES_score", "IFN_score", "COX_IS")
res <- NULL
for( i in 1:length(all_signatures)){
  file = paste( "./TLS_project/results/siganatures_ICB/" , all_signatures[ i ] , "/Overall/Response/", all_signatures[ i ], "_Log_Response.RData" , sep="" )
  if( file.exists( file ) ){
    load( file )
    res = rbind( res , c( meta_res , "sig" ) )
  }
}

res = as.data.frame( cbind( res , NA ) )
colnames(res) = c( "feature" , "coef" , "se_coef" , "CI95_low" , "CI95_high" , "Pval" , "I2" , "Pval_I2" , "type" , "Padj_BH" )

res$feature = as.character( res$feature )
res$type = as.character( res$type )
res$coef = round( as.numeric( as.character( res$coef ) ) , 2 )
res$se_coef = round( ifelse( as.numeric( as.character( res$se_coef ) ) >= 0.3 , 0.3 , as.numeric( as.character( res$se_coef ) ) ) , 2 )
res$CI95_low = round( as.numeric( as.character( res$CI95_low ) ) , 2 )
res$CI95_high = round( as.numeric( as.character( res$CI95_high ) ) , 2 )
res$Pval = as.numeric( as.character( res$Pval ) )

res$Padj_BH = round( p.adjust( res$Pval , method="BH" ) , 2)

response = res


res = NULL
for( i in 1:length(all_signatures)){
  file = paste( "./TLS_project/results/siganatures_ICB/" , all_signatures[ i ] , "/Overall/OS/", all_signatures[ i ], "_COX_OS.RData" , sep="" )
  if( file.exists( file ) ){
    load( file )
    res = rbind( res , c( meta_res , "sig" ) )
  }
}

res = as.data.frame( cbind( res , NA ) )
colnames(res) = c( "feature" , "HR" , "se_hr" , "CI95_low" , "CI95_high" , "Pval" , "I2" , "Pval_I2" , "type" , "Padj_BH" )

res$feature = as.character( res$feature )
res$type = as.character( res$type )
res$HR = round( as.numeric( as.character( res$HR ) ) , 2 )
res$se_hr = round( ifelse( as.numeric( as.character( res$se_hr ) ) >= 0.3 , 0.3 , as.numeric( as.character( res$se_hr ) ) ) , 2 )
res$CI95_low = round( as.numeric( as.character( res$CI95_low ) ) , 2 )
res$CI95_high = round( as.numeric( as.character( res$CI95_high ) ) , 2 )
res$Pval = as.numeric( as.character( res$Pval ) )

res$Padj_BH = round( p.adjust( res$Pval , method="BH" ) , 2)

os = res

data = cbind( response$coef , os$HR )
rownames(data) = response$feature
colnames(data) = c( "Response" , "OS" )

pval = cbind( response$Pval , os$Pval )
rownames(pval) = response$feature
colnames(pval) = c( "Response" , "OS" )
pval[ is.na(pval) ] = 1

padj = matrix( p.adjust( as.vector( as.matrix(pval) ) , method='fdr') , ncol=2 )
rownames(padj) = response$feature
colnames(padj) = c( "Response" , "OS" )
padj[ is.na(padj) ] = 1

annotation_col = data.frame(
                              OS_Sig = factor( ifelse( round( padj[ , 'OS' ] , 2 ) <= .05 , 'FDR' , ifelse( round( pval[ , "OS" ] , 2 ) <= 0.05 , "Pvalue.Only" ,  'NS' ) ) ) ,  
                              # PFS_Sig = factor( ifelse( round( padj[ , 'PFS' ] , 2 ) <= .05 , 'FDR' , ifelse( round( pval[ , "PFS" ] , 2 ) <= 0.05 , "Pvalue.Only" , 'NS' ) ) ) ,
                              Response_Sig = factor( ifelse( round( padj[ , 'Response' ] , 2 ) <= .05 , 'FDR' , ifelse( round( pval[ , "Response" ] , 2 ) <= 0.05 , "Pvalue.Only" , 'NS' ) ) )
)

rownames(annotation_col) = rownames(data)

# Signature_Type = c( sensitive = "#4caf50" , resistance = "#8e44ad" ) ,
ann_colors = list(
                    Response_Sig = c( FDR = "#222021", Pvalue.Only = "#828282", NS = "white" ) ,
                    # PFS_Sig = c( FDR = "#222021", Pvalue.Only = "#828282", NS = "white" ),
                    OS_Sig = c( FDR = "#222021", Pvalue.Only = "#828282", NS = "white" ) )

neg = seq( round( min( data , na.rm=TRUE ) , 1 ) , 0 , by=.05 )
neg = neg[ -length(neg)]
pos = seq( 0 , round( max( data , na.rm=TRUE ) , 1 ) , by=.05 )

col = c( colorRampPalette( rev( brewer.pal(4, "Blues") ) )( length(neg) ) ,
         colorRampPalette( brewer.pal(4, "OrRd") )( length(pos) ) 
)

pheatmap( t( data[ order( rowSums(data)) , ] ) , cluster_rows=FALSE , cluster_cols=FALSE , scale="none" , annotation_col = annotation_col, annotation_colors = ann_colors, 
          col = col , breaks = c( neg , pos ) , na_col="white" , border_color="#424242", angle_col = 315,  fontfamily="serif" ,
          #display_numbers = data , number_color="black" , fontsize_number= 10 , 
          filename="./TLS_project/results/summary_figures/heatmap/Heatmap_Summary_ICB.pdf", height=2.8, width=9.5 )  

pheatmap( t( data[ order( rowSums(data)) , ] ) , cluster_rows=FALSE , cluster_cols=FALSE , annotation_col = annotation_col, annotation_colors = ann_colors, 
          col = col , breaks = c( neg , pos ) , na_col="white" , border_color="#424242", angle_col = 315,  fontfamily="serif" ,
          #display_numbers = data , number_color="black" , fontsize_number= 10 , 
          filename="./TLS_project/results/summary_figures/heatmap/Heatmap_Summary_legend_ICB.pdf", height=4, width=9.5 )  
