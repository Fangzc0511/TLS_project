suppressMessages(library(meta))
suppressMessages(library(metafor))
suppressMessages(library(forestplot))
suppressMessages(library(dmetar))

Get_Cox_Forestplot = function( data , prefix , signatureID){
  data$study <- as.character( data$study )
  data$HR <- as.numeric(as.character( data$HR ))
  data$SE <- as.numeric(as.character( data$SE ))
  data$Pval <- as.numeric(as.character( data$Pval )) 
  
  data = data[ order( data$HR ) , ]
  
  m <- c( min( c( 0 , data$HR ) , na.rm=TRUE) - .5 , ( max( c( 0 , abs(data$HR) ) , na.rm=TRUE) ) + .5 )
  
  meta <- metagen( TE = HR,
                   seTE = SE ,
                   data = data ,
                   studlab = study ,
                   fixed = FALSE ,
                   random = TRUE ,
                   control = list( maxiter = 10000 , stepadj=0.5 ) )
  
  ######################################################################
  ######################################################################
  meta_res <- c( signatureID , 
                 meta$TE.random ,  
                 meta$seTE.random ,   
                 meta$lower.random ,   
                 meta$upper.random ,
                 meta$pval.random , 
                 meta$I2 ,
                 meta$pval.Q )
  names(meta_res) <- c( "study" , "logHR" , "se_logHR" , "CI95_low" , "CI95_high" , "Pval" , "I2" , "Pval_I2" )
  save( meta_res , file=paste( prefix, "Overall/OS/", signatureID, "_COX_OS.RData", sep="" ) )
  ######################################################################
  ######################################################################
  
  pdf( paste( prefix, "Overall/OS/", signatureID, "_COX_OS.pdf", sep="" ), height= 6, width= 10 , bg="transparent" , onefile=FALSE , family="Times" )
  forest( meta , 
          leftcols = c("studlab", "effect.ci" , "Pval" ),
          leftlabs= c( "Study" , "logHR [95%CI]" , "P-value" ) , 
          xlab = "logHR estimate",
          digits.se = 2 ,
          colgap.forest=unit(10, "mm") ,
          plotwidth = unit( 30 , "mm") , 
          pooled.totals = TRUE,
          smlab = " ",
          comb.random =TRUE,
          comb.fixed = FALSE,
          text.fixed.w = FALSE,
          layout = "JAMA",
          print.I2.ci = TRUE,
          print.Q = FALSE,
          print.pval.Q = TRUE,
          print.I2 = TRUE,
          print.tau2 = FALSE,
          resid.hetstat = FALSE,
          test.overall.random = TRUE,
          test.overall.fixed = FALSE,
          xlim = m ,  
          col.square= "black" ,  
          col.study= "black" ,  
          col.square.lines = "black" ,
          col.diamond.random  = "#1565c0"  ,
          col.diamond.lines.random  ="#1565c0" ,
          col.by = "#1565c0",
          addrow.subgroups=TRUE 
  )
  dev.off()
  
}


Get_LogReg_Forestplot = function( data , prefix  , signatureID ){
  data$study = as.character( data$study )
  data$coef = as.numeric(as.character( data$coef ))
  data$SE = as.numeric(as.character( data$SE ))
  data$Pval = as.numeric(as.character( data$Pval )) 
  
  data = data[ order( data$coef ) , ]

  
  m <- c( min( c( data$coef , 0 ) , na.rm=TRUE) - .5 , ( max( c( 0 , abs(data$coef) ) , na.rm=TRUE) ) + .5 )
  
  
  meta <- metagen( TE = coef ,
                   seTE = SE ,
                   data = data ,
                   studlab = study ,
                   fixed = FALSE ,
                   random = TRUE ,
                   control = list( maxiter = 10000 , stepadj=0.5 ) )
  
  ######################################################################
  ######################################################################
  ## Save the merged coef and Pvalue
  meta_res = c( signatureID , 
                meta$TE.random ,  
                meta$seTE.random ,   
                meta$lower.random ,   
                meta$upper.random ,
                meta$pval.random , 
                meta$I2 ,
                meta$pval.Q )
  names(meta_res) = c( "study" , "coef" , "se_coef" , "95CI_low" , "95CI_high" , "Pval" , "I2" , "Pval_I2" )
  save( meta_res , file=paste( prefix, "Overall/Response/", signatureID, "_Log_Response.RData", sep="" ) )
  ######################################################################
  ######################################################################
  
  pdf( paste( prefix, "Overall/Response/", signatureID, "_Log_Response.pdf", sep="" ), height= 6, width= 10 , bg="transparent" , onefile=FALSE , family="Times" )
  forest( meta , 
          leftcols = c("studlab", "effect.ci" , "Pval"),
          leftlabs= c( "Study" , "logOR [95%CI]" , "P-value" ) , 
          xlab = "Estimated logOR",
          digits.se = 2,
          colgap.forest=unit(10, "mm") ,
          plotwidth = unit( 30 , "mm") ,
          pooled.totals = TRUE,
          smlab = " ",
          comb.random =TRUE,
          comb.fixed = FALSE,
          text.fixed.w = FALSE,
          layout = "JAMA",
          print.I2.ci = TRUE,
          print.Q = FALSE,
          print.pval.Q = TRUE,
          print.I2 = TRUE,
          print.tau2 = FALSE,
          resid.hetstat = FALSE,
          test.overall.random = TRUE,
          test.overall.fixed = FALSE,
          xlim = m , 
          col.square= "black" ,  
          col.study= "black" ,  
          col.square.lines = "black" ,
          col.diamond.random  = "#1565c0"  ,
          col.diamond.lines.random  ="#1565c0" ,
          col.by = "#1565c0" ,
          addrow.subgroups=TRUE
  )
  
  dev.off()
  
}

Get_Cox_Forestplot_Pfi = function( data , prefix , signatureID ){
  data$study <- as.character( data$study )
  data$HR <- as.numeric(as.character( data$HR ))
  data$SE <- as.numeric(as.character( data$SE ))
  data$Pval <- as.numeric(as.character( data$Pval )) 
  
  data = data[ order( data$HR ) , ]
  # cancer = cancer[ order( cancer$HR ) , ]
  # seq = seq[ order( seq$HR ) , ]
  
  # Get xlim
  m <- c( min( c( 0 , data$HR ) , na.rm=TRUE) - .5 , ( max( c( 0 , abs(data$HR) ) , na.rm=TRUE) ) + .5 )
  
  meta <- metagen( TE = HR,
                   seTE = SE ,
                   data = data ,
                   studlab = study ,
                   fixed = FALSE ,
                   random = TRUE ,
                   control = list( maxiter = 10000 , stepadj=0.5 ) )
  
  ######################################################################
  ######################################################################
  meta_res <- c( signatureID , 
                 meta$TE.random ,  
                 meta$seTE.random ,   
                 meta$lower.random ,   
                 meta$upper.random ,
                 meta$pval.random , 
                 meta$I2 ,
                 meta$pval.Q )
  names(meta_res) <- c( "study" , "logHR" , "se_logHR" , "CI95_low" , "CI95_high" , "Pval" , "I2" , "Pval_I2" )
  save( meta_res , file=paste( prefix, "Overall/PFI/", signatureID, "_COX_PFI.RData", sep="" ) )
  ######################################################################
  ######################################################################
  
  pdf( paste( prefix, "Overall/PFI/", signatureID, "_COX_PFI.pdf", sep="" ), height= 10, width= 6 , bg="transparent" , onefile=FALSE , family="Times" )
  forest( meta , 
          leftcols = c("studlab", "effect.ci" , "Pval" ),
          leftlabs= c( "Study" , "logHR [95%CI]" , "P-value" ) , 
          xlab = "logHR estimate",
          digits.se = 2 ,
          colgap.forest=unit(10, "mm") ,
          plotwidth = unit( 30 , "mm") , 
          pooled.totals = TRUE,
          smlab = " ",
          comb.random =TRUE,
          comb.fixed = FALSE,
          text.fixed.w = FALSE,
          layout = "JAMA",
          print.I2.ci = TRUE,
          print.Q = FALSE,
          print.pval.Q = TRUE,
          print.I2 = TRUE,
          print.tau2 = FALSE,
          resid.hetstat = FALSE,
          test.overall.random = TRUE,
          test.overall.fixed = FALSE,
          xlim = m ,  
          col.square= "black" ,  
          col.study= "black" ,  
          col.square.lines = "black" ,
          col.diamond.random  = "#1565c0"  ,
          col.diamond.lines.random  ="#1565c0" ,
          col.by = "#1565c0",
          addrow.subgroups=TRUE 
  )
  dev.off()
  
}

