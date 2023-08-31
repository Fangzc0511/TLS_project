suppressMessages(library(survival))

Get_HR_continous = function( status , time , score , data ){
  status <- as.numeric(status)
  time <- as.numeric(time)
  score <- as.numeric(score)
  
  data = data.frame( status=status , time=time , score=score )
  data$time = as.numeric(as.character(data$time))
  data$variable = as.numeric( as.character(data$score) )
  
  cox = coxph( formula= Surv( time , status ) ~ score , data=data )
  mean <- round( coef( cox )[1] , 2 )
  low <- round( confint( cox , level=.95 )[ 1 , 1 ] , 2 )
  up <- round( confint( cox , level=.95 )[ 1 , 2 ] , 2 )
  pval <- summary(cox)$coef[1,5]
  
  c( mean , round( summary(cox)$coefficients[3] , 2 ) , low , up , pval )
}

Get_coef_continous = function( score , response ){
  score = as.numeric(score)
  response = as.numeric(response)
  
  fit = glm( response ~ score , family=binomial( link="logit" ) )
  coef <- round( summary(fit)$coefficients[ 2 ,  1  ] , 2 ) 
  se <- round( summary(fit)$coefficients[ 2 , 2 ] , 2 )
  low <- round( confint(fit)[ 2 , ] , 2 ) [ 1 ]
  up <- round( confint(fit)[ 2 , ] , 2 ) [ 2 ]
  pval <- summary(fit)$coefficients[ 2 , 4 ] 
  
  c( coef , se , low , up , pval )
}