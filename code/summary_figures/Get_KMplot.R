suppressMessages(library(survcomp))
suppressMessages(library(tidyverse))


Get_KMplot = function(cancer_type , status , time , score , data , dir ){
  
  status <- as.numeric(status)
  time <- as.numeric(time)
  score <- as.numeric(score)
  
  cutoff <- median(as.numeric(score))
  flag <- lapply(score, function(x){ifelse(x >= cutoff, 1, 0)}) %>% as.numeric()
  pdf( dir,bg="transparent", onefile=FALSE, family="Times")
  km.coxph.plot( formula.s=Surv( time, status ) ~ flag ,data.s=data, x.label="Time (Days)", y.label="Overall Survival", main.title=paste( cancer_type , "\n(cutoff=" , round( cutoff , 2 ) , ")" , sep="" ) ,sub.title="", 
                 leg.text=c( "Low" , "High"), 
                 leg.pos="topright", .col=c( "black","#e53935"),  show.n.risk=TRUE, n.risk.step=400, n.risk.cex=0.85, ylim=c(0,1), leg.inset=0,.lwd=3 , verbose=FALSE )
  dev.off()
}
