plot.subject = function(model, subject, value.vars ){
  sr = sqrt(length(value.vars))
 # pdf( 'traj_plot.pdf')
  #par(mfrow = c(2,2))
  
  for (value.var in value.vars){
    ylim = c(min(data[[value.var]],na.rm = TRUE), max(data[[value.var]],na.rm = TRUE))
    
    idx = (model$data_long[["idx"]][model$slices[[value.var]]] == model$subjects[subject])
    
    plot(model$grid, model$functions[[value.var]][subject,],t='l',ylim=ylim, col=3,xlab = "age",ylab = "value")
    lines(model$mean_curves[[value.var]],lwd=2)
    points(model$data_long[[time.var]][model$slices[[value.var]]][idx], model$data_long[["value"]][model$slices[[value.var]]][idx], col="grey50")
    title(value.var)
  }
  #par(mfrow = c(1, 1))
  #title(paste("Subject", subject), line = -1, outer = TRUE)
  
  #dev.off()
}

#### select protein and GMV
prot.vars = c ( grep( 'chem', value.vars ), grep( 'cyt', value.vars), grep( 'GMV_TIV', value.vars ) )
prot.vars <- value.vars[ prot.vars ]


pdf( 'GMV_TIV_out.pdf')
  par( mfrow = c( 2,2 ))
for ( i in 1:8 ){  
  plot.subject(model,i, value.vars = 'GMV_TIV_out' )
}
dev.off()


#extract trajectory for GMV_TIV_out
k <- which( value.vars == 'GMV_TIV_out' )
GMV_beta <- model$W[, c( 2*k-1, 2*k) ]




plot.subject(model,100, value.vars = prot.vars)
plot.subject(model,3)
plot.subject(model,4)
plot.subject(model,5)
plot.subject(model,6)
plot.subject(model,7)
