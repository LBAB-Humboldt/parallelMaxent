ThresholdBRT<-function(raw.threshold, brt.obj, sp.covs, map){
  #Use all default maxent thresholds
  
  if(is.numeric(raw.threshold)){
    preds<-predict(brt.obj, as.data.frame(sp.covs), n.trees=brt.obj$gbm.call$best.trees, type="response")
    thresholds <- quantile(preds, raw.threshold / 100,na.rm=T)
    tsuffix <- as.character(raw.threshold)
  }
  
  out.name <- paste("map_", tsuffix, sep="")
  assign(out.name, (map >= thresholds))
  return(get(out.name))
}