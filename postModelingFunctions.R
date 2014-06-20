Threshold2 <- function(raw.threshold, mxnt.obj, map, sp.occs){
  #Use all default maxent thresholds
  tnames.long <- c("Minimum.training.presence.logistic.threshold",
                   "X10.percentile.training.presence.logistic.threshold",
                   "Equal.training.sensitivity.and.specificity.logistic.threshold",
                   "Maximum.training.sensitivity.plus.specificity.logistic.threshold")
  tnames <- c("min","10p","ess","mss") #thresholds: Minimum training presence,10 percentile training presence,equal specificity and sensitivity,maximum specificity and sensitivity
  
  if(is.numeric(raw.threshold)){
    preds <- predict(mxnt.obj, mxnt.obj@presence)
    thresholds <- quantile(preds, raw.threshold / 100)
    tsuffix <- as.character(raw.threshold)
  }
  
  if(is.character(raw.threshold)){
    thres.idx <- match(raw.threshold, tnames)
    thresholds <- mxnt.obj@results[tnames.long[thres.idx], ]
    tsuffix <- tnames[thres.idx]
  }
  
  out.name <- paste("map_", tsuffix, sep="")
  assign(out.name, (map >= thresholds))
  return(get(out.name))
}
  
CutModel2 <- function(map, sp.points){
  tmp.mask <- map >= 0
  map[map==0] <- NA
  map.patch <- ConnCompLabel(map)
  pts.patch <- extract(map.patch, sp.points)
  pts.patch <- unique(pts.patch)
  pts.patch <- pts.patch[which(!is.na(pts.patch))]
  map.cut <- map.patch %in% pts.patch
  map.cut[is.na(tmp.mask)] <- NA
  return(map.cut)
}
