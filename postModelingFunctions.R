# postModelingFunctions.R
# Set of functions to process species distribution models after a model has been generated
# Author: Jorge Velásquez

#Threshold2.R
#Function to threshold continuous species distribution models
#Arguments
## raw.threshold(character or numeric vector): vector of thresholds to be applied to distribution models. 
#                 It can either be a character (min, 10p, ess, mss) or numeric,
#                 in which case it represents the percentile of training presence
#                 probabilities.
## mxn.obj(MaxEnt): maxent object
## map(raster): raster object corresponding to the projection of mxnt.obj into environmental
##      space with continuous values from 0 to 1.
#Returns:
#  A raster object for each threshold used.

Threshold2 <- function(raw.threshold, mxnt.obj, map){
  #Use all default maxent thresholds
  tnames.long <- c("Minimum.training.presence.logistic.threshold",
                   "X10.percentile.training.presence.logistic.threshold",
                   "Equal.training.sensitivity.and.specificity.logistic.threshold",
                   "Maximum.training.sensitivity.plus.specificity.logistic.threshold")
  tnames <- c("min","10p","ess","mss") #thresholds: Minimum training presence,10 percentile training presence,equal specificity and sensitivity,maximum specificity and sensitivity
  
  if(is.numeric(raw.threshold)){
    preds <- predict(mxnt.obj, mxnt.obj@presence)
    thresholds <- quantile(preds, raw.threshold / 100,na.rm=T)
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

#CutModel2
#This functions allows implementation of a patch rule to avoid overprediction in 
#thresholded species distribution models. For a given thresholded distribution model
#this funcion will return a model in which only distribution patches with evidence
#of being occupied are selected.

#Arguments:
## map(raster): raster object of presence/absence species distribution model
## sp.points(data frame or matrix): two-column matrix or data.frame, or SpatialPoints with locations of species
#             occurrence.
#Returns:
##  A raster object with distribution patches without evidence of occurrence deleted.

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

##ThresholdBRT
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
