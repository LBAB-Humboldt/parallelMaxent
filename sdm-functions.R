#AUXILIARY FUNCTIONS

#LOAD PACKAGES
loadLibraries<-function(pkg){
  loadPkg<-require(pkg,,character.only=TRUE)
  if(loadPkg){
    print(paste0(pkg," is loaded correctly"))
  } else {
    print(paste0("trying to install ",pkg))
    install.packages(pkg)
    if(require(pkg,character.only=TRUE)){
      print(paste0(pkg," installed and loaded"))
    } else {
      stop(paste0("could not install ",pkg))
    }
  }
}

mxModel<-function(occ,bg,mxntArgs,path,doSave=TRUE,doEval=FALSE,filename){
  mxnt.obj<-maxent(x=as.data.frame(rbind(occ,bg)),p=c(rep(1,nrow(occ)),rep(0,nrow(bg))),
                   removeDuplicates=FALSE,path=path,args=mxntArgs)
  ev.table=FALSE
  if(doSave){
    save(mxnt.obj,file=paste(path,"/",filename,".RData",sep=""))
    write.table(mxnt.obj@results,paste(path,"/",filename,"_results.csv",sep=""),sep=",")
  }
  if(doEval){
    folds=4
    Sample=kfold(occ, k=folds)
    ev.names<-c("np","na","auc","pauc")
    ev.table<-as.data.frame(matrix(NA,nrow=folds,ncol=4))
    colnames(ev.table)<-ev.names
    for (i in 1:folds){
      occtest=occ[Sample==i,]	# para evaluar el modelo y
      occtrain=occ[Sample!=i,]# para el aprendizaje o 'training'
      if(is.null(nrow(occtest))) occtest<-t(as.matrix(occtest)) #cuando solo existe un valor de testing forzar a que permanezca como matriz y no vector
      mxnt.obj_k<-maxent(x=as.data.frame(rbind(occtrain,bg)),p=c(rep(1,nrow(occtrain)),rep(0,nrow(bg))),removeDuplicates=FALSE,path=path,args=c("product=false","threshold=false"))
      ev=evaluate(model=mxnt.obj_k,p=occtest,a=bg)
      for(j in 1:length(ev.names)){
        evValue<-slot(ev,ev.names[j])
        if(length(evValue)==0){evValue=NA}
        ev.table[i,j]<-evValue
      }
    }
    write.table(ev.table,paste(path,"/",filename,"_kfoldEval.csv",sep=""),sep=",",row.names=FALSE)
  }
  return(list(ev.table,mxnt.obj))
}

#Add Function to compute other thresholds?

#Function to predict. Explore whether it makes sense to use another export format
mxPredict<-function(predictors,thresholds=NULL,choiceThres,doCut=FALSE,doWrite1=FALSE,doWrite2=FALSE,
                    spPoints,filePath,rootname){
  load(filePath)
  map<-predict(mxnt.obj,predictors,args=c("outputformat=logistic"),progress="text") #Argumentos features
  writeRaster(map,rootname,format="raster",overwrite=TRUE)
  if(thresholds){#Use all default maxent thresholds
    tnames_long<-c("Minimum.training.presence.logistic.threshold",
                   "X10.percentile.training.presence.logistic.threshold",
                   "Equal.training.sensitivity.and.specificity.logistic.threshold",
                   "Maximum.training.sensitivity.plus.specificity.logistic.threshold")
    tnames<-c("min","10p","ess","mss") #thresholds: Minimum training presence,10 percentile training presence,equal specificity and sensitivity,maximum specificity and sensitivity
    
    if(is.numeric(choiceThres)){
      preds<-predict(mxnt.obj,mxnt.obj@presence)
      thresholds<-quantile(preds,choiceThres / 100)
      tnames<-as.character(choiceThres)
    }
    if(is.character(choiceThres)){
      thresholds<-mxnt.obj@results[tnames_long[as.numeric(choiceThres)], ]
      tnames<-tnames[as.numeric(choiceThres)]
    }
    
    for(i in 1:length(thresholds)){
      outname<-paste("map_",tnames[i],sep="")
      assign(outname,(map>=thresholds[i]))
      if(doWrite1) writeRaster(get(outname),paste(rootname,tnames[i],sep="_"),format="raster",overwrite=TRUE)
    }
    
    #Cut models by patches with records
    if(doCut){
      for(j in 1:length(thresholds)){
        inName<-paste("map_",tnames[j],sep="")
        cutModel(get(inName),spPoints,doWrite=doWrite2,paste(rootname,tnames[j],"cut",sep="_"))
      }
    }
  }
}

cutModel<-function(map,spPoints,doWrite=FALSE,filename){
  map[map==0]<-NA
  map_patch<-ConnCompLabel(map)
  pts_patch<-extract(map_patch,spPoints)
  pts_patch<-unique(pts_patch)
  pts_patch<-pts_patch[which(!is.na(pts_patch))]
  map_cut<-map_patch %in% pts_patch
  if(doWrite) writeRaster(map_cut,filename,format="raster",overwrite=TRUE)
}

countUnique<-function(occData,maskRaster){
  indNotNA<-Which((!is.na(maskRaster)),cells=TRUE)
  maskRaster[indNotNA]<-1:length(indNotNA)
  cellValues<-extract(maskRaster,occData[,-1])
  uniqueRecs<-unique(data.frame(occData,cellid=cellValues))
  countRecs<-table(uniqueRecs[, 1])
  return(countRecs)
}

#This function creates a background in raster format by
#selecting polygons on a shapefile that intersect species
#records and clipping by an area of interest defined by 
#a raster layer.
#Arguments:
#recs: data frame with Longitude and Latitude fields for species records
#regions: shapefile with regions to use to define background (e.g ecoregions)
#field: shapefile field that identifies regions
#aoi: raster object that defines area of interest

#Note:
#Projection information is taken from AOI

createBkg<-function(recs,method,regions,field,aoi){
  inPts=SpatialPoints(cbind(recs$Longitude,recs$Latitude),proj4string=CRS(projection(aoi)))
  if(method=="regions"){
    proj4string(regions)=CRS(projection(aoi))
    units<-na.omit(unique(over(inPts,regions)[,field]))
    ind<-NULL
    for(unit in units){
      ind<-c(ind,which(regions@data[,field]==unit))
    }
    bkgShp<-regions[ind,]
    bkg<-rasterize(bkgShp,aoi,field=rep(1,length(ind)))
    return(bkg)
  } 
  if(method=="ch"){
    chObj<-convHull(inPts)
    bkg<-predict(chObj,aoi)
    bkg<-bkg*(!is.na(aoi))
    bkg[bkg==0]<-NA
    bkg[bkg==-1]<-1
    return(bkg)
  }
}