mxParallel<-function(){
  source("sdm-functions.R") #Load auxilliary functions
  pkgs<-c("dismo","maptools","sp","rJava","rgdal","spatstat","reshape2",
          "SDMTools","snowfall","raster","svDialogs")
  lapply(pkgs,loadLibraries) #Load required libraries
  
  #Determine working directory
  wd=dlgDir(title="Seleccione carpeta archivos salida")$res
  
  #Extract environmental data
  ruta=dlgDir(title="Seleccione carpeta con variables continuas")$res #Seleccione ruta de variables ambientales
  envVars<-dlgList(list.files(ruta,pattern="*.asc$"),multiple=TRUE,
                   title="Seleccione variables continuas")$res
  
  envVars<-stack(paste(ruta,envVars,sep="/"))
  if(projection(envVars)=="NA"){
    dlgMessage("Undefined environmental variables projection
               Setting projection to geographic",  type = "ok")
    projection(envVars)<-"+proj=longlat +ellps=WGS84 +datum=WGS84"
  }
  
  #Load species records
  spFile<-dlgOpen(title = "Select CSV species file", filters = "*.csv")$res
  occs=read.csv(spFile,h=T)
  
  #Count unique records and eliminate species with less than 10 records
  maskRaster<-calc(envVars, fun=function(x){
    if(sum(is.na(x))==0){return(1)} 
    else {return(NA)}})
  
  countOccs<-countUnique(occs,maskRaster)
  
  exclSp<-names(which(countOccs < 10))
  dlgMessage(paste0("Eliminating ", length(exclSp)," species with < 10 unique records out of ",
                    length(countOccs), " species"),type = "ok")
  
  occs <- occs[-(unlist(sapply(exclSp,function(x) which(x==occs[,1])))),]
  
  #Select species to model
  spList=dlgList(sort(unique(occs[,1])),multiple=TRUE)$res
  
  #Feature selection options
  features=list("autofeature","linear","quadratic","hinge","threshold","product")
  selFeats<-dlgList(features,multiple=TRUE,title="Select Maxent features")$res
  logFeats<-!is.na(match(features, selFeats))
  mxntArgs<-paste(features,as.character(logFeats),sep="=")
  
  #Extrapolation options
  projPred=c("none","extrapolate","doclamp")
  extOpts=list("extrapolate","doclamp")
  selExtOpts<-dlgList(projPred,multiple=TRUE,title="Select extrapolation options")$res
  if(!("none"%in%selExtOpts)){
    logExtOpts<-!is.na(match(extOpts, selExtOpts))
    mxntArgs<-c(mxntArgs,paste(extOpts,as.character(logExtOpts),sep="="))
  } else {
    mxntArgs<-c(mxntArgs,paste(extOpts,as.character(FALSE),sep="="))
  }
  
  #Threshold options
  thresOpts<-c("Minimum training presence",
               "Ten percentile training presence",
               "Maximum Sensitivity + Specificity",
               "Equal Sensitivity and Specificity",
               "Custom percentile(s)")
  
  selThres<-dlgList(thresOpts,multiple=TRUE,title="Select threshold(s)")$res
  if(("Custom percentile(s)" %in% selThres)){
    choiceThres<-dlgInput(message = "Enter threshold(s) based on percentiles (0-100) from training presence data, separated by commas")$res
    choiceThres<-as.numeric(strsplit(choiceThres,",")[[1]])
  } else {
    choiceThres<-as.character(match(selThres,thresOpts[1:4]))
  }
  
  
  #Do model evaluations?
  selEval=dlgMessage("Do model evaluation (will increase computation time)?",  type = c("yesno"))$res
  selEval <- selEval == "yes"
  
  #Cut models by patch rule?
  selCut <- dlgMessage("Cut models by patch rule?",  type = c("yesno"))$res
  selCut <- selCut == "yes"
  
  #Select area of interest. This will be used for model building
  shapeFile=NULL
  aoi<-dlgList(c("Raster Extent","Convex Hull","Regions"),multiple=FALSE,
               title="Set the area of interest")$res
  
  #Choose region file if AOI is defined by a shapefile
  shapeFile=NULL
  if(aoi=="Regions"){
    shapeFile<-dlgOpen(title = "Select AOI polygon", filters = "*.shp")$res
    inShape<-readShapePoly(shapeFile)
    fieldID=dlgList(colnames(inShape@data),
                    multiple=FALSE,title="Select field that defines regions")$res
  }
  
  #Choose background selection method
  bkgMethod<-dlgList(c("Random","Samples"),multiple=FALSE,
                     title="Choose background selection method")$res
  
  #Choose samples file if samples are defined by a file
  samples=NULL
  if(bkgMethod=="Samples"){
    samples=read.csv(dlgOpen(title = "Select CSV samples file", filters = "*.csv")$res,h=T)
  }
  
  #Get environmental data only if random background
  if(aoi=="Raster Extent"){
    occCovs_ran<-extract(envVars,occs[,2:3]) #Extract all environmental info
    if(bkgMethod=="Random"){
      bg_ran_xy<-randomPoints(maskRaster,10000,p=cbind(occs[,2],occs[,3])) #Select background XY that don't fall on sampled cells. Those samples will be added at the modeling step.
      bg_ran_covs<-extract(envVars,bg_ran_xy)#Extract background covariates
    } else {
      bg_ran_covs<-na.omit(extract(envVars,samples))#Extract background covariates
    }
  }
  
  #Correr funciones de maxent en paralelo
  
  nCPU=8 #Para un computador con 4 cores
  sfInit(parallel=T,cpus=nCPU)#Initialize nodes
  sfExportAll() #Export vars to all the nodes
  sfLibrary(dismo)
  sfLibrary(raster)
  sfLibrary(rJava)
  sfLibrary(rgdal)
  sfLibrary(SDMTools)
  sfClusterSetupRNG()
  
  #Run MAXENT models in parallel
  
  beg<-seq(1,length(spList),by=nCPU)
  if(length(spList)<nCPU){
    fin<-length(spList)
  } else {
    fin<-seq(nCPU,length(spList),by=nCPU)
    if (length(beg)!=length(fin)){
      fin<-c(fin,beg[length(beg)])
    }
    fin[length(fin)]<-length(spList)
  }

  for (iter in 1:length(beg)){
    print(iter)
    maxent.pc=sfLapply(beg[iter]:fin[iter],function(j){
      ind<-which(occs[,1]==spList[j])
      
      #Get modeling data for background from the entire study area
      if(aoi=="Raster Extent"){
        occCovs<- occCovs_ran[ind,]#Get occurrence environmental info
        bg_covs<- bg_ran_covs#Get background covariates
      }
      
      #Get modeling data for background from regions defined by a shapefile
      if(aoi=="Regions"){
        bkg<-createBkg(occ,method="regions",inShape,fieldID,envVars[[1]])
        tmpVars<-stack(bkg,envVars)
        occCovs<-extract(tmpVars,occs[ind,2:3]) #Extract all environmental info
        if(bkgMethod=="Random"){
          bg_ran_xy<-randomPoints(maskRaster,10000,p=occ[,2:3]) #Select background XY that don't fall on sampled cells. Those samples will be added at the modeling step.
          bg_covs<-extract(tmpVars,bg_ran_xy)#Extract background covariates
        } else {
          bg_covs<-na.omit(extract(tmpVars,samples))#Extract background covariates
        }
      } 
      
      if(aoi=="Convex Hull"){
        bkg<-createBkg(occ,method="ch",aoi=envVars[[1]])
        tmpVars<-stack(bkg,envVars)
        occCovs<-extract(tmpVars,occs[,2:3]) #Extract all environmental info
        if(bkgMethod=="Random"){
          bg_ran_xy<-randomPoints(maskRaster,10000,p=occ[,2:3]) #Select background XY that don't fall on sampled cells. Those samples will be added at the modeling step.
          bg_covs<-extract(tmpVars,bg_ran_xy)#Extract background covariates
        } else {
          bg_covs<-na.omit(extract(tmpVars,samples))#Extract background covariates
        }
      } 
      
      dir.create(paste(wd,"/core",j,sep=""))
      
      mxModel(occCovs,bg_covs,mxntArgs,paste(wd,"/core",j,sep=""),doSave=TRUE,doEval=selEval,spList[j])
      mxPredict(predictors=envVars,thresholds=TRUE,choiceThres,doCut=selCut,doWrite1=TRUE,doWrite2=TRUE,
                spPoints=occs[ind,2:3],filePath=paste0(wd,"/core",j,"/",spList[j],".RData"),rootname=paste0(wd,"/",spList[j]))
    })
  }
  
  sfStop()
  return(list(outputDirectory=wd,rastersPath=ruta,rasters=names(envVars),speciesFile=spFile,
              speciesModeled=spList,features=selFeats,extrapolate=selExtOpts,thresholds=selThres,
              modelEval=selEval,modelCut=selCut,areaOfInterest=aoi,aoiShape=shapeFile,
              backgroundMethod=bkgMethod,bkgSampleFile=samples))
}

