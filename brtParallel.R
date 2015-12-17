
brtParallel<-function(occ.file,env.dir,env.files,dist=1000,bkg.aoi,bkg.type,
                      n.bkg,sample.bkg,wd,folds,do.threshold,raw.threshold,do.cut,brt.params,n.cpu,do.eval){
  #Create log file
  sink(paste0(wd,"/log.txt"))
  on.exit(sink())
  
  #Load Functions
  library(devtools)
  source_url("https://raw.githubusercontent.com/LBAB-Humboldt/parallelMaxent/master/preModelingFunctions.R")
  source_url("https://raw.githubusercontent.com/LBAB-Humboldt/parallelMaxent/master/evaluationFunctions.R")
  source_url("https://raw.githubusercontent.com/LBAB-Humboldt/parallelMaxent/master/postModelingFunctions.R")
  LoadLibraries()
  cat(paste(Sys.time(), "Functions and libraries loaded\n"))
  
  #Load and clean data
  occs <- LoadOccs(occ.file)
  current.spp <- length(unique(occs$species))
  current.recs <- nrow(occs)
  cat(paste(Sys.time(), "Loaded",current.recs,"records corresponding to",
            current.spp, "species from file", occ.file,"\n"))
  
  env.vars <- stack(paste0(env.dir,"/",env.files))
  cat(paste(Sys.time(), "Loaded environmental layers", paste(as.character(env.files), collapse=","), "from directory", env.dir,"\n"))
  
  if(is.na(projection(env.vars))|projection(env.vars)=="NA"){
    cat(paste(Sys.time(), "WARNING: Undefined environmental variables projection\n"))
    cat(paste(Sys.time(), "WARNING: Setting projection to geographic\n"))
    projection(env.vars)<-"+proj=longlat +ellps=WGS84 +datum=WGS84"
  }
  
  ## Remove records within radius defined by variable "dist"
  occs <- ddply(occs,.(species),IdNeighbors,dist=dist)
  current.spp <- length(unique(occs$species))
  current.recs <- nrow(occs)
  cat(paste(Sys.time(), "After removing  points within",dist, "meters of each other, ",
            current.recs, "records corresponding to", current.spp, "species remain\n"))
  
  
  
  #Extract covariate data for presences (and background if bkg.aoi="extent")
  occs.covs <- extract(env.vars, cbind(occs$lon,occs$lat))
  nna.rows <- which(apply(!is.na(occs.covs), 1, any))
  occs.covs <- occs.covs[nna.rows, ]
  occs <- occs[nna.rows, ]
  
  if (bkg.aoi == "extent"){
    train.bkg <- GenerateBkg(n.bkg, env.vars, bkg.type, sample.bkg)
    test.bkg <- GenerateBkg(n.bkg, env.vars, bkg.type, sample.bkg)
    cat(paste(Sys.time(), "Background generated for raster extent using",bkg.type, "sampling \n"))
  }
  
  ## Define list of species with more than 10 records
  sp.list <- FilterSpeciesByRecords(occs, 10)
  if(length(sp.list)==0){
    return()
  }
  current.spp <- length(sp.list)
  cat(paste(Sys.time(), "After removing species with less than 10 unique records",
            current.spp, "species remain \n"))
  
  cat(paste(Sys.time(), "Began parallel loop using", n.cpu, "cores \n"))
  sink()
  
  #Begin modeling loop  
  sfInit(parallel=T,cpus=n.cpu)#Initialize nodes
  sfExportAll() #Export vars to all the nodes
  sfClusterSetupRNG()
  sfClusterApplyLB(1:length(sp.list),function(i){
    tmp.dir <- tempdir()
    sink(paste0(wd,"/log.txt"), append=TRUE)
    on.exit(sink())
    
    LoadLibraries()
    #Get species data
    sp.name <- sp.list[i]
    sp.idx <- which(occs$species == sp.list[i])
    sp.occs <- occs[sp.idx, ]
    
    #Generate covariate data for background (when bkg.type!="extent")
    if(!exists("train.bkg")&!exists("test.bkg")){
      train.bkg <- GenerateSpBkg(sp.occs, n.bkg, env.vars, bkg.type, bkg.aoi, 
                                 regions, field, sample.bkg, buffer)
      test.bkg <- GenerateSpBkg(sp.occs, n.bkg, env.vars, bkg.type, bkg.aoi, 
                                regions, field, sample.bkg, buffer)
      cat(paste(Sys.time(), "Background generated for species", sp.name, 
                "using area defined by", bkg.aoi, "and", bkg.type, "sampling \n"))
    }
    
    
    #Do model evaluation
    if(do.eval){
      sp.eval <- EvaluateBRTModel(folds, covs.pres=occs.covs[sp.idx, ], covs.bkg.train=train.bkg, covs.bkg.test=test.bkg, brt.params)
      write.csv(sp.eval, paste0(wd, "/", sp.list[i],"_evaluation_brt.csv"), row.names=F)
      cat(paste(Sys.time(), "Performed model evaluation for", sp.name, "\n")) 
    }
    
    #Start modeling
    df <- data.frame(c(rep(1,length(sp.idx)),rep(0,nrow(train.bkg))), rbind(occs.covs[sp.idx,],train.bkg))
    
    brt.obj <- gbm.step(data=df, gbm.x = 2:ncol(df), gbm.y = 1,
                        family = "bernoulli", tree.complexity = brt.params[1], 
                        learning.rate = min(sp.eval$lr), bag.fraction = brt.params[3],prev.stratify=FALSE,
                        site.weights=c(rep(1,length(sp.idx)), rep(length(sp.idx)/nrow(train.bkg), nrow(train.bkg))))
    
    save(brt.obj, file=paste0(wd, "/", sp.list[i], "_brt.RData"))
    cat(paste(Sys.time(), "Generated BRT distribution model for", sp.name, "\n"))
    map <- predict(env.vars, brt.obj, n.trees=brt.obj$gbm.call$best.trees, type="response")
    
    writeRaster(map, paste0(wd, "/", sp.list[i], "_brt.tif"), format="GTiff",
                overwrite=TRUE, NAflag=-9999)
    cat(paste(Sys.time(), "Generated prediction of BRT distribution model for", sp.name, "\n"))
    
    write.csv(occs[sp.idx, ], paste0(wd, "/", sp.list[i], "_brt.csv"), row.names=FALSE)
    
    #Post-processing: threshold & cut
    if(do.threshold){
      thres.maps <- sapply(raw.threshold, FUN=ThresholdBRT, brt.obj=brt.obj, sp.covs=occs.covs[sp.idx, ],
                           map=map)
      for(j in 1:length(raw.threshold)){
        writeRaster(thres.maps[[j]],filename=paste0(wd, "/", sp.name,"_", raw.threshold[j], "_brt.tif"), 
                    format="GTiff",overwrite=TRUE, NAflag=-9999)
      }
      cat(paste(Sys.time(), "Generated thresholded prediction of BRT distribution model
                  using thresholds ", paste(raw.threshold,collapse=", "), "for", sp.name, "\n"))
      if(do.cut){
        cut.maps <- sapply(thres.maps, FUN=CutModel2, sp.points=cbind(sp.occs$lon,sp.occs$lat))
        for(j in 1:length(raw.threshold)){
          writeRaster(cut.maps[[j]],filename=paste0(wd, "/", sp.name,"_",raw.threshold[j], "_cut_brt.tif"), 
                      format="GTiff",overwrite=TRUE, NAflag=-9999)
        }
        cat(paste(Sys.time(), "Cut thresholded prediction(s) of BRT distribution model for", sp.name, "\n"))
      }
    }
    
    #Remove temporary files
    removeTmpFiles(2)
  })
  sfStop()
}