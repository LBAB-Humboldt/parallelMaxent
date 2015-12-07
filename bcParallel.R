#bcParallel
# Prepares data for species distribution modeling in Bioclim, and runs
# distribution models in parallel using Bioclim

# Arguments:
## Data input, cleaning and output (mandatory):
### occ.file(string):  Full path to species occurrence file. 
###                    Occurrence file must be comma separated and must contain fields
###                    id, species, lat, lon (lowercase, order not important)
### env.dir(string):   Path to directory that contains environmental layers
### env.files(string): File names of environmental layers to be used in distribution models.
###                    Must contain extension, e.g. bio_1.tif.
### wd(string):        Path to directory where output files will be saved
### dist(numeric):     Distance (in meters) within which records are considered duplicates.
###                    Only one record is kept within dist.Default is 1000.
##
## Background generation options:
### bkg.aoi(string):   Keyword that defines where background will be sampled from. 
###                      extent (default): background will be sampled from raster extent
###                      regions: background will be species specific, and it will correspond
###                               to the polygons of a shapefile that completely contain the
###                               species records.
###                      ch: background will be species specific and it will correspond to the
###                          convex hull of the species' records.
### bkg.type(string):  Keyword that defines how the background will be sampled from bkg.aoi.
###                    random (default): background will be sampled randomly from bkg.aoi
###                    samples: get background samples from a file.
### Optional arguments:
###   n.bkg(numeric):           number of background samples. 
###                             Used when bkg.type="random". Default is 10000.
###   sample.bkg(string):       Path to comma separated file containing background samples.
###                             Must include the fields lat and lon (lowercase, order doesn't matter).
###                             Used only when bkg.type="samples"
###   regions(SpatialPolygons): SpatialPolygons object with the regions that will be used to
###                             define species background.
###                             Used only when bkg.aoi="regions"
###   field(string):            field (column name) that defines the regions.
###                             Used only when bkg.aoi="regions"
###   buffer(numeric):          Buffer in meters to be applied to convex polygons.
###                             Used only when bkg.aoi="ch".
## Evaluation aoptions:
### do.eval(logical):   Do model evaluation? Default TRUE.
### folds(numeric):     Number of folds for k-fold partitioning used in evaluation 
###                     and regularization optimization. Default = 10.
## Modeling options:
### n.cpu:               Number of cores to uses for parallel processing
## Post-processing options:
### do.threshold(logical): Threshold distribution models?
### raw.threshold(vector): numeric or character vector. If numeric, this will
###                        specify the percentiles (from 0 to 100) at which
###                        models should be thresholded according to the 
###                        "probability of occurrence" at training sites.
###                        If character, this should be any combination of the
###                        following keywords: min, 10p, ess, mss.
### do.cut(logical):       Select distribution patches with evidence of occurrence
###                        from thresholded distribution models?
# 
#                         
# Example:
#
# bcParallel(occ.file="C:/Modelos/Primates/Primates_05052014.csv",
#   env.dir="D:/Datos/worldclim/aoi",
#   env.files=c(paste0("bio_",c(1,2,3,4,12,15,18),".asc")),
#   dist=1000,
#   bkg.aoi = "extent",
#   bkg.type="random",
#   n.bkg = 10000,
#   sample.bkg = NULL,
#   folds=5,
#   wd="C:/Workspace",
#   do.eval=TRUE,
#   n.cpu=4,
#   do.threshold=TRUE,
#   raw.threshold=c(0,10,20,30),
#   do.cut=TRUE)

bcParallel <- function(occ.file, env.dir, env.files, dist, bkg.aoi, bkg.type,
                       n.bkg, sample.bkg, folds, wd, do.eval, n.cpu, do.threshold, raw.threshold,
                       do.cut){

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

  #Extract covariate data for presences
  occs.covs <- extract(env.vars, cbind(occs$lon,occs$lat))
  nna.rows <- which(apply(!is.na(occs.covs), 1, any))
  occs.covs <- occs.covs[nna.rows, ]
  occs <- occs[nna.rows, ]

  if (bkg.aoi == "extent"){
    train.bkg <- GenerateBkg(n.bkg, env.vars, bkg.type, sample.bkg)
    test.bkg <- GenerateBkg(n.bkg, env.vars, bkg.type, sample.bkg)
    cat(paste(Sys.time(), "Background generated for raster extent using",bkg.type, "sampling \n"))
  }

  ## Define list of species with 5 or more records
  sp.list <- FilterSpeciesByRecords(occs, 5)
  if(length(sp.list)==0){
    return()
  }
  
  current.spp <- length(sp.list)
  cat(paste(Sys.time(), "After removing species with less than 5 unique records",
            current.spp, "species remain \n"))
  
  cat(paste(Sys.time(), "Began parallel loop using", n.cpu, "cores \n"))
  sink()
  
  #Begin modeling loop  
  sfInit(parallel=T,cpus=n.cpu)#Initialize nodes
  sfExportAll() #Export vars to all the nodes
  sfClusterSetupRNG()
  sfClusterApplyLB(1:length(sp.list),function(i){
    print(sp.list[i])
    tmp.dir <- tempdir()
    sink(paste0(wd,"/log.txt"), append=TRUE)
    on.exit(sink())
    LoadLibraries()
    
    #Get species data
    sp.name <- sp.list[i]
    sp.idx <- which(occs$species == sp.list[i])
    sp.occs <- occs[sp.idx, ]
    
    #Do model evaluation
    if(do.eval){
      sp.eval <- EvaluateBioclimModel(folds, covs.pres=occs.covs[sp.idx, ], covs.bkg.train=train.bkg, covs.bkg.test=test.bkg)
      write.csv(sp.eval, paste0(wd, "/", sp.list[i],"_evaluation_bc.csv"), row.names=F)
      cat(paste(Sys.time(), "Performed model evaluation for", sp.name, "\n")) 
    }
    
    #Start modeling
    bc.obj <- bioclim(occs.covs[sp.idx, ])
      
    save(bc.obj, file=paste0(wd, "/", sp.list[i], "_bc.RData"))
    cat(paste(Sys.time(), "Generated Bioclim distribution model for", sp.name, "\n"))
    map <- predict(env.vars, bc.obj)
    
    writeRaster(map, paste0(wd, "/", sp.list[i], "_bc.tif"), format="GTiff",
                overwrite=TRUE, NAflag=-9999)
    cat(paste(Sys.time(), "Generated prediction of maxent distribution model for", sp.name, "\n"))
    
    write.csv(occs[sp.idx, ], paste0(wd, "/", sp.list[i], "_bc.csv"), row.names=FALSE)
    
    #Post-processing: threshold & cut
    if(do.threshold){
      thres.maps <- sapply(raw.threshold, FUN=ThresholdBRT, brt.obj=bc.obj, sp.covs=occs.covs[sp.idx, ],
                           map=map)
      for(j in 1:length(raw.threshold)){
        writeRaster(thres.maps[[j]],filename=paste0(wd, "/", sp.name,"_", raw.threshold[j], "_bc.tif"), 
                    format="GTiff",overwrite=TRUE, NAflag=-9999)
      }
      cat(paste(Sys.time(), "Generated thresholded prediction of Bioclim distribution model
                using thresholds ", paste(raw.threshold,collapse=", "), "for", sp.name, "\n"))
      if(do.cut){
        cut.maps <- sapply(thres.maps, FUN=CutModel2, sp.points=cbind(sp.occs$lon,sp.occs$lat))
        for(j in 1:length(raw.threshold)){
          writeRaster(cut.maps[[j]],filename=paste0(wd, "/", sp.name,"_",raw.threshold[j], "_cut_bc.tif"), 
                      format="GTiff",overwrite=TRUE, NAflag=-9999)
        }
        cat(paste(Sys.time(), "Cut thresholded prediction(s) of Bioclim distribution model for", sp.name, "\n"))
      }
    }
  
    #Remove temporary files
    removeTmpFiles(2)
  })
  sfStop()
}