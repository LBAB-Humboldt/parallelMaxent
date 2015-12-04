#preModelingFunctions.R
#Set of functions to prepare data and environment to run species distribution models
#Author: Jorge Velásquez

#LoadLibraries
#Loads all libraries needed to run functions called by mxParallel
#Arguments: 
# memory(string): string specifying how much memmory allocate to rJava, e.g. "512m"
#  for 512 megabytes, "2g" for 2 gigabytes.

LoadLibraries<-function(memory="4g"){
  options(java.parameters = paste0("-Xmx",memory))
  require("gbm")
  require("devtools")
  require("dismo")
  require("maptools")
  require("plyr")
  require("raster")
  require("reshape2")
  require("rJava")
  require("rgdal")
  require("rgeos")
  require("SDMTools")
  require("sp")
  require("spatstat")
  require("snowfall")
  require("rlecuyer")
}

#LoadOccs
#Loads and verifies species occurrence file
##Arguments:
# occ.file(string): path to comma separated file of occurrences. Must have fields id, species, lon, lat
##Returns:
# data frame of occurrences.

LoadOccs<-function(occ.file){
  if(is.data.frame(occ.file)){
    occs <- occ.file
  } else {
    occs <- read.csv(occ.file,as.is=T)
  }

  with(occs, if(nrow(occs)==0){
    stop("Occurrence file has 0 rows")
  })
  
  with(occs, if(!exists("id")){
    stop("Variable id is missing from occurrence file")
  })
  
  with(occs, if(!exists("species")){
    stop("Variable species is missing from occurrence file")
  })
  
  with(occs, if(!exists("lon")){
    stop("Variable lon is missing from occurrence file")
  })
  
  with(occs, if(!exists("lat")){
    stop("Variable lat is missing from occurrence file")
  })

  lon.errors <- with(occs, which((lon > 180)|(lon < -180)))
  lat.errors <- with(occs, which((lat > 90)|(lat < -90)))
  lon.na <- with(occs, which(is.na(lon)))
  lat.na <- with(occs, which(is.na(lon)))
  rem.idx <- unique(c(lon.errors,lat.errors,lon.na,lat.na))
  if(length(rem.idx) > 0){
    occs<-occs[-rem.idx, ]
    message(paste0("Removed ",length(lon.errors)," longitude values >180 or < -180\n",
                   "Removed ",length(lat.errors)," latitude values >90 or < -90\n",
                   "Removed ",length(lon.na)," NA longitude values\n",
                   "Removed ",length(lat.na)," NA latitude values"))
  }
  if(nrow(occs)==0){
    stop("Occurrence file has 0 rows")
  }
  return(occs)
}

#CleanOccs
#Extracts environmental values associated with occurrences after removing duplicates
#and eliminating records at a particular distance
##Arguments:
##  occs(data frame): data frame object of occurrences
##  env.vars(raster): raster or stack of environmental variables from which to extract
##   environmental values. 
##  dist(numeric): distance below which two coordinates are considered a duplicate.

CleanOccs<-function(occs,env.vars,dist){
  occs <- ddply(occs,.(species),IdNeighbors,dist=1000) #Apply the IdNeighbors function to each species
  occs.covs <- extract(env.vars, cbind(occs$lon,occs$lat))
  return(list(occs=occs, occs.covs=occs.covs))
}

#IdNeighbors
#Identifies records below a specified threshold distance and deletes them from occurrence
#file.
#Arguments:
##  occs(data frame): data frame object of occurrences
##  dist(numeric): distance below which two coordinates are considered a duplicate.
##  longlat(logical): Are coordinates geographic?
#Returns:
##  data frame object of occurrences without duplicate coordinates.

IdNeighbors<-function(occs,dist,longlat=TRUE){
  if(nrow(occs)<2){
    return(occs)
  }
  coords <- cbind(occs$lon,occs$lat)
  dst <- pointDistance(coords,longlat=longlat)
  diag(dst) <- NA
  rmv.idx <- which(dst < dist,arr.ind=T)
  if(nrow(rmv.idx)==0){
    return(occs)
  } else {
    occs <- occs[-rmv.idx[, 1], ]
    return(occs)
  }
}

#FilterSpeciesByRecords
#Create list of species with more than min.recs records
#Arguments:
##  occs(data frame): data frame object of occurrences
##  min.recs(numeric): minimum number of records to be included in list.
#Returns:
##  character vector of species with more than min.recs records

FilterSpeciesByRecords <- function(occs, min.recs){
  df <- ddply(occs,"species",summarise,N=length(species))
  sp.list <- df$species[which(df$N >= min.recs)]
  if(length(sp.list) == 0){
    stop(paste0("None of the species in occurrence file has more than ", min.recs, " records"))
  } else {
    return(sp.list)
  }
}

#GenerateBkg
#Generates background to use in species distribution models
#Arguments:
##  n (numeric): size of background dataset
##  env.vars(raster): raster or stack of environmental variables from which background
##   will be extracted.
##  bkg.type(string): keyword that defines how the background will be sampled from bkg.aoi.
##    random (default): background will be sampled randomly from bkg.aoi 
##    samples: get background samples from a file.
##  sample.bkg(string): data frame (should better be a csv??) with coordinates lon lat for
##   each record.Background is defined by this data frame.
#Returns:
##  data frame with environmental conditions associated with background.

GenerateBkg <- function(n, env.vars, bkg.type="random", sample.bkg=NULL){
  if(bkg.type == "random"){
    bkg.covs <- sampleRaster(env.vars, n)
  } 
  if(bkg.type == "samples"){
    sample.coords <- cbind(sample.bkg$lon, sample.bkg$lat)
    bkg.covs <- extract(env.vars, sample.coords)
  }
  return(as.data.frame(bkg.covs))
}

#GenerateSpBkg
#Generates species-specific psudoabsences or background
#Arguments:
##  occs(data frame or matrix): 2-column matrix or data frame of species' occurrences.
##  n(numeric): size of background dataset
##  env.vars(raster): raster or stack of environmental variables from which background
##   will be extracted.
##  bkg.type(string): keyword that defines how the background will be sampled from bkg.aoi.
##    random (default): background will be sampled randomly from bkg.aoi 
##    samples: get background samples from a file.
##  bkg.aoi(string): keyword that defines where background will be sampled from. 
##    regions: background will be species specific, and it will correspond
##            to the polygons of a shapefile that completely contain the
##            species records.
##    ch: background will be species specific and it will correspond to the
##       convex hull of the species' records.
##  regions(SpatialPolygons): SpatialPolygons object with the regions that will be used to
##                            define species background.
##  field(string):            field (column name) that defines the regions.
##                            Used only when bkg.aoi="regions"
##  sample.bkg(data frame): data frame (should better be a csv??) with coordinates lon lat for
##              each record.Background is defined by this data frame.
##  buffer(numeric): Buffer in meters to be applied to convex polygons.
##                   Used only when bkg.aoi="ch".
#Returns:
##  data frame with environmental conditions associated with background.

GenerateSpBkg <- function(occs, n, env.vars, bkg.type="random", bkg.aoi, 
                          regions, field, sample.bkg=NULL, buffer=NULL){
  bkg <- CreateAOI(occs, method=bkg.aoi, env.vars, regions, field, buffer)
  tmp.stack <- stack(bkg, env.vars)
  if(bkg.type == "random"){
    bkg.covs <- sampleRaster(tmp.stack, n)[, 2:(nlayers(env.vars)+1)]
  } 
  if(bkg.type == "samples") {
    if(is.null(sample.bkg)){
      stop("Missing target background samples file")
    }
    with(sample.bkg, if(!exists("lon")){
      stop("Variable lon is missing from occurrence file")
    })
    
    with(sample.bkg, if(!exists("lat")){
      stop("Variable lat is missing from occurrence file")
    })
    bkg.covs <- na.omit(extract(tmp.stack, cbind(sample.bkg$lon,sample.bkg$lat)))[, 2:(nlayers(env.vars)+1)]
  }
  return(list(bkg.aoi=bkg,bkg.covs=as.data.frame(bkg.covs)))
}

#CreateAOI
#Create raster of area of interest for modeling
#Arguments:
##  occs(matrix or data frame): 2-column matrix or data frame of species' occurrences.
##  method(string): either regions or ch depending on whether the area of interest is 
##   defined by the polygons of a shapefile that contain species occurrences or by a
##   convex hull from occurrences
##  aoi(raster): a raster object with the extent, resolution and projection of the area of interest.
##  regions(SpatialPolygons): SpatialPolygons object with the regions that will be used to
##                            define species background.
##  field(string): field (column name) that defines the regions.
##  buffer(numeric): buffer in meters to be applied to convex polygons.
#Returns:
# A raster object with area of interest for modeling.

CreateAOI<-function(occs, method, aoi, regions, field, buffer){
  in.pts <- SpatialPoints(cbind(occs$lon, occs$lat), proj4string = CRS(projection(aoi)))
  if(method == "regions"){
    if(missing(regions)){
      stop("Missing regions argument")
    }
    if(missing(field)){
      stop("Missing field argument")
    }
    proj4string(regions) <- CRS(projection(aoi))
    units <- na.omit(unique(over(in.pts, regions)[, field]))
    ind <- which(regions@data[,field] %in% units)
    bkg.shp <- regions[ind,]
    bkg <- rasterize(bkg.shp, aoi, field=1)
    return(bkg)
  } 
  if(method=="ch"){
    ch.shp <- convHull(in.pts)@polygons
    if(!is.null(buffer)){
      if(buffer>0){
        ch.shp <- gBuffer(ch.shp, width=buffer)
      }
    } 
    bkg <- rasterize(ch.shp, aoi, field=1)
    bkg <- bkg * (!is.na(aoi[[1]]))
    return(bkg)
  }
}

#sampleRaster
#Function to sample n points randomly from a raster object
#Arguments:
## raster.obj(raster): raster object to sample coordinates pairs from
## n(numeric): number of coordinate pairs to sample
#Returns:
## data frame of sampled coordinates

sampleRaster<-function(raster.obj,n){
  if(nlayers(raster.obj)>1){
    mask <- sum(raster.obj)
    cells <- Which(!is.na(mask),cells=T)
  } else {
    cells <- Which(!is.na(raster.obj),cells=T)
  }
  if(length(cells)<n){
    n <- length(cells)
    warning("n value exceeds the number of cells with data")
  }
  sel.cells <- sample(cells, n)
  output <- raster.obj[sel.cells]
  return(output)
}
