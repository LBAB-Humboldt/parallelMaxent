#Run Bioclim models

#Run bioclim
library(dismo)
library(devtools)

# occs<-read.csv("C:/Users/GIC 14/Documents/Humboldt/Talleres/Zamias/Datos/ZM_BD.csv",as.is=T)
# covs<-stack(paste0("D:/Datos/worldclim/aoi/bio_",1:19,".asc"))
# 
# source_url("https://raw.githubusercontent.com/LBAB-Humboldt/parallelMaxent/master/preModelingFunctions.R")
# source_url("https://raw.githubusercontent.com/LBAB-Humboldt/parallelMaxent/master/evaluationFunctions.R")
# source_url("https://raw.githubusercontent.com/LBAB-Humboldt/parallelMaxent/master/postModelingFunctions.R")
# 
# s<-sampleRaster(covs,10000)
# 
# reduceCorrelation <- function(s, cor.threshold){
#   cor.s <- cor(s)
#   diag(cor.s) <- NA
#   i<-1
#   while (i <= nrow(cor.s)){
#     ind <- which(abs(cor.s[,i]) > cor.threshold)
#     if(length(ind) == 0){
#       i <- i+1
#     } else {
#       cor.s <- cor.s[-ind, -ind]
#       i <- i+1
#     }
#   }
#   return(list(ind = match(colnames(cor.s),colnames(s)), cor.s = cor.s))
# }
# 
# cor.s <- reduceCorrelation(s, 0.7)
# covs.red <- covs[[cor.s$ind]]
# rm(list=ls())
###
occ.file="C:/Users/GIC 14/Documents/Humboldt/Talleres/Zamias/Datos/ZM_BD.csv"
env.dir="D:/Datos/worldclim/aoi"
env.files=c(paste0("bio_",c(1,2,3,4,12,15,18),".asc"))
dist=1000*sqrt(2)
bkg.aoi = "extent"
bkg.type="random"
n.bkg = 10000 
sample.bkg = NULL
wd="C:/Users/GIC 14/Documents/Humboldt/Talleres/Zamias/Modelos/Bioclim"
do.eval=TRUE
do.threshold=TRUE
raw.threshold=c(0,10,20,30)
do.cut=TRUE



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

## Define list of species with more than 10 records
sp.list <- FilterSpeciesByRecords(occs, 5)
if(length(sp.list)==0){
  return()
}

current.spp <- length(sp.list)
cat(paste(Sys.time(), "After removing species with less than 10 unique records",
          current.spp, "species remain \n"))

 