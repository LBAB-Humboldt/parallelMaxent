#convert2PNG.R
#
#This function creates a KMZ file, georeferenced PNG and a thumb PNG for any given tif file.
#Projection, extent and color scheme have been optimized for BioModelos
#This function can be wrapped in a for loop or apply (sapply or sfClusterApplyLB) to
#convert all tif files on a list. See example below
#
#Args:
# sp.raster = character string with tif filename (including extension)
# in.folder = folder that contains the .tif file to convert
# fill = color to be used as fill of binary files
# add.trans = add trasparent color to palette? Usually you would use TRUE when the tif file
# contains NA, 0 and 1 values, whereas you would use FALSE when the tif only
# has NA and 1 values.
#
#Usage:
#   in.folder = "~/Modelos/librorojo2"
#   fill = rgb(193,140,40,maxColorValue=255)
#   sp.raster = "Anas_bahamensis_0.tif"
#   convert2PNG(sp.raster, in.folder, fill, TRUE)
#
#Example on parallel loop
# require(snowfall)
# sfInit(parallel=T,cpus=16)#Initialize nodes
# sfExportAll() #Export vars to all the nodes
# sfClusterSetupRNG()
# sfClusterApplyLB(sp.list, convert2PNG, in.folder=in.folder, fill=fill, add.trans=TRUE)
# sfStop()
#
#Author: Jorge Velásquez
#Date: 05-09-2014

convert2PNG<-function(sp.raster, in.folder, fill, add.trans){
  require(raster)
  require(sp)
  require(rgdal)
  
  #Set transparent color
  tr <- rgb(255, 255, 255, 0, maxColorValue=255)
  if(add.trans){
    col.pal <- c(tr, fill)
  } else {
    col.pal <- c(fill)
  }
  
  #Load layers
  if(!file.exists("baseLayers.RData")){
    download.file("https://github.com/LBAB-Humboldt/parallelMaxent/raw/master/baseLayers.RData",
                  "baseLayers.RData", method="wget")
  }
  load("baseLayers.RData")
  
  #Create dirs
  dir.create(paste0(in.folder,"/PNG"), recursive=T)
  dir.create(paste0(in.folder,"/KMZ"), recursive=T)
  dir.create(paste0(in.folder,"/thumb"), recursive=T)
  
  #Plots for geovisor
  in.raster <- raster(paste0(in.folder, "/", sp.raster))
  if(is.na(projection(in.raster))){
    projection(inRaster)<-"+proj=longlat +ellps=WGS84 +datum=WGS84"
  }
  name <- strsplit(sp.raster,"[.]")[[1]][1]
  KML(in.raster, filename=paste0(in.folder,"/KMZ/",name,".kmz"),
      maxpixels=ncell(in.raster), col=col.pal, overwrite=T)
  unzip(paste0(in.folder, "/KMZ/", name, ".kmz"), exdir=paste0(in.folder,"/PNG"))
  file.remove(paste0(in.folder, "/PNG/", name,".kml"))
  
  #Generate thumbnails
  in.raster.co <- crop(in.raster, thumb.aoi) * mask.co
  png(paste0(in.folder, "/thumb/", name, "_thumb.png"),
      width=145, height=205, units="px", type="cairo")
  op <- par(mar = rep(0, 4), bg=NA)
  image(dem1000.co, axes=F, xlab="", ylab="", col=c(tr, "grey90"))
  image(in.raster.co, col=col.pal, axes=FALSE, add=T)
  plot(colombia, add=T, lwd=1, border="grey40",)
  dev.off()
  unlink(list.files(tempdir()),recursive=T)
}
