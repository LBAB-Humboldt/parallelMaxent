#convert2PNG.R
#
#Converts all the tif files in a folder to KMZ and PNG
#as long as their filename matches any of the strings in the argument include. The color
#scheme is thought only for the binary color scheme used in BioModelos.
#
#Args:
# inFolder = folder that contains .tif files to convert
# files2process = character vector with the filenames of files to process

#
#Usage:
#  inFolder<-"W:/Modelos/20131122"
#  convert2PNG(inFolder,files2process)
#
#Author: Jorge Velásquez
#Date: 26-02-2014

convert2PNG<-function(inFolder,files2process){
  require(raster)
  require(sp)
  require(rgdal)

  #Plotting parameters
  tr=rgb(255,255,255,0,maxColorValue=255)
  fill=rgb(193,140,40,maxColorValue=255)
  
  #Create base rasters
  colombia<-readOGR(dsn="C:/Users/aves/Google Drive H/Google Drive/SDM_BaseFiles/tmp","COL_adm0")
  thumb_aoi<-readOGR(dsn="C:/Users/aves/Google Drive H/Google Drive/SDM_BaseFiles/tmp","thumbnail_area")
  dem<-raster("C:/Users/aves/Google Drive H/Google Drive/SDM_BaseFiles/tmp/alt.asc")
  dem1000_co <- crop((dem > 1000),thumb_aoi)
  
    
  mask.co <- rasterize(colombia,dem1000_co,field=1)
  dem1000_co <- dem1000_co * mask.co
  
  inRaster<-raster(paste0(inFolder,"/",files2process[1])) #Assumes all layers in inFolder have the same extent and res
  mask.co <- rasterize(colombia,inRaster,field=1)
  #Create dirs
  dir.create(paste0(inFolder,"/PNG"),recursive=T)
  dir.create(paste0(inFolder,"/KMZ"),recursive=T)
  dir.create(paste0(inFolder,"/thumb"),recursive=T)
  
  #Start plotting loop
  for(j in files2process){
    #Plots for geovisor
    inRaster<-raster(paste0(inFolder,"/",j))*1
    if(is.na(projection(inRaster))){
      projection(inRaster)<-"+proj=longlat +ellps=WGS84 +datum=WGS84"
    }
    name=strsplit(j,"[.]")[[1]][1]
    print(name)
    KML(inRaster,filename=paste0(inFolder,"/KMZ/",name,".kmz"),
        maxpixels=ncell(inRaster),col=c(tr,fill),overwrite=T)
    unzip(paste0(inFolder,"/KMZ/",name,".kmz"),exdir=paste0(inFolder,"/PNG"))
    file.remove(paste0(inFolder,"/PNG/",name,".kml"))
    
    #Generate thumbnails
    in.raster_co <- crop(inRaster,thumb_aoi)*mask.co
    png(paste0(inFolder,"/thumb/",name,"_thumb.png"),width=145,height=205,units="px",type="cairo")
    op <- par(mar = rep(0, 4),bg=NA)
    image(dem1000_co,axes=F,xlab="",ylab="",col=c(tr,"grey90"))
    image(in.raster_co,col=c(tr,fill),axes=FALSE,add=T)
    plot(colombia,add=T,lwd=1,border="grey40",)
    dev.off()
    unlink(list.files(tempdir()),recursive=T)
  }
}