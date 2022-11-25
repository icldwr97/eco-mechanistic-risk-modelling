
library(maptools)
data(wrld_simpl)
library(rgdal)
library(raster)


dats<-read.csv("E:\\Dropbox\\ebola\\wrld_simpl_data_child_mortality_GNI2.csv")

dats2<-dats[match(wrld_simpl@data$ISO3,dats$ISO3),]

ws2<-wrld_simpl
ws2@data<-dats2

template<-raster(nrow=3432, ncol=8640,ext=extent(-180, 180, -58, 85),crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

f<-names(dats2)[13:21]

for (i in 1:length(f)){
  rasterx<-rasterize(ws2,template,field=f[i],fun="last")
  names(rasterx)<-f[i]
  plot(rasterx)
  writeRaster(rasterx,format="GTiff",paste("E:\\Dropbox\\ebola\\",f[i],".tiff",sep=""))
}
