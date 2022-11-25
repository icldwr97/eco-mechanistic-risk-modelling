library(raster)
library(rgdal)
library(ncdf4)

nc2<- stack(list.files("C:\\LUH2_v2f_beta\\LUH2_v2f_beta\\",pattern="states",recursive=TRUE,full.names=TRUE)[1])

landcov<-NULL
for (j in 1:length(list.files("C:\\LUH2_v2f_beta\\LUH2_v2f_beta\\",pattern="states",recursive=TRUE,full.names=TRUE))){

nc<- nc_open(list.files("C:\\LUH2_v2f_beta\\LUH2_v2f_beta\\",pattern="states",recursive=TRUE,full.names=TRUE)[j])

nam1<-list.files("C:\\LUH2_v2f_beta\\LUH2_v2f_beta\\",pattern="states",recursive=TRUE,full.names=FALSE)[j]
nam1<-gsub("LUH2_v2f_beta_","",nam1)
nam1<-gsub("/states.nc","",nam1)

#print(nc)
#attributes(nc)$names
#attributes(nc$var)$names

#for (i in 1:length(attributes(nc$var)$names)){
for (i in c(1:2,5:12)){
#print(ncatt_get(nc, attributes(nc$var)$names[i])$long_name)
#print(i)}

# Retrieve a matrix of the chlorophyll data using the ncvar_get function:
chla_mean <- ncvar_get(nc, attributes(nc$var)$names[i])

layer_name<-ncatt_get(nc, attributes(nc$var)$names[i])$long_name

# Print the data's dimensions
dims1<-dim(chla_mean)

test1<-brick(nrow=dims1[2],ncol=dims1[1],nl=dims1[3])

test2<-setValues(test1,chla_mean)

test3<-subset(test2,c(1,16,36,56,86))

names(test3)<-paste(nam1,layer_name,gsub("X","",names(nc2))[c(1,16,36,56,86)],sep="_")

if(is.null(landcov)){landcov<-test3} else
 {landcov<-stack(landcov,test3)}

print(i)

}##end i loop

print(j)

}##end of j loop
gc()

landcov2<-disaggregate(landcov,fact=6, method='')
landcov3<-crop(landcov2,template)
lcSSP<-values(landcov3)
lcSSP3<-subset(landcov3,101:150)

require(foreach)
require(doParallel)
require(tcltk)

    #Determine optimal block size for loading in MODIS stack data
    block_width = 15
    nrows = dim(landcov3)[1]
    nblocks <- nrows%/%block_width
    bs_rows <- seq(1,nblocks*block_width+1,block_width)
    bs_nrows <- rbind(matrix(block_width,length(bs_rows)-1,1),nrows-bs_rows[length(bs_rows)]+1)
    print('Working on the following rows')
    print(paste(bs_rows))

    #Register the parallel backend
    cl <- makeCluster(4)
    registerDoParallel(cl)
    n <- length(bs_rows)
    clusterExport(cl, c("n")) # Export max number of iteration to workers

    result <- foreach(i = 1:n, .combine = rbind) %dopar% {
     if(!exists("pb")) pb <- tcltk::tkProgressBar("Parallel task", min=1, max=n)
    tcltk::setTkProgressBar(pb, i)
    stack_values = raster::getValues(landcov3, bs_rows[i], bs_nrows[i])
      return((stack_values))
     }

  #stopImplicitCluster()
  stopCluster(cl)

save(lcSSP,file="E:\\Dropbox\\legion_script\\all_RCP_SSP_landcover.r")
#write.csv(lcSSP,file="E:\\Dropbox\\legion_script\\all_RCP_SSP_landcover.csv")

#rrr<-xyFromCell(template,coordinates(template))
#landcovE<-data.frame(cell.id=1:ncell(template))
#for (x in 1:nlayers(landcov)){
#	fff<-subset(landcov,x)
#	landcovE[,names(fff)]<-extract(fff,rrr,small=TRUE)
#	gc()
#	print(x)
#	}

nc2<- stack(list.files("C:\\Spatial_population_scenarios_GeoTIFF\\GeoTIFF\\",pattern="0.tif",recursive=TRUE,full.names=TRUE)[seq(from=1,to=600,by=4)])

test<-resample(nc2$ssp1rur2010,template,method="ngb")

sum(as.numeric(values(nc2$ssp1rur2010)),na.rm=TRUE)

sum(values(test),na.rm=TRUE)






template<-raster(nrow=3432, ncol=8640,ext=extent(-180, 180, -58, 85),crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
load("./WORLD_MASK.r")
#load(file="./livestock_future2.r")
load(file="./landcover_future3.r")
load(file="./pop_gdp2.r")

