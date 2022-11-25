library(raster)
library(dismo)
#install.packages(c("maps","GISTools"))
library(maps)
library(GISTools)
###from here

brks <- seq(0, 1, by=0.1) 
nb <- length(brks)-1 
cols <- rev(terrain.colors(nb))

bats1<-read.csv("E:\\Dropbox\\ebola\\all_fruit_bats_gbif_africa2.csv",stringsAsFactors=F)
specs<-unique(bats1$species)

predictors<-raster::stack("E:\\Dropbox\\ebola\\data\\all_fruit_bats_predictors.tif")
names1<-read.csv(file="E:\\Dropbox\\ebola\\data\\all_fruit_bats_predictors_names.csv",stringsAsFactors=F)
names(predictors)<-names1$x

#y=predictors$landcover
#lunq<-function(...) length(unique(na.omit(...)))
#newvals<-extract(y,xyFromCell(y,1:ncell(y),spatial=T),fun=lunq,buffer=10000)

##africa mask
afr<-raster("E:\\Dropbox\\Public\\Africa.tif")
afr2<-crop(afr,predictors)

#mines<-raster("E:\\Dropbox\\Public\\mines.tif")
#mines[mines==0]<-NA
#mines2<-distance(mines)
#mines3<-mask(mines2,afr2)
#writeRaster(mines3,file="E:\\Dropbox\\Public\\mines_distance.tif")

###future climate
test<-stack(unzip("E:\\Dropbox\\Public\\bio2070_hd60.zip"))
test<-stack(unzip("E:\\Dropbox\\Public\\bio2070_hd85.zip"))
test<-stack(unzip("E:\\Dropbox\\Public\\bio2070_hd26.zip"))

#landf<-raster("E:\\Dropbox\\Disease_analyses\\modis\\modis2070_high8915546.tif")
landf<-raster("E:\\Dropbox\\Disease_analyses\\modis\\land2050_high_modis_4778046.tif")
#landf<-raster("E:\\Dropbox\\Disease_analyses\\modis\\modis2070_low6485875.tif")
landf2b<-crop(landf,predictors)
landf2<-resample(landf2b,predictors,method="ngb")
landf3<-mask(landf2,afr2)
landf3[landf3==15]<-16

#test2<-raster(unzip("E:\\Dropbox\\Public\\modis2070min.zip"))
fbios<-crop(test,predictors)
names(fbios)<-c("bio10","bio11","bio12","bio13","bio14","bio2","bio5","bio6","bio7")
#fmod<-crop(test2,predictors)
#fmod2<-resample(fmod,predictors,method="ngb")
#f2070<-as.data.frame(table(values(fmod2)))[-c(4,16),]
#p2005<-as.data.frame(table(values(predictors$landcover)))
#chg<-f2070$Freq-p2005$Freq

biosx<-read.csv("E:\\Dropbox\\Public\\bioclim_trans.csv",stringsAsFactors=F)
biosx2<-biosx[biosx$abbr1 %in% names(fbios),]
biosx3<-biosx2[order(match(biosx2$abbr1, names(fbios))),]
names(fbios)<-biosx3$full4

names(predictors)[4:12]<-names(fbios)


for (v in 1:10){
#dev.off();dev.off()
###########################
  for(w in c(2,4,19,20)){

    id<-sample(1:1000000,1)
    e<-NULL;auc<-NULL;test2<-NULL;test3<-NULL;p<-NULL;pb<-NULL;e2<-NULL
    rm(sdmdata,testdata,pres_train,datax)
    datax<-bats1[bats1$species==specs[w],]

    dtx<-datax
    coordinates(dtx)<-~lon+lat

    templ<-raster(predictors)
    templ<-aggregate(templ,fact=5)
    datax$cellno<-cellFromXY(templ,dtx)

    dtx<-dtx[!duplicated(datax$cellno),]

    #plot(afr2)
    #points(dtx)

    ##sub-species rosettus and gambias
    if(w==2){dtx<-dtx[dtx$lat>(-20)&dtx$lon<(22),]}
    if(w==4){dtx<-dtx[dtx$lat>(0)&dtx$lon<(22),]}


    qtx<-bats1[!bats1$genus==strsplit(specs[w]," ")[[1]][1],]
    coordinates(qtx)<-~lon+lat
    qtx$cellno<-cellFromXY(templ,qtx)
    qtx<-qtx[!qtx$species==specs[w],]

    qtx<-bats1[!bats1$species==specs[w],]
    coordinates(qtx)<-~lon+lat
    qtx$cellno<-cellFromXY(templ,qtx)

    #print(specs[w])
    #print(nrow(dtx))
    #}

    #col1=c("pb","landcover","altitude","Mean_Temperature_of_Warmest_Quarter","Mean_Temperature_of_Coldest_Quarter","Annual_Precipitation","Precipitation_of_Wettest_Month","Precipitation_of_Driest_Month","Mean_Diurnal_Range","Max_Temperature_of_Warmest_Month","Min_Temperature_of_Coldest_Month","Temperature_Annual_Range")

    col1=c("pb","altitude","Mean_Temperature_of_Warmest_Quarter","Mean_Temperature_of_Coldest_Quarter","Annual_Precipitation","Precipitation_of_Wettest_Month","Precipitation_of_Driest_Month","Mean_Diurnal_Range","Max_Temperature_of_Warmest_Month","Min_Temperature_of_Coldest_Month","Temperature_Annual_Range")

    auc<-NULL;e2<-NULL;test3<-NULL
    for (i in 1:1){

      try1<-applyGBM(presence=dtx,absence=NULL,num1=250,predictors,col1)

      if(i==1){auc<-try1$e@auc;test3<-try1$test2;e2<-try1$e}
    if(auc<try1$e@auc){auc<-try1$e@auc;test3<-try1$test2;e2<-try1$e}

      }###end of loop

    sdmdata<-try1$sdmdata

    save(test3,file=paste("E:\\Dropbox\\ebola\\batBRM\\",specs[w],id,"_modelOLD.r",sep=""))

    write.csv(data.frame(auc=e2@auc,cor=e2@cor,pcor=e2@pcor),paste("E:\\Dropbox\\ebola\\batBRM\\",specs[w],id,"_prediction_evaluationOLD.csv",sep=""))

    tiff(file=paste("E:\\Dropbox\\ebola\\batBRM\\",specs[w],id,"_prediction_response2OLD.tiff",sep=""), width = 1080, height = 1080, units = "px", pointsize = 30,   compression = "none")
    gbm.plot(test3)
    dev.off()

    ##present
    predictors2<-subset(predictors,names(sdmdata[,2:ncol(sdmdata)]))
    #predictors2$gdp_growth[!is.na(values(predictors2$gdp_growth))]<-0                    
    #predictors2$government_effectiveness[!is.na(values(predictors2$government_effectiveness))]<-0
    #predictors2$GDP_NL2[!is.na(values(predictors2$GDP_NL2))]<-0
    predictors2<-mask(predictors2,afr2)

    p <- predict(predictors2, test3, n.trees=test3$gbm.call$best.trees, type="response")

    new1<-stack(subset(predictors,"landcover"),p)

    auc<-NULL;e2<-NULL;test3<-NULL
    for (i in 1:1){

      try1<-applyGBM(presence=dtx,absence=NULL,num1=250,predictors=new1,col1=NULL)

      #gbm.plot(try1$test)

      if(i==1){auc<-try1$e@auc;test3<-try1$test2;e2<-try1$e}
      if(auc<try1$e@auc){auc<-try1$e@auc;test3<-try1$test2;e2<-try1$e}

    }###end of loop

  sdmdata<-try1$sdmdata

  save(test3,file=paste("E:\\Dropbox\\ebola\\batBRM\\",specs[w],id,"_model.r",sep=""))

  write.csv(data.frame(auc=e2@auc,cor=e2@cor,pcor=e2@pcor),paste("E:\\Dropbox\\ebola\\batBRM\\",specs[w],id,"_prediction_evaluation.csv",sep=""))

  tiff(file=paste("E:\\Dropbox\\ebola\\batBRM\\",specs[w],id,"_prediction_response2.tiff",sep=""), width = 1080, height = 1080, units = "px", pointsize = 30,   compression = "none")
  gbm.plot(test3)
  dev.off()


  ##present
  predictors2<-subset(predictors,names(sdmdata[,2:ncol(sdmdata)]))
  #predictors2$gdp_growth[!is.na(values(predictors2$gdp_growth))]<-0                    
  #predictors2$government_effectiveness[!is.na(values(predictors2$government_effectiveness))]<-0
  #predictors2$GDP_NL2[!is.na(values(predictors2$GDP_NL2))]<-0
  predictors2<-mask(predictors2,afr2)

  p <- predict(predictors2, test3, n.trees=test3$gbm.call$best.trees, type="response")


  tiff(file=paste("E:\\Dropbox\\ebola\\batBRM\\",specs[w],id,"_prediction_image_present.tiff",sep=""), width = 1080, height = 1080, units = "px", pointsize = 30,   compression = "none")
  plot(p,cex.axis=2,cex.lab=2,breaks=brks, col=cols, lab.breaks=brks)
  points(dtx,pch=20,cex=0.8)
  maps::map.scale(x=-15.5, y=-32.75, ratio=FALSE, relwidth=0.25)
  north.arrow(xb=50, yb=-35, len=1, lab="N",col="black")
  dev.off()

  writeRaster(p,format="GTiff",file=paste("E:\\Dropbox\\ebola\\batBRM\\",specs[w],id,"_prediction_present.tif",sep=""))

  ##climate
  predictors2<-subset(predictors,names(sdmdata[,2:ncol(sdmdata)]))
  #predictors2$gdp_growth[!is.na(values(predictors2$gdp_growth))]<-0                    
  #predictors2$government_effectiveness[!is.na(values(predictors2$government_effectiveness))]<-0
  #predictors2$GDP_NL2[!is.na(values(predictors2$GDP_NL2))]<-0
  predictors2<-subset(predictors2,c(1:2))###remove old climate
  predictors2<-stack(predictors2,fbios)
  predictors2<-mask(predictors2,afr2)

  p <- predict(predictors2, test3, n.trees=test3$gbm.call$best.trees, type="response")

  tiff(file=paste("E:\\Dropbox\\ebola\\batBRM\\",specs[w],id,"_prediction_image_future_climate.tiff",sep=""), width = 1080, height = 1080, units = "px", pointsize = 30,   compression = "none")
  plot(p,cex.axis=2,cex.lab=2,breaks=brks, col=cols, lab.breaks=brks)
  maps::map.scale(x=-15.5, y=-32.75, ratio=FALSE, relwidth=0.25)
  north.arrow(xb=50, yb=-35, len=1, lab="N",col="black")
  dev.off()

  writeRaster(p,format="GTiff",file=paste("E:\\Dropbox\\ebola\\batBRM\\",specs[w],id,"_prediction_future_climate.tif",sep=""))

  ##both
  predictors2<-subset(predictors,names(sdmdata[,2:ncol(sdmdata)]))
  #predictors2$gdp_growth[!is.na(values(predictors2$gdp_growth))]<-0                    
  #predictors2$government_effectiveness[!is.na(values(predictors2$government_effectiveness))]<-0
  #predictors2$GDP_NL2[!is.na(values(predictors2$GDP_NL2))]<-0
  predictors2<-subset(predictors2,c(2))###remove old climate
  predictors2<-stack(landf3,predictors2,fbios)
  predictors2<-mask(predictors2,afr2)
  names(predictors2)[1]<-"landcover"

  p <- predict(predictors2, test3, n.trees=test3$gbm.call$best.trees, type="response")

  tiff(file=paste("E:\\Dropbox\\ebola\\batBRM\\",specs[w],id,"_prediction_image_future2_both.tiff",sep=""), width = 1080, height = 1080, units = "px", pointsize = 30,   compression = "none")
  plot(p,cex.axis=2,cex.lab=2,breaks=brks, col=cols, lab.breaks=brks)
  maps::map.scale(x=-15.5, y=-32.75, ratio=FALSE, relwidth=0.25)
  north.arrow(xb=50, yb=-35, len=1, lab="N",col="black")
  dev.off()

  writeRaster(p,format="GTiff",file=paste("E:\\Dropbox\\ebola\\batBRM\\",specs[w],id,"_prediction_future2_both.tif"))


  ##just land
  predictors2<-subset(predictors,names(sdmdata[,2:ncol(sdmdata)]))
  #predictors2$gdp_growth[!is.na(values(predictors2$gdp_growth))]<-0                    
  #predictors2$government_effectiveness[!is.na(values(predictors2$government_effectiveness))]<-0
  #predictors2$GDP_NL2[!is.na(values(predictors2$GDP_NL2))]<-0
  predictors2<-subset(predictors2,c(2:11))###remove old climate
  predictors2<-stack(landf3,predictors2)
  predictors2<-mask(predictors2,afr2)
  names(predictors2)[1]<-"landcover"

  p <- predict(predictors2, test3, n.trees=test3$gbm.call$best.trees, type="response")

  tiff(file=paste("E:\\Dropbox\\ebola\\batBRM\\",specs[w],id,"_prediction_image_future2_just_landcover.tiff",sep=""), width = 1080, height = 1080, units = "px", pointsize = 30,   compression = "none")
  plot(p,cex.axis=2,cex.lab=2,breaks=brks, col=cols, lab.breaks=brks)
  maps::map.scale(x=-15.5, y=-32.75, ratio=FALSE, relwidth=0.25)
  north.arrow(xb=50, yb=-35, len=1, lab="N",col="black")
  dev.off()

  writeRaster(p,format="GTiff",file=paste("E:\\Dropbox\\ebola\\batBRM\\",specs[w],id,"_prediction_future2_just_landcover.tif"))


  }### end of w bats loop

}#end of v












install.packages("RColorBrewer")
library(RColorBrewer)
par(mar = c(0, 4, 0, 0))
display.brewer.all()

fi<-list.files("E:\\Dropbox\\ebola\\batBRM\\",pattern="2.tif",full.names=T)[-c(3,6,9,12)]
fum<-list.files("E:\\Dropbox\\ebola\\batBRM\\",pattern="2.tif")[-c(3,6,9,12)]

for (w in 1:length(fi)){

brks <- seq(0, 1, by=0.1) 
nb <- length(brks)-1 
cols <- rev(terrain.colors(nb))


fo<-raster(fi[w])
tiff(file=paste("E:\\Dropbox\\ebola\\batBRM\\","Image_of_",fum[w],"f",sep=""), width = 1080, height = 1080, units = "px", pointsize = 30,   compression = "none")
plot(fo,cex.axis=2,cex.lab=2,breaks=brks, col=cols, lab.breaks=brks)#,col=rev(brewer.pal(9, "RdYlBu")))
maps::map.scale(x=-15.5, y=-32.75, ratio=FALSE, relwidth=0.25)
north.arrow(xb=50, yb=-35, len=1, lab="N",col="black")
dev.off()

print(w)

}
