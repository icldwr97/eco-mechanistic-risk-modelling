

mortsX<-read.csv("C:\\Users\\kjones\\Downloads\\sh.dyn.mort_Indicator_en_csv_v2\\sh.dyn.mort_Indicator_en_csv_v2.csv",stringsAsFactors=F)

for(i in 1:nrow(mortsX)){


tt<-as.data.frame(t(mortsX[i,5:ncol(mortsX)]))
tt$year=1960:2013
names(tt)[1]<-"mort"
tt<-tt[!is.na(tt$mort),]

fit1<-lm(mort~year,data=tt)
summary(fit1)

tt3<-data.frame(year=1960:2070,mort=NA)
plot(tt$year,(tt$mort),xlim=c(1960,max(tt3$year)),ylim=c(0,max(tt$mort)))

mod <- nls(mort ~   exp(b-a * year), data = tt, start = list(a = 0, b=0))
summary(mod)

# add fitted curve
lines(1960:2070, predict(mod, tt3))

tt2<-tt[1:2,]
tt2$year<-c(2050,2070)
tt2$mort<-NA
res1<-data.frame(country=mortsX$Country.Name[i],years=c(2050,2070),est=predict(mod,tt2))

if(i==1){res2<-res1;next}

res2<-rbind(res2,res1)

} #end of loop

res2$est2<-round(res2$est,3)
write.csv(res2,file="C:\\Users\\kjones\\Dropbox\\ebola\\child_mort_prediction.csv")

ests<-read.csv(file="C:\\Users\\kjones\\Dropbox\\ebola\\child_mort_prediction.csv",stringsAsFactors=F)

library(maptools)
data(wrld_simpl)
library(rgdal)
library(raster)


dats<-read.csv("C:\\Users\\kjones\\Dropbox\\ebola\\wrld_simpl_data_child_mortality.csv",stringsAsFactors=F)

dats2<-dats[match(wrld_simpl@data$ISO3,dats$ISO3),]

ws2<-wrld_simpl

ws2@data<-dats2

template<-raster("C:\\Users\\kjones\\Dropbox\\Public\\world_pop.tif")

ws3<-crop(ws2,template)

mort2013<-rasterize(ws3,template,field="X2013",fun="last")
mort2050<-rasterize(ws3,template,field="X2050",fun="last")
mort2070<-rasterize(ws3,template,field="X2070",fun="last")

writeRaster(mort2013,format="GTiff",file="C:\\Users\\kjones\\Dropbox\\Public\\cild_mortality13.tif")
writeRaster(mort2050,format="GTiff",file="C:\\Users\\kjones\\Dropbox\\Public\\cild_mortality50.tif")
writeRaster(mort2070,format="GTiff",file="C:\\Users\\kjones\\Dropbox\\Public\\cild_mortality70.tif")



