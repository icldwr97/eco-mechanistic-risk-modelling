
library(maps)
library(spatstat)
library(maptools)
library(raster)
library(dismo)
library(fields)
library(GISTools)


data(world.cities)
wc<-world.cities[world.cities$pop>1000000,]
coordinates(wc)<-~long+lat

template<-raster(nrow=3432, ncol=8640,ext=extent(-180, 180, -58, 85),crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

##africa mask
all_afr<-raster("E:\\Dropbox\\ebola\\batBRM\\final3\\high_climate.tif")
all_afr[!is.na(all_afr)]<-1
afr<-raster("E:\\Dropbox\\Public\\conts\\Africa.tif")
afr2b<-crop(afr,all_afr)

#alt1<-subset(all_afr,1)
#alt1<-mask(alt1,afr2)

csx<-read.csv("E:\\Dropbox\\ebola\\old_workings\\ebola_cases_africa4.csv",stringsAsFactors=F)
coordinates(csx)<-~LONG+LAT


load("E:\\Dropbox\\ebola\\old_workings\\all_points_all_ebola.r")

save3$rcpssp<-paste(save3$RCP,save3$SSP,sep="_")

save3$group2<-save3$group
save3$group2[save3$group2=="verylarge"]<-"large"
save4<-save3[!save3$group=="cases",]
save4$group2<-factor(save4$group2)
save4$group<-factor(save4$group)


sv4<-xtabs(dummy~group2+rcpssp,data=save4)
ch4<-chisq.test(sv4)
ch4$residuals

sv3<-xtabs(dummy~group+rcpssp,data=save4)
ch3<-chisq.test(sv3)
ch3$residuals





#cols<-c("#5b0004","#8e0401","#a40000","#ee8609","#f7d207","#7dc000","#3ea280","#1b758f","#0c385b","#020913")#GLOBAIA1
cols<-c("#5b0004","#8e0401","#a40000","#ee8609","#f7d207","#7dc000","#3ea280","#0c385b","#020913")#GLOBAIA2
cols<-c("#5b0004","#8e0401","#a40000","#ee8609","#f7d207","#7dc000","#3ea280","white")#GLOBAIA3
cols<-c("#5b0004","#8e0401","#a40000","#ee8609","#f7d207","#7dc000","#3ea280","#0c385b")#GLOBAIA4
cols<-c("#e03730","#f2a490","#fce9af","#92b7ca")#piggot
cols<-c("#ffffb2","#fed976","#feb24c","#fd8d3c","#f03b20","#bd0026")#6-class YlOrRd
cols<-c("#ffffcc","#c7e9b4","#7fcdbb","#41b6c4","#2c7fb8","#253494")#6-class YlGnBu
cols<-c("#edf8fb","#bfd3e6","#9ebcda","#8c96c6","#8856a7","#810f7c")#6-class BuPu
cols<-c("#e6f598","#fee08b","#fdae61","#f46d43","#d53e4f","#9e0142")#spatial.ly
cols<-c("#5b0004","#8e0401","#a40000","#ee8609","#f7d207","#7dc000","#3ea280","#0c385b")#GLOBAIA4
#colXX<-colorRampPalette(cols,space="rgb",interpolate ="spline",bias = 0.8)
#colXX<-colorRampPalette(rev(cols),space="rgb",interpolate ="spline")

## divide up cases
	cs<-csx[csx$Virus=="Ebola Zaire",]
	epi1<-csx[csx$UNIQ_ID=="x6" & csx$CASE_TYPE=="index",]
	css<-cs[cs$CASE_TYPE=="secondary",]
	cs2<-cs[cs$CASE_TYPE=="index",]
	csi<-cs[cs$CASE_TYPE=="import",]
	csud<-csx[csx$Virus=="Ebola Sudan" & csx$CASE_TYPE=="index",]
	csuds<-csx[csx$Virus=="Ebola Sudan" & csx$CASE_TYPE=="secondary",]
	cstai<-csx[csx$Virus=="Ebola Ivory Coast"& csx$CASE_TYPE=="index",]
	csbb<-csx[csx$Virus=="Ebola Bundibugyo"& csx$CASE_TYPE=="index",]

	allind<-csx[csx$CASE_TYPE=="index",]
	
	

#plotRGB(ocean2)
#plot(afr2,col="#edf8fb",add=T)
#plot(afr2,col="#0c385b",add=T)
#plot(afr2,col="#e6f598",add=T)
#plot(afr2,col="#020913",add=T)

save3$uni<-paste(save3$group,save3$RCP,save3$SSP,sep="_")

		library(maptools)
	data(wrld_simpl)

	library(viridis)

ttt<-unique(save3$uni)
res2<-NULL
nums<-c(seq(1,length(ttt),by=4),seq(2,length(ttt),by=4),seq(3,length(ttt),by=4))

#for (y in seq(1,length(ttt),by=4)){
#for (y in seq(3,length(ttt),by=4)){
#for (y in seq(2,length(ttt),by=4)){

for (y in nums){
	print(ttt[y])
	pointsX<-save3[save3$uni==ttt[y],]
	pointsX$xy<-paste(pointsX$x,pointsX$y)
	pointsX<-pointsX[!duplicated(pointsX$xy),]
	coordinates(pointsX)<-~x+y

	xxx<-raster(all_afr)
	res(xxx)<-0.5
	afr2<-resample(afr2b,xxx)
	xxx<-crop(xxx,afr2)
	pres<-rasterize(pointsX,xxx,field="dummy",fun=sum)
	pres[is.na(pres)]<-0
	pres<-raster::mask(pres,afr2)
	#table(values(pres))
	#pres[pres>20]<-20
	#plot(pres,col=c("light blue",rev((heat.colors(50)))),colNA="grey")

	pres2<-disaggregate(pres,fact=12, method= 'bilinear')
	pres2<-mask(pres2,all_afr)
	pres2[is.na(pres2)]<-0
	p<-focal(pres2, w=matrix(1,25,25), fun=mean)
	p<-focal(p, w=matrix(1,25,25), fun=mean)
	#p<-focal(p, w=matrix(1,25,25), fun=mean)
	print(p)
	p[p<0]<-0
	p<-mask(p,all_afr)
	#plot(p,col=c("light grey",rev((heat.colors(50)))),colNA="white")
	#if(y %in% seq(1,length(ttt),by=4)){p[1,1]<-120}
	p<-p/365
	
	##set up difference
	if(y==1){cases_pres<-p}
	if(y==2){index_pres<-p}
	if(y==3){larg_pres<-p}
	
	
		#cols<-c("#5b0004","#8e0401","#a40000","#ee8609","#f7d207","#7dc000","#3ea280","#0c385b")#GLOBAIA4
	#colXX<-colorRampPalette(cols,space="rgb",interpolate ="spline")

	#if(y %in%  seq(1,length(ttt),by=4)){ccc<-css} else {ccc<-cs2}
	#p2<-p/max(values(p),na.rm=TRUE)
	#p2[p2<=quantile(values(p2)[values(p2)>0 & !is.na(values(p2))],0.05)]<-0
	#p2<-crop(p2,extent(c(-18,53,-20,20)))
	#p3<-p
	#p3[p3>quantile(values(p2)[values(p2)>0 & !is.na(values(p2))],0.05)]<-1
	#p3[p3<=quantile(values(p2)[values(p2)>0 & !is.na(values(p2))],0.05)]<-0

	#p2b<-p2
	#p2b<-crop(p2,c(-18,53,-10,10))###
	#p2b<-crop(p2,c(-18,53,-15,15))###
	#p2b<-afr2b
	
	#cases_pres<-p
	
	pdf(file=paste("E:\\Dropbox\\ebola\\",ttt[y],"_",sample(1:1000,1),"NNN.pdf",sep=""),width=10,height=10)

	#main plot
	#par(bg=NA) 
	#plot(p,col=rev(colXX(100)),legend=F,bg=NA,cex.axis=2.5)
	if(y %in% seq(5,24,by=4)){dif<-p-cases_pres}
	if(y %in% seq(6,24,by=4)){dif<-p-index_pres}
	if(y %in% seq(7,24,by=4)){dif<-p-larg_pres}
	
	dif[dif>abs(min(values(dif),na.rm=TRUE))]<-abs(min(values(dif),na.rm=TRUE))
	plot(dif,colNA="white",col=c(rev(viridis(30)[1:20]),viridis(10,option="A")[2],(viridis(30,option="A"))),cex.axis=2.5)
	plot(p,colNA="white",col=c("light grey",(viridis(30,option="D"))),cex.axis=2.5)
	
		#plotRGB(ocean2,add=T,bg=NA)
	#plot(p,col=rev(colXX(100)),add=T,maxpixels=10000000,bg=NA,cex.axis=2.5,legend=T)
	#legend(x = 10, y = 10, legend = seq(p@data@min,p@data@max,by=0.05), fFill = rev(colXX(100)))
	#scalebar(d=1000,type="bar",xy=c(41,-38),cex=2.5)
	#north.arrow(len=1,yb=-31,xb=45,col="black",lab="N",cex=2.5)
	
	#borders
	ws2<-wrld_simpl[wrld_simpl$REGION %in% c(2,142)  & wrld_simpl$AREA>1000 |wrld_simpl$NAME=="Western Sahara" ,]
	plot(ws2,border="#FFFFFF60",lwd=0.25,add=T)#,lty=3)
	
	
	if(y==1){
		#add disease points
		leg<-data.frame(x=(-10),y=c(-23,-30,-34))
		leg2<-data.frame(x=(2),y=c(-21.5,-29,-33.5))
		text(leg2,c("10000 cases","500 cases","10 cases"),cex=1.5)		
		points(leg,col="black",pch='+',lwd=3,cex=c(log(10000),log(500),log(10)))
		points(css,col="#FFFFFF99",pch='+',lwd=3,cex=log(css$ADJ_cases))
		points(csuds,col="#FFFFFF99",pch='+',lwd=3,cex=log(csuds$ADJ_cases))
		#points(cs2,col="light grey",cex=3,pch=19)
		#points(csud,col="light grey",cex=3,pch=15)
		#points(cstai,col="light grey",cex=3,pch=24,bg="white")
		#points(csbb,col="light grey",cex=3,pch=25,bg="white")
		
		}
	if(y==3){
		#add disease points
		#points(css,col="#000000",pch=1,lwd=0.5,cex=log(css$ADJ_cases))
		#points(csuds,col="black",pch=1,lwd=2,cex=log(csuds$ADJ_cases)*3)
		points(epi1,col="#FFFFFF",cex=4,lwd=2,pch=4)
	  points(css[css$STR_YEAR==2018,],col="#FFFFFF",pch=8,lwd=2,cex=log(css$ADJ_cases)+1)
	  
		#points(csud,col="white",cex=2,pch=15)
		#points(cstai,col="white",cex=2,pch=24,bg="black")
		#points(csbb,col="white",cex=2,pch=25,bg="black")
		#leg<-data.frame(x=(-10),y=c(-20,-26,-30))
		#points(leg,col="#000000",pch=1,lwd=2,cex=c(log(10000),log(1000),log(10)))
		}
	if(y==2){
		#add disease points
		#points(css,col="#000000",pch=1,lwd=0.5,cex=log(css$ADJ_cases))
		#points(csuds,col="black",pch=1,lwd=2,cex=log(csuds$ADJ_cases)*3)
		points(cs2,col="white",cex=2,lwd=3,pch=0)
		points(csud,col="white",cex=2,lwd=3,pch=1)
		points(cstai,col="white",cex=2,lwd=3,pch=2)
		points(csbb,col="white",cex=2,lwd=3,pch=5)
		#leg<-data.frame(x=(-10),y=c(-20,-26,-30))
		#points(leg,col="#000000",pch=1,lwd=2,cex=c(log(10000),log(1000),log(10)))
		}
	
	
	dev.off()

	}### end of y loop
#}############end of h loop


library(ggplot2)
install.packages("corrplot")
library(corrplot)
#save(results3,file="E:\\Dropbox\\ebola\\all_ebola_results_by_IDb.r")
load(file="E:\\Dropbox\\ebola\\all_ebola_results_by_ID.r")
load(file="E:\\Dropbox\\ebola\\all_ebola_results_by_IDb.r")
results3b<-results3
load(file="E:\\Dropbox\\ebola\\all_ebola_results_by_IDc.r")


results3$outbreak=results3$dead+results3$recovered

results5<-results3


#results5$outbreak[results5$outbreak<500]<-500
#results5<-results5[results5$outbreak<0,]
results5$large<-"less than 3"
results5$large[results5$outbreak>3 & results5$outbreak<1500]<-"3 to 1500"
results5$large[results5$outbreak>1500]<-"1500+"
results5$condition3<-gsub("_"," ",results5$condition3)
zz<-xtabs(~condition3+large,data=results5)
zz2<-chisq.test(zz)	
#zz2$expected
#zz2$observed
#zz2$residuals
#zz2$stdres

png(height=10000,width=10000,file="E:\\Dropbox\\ebola\\corrplot3.png",antialias = "cleartype" )
corrplot(zz2$stdres,is.cor=FALSE,tl.col="black", tl.srt=45,tl.cex=20,cl.cex=20,cl.pos="b",cl.length=5,cl.lim=c(-20,20),col=rev(col4(50)))
dev.off()

col1 <- colorRampPalette(c("#7F0000","red","#FF7F00","yellow","white", 
        "cyan", "#007FFF", "blue","#00007F"))
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
            "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))  
col3 <- colorRampPalette(c("red", "white", "blue")) 
col4 <- colorRampPalette(c("#7F0000","red","#FF7F00","yellow","#7FFF7F", 
            "cyan", "#007FFF", "blue","#00007F"))   
wb <- c("white","black")



results3$dummy<-1
rrr1<-aggregate(results3$dummy,by=list(results3$type),FUN=sum)
result3

# Error bars represent standard error of the mean
ggplot(results5, aes(x=condition3, y=log(outbreak))) + 
  geom_violin(fill = "grey",colour="grey",draw_quantiles = c(0.25, 0.5, 0.75))+ 
stat_summary(fun.data = "mean_cl_boot", geom = "errorbar",
                         colour = "black", width = 0.15)   + 
stat_summary(fun.data = "mean_cl_boot", geom = "point",
                         colour = "black", width = 0.3)  +
xlab("RCP Scenario") +
ylab("Log Mean EVD Cases per Year") +
ylim(1,8)+
 theme_bw() +
    theme(axis.title.y = element_text(size = rel(1.5), angle = 90)) +
   theme(axis.title.x = element_text(size = rel(1.8), angle = 0)) +
 theme(axis.text.x = element_text(size = rel(1.8)))+
 theme(axis.text.y = element_text(size = rel(1.8)))+
 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    panel.background = element_blank(),   panel.border = element_blank(),axis.line = element_line(colour = "black"))



p<-xx
r.range <- c(minValue(p), maxValue(p))


png(file="E:\\Dropbox\\ebola\\legend3.png",width=500,height=600,bg="transparent")
plot(p, legend.only=TRUE, col=rev(terrain.colors(50)),bg="transparent",
     legend.width=1, legend.shrink=20,
     axis.args=list(at=seq(r.range[1], r.range[2], 0.2),
                    labels=round(seq(r.range[1], r.range[2], 0.2),1), 
                    cex.axis=2),
     legend.args=list(text='Probability of presence', side=4, font=3, line=5, cex=2),
smallplot=c(.15, .17, .5, .85))

dev.off()

plot(cs2)

png(file="E:\\Dropbox\\ebola\\legend2.png",width=500,height=600,bg="transparent")
plot(p2, legend.only=TRUE, col=rev(colXX(100)),bg="transparent",
     legend.width=1, legend.shrink=20,
     axis.args=list(at=seq(r.range[1], r.range[2], 0.001),
                    labels=seq(r.range[1], r.range[2], 0.001), 
                    cex.axis=2),
     legend.args=list(text='Cell probability per run', side=4, font=3, line=8, cex=2))#,
#smallplot=c(.15, .17, .5, .85))
points(10,10,cex=100)
dev.off()

image(p,  col=rev(terrain.colors(255))) 
plot(p, horizontal=FALSE, smallplot=c(.15, .5, .84, .86), 
legend.only=TRUE) 

image(p,  col=rev(terrain.colors(255))) 
plot(p, smallplot=c(.15, .17, .5, .85), legend.only=TRUE,legend.width=10) 

####present and future host

xx<-raster("E:\\Dropbox\\ebola\\batBRM\\final3\\med_both.tif")
xx<-raster("E:\\Dropbox\\ebola\\batBRM\\final3\\med_present.tif")
xx<-xx/max(values(xx),na.rm=T)
tiff(file="E:\\Dropbox\\ebola\\host_present.tiff",width=1018,height=1080,bg=NA)
plot(xx,cex.axis=2.5)
scalebar(d=1000,type="bar",xy=c(41,-38),cex=2.5)
north.arrow(len=1,yb=-31,xb=45,col="black",lab="N",cex=2.5)
dev.off()


