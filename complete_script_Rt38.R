
args <- commandArgs(TRUE)
input_data <- args[1]

### LOAD DATASET ####
diseasesxf <- input_data
diseasesx <- basename(input_data)
load(diseasesxf)

setwd('/scratch/scratch/ucbtdw0/diseases')

suppressMessages(library(foreach))
suppressMessages(library(doParallel))
suppressMessages(library(methods))
suppressMessages(library(raster))

#setwd('E:\\Dropbox\\legion_script\\')
source('./functions23.R')

##template
template<-raster(nrow=3432, ncol=8640,ext=extent(-180, 180, -58, 85),crs='+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')

##### AIRLINE DATA #############
airlines1<-read.csv('./air_network.csv',header=T)
flights_years<-read.csv('./wrld_simpl_data.csv',header=T)
airlines1$cell_id_start<-raster::cellFromXY(template,cbind(airlines1$sLongitude,airlines1$sLatitude))
airlines1$cell_id_dest<-raster::cellFromXY(template,cbind(airlines1$dLongitude,airlines1$dLatitude))
flights<-merge(airlines1,flights_years,by.x='Country',by.y='NAME',all.x=T)

##### ROADS SPATIAL NETWORK DATA #######################
load(file='./global_roads_rasterized_extracted.r')
load(file='./cell_ids_probs_roads.r')
load(file='./pop_gdp2.r')

#load('/home/ucbtdw0/Scratch/diseases/datasetsz23/209_dataset_ebola_all bats apes_high_both_both_3225.r')
#load('E:\\Documents\\complete_maxent_runs\\ebolas2\\209_dataset_ebola_all bats apes_med_clim_clim_2148.r')
#diseasesx <- '209_dataset_ebola_all bats apes_high_both_both_3225.r'

dataset<-merge(dataset,dd2,all.x=TRUE,all.y=FALSE,by.x="cell.id",by.y="cell.id")
rm(dd2)

## get disease data
trans<-read.csv('./diseasedata12.csv',stringsAsFactors=F)

###define some constants
e <- simpleError('test error')

###get details
details<-strsplit(diseasesx,'_')
group<-details[[1]][1]
dis<-details[[1]][3]
state<-details[[1]][5]
factors<-paste(details[[1]][6],details[[1]][7],sep='_')
spec<-details[[1]][4]

###get & permute data
transx2<-trans[trans$name==dis,]
if(nrow(transx2)==0){writeLines(text=paste(dis,' not in data file',sep=''),con=paste('./results/',dis,sample(1:10000,1),'error2.txt',sep=''))}
####add in roadn column
distsr2<-data.frame(cell.id=1:length(template),roadn=distsr)
distsr3<-distsr2[distsr2$cell.id %in% dataset$cell.id,]
distsr3<-distsr3[match(distsr3$cell.id , dataset$cell.id),]
dataset$roadn<-distsr3$roadn
rm(distsr,distsr2,distsr3)

###rid of NA's in pop
dataset$pop2015[is.na(dataset$pop2015)]<-(dataset$pop2005[is.na(dataset$pop2015)])*30
dataset$pop2015[is.na(dataset$pop2015)]<-0
dataset$pop2015[dataset$pop2015>quantile(dataset$pop2015,0.999)]<-quantile(dataset$pop2015,0.999)
dx<-dataset

### names of modis
xx<-c('Water','Evergreen.Needleleaf.forest','Evergreen.Broadleaf.forest','Deciduous.Needleleaf.forest','Deciduous.Broadleaf.forest','Mixed.forest','Closed.shrublands','Open.shrublands','Woody.savannas','Savannas','Grasslands','Permanent.wetlands','Croplands','Urban.and.built.up','Cropland.Natural.vegetation.mosaic','Snow.and.ice','Barren.or.sparsely.vegetated')
yy<-c('vWater','vEvergreen.Needleleaf.forest','vEvergreen.Broadleaf.forest','vDeciduous.Needleleaf.forest','vDeciduous.Broadleaf.forest','vMixed.forest','vClosed.shrublands','vOpen.shrublands','vWoody.savannas','vSavannas','vGrasslands','vPermanent.wetlands','vCroplands','vUrban.and.built.up','vCropland.Natural.vegetation.mosaic','vSnow.and.ice','vBarren.or.sparsely.vegetated')

######start of loop
#foreach(i=1:8) %dopar% {
#foreach(i=1:50) %do% {
for(i in 1:50){

dataset<-dx
##sort out future populations
dataset$GNI2015<-dataset$GDP2015
state2<-sample(1:3,1)
if(state2==3){dataset$GNI2050<-dataset$SSP3_GDP_2050;dataset$GNI2070<-dataset$SSP3_GDP_2070}
if(state2==2){dataset$GNI2050<-dataset$SSP2_GDP_2050;dataset$GNI2070<-dataset$SSP2_GDP_2070}
if(state2==1){dataset$GNI2050<-dataset$SSP1_GDP_2050;dataset$GNI2070<-dataset$SSP1_GDP_2070}
if(state2==3){dataset$pop2050<-dataset$SSP3_POP_2050;dataset$pop2070<-dataset$SSP3_POP_2070}
if(state2==2){dataset$pop2050<-dataset$SSP2_POP_2050;dataset$pop2070<-dataset$SSP2_POP_2070}
if(state2==1){dataset$pop2050<-dataset$SSP1_POP_2050;dataset$pop2070<-dataset$SSP1_POP_2070}
dataset$pop2050[is.na(dataset$pop2050)]<-median(dataset$pop2050,na.rm=TRUE)
dataset$pop2070[is.na(dataset$pop2070)]<-median(dataset$pop2070,na.rm=TRUE)
dataset$pop2050[dataset$pop2050>quantile(dataset$pop2050,0.999)]<-quantile(dataset$pop2050,0.999)
dataset$pop2070[dataset$pop2070>quantile(dataset$pop2070,0.999)]<-quantile(dataset$pop2070,0.999)
dataset$pop2005<-round(dataset$pop2015,0)
dataset$pop2050<-round(dataset$pop2050,0)
dataset$pop2070<-round(dataset$pop2070,0)

### flight probability change.
dataset$flight2015<-exp((1.183513*log(dataset$GNI2015))-12.96993)
dataset$flight2050<-exp((1.183513*log(dataset$GNI2050))-12.96993)
dataset$flight2070<-exp((1.183513*log(dataset$GNI2070))-12.96993)

##get rid of unlikely
if(state=="high" & state2==1){next}
if(state=="high" & state2==2){next}
if(state=="med" & state2==1){next}
if(state=="low" & state2==3){next}

##sort out names
if(state=="high"){statex="RCP8.5"}
if(state=="med"){statex="RCP6"}
if(state=="low"){statex="RCP4.5"}
stateX<-paste(statex,"_SSP",state2,sep="")

###rid of NA's in host
dataset$X2005[is.na(dataset$X2005)]<-0
dataset$X2050[is.na(dataset$X2050)]<-0
dataset$X2070[is.na(dataset$X2070)]<-0

#cl<-makeCluster(8)
#registerDoParallel(cl)
#registerDoParallel(8)
breaks1<-seq(from=25, to=50000,by=25)
dx2<-dataset

##trans
num<-paste(sample(c(0:9, letters, LETTERS),size=12, replace=TRUE),collapse='')
transx<-suppressWarnings(create_permute_trans(transx2,variation=2))
transx$w1<-(-1)*abs(sample(qnorm(c(0.05,0.25,0.50,0.75,0.95),mean=(0.025),sd=(0.025)),1,prob = c(0.05,0.25,0.50,0.25,0.05)))
transx$w2<-(-1)*sample(qnorm(c(0.05,0.25,0.50,0.75,0.95),mean=(0.82),sd=(0.025)),1,prob = c(0.05,0.25,0.50,0.25,0.05))
write.table(transx, row.names = FALSE,sep=',',file=paste('./results/A_',num,'_transx_XXX.csv',sep=''))

#zzz<-sample(2:6,5,replace=T)
###set id tag qq and k
for (z in c(1,3)){
dataset<-dx2

dataset$CM2013b<-exp(4.73106997+ dataset$GNI2015* -0.28795845 +  0.01056998*I(dataset$GNI2015^2))
dataset$CM2050b<-exp(4.73106997+ dataset$GNI2050* -0.28795845 +  0.01056998*I(dataset$GNI2050^2))
dataset$CM2050b[dataset$CM2050b>max(dataset$CM2013)]<-max(dataset$CM2013)
dataset$CM2070b<-exp(4.73106997+ dataset$GNI2070* -0.28795845 +  0.01056998*I(dataset$GNI2070^2))
dataset$CM2070b[dataset$CM2070b>max(dataset$CM2013)]<-max(dataset$CM2013)
  
### alter focal compartments to reflect year
### 2010
if(z==1){dataset$susceptible<-dataset$pop2005;dataset$host_densityx<-dataset$X2005;qq<-'2015_full';dataset$vectorx<-1;dataset$weighting<-exp(transx$w1*dataset$GNI2015);dataset$landc<-dataset$landc;  
 dataset$Rtgrad<-1.0475*(dataset$CM2013b^(transx$w2));dataset$daily_walking_distance<-dataset$WD2014;dataset$landc<-dataset$landc;dataset$flightsPersonYear<-dataset$flight2015}
### 2050
#if(z==2){dataset$susceptible<-dataset$pop2050;dataset$host_densityx<-dataset$X2050; qq<-'2050_full';dataset$vectorx<-1;dataset$weighting<-exp(transx$w1*dataset$GNI2050);dataset$landc<-dataset$fland2050;
#  dataset$Rtgrad<-1.0475*(dataset$CM2050b^(transx$w2));dataset$daily_walking_distance<-dataset$WD2050;dataset$landc<-dataset$fland2050;dataset$flightsPersonYear<-dataset$flight2050}
### 2070
if(z==3){dataset$susceptible<-dataset$pop2070;dataset$host_densityx<-dataset$X2070;qq<-'2070_full';dataset$vectorx<-1;dataset$weighting<-exp(transx$w1*dataset$GNI2070);dataset$landc<-dataset$fland2070;
  dataset$Rtgrad<-1.0475*(dataset$CM2070b^(transx$w2));dataset$daily_walking_distance<-dataset$WD2070;dataset$landc<-dataset$fland2070;dataset$flightsPersonYear<-dataset$flight2070}

### reduce to where there are landcover desingations
dataset<-dataset[!is.na(dataset$landc),]
if(nrow(dataset)==0){break}
dataset[dataset$susceptible>mean(dataset$susceptible[dataset$landc==13]),'landc']<-13
dataset[dataset$susceptible>quantile(dataset$susceptible,0.95),'landc']<-13

###give landcover a number as per raster
landmat<-data.frame(landcover=0:16,probs=1)#(as.numeric(transx2[1,xx])))
#landmatv<-data.frame(landcover=0:16,probs=(as.numeric(transx2[1,yy])))

##make unlikely but not impossible if landcover is zero
landmat[!is.na(landmat$probs==0)&landmat$probs==0,2]<-0.01
#landmatv[!is.na(landmatv$probs==0)&landmatv$probs==0,2]<-0.01

#make very urban place less likely
landmat[landmat$landcover==13,'probs']<-0.1

###set up host parameters
	dataset$probs<-as.numeric(as.character(landmat$probs[match(dataset$landc,landmat$landcover)]))
	###remember spread populations over the whole area
      dataset$host_density<-dataset$host_densityx*dataset$probs*transx$density #host density corrected by suitability

dataset$host_density[is.na(dataset$host_density)]<-0
dataset$host_density<-round(dataset$host_density,0)

if(z==1){
	##calculate global average illness interactions fully naive population
	### work out spill over need total contact i.e. time full pop times year
	datasetx<-dataset[dataset$susceptible>0 & dataset$host_density>0,]
	if(nrow(datasetx)==0){print('no rows datasetx');break}
	datasetx$susceptible[datasetx$susceptible>quantile(datasetx$susceptible,0.95)]<-ceiling(quantile(datasetx$susceptible,0.95))### random ceiling
	
	### over 2 because to estimates historic populations
	total_yearly_host_human_interactions<-sum(365*((((datasetx$susceptible/2))*0.67*0.8)*(datasetx$host_density)*transx$d*sqrt(transx$host_distance^2+mean(datasetx$daily_walking_distance)^2)),na.rm=T)
	#total_yearly_host_human_interactions<-sum(365*(log(datasetx$susceptible)*(datasetx$host_density)),na.rm=T)
	
	##per contact spill-over rate weight by present day cases
	##computational convenience
	mult=1.5
	transx$SOR<-((transx$cases_per_year*mult)/(total_yearly_host_human_interactions))

	### work out R0 need per person total i.e.e r0/ (daily times illness) but remeber recalculating Rt ## no need Rt is R0 when t=1
	transx$average_per_person_total_interactions<-mean(1*((datasetx$susceptible-1)*0.67*0.8)*0.005*5.6*sqrt(mean(datasetx$daily_walking_distance)^2+mean(datasetx$daily_walking_distance)^2),na.rm=T)## 0.67 due to groups  # 0.8 due to babies and old people not moving
}##end of z loop to construct present SOR

#if(z==1 & length(zzz)==2){next}

### clean up ####
dataset$Rtgrad[dataset$Rtgrad>(transx$R0/3)]<-transx$R0/3## give at least 3 days to react
dataset$weighting<-dataset$weighting*(nrow(dataset)/sum(dataset$weighting))
dataset$CFR<-transx$CFR*dataset$weighting
dataset$CFR[dataset$CFR<=0.01]<-0.01  
dataset$CFR[dataset$CFR>1]<-1
dataset$Rtgrad[is.na(dataset$Rtgrad)& dataset$susceptible>0]<-median(dataset$Rtgrad,na.rm=T)
dataset$CFR[is.na(dataset$CFR) & dataset$susceptible>0]<-median(dataset$CFR,na.rm=T)
dataset$host_density[is.na(dataset$host_density)]<-0
dataset$Rt<-transx$R0
dataset$SOR<-transx$SOR*dataset$weighting

#dx<-dataset
#dataset<-dx
#dataset=dataset[,c('cell.id','flightsPersonYear','x','y','airport2','susceptible','exposed','infectedA','infectedS','convalesence','spill','dead','recovered','host_density','daily_walking_distance','Rt','Rtgrad','CFR','roadn','SOR')]
#list2<-c('dx','change_compart','create_permute_trans','dataset','dis','e','flights','human_human','human_movementr','human_movementa','no_pop','rm4','run_script','SpatialLinesNetwork','spillover','stateX','template','transx','qq','prob_curve','input_data','factors')
#rm(list=ls()[!ls() %in% list2])

###last changes
dataset[is.na(dataset$airport2),'airport2']<-(-999)
dataset$ID<-NA
dataset$export<-NA
breaks1<-seq(from=25, to=50000,by=25)


       for (y in 1:365){      
       #for (y in c(1:5,25)){   
        dataset<-tryCatch(run_script(dataset=dataset,flights=flights,rm4=rm4,transx=transx,template=template),error = function(e) e)
		#print(dataset[dataset$exposed>0,c('susceptible','exposed','Rt')]) 
		#points(dataset$x[dataset$exposed>0],dataset$y[dataset$exposed>0],pch=20,cex=1,col='red')  
		print(colSums(dataset[,c('exposed','infectedS','spill','dead','recovered')] ))
	     print(table(dataset$ID));print(y)
         if(class(dataset)[1]=='simpleError'){writeLines(text=as.character(head(dataset)),con=paste('./results/',dis,'_',y,'_','_',num,'_',qq,'_','error_XXX2.txt',sep=''));break}
               if(y %in% breaks1){
			   resq<-dataset[,c('cell.id','spill','dead','recovered','Rt','SOR','ID','export'),drop=FALSE]
			   resq2<-resq[rowSums(resq[,2:(ncol(resq)-2)])>0,]
			    if(nrow(resq2)>0){
			    	save(resq2, compress='xz',file=paste('./results/',dis,'_',y,'_',stateX,'_',factors,'_',num,'_',qq,'_res','_XXX2.r',sep=''))
			    	write.csv(colSums(resq2[,2:(ncol(resq)-2),],na.rm=T),file=paste('./results/',dis,'_',y,'_',stateX,'_',factors,'_',num,'_',qq,'_summary','_XXX2.csv',sep=''))
				} else {
				resq2<-data.frame(cell.id=0,spill=0,dead=0,recovered=0,Rt=0,SOR=0)
				write.csv(resq2,file=paste('./results/',dis,'_',y,'_',stateX,'_',factors,'_',num,'_',qq,'_summary','_XXX2.csv',sep=''))
				}#resq
         		  }##### finish      
         }#end of y
	

} ### end of loop

}## end of dopar

 #library(maptools);data(wrld_simpl)
 #ttt<-template
 #values(ttt)<-NA
 #values(ttt)[resq$cell.id]<-resq$dead+resq$recovered
 #values(ttt)[dataset$cell.id]<-dataset$dead+dataset$recovered
 #values(ttt)[dataset$cell.id]<-dataset$host_density
 #values(ttt)[dataset$cell.id]<-dataset$landc
 #values(ttt)[dataset$cell.id]<-dataset$contactsx
 #values(ttt)[dataset$cell.id]<-dataset$Rtgrad
 #values(ttt)[dataset$cell.id]<-dataset$weighting
 #values(ttt)[dataset$cell.id]<-dataset$SOR
 #ttt2<-crop(ttt,extent(min(dataset$x),max(dataset$x),min(dataset$y),max(dataset$y)))
 #plot(ttt2,colNA='light blue')
 #plot(wrld_simpl,add=T)
 #points(-11.077209,7.994357)	
#points(dataset$x[dataset$airport2>0],dataset$y[dataset$airport2>0],pch=20,cex=1,col='green')  
		#