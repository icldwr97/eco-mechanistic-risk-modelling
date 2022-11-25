
#install.packages("ggplot2")
library(ggplot2)

#install.packages("psych")
library(psych)

source("E:\\Dropbox\\R_scripts\\multiplot_ggplot.r")
source("E:\\Dropbox\\R_scripts\\summarySE.r")

trmean<-function(x,y) mean(x,na.rm=T,trim=0.1)
trse<-function(x,y) (sd((x[order(x)])[round(length(x)*y,0):round(length(x)*(1-y),0)]))/sqrt(round(length(x)*(1-y),0)-round(length(x)*y,0))

res1y<-list.files("C:\\ResultsX\\results\\",pattern="summary_XXX.csv",full.names=T)
res1ys<-list.files("C:\\ResultsX\\results\\",pattern="summary_XXX.csv")

res1y2<-list.files("C:\\ResultsX\\results2\\",pattern="summary_XXX.csv",full.names=T)
res1ys2<-list.files("C:\\ResultsX\\results2\\",pattern="summary_XXX.csv")

res1y<-c(res1y,res1y2)
res1ys<-c(res1ys,res1ys2)
length(res1y)

res1ys<-gsub("_summary.csv","",res1ys,fixed=TRUE)
res1ys<-gsub("_full","_full_full",res1ys,fixed=TRUE)
res1ys<-gsub("no_pop_no_pov","nopoppov_nopoppov",res1ys,fixed=TRUE)
res1ys<-gsub("poponly","poponly_poponly",res1ys,fixed=TRUE)
res1ys<-gsub("_only_","_only_only_",res1ys,fixed=TRUE)

details<-read.table(text=res1ys,sep="_",stringsAsFactors=F)
names(details)<-c("disease","day","state2","clim_land","clim_land2","id","year2","condition","condition2")
details$scrap<-NA
details$scrap2<-NA

###########################

for (i in 1:length(res1y)){

cells<-read.table(res1y[i],header=T,sep=",",stringsAsFactors=F,blank.lines.skip =F)
cells2<-as.data.frame(t(cells[,2]))
names(cells2)<-cells$X
cells2<-cbind(cells2,details[i,])
if(ncol(cells2)==10){cells2<-cbind(data.frame(cell.id=0,spill=0,dead=0,recovered=0,Rt=0,SOR=0),cells2[,2:ncol(cells2)])}       
if(ncol(cells2)==14){cells2<-cells2[,-6]}
if(ncol(cells2)==15){cells2<-cells2[,-6];cells2<-cells2[,-6]}

if(i==1){resxx2<-cells2;next}
resxx2<-rbind(resxx2,cells2)
print(i)}

numdone<-i/length(res1y)
resxx3<-resxx2

#### other stuff
res1x<-list.files("C:\\ResultsX\\results\\",pattern="transx_XXX.csv",full.names=T)
res1xs<-list.files("C:\\ResultsX\\results\\",pattern="transx_XXX.csv",)

detailsx<-read.table(text=res1xs,sep="_",stringsAsFactors=F)
names(detailsx)<-c("scrap1","id","scrap2")
detailsx$scrap3<-NA

#### other stuff
res1x2<-list.files("C:\\ResultsX\\results2\\",pattern="transx_XXX.csv",full.names=T)
res1xs2<-list.files("C:\\ResultsX\\results2\\",pattern="transx_XXX.csv",)

detailsx2<-read.table(text=res1xs2,sep="_",stringsAsFactors=F)
names(detailsx2)<-c("scrap1","id","scrap2")
detailsx2$scrap3<-NA

res1x3<-c(res1x,res1x2)
details3<-rbind(detailsx,detailsx2)

##merge runs with esimtates & parameters
dt3<-NULL
for (i in 1:length(res1x3)){

cellsx<-read.table(res1x3[i],header=T,sep=",",stringsAsFactors=F,blank.lines.skip =F)
cellsx$id=details3$id[i]
#cellsx<-cellsx[,c("illness_length","R0","incubation","CFR","immunity","spillover_rate","cases_per_year","id")]

dt1<-resxx3[resxx3$id==cellsx$id,c("spill","dead","recovered","day","clim_land2","state2","year2","condition2") ]
if(nrow(dt1)==0){print("blah");next}#dt1<-data.frame(spill=0,dead=0,recovered=0,day=999,state2="xxx",year2=9999,condition2="xxx")}

if(i>length(res1x)){dt1$run=2}else{dt1$run=1}

dt2<-cbind(dt1,cellsx)

if(is.null(dt3)){dt3<-dt2;next}

dt3<-rbind(dt3,dt2)
print(i)}

numdone2<-i/length(res1x3)

dt3$type<-paste(dt3$id,dt3$state2,dt3$year2,dt3$condition2,dt3$clim_land2,dt3$run,sep="_")
xxx<-unique(dt3$type)
dt3$outbreak<-dt3$recovered+dt3$dead

e <- simpleError('test error')

###estimate future amount if not 300
dt5<-NULL
for(j in 1:length(xxx)){

dt4<-dt3[dt3$type==xxx[j],]

dt4max<-dt4[dt4$day==max(dt4$day),]

	if(nrow(dt4)<4){print("problem");next}
		if(max(dt4$day)<200){
			smooth1<-tryCatch(smooth.spline(dt4$day,dt4$outbreak,spar=0.25),error = function(e) e)
			if((class(smooth1)[1]=="simpleError")==TRUE){break}
			plot(dt4$day,dt4$outbreak)
			lines(smooth1)
			dt4max$outbreak<-predict(smooth1,350)$y
		}
	if(is.null(dt5)){dt5<-dt4max}else{dt5<-rbind(dt5,dt4max)}

print(j)
}
numdone3<-j/length(xxx)

save(dtx, file="all_ebola_outbreaks.r")
#dtx<-dt5
#dt5$outbreak[dt5$outbreak>max(dt3$outbreak)]<-max(dt3$outbreak)
dt5$outbreak[dt5$outbreak<0]<-0

dt5<-dt3[dt3$day==250 & dt5$outbreak<exp(17.5),]

#dt6<-dt5
dt5[dt5$year2==2015 & dt5$condition2 %in% c("nopoppov","pov"),c("clim_land2","state2","condition2")]<-"nopov"
dt5[dt5$year2==2015 & !dt5$condition2=="nopov",c("clim_land2","state2","condition2")]<-"all"

#dt5$condition<-paste(dt5$state2,dt5$year2,dt5$condition2,dt5$clim_land2,sep="_")
dt5$condition<-paste(dt5$year2,dt5$condition2,dt5$clim_land2,sep="_")
ttt<-unique(dt5$condition)

table(dt5$condition)


resxxd2<-dt5[dt5$condition %in% c("2015_all_all","2070_full_both","2070_pov_both","2015_nopov_nopov"),]
#resxxd2<-dt5[dt5$condition %in% c("2015_all_all","2070_only_land","2070_full_both","2070_pov_both","2070_only_clim","2070_poponly_both","2015_nopov_nopov"),]
dim(resxxd2)
resxxdX<-resxxd[resxxd2$outbreak>10000000,]
resxxd<-resxxd2#[resxxd2$outbreak<3000,]
tall<-as.data.frame(table(resxxd2$condition))
tlarge<-as.data.frame(table(resxxdX$condition))
tlarge$prop<-tlarge$Freq/tall$Freq
tlarge

resxxd$condition2<-as.factor(resxxd$condition)
#levels(resxxd$condition2)<-c("Present","Present w/o GDP","Future","Climate Only","Land Use Only","Demography Only","Future w/o GDP")
#resxxd$condition3 <- factor(resxxd$condition2, levels = c("Present","Future","Present w/o GDP","Future w/o GDP","Climate Only","Land Use Only","Demography Only"))
levels(resxxd$condition2)<-c("Present","Present w/o GDP","Future","Future w/o GDP")
resxxd$condition3 <- factor(resxxd$condition2, levels = c("Present","Future","Present w/o GDP","Future w/o GDP"))

#resxxd<-resxxd[resxxd$condition3 %in% c("Present","Future"),]

resxxd$state3<-resxxd$state2
resxxd$state3[resxxd$state3=="med"]<-"RCP4.5"
resxxd$state3[resxxd$state3=="low"]<-"RCP6"
resxxd$state3[resxxd$state3=="high"]<-"RCP8.5"
resxxd$state3[resxxd$state3=="all"]<-"Present"
resxxd$state3[resxxd$state3=="nopov"]<-"Present"

resxxd$state3 <- factor(resxxd$state3, levels = c("Present","RCP4.5","RCP6","RCP8.5"))


# Error bars represent standard error of the mean
ggplot(resxxd, aes(x=state3, y=log(outbreak+1)),colour=condition3) + 
  geom_violin(fill = "grey",colour="grey",draw_quantiles = c(0.25, 0.5, 0.75))+ 
stat_summary(fun.data = "mean_cl_boot", geom = "errorbar",
                         colour = "black", width = 0.15)   + 
stat_summary(fun.data = "mean_cl_boot", geom = "point",
                         colour = "black", width = 0.3)  +
xlab("RCP Scenario") +
ylab("Log Mean EVD Cases per Year") +
 theme_bw() +
    theme(axis.title.y = element_text(size = rel(1.5), angle = 90)) +
   theme(axis.title.x = element_text(size = rel(1.8), angle = 0)) +
 theme(axis.text.x = element_text(size = rel(1.8)))+
 theme(axis.text.y = element_text(size = rel(1.8)))+
 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    panel.background = element_blank(),   panel.border = element_blank(),axis.line = element_line(colour = "black"))

scale_y_continuous(breaks=c(10,11,12,13,14,15,16,17,18,19,20,21), labels=c(round(exp(10),0)-1,round(exp(11),0)-1,round(exp(12),0)-1, round(exp(13),0)-1, round(exp(14),0)-1, round(exp(15),0)-1, round(exp(16),0), round(exp(17),0)-1, round(exp(18),0)-1, round(exp(19),0)-1, round(exp(20),0)-1, round(exp(21),0)-1))+

###MAIN FIGURE 3

resxxd$outpp<-resxxd$outbreak#/resxxd$sum.dataset.pop1.
trmean<-function(x,y,z) mean(c(x,rep(0,z)),na.rm=T,trim=y)
#trse<-function(x,y,z) (sd((c(x,rep(0,z))[order(c(x,rep(0,z)))])[round(length(c(x,rep(0,z)))*y,0):round(length(c(x,rep(0,z)))*(1-y),0)]))/sqrt(round(length(c(x,rep(0,z)))*(1-y),0)-round(length(c(x,rep(0,z)))*y,0))

resxxd2<-resxxd#[resxxd2$outpp<10000,]

#summary stats
trim1<-0.2#.01#.01
zz<-tapply(resxxd2$outpp, resxxd2$condition,median)#,trmean,y=trim1,z=nrow(resxxd3))
zz2<-tapply(resxxd2$outpp, resxxd2$condition,sd)
zz3<-tapply(resxxd2$outpp, resxxd2$condition,length)
dfc<-data.frame(condition=names(zz),N=zz3,subject=zz,se=zz2)


#dfc$condition<-factor(c("2015 with GDP","2015 Present","2070 Climate","2070 with GDP","2070 Land Use","2070 no Pop","2070 Future"),levels=c("2015 Present","2015 with GDP","2070 Future","2070 with GDP","2070 Climate","2070 Land Use","2070 no Pop"))
#dfc2<-dfc#[dfc$condition %in% c("2015 with GDP","2070 with GDP"),]
#dfc2<-dfc[!dfc$condition %in% c("2015 with GDP","2070 with GDP"),]


# Error bars represent standard error of the mean #p<-
ggplot(dfc, aes(x=condition, y=subject,color=condition)) + 
       geom_bar(position=position_dodge(), stat="identity",size=1.2,show_guide=FALSE,fill="white",color="black") +
    geom_errorbar(aes(ymin=(subject)-se, ymax=(subject)+se),
                  size=.2,    # Thinner lines
                  width=.2,
                  position=position_dodge(.9)) +
 xlab(NULL)+
       ylab("Mean EVD cases per year")+
        theme_bw() +
    theme(axis.title.y = element_text(size = rel(1.5), angle = 90)) +
   theme(axis.title.x = element_text(size = rel(1.8), angle = 90)) +
 theme(axis.text.x = element_text(size = rel(1.8)))+
 theme(axis.text.y = element_text(size = rel(1.8)))+
 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    panel.background = element_blank(),   panel.border = element_blank(),axis.line = element_line(colour = "black"))


### STATES figure 2c

resxxd<-resxx[resxx$day==365 & resxx$condition=="both" ,]

multp<-365/max(resxxd$day)

resxxd$outpp<-resxxd$outbreak#/resxxd$sum.dataset.pop1.

# summary stats
trim1<-0.01
zz<-tapply(resxxd$outpp, resxxd$state,trmean,y=trim1)
zz2<-tapply(resxxd$outpp, resxxd$state,trse,y=trim1)
zz3<-tapply(resxxd$outpp, resxxd$state,length)
dfc2<-data.frame(state=names(zz),N=zz3,subject=zz*multp,se=zz2)

dfc2$state<-factor(c("High Change","Low Change","Business as Usual"),levels=c("High Change","Business as Usual","Low Change"))


# Error bars represent standard error of the mean #p<-
ggplot(dfc2, aes(x=state, y=subject)) + 
       geom_bar(position=position_dodge(), width=0.5,stat="identity",show_guide=FALSE,fill="white",color="black") +
    geom_errorbar(aes(ymin=(subject)-se, ymax=(subject)+se),
                  size=.3,    # Thinner lines
                  width=.2,
                  position=position_dodge(.9)) +
 xlab(NULL)+
       ylab("Mean EVD cases per year")+
        theme_bw() +
    theme(axis.title.y = element_text(size = rel(1.8), angle = 90)) +
   theme(axis.title.x = element_text(size = rel(1.8), angle = 0)) +
 theme(axis.text.x = element_text(size = rel(1.8)))+
 theme(axis.text.y = element_text(size = rel(1.8)))+
 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    panel.background = element_blank(),   panel.border = element_blank(),axis.line = element_line(colour = "black"))+
geom_hline(yintercept=dfc$subject[6],linetype=2)+
geom_hline(yintercept=dfc$subject[6],linetype=1,size=15,alpha=0.5,colour="grey")

### STATES figure 2b

resxxd<-resxx[resxx$day==365 & resxx$condition=="both" ,]

resxxc<-mean((resxx[resxx$day==365 & resxx$condition=="present" ,"spill_overs"])/max(resxxd$day),trim=0.05)
resxxc2<-trse((resxx[resxx$day==365 & resxx$condition=="present" ,"spill_overs"])/max(resxxd$day),y=0.05)

multp<-365/max(resxxd$day)

# summary stats
trim1<-0.05
zz<-tapply(resxxd$spill_overs/resxxd$day, resxxd$state,trmean,y=trim1)
zz2<-tapply(resxxd$spill_overs/resxxd$day, resxxd$state,trse,y=trim1)
zz3<-tapply(resxxd$spill_overs/resxxd$day, resxxd$state,length)
dfc2<-data.frame(state=names(zz),N=zz3,subject=zz*multp,se=zz2)

dfc2$state<-factor(c("High Change","Low Change","Business as Usual"),levels=c("High Change","Business as Usual","Low Change"))



# Error bars represent standard error of the mean #p<-
ggplot(dfc2, aes(x=state, y=subject)) + 
       geom_bar(position=position_dodge(), width=0.5,stat="identity",show_guide=FALSE,fill="white",color="black") +
    geom_errorbar(aes(ymin=(subject)-se, ymax=(subject)+se),
                  size=.3,    # Thinner lines
                  width=.2,
                  position=position_dodge(.9)) +
       xlab(NULL)+
       ylab(" Mean index cases per year")+
ylim(0,3)+
        theme_bw() +
    theme(axis.title.y = element_text(size = rel(1.8), angle = 90)) +
   theme(axis.title.x = element_text(size = rel(1.8), angle = 0)) +
 theme(axis.text.x = element_text(size = rel(1.8)))+
 theme(axis.text.y = element_text(size = rel(1.8)))+
 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    panel.background = element_blank(),   panel.border = element_blank(),axis.line = element_line(colour = "black"))+
geom_hline(yintercept=resxxc,linetype=2)+
geom_hline(yintercept=resxxc,linetype=1,size=resxxc2*100,alpha=0.5,colour="grey")


###final figure 2a

load(file="ebola_cell_numbers2.r")

resxxd<-res2[res2$day==365 & res2$condition=="both" & res2$cells>0 ,]
#resxxd<-res2[res2$day==225 & res2$condition=="present" & res2$cells>0,]

resxxd$fact1<-1
resxxd$fact1[resxxd$names=="Dead" | resxxd$names=="Recovered"]<-"yes"
resxxd$fact1[resxxd$names=="Infectious" | resxxd$names=="Index"]<-"no"

resxxd<-resxxd[resxxd$fact1=="yes",]
#fut1<-resxxd$cells

dfc <- summarySEwithin(data=resxxd,withinvars="year", measurevar="cells",na.rm=TRUE)  
names(dfc)[3]<-"subject"
dfc$subject<-dfc$subject*(5.6^2)
dfc$year<-factor(c("High Change","Low Change","Business as Usual"),levels=c("High Change","Business as Usual","Low Change"))

resxxc<-995;mean((resxxd[ ,"cells"])*(5.6^2),trim=0.01)
resxxc2<-72;trse((resxxd[ ,"cells"])*(5.6^2),y=0.01)

# Error bars represent standard error of the mean #p<-
ggplot(dfc, aes(x=year, y=subject)) + 
       geom_bar(position=position_dodge(), width=0.5, stat="identity",show_guide=FALSE,fill="white",color="black") +
    geom_errorbar(aes(ymin=(subject)-sd, ymax=(subject)+sd),
                  size=.3,    # Thinner lines
                  width=.2,
                  position=position_dodge(.9)) +
       xlab(NULL)+
       ylab("Mean outbreak area per year (km2)")+
       ylim(0,2150)+
        theme_bw() +
    theme(axis.title.y = element_text(size = rel(1.5), angle = 90)) +
   theme(axis.title.x = element_text(size = rel(1.8), angle = 0)) +
 theme(axis.text.x = element_text(size = rel(1.8)))+
 theme(axis.text.y = element_text(size = rel(1.8)))+
 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    panel.background = element_blank(),   panel.border = element_blank(),axis.line = element_line(colour = "black"))+

geom_hline(yintercept=resxxc,linetype=2)+
geom_hline(yintercept=resxxc,linetype=1,size=resxxc2/10,alpha=0.5,colour="grey")



