
####################
#require(raster)
#require(igraph)
#require(foreach)
#require(doParallel)
#registerDoParallel(6)
#require(compiler)

suppressMessages(library(sp))
suppressMessages(library(raster))
suppressMessages(library(foreach))


################################################
############### CREATE PERMUTE TRANS ############
create_permute_trans<-function(transx2,permute=FALSE,quantiles=TRUE,variation=5){
    if(nrow(transx2)==2){transx2<-transx2[sample(c(1,2),1,prob=transx2$spillover_rate..low.to.high.),]}
    res1=data.frame(
        illness_length=mean(transx2$acute_low,transx2$acute_high,na.rm=T),
        R0=transx2$R0,
        chronic.morbidity=transx2$chronic.morbidity,
        convalescence_length=transx2$Convalescence_length,
        chronic_length=mean(transx2$Chronicduration_low,transx2$Chronicduration_high,na.rm=T),
        incubation=mean(transx2$Incubation.low,transx2$Incubation.high,na.rm=T),
        CFR=mean(transx2$CFR.low,transx2$CFR.high,na.rm=T),
        immunity=transx2$proportion.of.recovered.that.are.immmune,
        spillover_rate=transx2$spillover_rate_outbreaks,
        chronic.carrier=transx2$X..chronic.carrier,
	  cases_per_year=transx2$cases_per_year,
	host_distance=transx2$daily_dist,
	d=transx2$d,
	density=transx2$density,
	vdensity=transx2$vdensity,
	vhost_distance=transx2$vdaily_dist,
	v_d=transx2$v_d,
	stringsAsFactors = F)#end of data.frame#
	is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
	cols<-(1:ncol(res1))[!is.na(res1)]
 
       if(permute==TRUE){
                for (y in cols){
                if(is.wholenumber(res1[,y])){res1[,y]<-round(rnorm(1,mean=res1[,y],sd=res1[,y]/variation),0)}else{
                        res1[,y]<-rnorm(1,mean=res1[,y],sd=res1[,y]/variation)}
                        if(res1[,y]<0){res1[,y]<-0}
                }#end of y loop
      	}
	if(quantiles==TRUE){
                for (y in  cols){
                if(is.wholenumber(res1[,y])){res1[,y]<-sample(round(qnorm(c(0.05,0.25,0.50,0.75,0.95),mean=res1[,y],sd=res1[,y]/variation ),0),1,prob = c(0.05,0.25,0.50,0.25,0.05))
				}else{
                  	res1[,y]<-sample(qnorm(c(0.05,0.25,0.50,0.75,0.95),mean=res1[,y],sd=res1[,y]/variation),1,prob = c(0.05,0.25,0.50,0.25,0.05))
				}
                if(res1[,y]<0){res1[,y]<-0}
                }#end of y loop
		}
        res1[res1<0]<-0
        if(res1$CFR>100){res1$CFR<-100}
	  if(res1$immunity>100){res1$immunity<-100}
        if(res1$chronic.carrier>100){res1$chronic.carrier<-100}
        if(res1$illness_length==0){res1$illness_length<-1}
	  if(res1$incubation==0){res1$incubation<-1}
	 if(res1$cases_per_year==0){res1$cases_per_year<-1}

res1$CFR=res1$CFR/100
res1$immunity=res1$immunity/100
res1$travel=1

return(res1)
}#end of function



####################################
########### SPILLOVER ###########

spillover<-function(dataset,R0,host_distance,d){
       namesx<-names(dataset)
       ## create full probs
	 dataset2<-dataset[dataset$susceptible>0 & dataset$host_density>0,]
	   ### gas model ### per person contact rate per day contact rate
	   dataset2$contactsx=round((dataset2$susceptible*0.67*0.8)*dataset2$host_density*d*sqrt(host_distance^2+(dataset2$daily_walking_distance)^2),0)#/dataset2$susceptible)
	   dataset2<-dataset2[dataset2$contactsx>0,]
	   dataset2$contactsx[dataset2$contactsx>quantile(dataset2$contactsx,0.95)]<-ceiling(quantile(dataset2$contactsx,0.95))### random ceiling

         if(nrow(dataset2)==0){return(dataset[,namesx])}## end of if

         #weight by density and vector efficiency # how to get to chance per susceptible
         no_infected<-function(x,susceptible,eff_cont_rate) sum(rbinom(n=susceptible[x],size=1,prob=eff_cont_rate[x]),na.rm=T) 
      
	   ### change to per day spillover#   nrow(dataset2)
	   #dataset2$spills<-sapply(1:nrow(dataset2), FUN=no_infected,contacts=dataset2$contacts,eff_cont_rate=spillover_rate[x],susceptible=dataset2$susceptible)
	   dataset2$spills2<-sapply(1:nrow(dataset2), FUN=no_infected,eff_cont_rate=dataset$SOR,susceptible=dataset2$contactsx)
	  #sum(dataset2$spills2)
        #dataset2$spills[dataset2$spills>dataset2$susceptible]<-dataset2$susceptible
	  if(sum(dataset2$spills2)==0){return(dataset[,namesx])}
	 dataset2<-dataset2[dataset2$spills2>0,]
	  sample_infected<-function(x,susceptible,infections) length(unique(sample(1:susceptible[x],infections[x],replace=T)))
        dataset2$spills<-sapply(1:nrow(dataset2), FUN=sample_infected,susceptible=dataset2$susceptible,infections= dataset2$spills2)
	 
       dataset$change<-0
       dataset[match(dataset2$cell.id,dataset$cell.id),"change"]<-dataset2$spills
       dataset[match(dataset2$cell.id,dataset$cell.id),"ID"]<-paste(dataset[match(dataset2$cell.id,dataset$cell.id),"ID"],paste(sample(c(0:9, letters, LETTERS),size=12, replace=TRUE),collapse=''),sep=";")
	dataset$ID<-gsub("NA;","",dataset$ID,fixed=TRUE)
       dataset$Rt[dataset$change>0 & dataset$Rt==0]<-R0 ### this should stop cell going back up to R0 for common diseases
       dataset$exposed<-dataset$exposed+dataset$change
       dataset$susceptible<-dataset$susceptible-dataset$change
	if(min(dataset$susceptible)<0){print("spillover");break}
      return(dataset[,namesx])
	   
}# end of spillover



##########################################
########### HUMAN TO HUMAN ###########
##############################################
##human-human within square
##morbidity=0.5;infected="exposed"

human_human<-function(dataset,infected,morbidity,average_per_person_total_interactions,illness_length){
       namesx<-names(dataset)
       dataset$sortx<-dataset[,infected]
       infect_sq<-dataset[dataset$sortx>0 & dataset$Rt>0,] ###all squares with infected people
       infect_sq<-infect_sq[infect_sq$susceptible>0,] ###all squares with people in ## what about 2050 # generalise
       if(nrow(infect_sq)==0){return(dataset[,namesx])}
	 
	### total pop of the square
        infect_sq$pop1<-rowSums(infect_sq[,c("susceptible","exposed","infectedA","infectedS","convalesence","spill","dead","recovered")])
  
       ### work effective contact rate for this group
       infect_sq$eff_cont_rate<-infect_sq$Rt/(average_per_person_total_interactions*illness_length) #### global r0 and then reduce to less sucepitbles
	 
         ### THIS IS DAILY EFF_CONTACT RATE so R0 over the chance that day than an effective contact will be made
	   ### DAILY EFF_CONTACT RATE = suc_infections/(tot_contacts when ill)/(average_per_person_per_day_interactions)
         ### DAILY EFF_CONTACT RATE = suc_infections/(daily_contacts*days)/days ##cancel days
	   ### DAILY EFF_CONTACT RATE = suc_infections/daily_contacts	
	   ### PER CONTACT EFF_CONTACT RATE = pp suc_infections/ pp tot_contacts when ill <-
	
	 infect_sq$eff_cont_rate[infect_sq$eff_cont_rate>1]<-1

	 #### cell specific contact rate - how many suceptibles will someone bump into fucked up
	 #infect_sq$per_day_interactions<-round(1*(infect_sq$susceptible*0.67*0.8)*0.0005*infect_sq$daily_walking_distance,0)## 0.67 due to groups  # 0.8 due to babies and old people not moving
	 infect_sq$per_day_interactions<-round(1*(infect_sq$susceptible)*0.67*0.8*0.005*sqrt((infect_sq$daily_walking_distance)^2+(infect_sq$daily_walking_distance)^2),0) 
       infect_sq$per_day_interactions[infect_sq$per_day_interactions>quantile(infect_sq$per_day_interactions,0.95)]<-ceiling(quantile(infect_sq$per_day_interactions,0.95))### random ceiling

       ### function people and contacts
       no_infected<-function(x,susceptible,eff_cont_rate,contacts) sum(rbinom(n=susceptible[x],size=contacts[x],prob=eff_cont_rate[x]),na.rm=T) 
       infect_sq$new_inf2<-sapply(1:nrow(infect_sq), FUN=no_infected,susceptible=infect_sq$sortx,eff_cont_rate=infect_sq$eff_cont_rate,contacts=infect_sq$per_day_interactions)
       if(sum(infect_sq$new_inf2)==0){return(dataset[,namesx])}

	 sample_infected<-function(x,susceptible,infections) length(unique(sample(1:susceptible[x],infections[x],replace=T)))
       infect_sq$new_inf<-sapply(1:nrow(infect_sq), FUN=sample_infected,susceptible=infect_sq$susceptible,infections=infect_sq$new_inf2)
    
 
	#change respective values
       infect_sq2<-infect_sq[infect_sq$new_inf>0,c("cell.id","new_inf")]
       dataset$new_inf<-0
       dataset[match(infect_sq2$cell.id,dataset$cell.id),"new_inf"]<-infect_sq2$new_inf

       #### always go to incubating even if for just one day
         dataset$susceptible<-dataset$susceptible-dataset$new_inf
	if(min(dataset$susceptible)<0){
			dd<-rep(0,length(dataset$susceptible))
			dd[dataset$susceptible<0]<-dataset$susceptible[dataset$susceptible<0]
			dataset$new_inf<-dataset$new_inf+dd
			dataset$susceptible[dataset$susceptible<0]<-0
					}
       dataset$exposed<-dataset$exposed+dataset$new_inf
	if(min(dataset$susceptible)<0){print("human_human");break}
       return(dataset[,namesx])
}# end of human2human

##########################################
########### CHANGE COMPARTMENT ###########

#### start
change_compart<-function(dataset=NULL,from=NULL,to=NULL,to2=NULL,compartment_duration=NULL,prob=0){
     namesx<-names(dataset)
     dataset$from=dataset[,from]
     dataset$prob<-prob
     dataset$compartment_duration<-compartment_duration   
     if(sum(dataset$from,na.rm=T)==0){return(dataset[,namesx])}
     change_comp<-function(x,from,prob) sum(rbinom(n=from[x],size=1,prob=prob[x]))  
     #change_comp<-cmpfun(change_comp)
     if(compartment_duration==Inf){dataset$decided<-dataset$from}else{
            dataset$decided<-sapply(1:nrow(dataset), FUN=change_comp, from=dataset$from,prob=dataset$compartment_duration)}
     if(sum(dataset$decided,na.rm=T)==0){return(dataset[,namesx])}
     dataset[,from]<-dataset[,from]-dataset$decided
     if(!is.null(to2)& max(prob,na.rm=T)>0){
          dataset$change<-sapply(1:nrow(dataset), FUN=change_comp, from=dataset$decided,prob=dataset$prob)
          dataset[,to2]<-dataset[,to2]+dataset$change
          if(sum(dataset$decided-dataset$change,na.rm=T)==0){return(dataset[,namesx])}
          dataset[,to]<-dataset[,to]+(dataset$decided-dataset$change)
     }else{
          dataset[,to]<-dataset[,to]+dataset$decided
     }##end of to2
	if(min(dataset$susceptible)<0){print("change_comp");break}
     return(dataset[,namesx])
}###end of change comp


##########################################
########### HUMAN MOVEMENT ###########

human_movementa<-function(dataset,infected,morbidity,flights,template,travel,R0){
     namesx<-names(dataset)
     if(morbidity==1){return(dataset)}
     dataset$sortx<-dataset[,infected]
     infect_sq<-dataset[dataset$sortx>0,] ###all squares with infected people
     if(nrow(infect_sq)==0){return(dataset[,namesx])}
    
     ### airlines ###
     #dataset$infectedS[sample(1:length(dataset$sortx),1000)]<-1
     results_a<-NULL
     ### run function to find out if any infected in squares
     infect_sq2<-infect_sq[!infect_sq$airport2==(-999),]
  	 if(!nrow(infect_sq2)==0){
          for(i in 1:nrow(infect_sq2)){
          #for(i in 1:nrow(infect_sq2)){
               infect_sq3<-infect_sq2[i,]
               infect_sq3$prob<-infect_sq3$flightsPersonYear *(infect_sq3$Rt/R0) ### weight by Rt such that airports closed down after infections in square
               ##infect_sq3$sortx<-round(infect_sq3$sortx*(1-morbidity),0)
               if(!nrow(infect_sq3[!infect_sq3$sortx==0,])==0) {
                   no_flights<-sum(rbinom(infect_sq3$sortx,1,prob=travel*infect_sq3$prob))
                   if(!no_flights==0){
                      flights2<-flights[flights$cell_id_start==infect_sq3$airport2,]
                      flights3<-flights[flights$cell_id_start %in% flights2$cell_id_dest,]
                      flights3<-flights3[!flights3$cell_id_dest==infect_sq3$cell.id,]
                      cell_dst<-sample(c(flights2$cell_id_dest,flights3$cell_id_dest),no_flights,replace=T)
                      cell_nos<-as.data.frame(table(cell_dst))#,stringsAsFactors=F)
					if(is.null(results_a)){results_a<-data.frame(cell_from=infect_sq3$cell.id,cell_to=levels(cell_nos$cell_dst),num=cell_nos$Freq,ID=infect_sq3$ID,stringsAsFactors=F)
					} else {
					 results_a<-rbind(results_a,data.frame(cell_from=infect_sq3$cell.id,cell_to=levels(cell_nos$cell_dst),ID=infect_sq3$ID,num=cell_nos$Freq,stringsAsFactors=F))
					}  
                    }
               }
           }
     }#end of aero loops

 
     if(!is.null(results_a)){
          res1<-results_a
          res2a<-rbind(data.frame(ID=res1$ID,cellx=as.numeric(res1[,2]),change=as.numeric(res1[,3])),data.frame(ID=res1$ID,cellx=as.numeric(res1[,1]),change=as.numeric(res1[,3])*(-1)))
          res2<-aggregate(res2a$change,by=list(res2a$cellx,res2a$ID),FUN='sum')
          res3<-res2[res2$Group.1 %in% dataset$cell.id,] ### res3 subsetting res2 by where cell.ids are in dataset
          res4<-res2[!res2$Group.1 %in% dataset$cell.id,] ### res3 subsetting res2 by where cell.ids are in dataset
		if(nrow(res4)>0){dataset[1,"export"]<-paste(dataset[1,"export"],paste(as.vector(res4$Group.1),collapse=";"),sep=";");dataset$export[1]<-gsub("NA;","",dataset$export[1],fixed=TRUE)}
          dataset$change<-0
          dataset[match(res3$Group.1,dataset$cell.id),"change"]<-res3$x
          dataset[match(res3$Group.1,dataset$cell.id),"ID"]<-as.vector(res3$Group.2)
	    dataset[,infected]<-(dataset[,infected])+dataset$change
		dataset[,infected][(dataset[,infected])<0]<-0
          #if(min(dataset[,infected])<0){print("air_movement");break}
     return(dataset[,namesx])
     }##end if null
     return(dataset[,namesx])
}#end of human movement

####roads
human_movementr<-function(dataset,infected,morbidity,rm4,template,travel,R0){
     namesx<-names(dataset)
     if(morbidity==1){return(dataset)}
     dataset$sortx<-dataset[,infected]
     infect_sq<-dataset[dataset$sortx>0,] ###all squares with infected people
     if(nrow(infect_sq)==0){return(dataset[,namesx])}
  

     ####only a thousand but cells with more infections have higher chance of being chosen
	results_a=NULL
	
	roadtr2<-function(x,from2,rm4,infected,travel){
		ce<-from2$cell.id[x]
		zzz<-rm4[rm4$from==ce,]
		return(data.frame(cell_from=ce,cell_to=(sample(strsplit(zzz$to,";")[[1]],from2[x,infected],replace=T,prob=travel*as.numeric(strsplit(zzz$probability,";")[[1]])))))
		}
 
 
	#######################
     ####### roads ########
     	infect_sq2<-infect_sq[!is.na(infect_sq$roadn),]
	infect_sq2<-infect_sq2[infect_sq2$cell.id %in% rm4$from,]

   if(nrow(infect_sq2)>0){
     if(nrow(infect_sq2)>10000){infect_sq2<-infect_sq2[sample(1:nrow(infect_sq2),10000,prob=infect_sq2$sortx,replace=F),]}
       
    #### distance to other nodes in road network
           moves<-sapply(1:nrow(infect_sq2), FUN=roadtr2,from2=infect_sq2,infected=infected,rm4=rm4,travel=travel,simplify=FALSE)
		dff <- do.call("rbind", moves)
		dff$cell_to<-suppressWarnings(as.numeric(as.vector(dff$cell_to)))
		dff$num=1
		results_a<-dff[!is.na(dff$cell_to),]
		if(nrow(results_a)==0){results_a=NULL}else{
		results_a$ID<-dataset[match(results_a$cell_from,dataset$cell.id),"ID"]}
    }
      
     if(!is.null(results_a)){
         res1<-results_a
          res2a<-rbind(data.frame(ID=res1$ID,cellx=as.numeric(res1[,2]),change=as.numeric(res1[,3])),data.frame(ID=res1$ID,cellx=as.numeric(res1[,1]),change=as.numeric(res1[,3])*(-1)))
          res2<-aggregate(res2a$change,by=list(res2a$cellx,res2a$ID),FUN='sum')
          res3<-res2[res2$Group.1 %in% dataset$cell.id,] ### res3 subsetting res2 by where cell.ids are in dataset
          res4<-res2[!res2$Group.1 %in% dataset$cell.id,] ### res3 subsetting res2 by where cell.ids are in dataset
		#print(res3)
		#if(nrow(res4)>0){dataset[1,"export"]<-paste(dataset[1,"export"],paste(res4$Group.1,collapse=";"),sep=";");dataset$export[1]<-gsub("NA;","",dataset$export[1],fixed=TRUE)}
          dataset$change<-0
         dataset[match(res3$Group.1,dataset$cell.id),"change"]<-res3$x
         dataset[match(res3$Group.1,dataset$cell.id),"ID"]<-as.vector(res3$Group.2)
	     dataset[,infected]<-(dataset[,infected])+dataset$change
		dataset[,infected][(dataset[,infected])<0]<-0
      #if(min(dataset[,infected])<0){print("road_movement");break}
     return(dataset[,namesx])
     }##end if null
     return(dataset[,namesx])
}#end of human movement


####walking
human_movementw<-function(dataset,infected,morbidity,template,travel,R0){
     namesx<-names(dataset)
     if(morbidity==1){return(dataset)}
     dataset$sortx<-dataset[,infected]
     infect_sq2<-dataset[dataset$sortx>0,] ###all squares with infected people
     if(nrow(infect_sq2)==0){return(dataset[,namesx])}
  	
	walk1<-function(x) sum(rbinom(x,1,prob=0.05))

     infect_sq2$travels<-sapply(infect_sq2$sortx,FUN=walk1)
	infect_sq2<-infect_sq2[infect_sq2$travels>0,] ###all squares with walkers people

     if(nrow(infect_sq2)==0){return(dataset[,namesx])}

     ####only a thousand but cells with more infections have higher chance of being chosen
	results_a=NULL
	
	walk2<-function(x,from2,template,infected,travel){
		ce<-from2$cell.id[x]
		zzz<-adjacent(template,ce)
		return(data.frame(cell_from=ce,cell_to=(sample(x=as.vector(zzz[,2]),from2[x,infected],replace=TRUE,rep(1/nrow(zzz),nrow(zzz))))))
		} 
	
     if(nrow(infect_sq2)>10000){infect_sq2<-infect_sq2[sample(1:nrow(infect_sq2),10000,prob=infect_sq2$sortx,replace=F),]}

           	moves<-sapply(1:nrow(infect_sq2), FUN=walk2,from2=infect_sq2,infected="travels",travel=travel,template=template,simplify=FALSE)
		dff <- do.call("rbind", moves)
		dff$cell_to<-suppressWarnings(as.numeric(as.vector(dff$cell_to)))
		dff$num=1
		results_a<-dff[!is.na(dff$cell_to),]
		if(nrow(results_a)==0){results_a=NULL}
		
     if(!is.null(results_a)){
          res1<-results_a
          res2a<-rbind(data.frame(ID=res1$ID,cellx=as.numeric(res1[,2]),change=as.numeric(res1[,3])),data.frame(ID=res1$ID,cellx=as.numeric(res1[,1]),change=as.numeric(res1[,3])*(-1)))
          res2<-aggregate(res2a$change,by=list(res2a$cellx,res2a$ID),FUN='sum')
          res3<-res2[res2$Group.1 %in% dataset$cell.id,] ### res3 subsetting res2 by where cell.ids are in dataset
          #res4<-res2[!res2$Group.1 %in% dataset$cell.id,] ### res3 subsetting res2 by where cell.ids are in dataset
		#if(nrow(res4)>0){dataset[1,"export"]<-paste(dataset[1,"export"],paste(res4$Group.1,collapse=";"),sep=";");dataset$export[1]<-gsub("NA;","",dataset$export[1],fixed=TRUE)}
          dataset$change<-0
          dataset[match(res3$Group.1,dataset$cell.id),"change"]<-res3$x
          dataset[match(res3$Group.1,dataset$cell.id),"ID"]<-as.vector(res3$Group.2)
	    dataset[,infected]<-(dataset[,infected])+dataset$change
    if(min(dataset[,infected])<0){print("walking_movement");break}
		return(dataset[,namesx])
     }else{return(dataset[,namesx])}
  
}#end of human movement

#table(dataset$ID)

##########################################
run_script<-function(dataset,rm4,flights,transx,template){
#for (i in 1:100){
	#dx2<-dataset
     #library(parallel)
     #cl <- makeCluster(4)
     #print(i)
     #SPILLOVER
     expos<-dataset$exposed
     dataset<-spillover(dataset=dataset,R0=transx$R0,host_distance=transx$host_distance,d=transx$d)
     dataset$spill<-dataset$spill+(dataset$exposed-expos)

     # EXPOSED
     if(sum(dataset$exposed)>0){
        if(!transx$R0==0){
   	 	dataset<-human_movementa(dataset=dataset,infected="exposed",morbidity=0,template=template,flights=flights,R0=transx$R0,travel=transx$travel)
      	 dataset<-human_movementr(dataset=dataset,infected="exposed",morbidity=0,template=template,rm4=rm4,R0=transx$R0,travel=transx$travel)}## end of R0=0
          if(transx$chronic.carrier>0){
           dataset<-change_compart(dataset=dataset,from="exposed",to="infectedS",to2="infectedA",compartment_duration=(1/transx$incubation),prob=transx$chronic.carrier/100)
       }else{ 
           dataset<-change_compart(dataset=dataset,from="exposed",to="infectedS",compartment_duration=1/transx$incubation)
       }##end of chronic
     }
#	print(colSums(dataset[,c("susceptible","exposed","infectedA", "infectedS","convalesence","spill","dead","recovered") ] ,na.rm=T))
      

     # INFECTED SYMPTOMATIC
     if(sum(dataset$infectedS)>0){
       if(!transx$R0==0){dataset<-human_human(dataset,infected="infectedS",morbidity=0.5,average_per_person_total_interactions=transx$average_per_person_total_interactions,illness_length=transx$illness_length)
       dataset<-human_movementa(dataset=dataset,infected="infectedS",morbidity=0.5,template=template,flights=flights,R0=transx$R0,travel=transx$travel)
       dataset<-human_movementr(dataset=dataset,infected="infectedS",morbidity=0.5,template=template,rm4=rm4,R0=transx$R0,travel=transx$travel)}## end of R0=0
       dataset<-change_compart(dataset=dataset,from="infectedS",to="convalesence",to2="dead",compartment_duration=1/transx$illness_length,prob=dataset$CFR)
     }
     
     # INFECTED CHRONIC/ASYMPTOMATIC
     if(sum(dataset$infectedA)>0){
       if(!transx$R0==0){dataset<-human_human(dataset,infected="infectedA",morbidity=transx$chronic.morbidity/100,average_per_person_total_interactions=transx$average_per_person_total_interactions,illness_length=transx$illness_length)
       dataset<-human_movementa(dataset=dataset,infected="infectedA",morbidity=0.5,template=template,flights=flights,R0=transx$R0,travel=transx$travel)
       dataset<-human_movementr(dataset=dataset,infected="infectedA",morbidity=0.5,template=template,rm4=rm4,R0=transx$R0,travel=transx$travel)}## end of R0=0
       dataset<-change_compart(dataset=dataset,from="infectedA",to="convalesence",compartment_duration=1/transx$chronic_length)
     }
     
     # CONVALESCENCE
     if(sum(dataset$convalesence)>0){
     dataset<-change_compart(dataset=dataset,from="convalesence",to="recovered",to2="susceptible",compartment_duration=1/transx$convalescence_length,prob=1-(transx$immunity))
     }
     
     # FUNERAL
     #if(sum(dataset$funeral)>0){
     #       dataset<-human_human(dataset,infected="funeral",morbidity=0,transx=transx)
     #       dataset<-change_compart(dataset=dataset,from="funeral",to="dead",compartment_duration=1/transx$funeral)
     #}
     dataset$Rt[dataset$infectedS>0|dataset$infectedA>0]<-dataset$Rt[dataset$infectedS>0|dataset$infectedA>0]-dataset[dataset$infectedS>0|dataset$infectedA>0,"Rtgrad"]
     dataset$Rt[dataset$Rt<0]<-0

	#print(colSums(dataset[,c("susceptible","exposed","infectedA", "infectedS","convalesence","spill","dead","recovered") ] ,na.rm=T))
     return(dataset)
}# end run script function

#stopCluster(cl)
