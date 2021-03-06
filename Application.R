setwd(dir = "C:/Users/Timoth�e/Documents/GitHub/CarpentryExercise/")
library("S4imul")

################ Random seed
set.seed(12)
################ Counter for the IDs
CID<-as.integer(1) 
################ Counter for the current year
YR<-0 
############### Creating an initial population with 10 individuals
pop<-replicate(n = 30,expr = new("Leprechaun"))
############### List of living individuals [their indices], this will save time later, because dead individuals are not looped over
ALIVE<-1:length(pop)
filename<-"pop.txt"

cat("t\tID\tz\tbvs\thunting\tbvh\tC\ts\tARS\tage\tp1\tp2\tphi",file=filename,append=FALSE)

MeanCamembert<-5000
SDCamembert<-10
SurvivalSelection1<-1
SurvivalSelection2<-0
fertilitySelection1<-0.1
fertilitySelection2<-0
camembertSelection<-0.5
lowBoundGrowth<-0.9
highBoundGrowth<-1.2

for(YR in 1:10){
  camembert<-abs(round(rnorm(n=1,mean=MeanCamembert,sd=SDCamembert),digits=0)) # resources for year YR
  
  print_info(YR,ALIVE,CID,camembert)
  
  #### Competition for resources
  HuntingAlive<-as.numeric(lapply(pop[ALIVE],Hunting))
  HunterQualities<- abs(HuntingAlive - mean(HuntingAlive)) # here the most original individuals have a strong advantage in the competition (f-dpd selection)
  if(sum(HunterQualities)==0){HunterQualities=rep(x=1,length(HunterQualities))}
  podium<-table(factor(sample(as.character(ALIVE),size=camembert,replace=T,prob=HunterQualities),levels=ALIVE))
  camams<-as.numeric(podium[match(ALIVE,names(podium))])
  
  pop[ALIVE]<-lapply(1:length(ALIVE), function(x) Food(pop[[ALIVE[x]]],camams[x]))
  
  #### Survival
  DEAD<-c()
  pop[ALIVE]<-lapply(pop[ALIVE],Surv)
  ALIVE<-ALIVE[!(ALIVE %in% DEAD)]
  
  if(length(ALIVE)==0){
    for(i in DEAD){
      cat("\n",YR,"\t",pop[[i]]@ID,"\t",pop[[i]]@size,"\t",pop[[i]]@bvs,"\t",pop[[i]]@hunting,"\t",pop[[i]]@bvh,"\t",pop[[i]]@camemberts,"\t",pop[[i]]@sex,"\t",pop[[i]]@ARS,"\t",pop[[i]]@age,"\t",pop[[i]]@pID[1],"\t",pop[[i]]@pID[2],"\t",0,file=filename,append=TRUE)
    }
    break
  }
  #### Age+1 and growth
  pop[ALIVE]<-lapply(pop[ALIVE],Age)
  pop[ALIVE]<-lapply(pop[ALIVE],Grow)
  
  #### Reproduction 
  
  ##########
  ### Part dedicated to retrieving the indices of all living males and of all living females
  ##########
  males<-lapply(pop,Sex)=="M" # Determine which individuals are males -- logicaly all other individuals should be females... However, this includes dead ones... Simply a list of T,T,F,F,T,F,F,....
  females<-which(!males) # Get the indices of the non-males (that is females..)
  males<-which(males)   # Get the indices of the males
  females<-intersect(females,ALIVE) # Retrieve the indices of the living(!) females
  males <-intersect(males,ALIVE) # Retrieve the indises of the living males
  #cat(females)
  
  ##########
  # Part dedicated to breeding.. 
  ##########
  
  # We take a female based approach: we determine for each females 
  from<-CID
  pop[females]<-lapply(pop[females],Num_off)
  for(i in females){
    Noffs<-pop[[i]]@ARS
    if(Noffs>0 & length(males>0)){
      #Determine the father
      fat<-sample(males,1)
      for(j in 1:Noffs){
        #cat(j,"\n")
        # Create the offspring
        pop<-c(pop,new("Leprechaun",parent1=i,parent2=fat))
      }
    }		
  }
  
  if(from!=CID){
    ALIVE<-c(ALIVE,(from):(CID-1))
  }
  dfwrite(pop=pop,ALIVE=ALIVE,DEAD=DEAD)
}


