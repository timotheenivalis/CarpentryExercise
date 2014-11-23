################ Basic initializations ################


################ Base population parameters
MeanBirthSize<-10
meanGrowth<-prod(runif(10000,0.9,1.2))^(1/10000)
MeanRepro<-2
survivalPenaltyForRepro<-0

SDZ<-1
SDH<-1
############### Genetic determinisms ###############
############### Genetic determinism Z
dominance<-1 # for additive effects only, must be 0
overdominance<-0 # non-null values generate overdominance

nbLoci<-10 #number of loci controling the trait phenotype
nbAlleles<-10 #number of existing alleles per loci

gvaluesZ<-array(data=NA,dim=c(nbAlleles,nbAlleles,nbLoci),dimnames=list(paste("A",1: nbAlleles,sep=""),paste("A",1: nbAlleles,sep=""),paste("L",1:nbLoci,sep=""))) # Initialising a matrix that will contain the genotypic effects on the/a trait

for(L in 1:nbLoci)
{
  # Setting the effects for the homozygotes [all loci]
  effect<-abs(rnorm(n=1,mean=0,sd=SDZ))# alter the locus importance in a realistic way (many small-effect loci, few major loci)
  diag(gvaluesZ[,,L])<-2*rnorm(n=dim(gvaluesZ)[1],mean=0,sd=effect)
  # Setting the effects for the heterozygotes
  for(A in 1:(nbAlleles-1))# loop for off-diagonal = heterozygotes (additive and dominance effects)
  {
    for (D in (A+1):nbAlleles)
    {
      d<-dominance*runif(n=1,min=-0.5-overdominance,max=0.5+overdominance)
      gvaluesZ[A,D,L]<-(0.5-d)*gvaluesZ[A,A,L]+(0.5+d)*gvaluesZ[D,D,L] # mean of additive effects + dominance, over diagonal
      gvaluesZ[D,A,L]<-(0.5-d)*gvaluesZ[A,A,L]+(0.5+d)*gvaluesZ[D,D,L] # the same below diagonal    
    }
  }
}
############### Genetic determinism H 
dominance<-1 # for additive effects only, must be 0
overdominance<-0 # non-null values generate overdominance

nbLoci<-10 #number of loci controling the trait phenotype
nbAlleles<-10 #number of existing alleles per loci

gvaluesH<-array(data=NA,dim=c(nbAlleles,nbAlleles,nbLoci),dimnames=list(paste("A",1: nbAlleles,sep=""),paste("A",1: nbAlleles,sep=""),paste("L",1:nbLoci,sep=""))) # Initialising a matrix that will contain the genotypic effects on the/a trait

for(L in 1:nbLoci)
{
  # Setting the effects for the homozygotes [all loci]
  effect<-abs(rnorm(n=1,mean=0,sd=SDH))# alter the locus importance in a realistic way (many small-effect loci, few major loci)
  diag(gvaluesH[,,L])<-2*rnorm(n=dim(gvaluesH)[1],mean=0,sd=effect)
  # Setting the effects for the heterozygotes
  for(A in 1:(nbAlleles-1))# loop for off-diagonal = heterozygotes (additive and dominance effects)
  {
    for (D in (A+1):nbAlleles)
    {
      d<-dominance*runif(n=1,min=-0.5-overdominance,max=0.5+overdominance)
      gvaluesH[A,D,L]<-(0.5-d)*gvaluesH[A,A,L]+(0.5+d)*gvaluesH[D,D,L] # mean of additive effects + dominance, over diagonal
      gvaluesH[D,A,L]<-(0.5-d)*gvaluesH[A,A,L]+(0.5+d)*gvaluesH[D,D,L] # the same below diagonal    
    }
  }
}


############### Definition of the class ############
setClass(
  Class="Leprechaun",
  representation=representation(
    ID = "integer",
    pID = "integer",
    age = "integer",
    Birth = "integer",
    alive = "logical",
    size = "numeric",
    hunting = "numeric",
    camemberts = "integer",
    sex = "character",
    DNAZ = "matrix",
    DNAH = "matrix",
    bvs = "numeric",
    bvh = "numeric",
    ARS = "integer"
  )  
)

############### Definition of the basic methods (for printing to the screen and initialisation)
setMethod("show","Leprechaun",
          function(object){
            cat(object@ID,"\t",object@size,"\t",object@camemberts,"\t",object@age,"\t",object@sex,"\t(",object@pID[1],",",object@pID[2],")\t",object@Birth,"\t",object@alive,"\n",sep="")
          }
)
setMethod("initialize","Leprechaun",function(.Object,parent1,parent2){
  .Object@DNAZ<-matrix(NA,nrow=2,ncol=nbLoci)
  .Object@DNAH<-matrix(NA,nrow=2,ncol=nbLoci)
  if(missing(parent1)){
    parent1<-NA
    .Object@DNAZ[1,]<-floor(runif(nbLoci,min=1,max=nbAlleles+1))
    .Object@DNAH[1,]<-floor(runif(nbLoci,min=1,max=nbAlleles+1))
  }else{
    #weight1<-pop[[parent1]]@size
    .Object@DNAZ[1,]<-pop[[parent1]]@DNAZ[cbind(floor(runif(n=nbLoci,min=1,max=3)),1:nbLoci)]
    .Object@DNAH[1,]<-pop[[parent1]]@DNAH[cbind(floor(runif(n=nbLoci,min=1,max=3)),1:nbLoci)]
  }
  if(missing(parent2)){
    parent2<-NA
    .Object@DNAZ[2,]<-floor(runif(nbLoci,min=1,max=nbAlleles+1))
    .Object@DNAH[2,]<-floor(runif(nbLoci,min=1,max=nbAlleles+1))
  }else{
    .Object@DNAZ[2,]<-pop[[parent2]]@DNAZ[cbind(floor(runif(n=nbLoci,min=1,max=3)),1:nbLoci)]
    .Object@DNAH[2,]<-pop[[parent2]]@DNAH[cbind(floor(runif(n=nbLoci,min=1,max=3)),1:nbLoci)]
  }
  .Object@age<-as.integer(0)
  .Object@ID<-CID
  .Object@pID<-c(as.integer(parent1),as.integer(parent2))
  .Object@Birth<-as.integer(YR)
  .Object@alive<-TRUE
  BreedingValueSize<-0
  for (Locus in 1:nbLoci)#take the mean of genetic values
  {
    BreedingValueSize<-BreedingValueSize+(gvaluesZ[ .Object@DNAZ[1,Locus], .Object@DNAZ[2,Locus], Locus]/nbLoci)
  }
  .Object@bvs<-BreedingValueSize
  size<-MeanBirthSize+BreedingValueSize
  .Object@size<-rnorm(n=1,mean=size,sd=0) # sd plasticity birth size
  
  BreedingValueHunting<-0
  for (Locus in 1:nbLoci)#take the mean of genetic values
  {
    BreedingValueHunting<-BreedingValueHunting+(gvaluesH[ .Object@DNAH[1,Locus], .Object@DNAH[2,Locus], Locus]/nbLoci)
  }
  .Object@bvh<-BreedingValueHunting
  .Object@hunting<-rnorm(n=1,mean=.Object@bvh,sd=0.5) # sd plasticity hunting quality
  
  .Object@camemberts<-as.integer(0)
  
  .Object@ARS<-as.integer(0)#annual reproductive success
  
  if(runif(1)>0.5){.Object@sex<-'F'}else{.Object@sex<-'M'}
  
  CID<<-as.integer(CID+1)
  
  return(.Object)
})

############## Definition of more biologically relevant methods ##############
# Implementing the famous bathtub, ages 1 to 20
bathtub<-function(age){
  p<-0.6*exp(-age/4)+(-1+exp(age*log(2)/20))
  p[p>1]<-1
  return(p)
}
# Incorporate the effect of size to the survival function
sizeSurvival<-function(age,size,camemberts){
  sizedeviation<-size-(MeanBirthSize*meanGrowth^age)
  p<-bathtub(age)
  if(p<1)#because size does not prevent animals of maximal age to die out
  {
    plogit<-log(p/(1-p))
    Philogit<-plogit-SurvivalSelection1*sizedeviation-SurvivalSelection2*sizedeviation^2
    p<-exp(Philogit)/(1+exp(Philogit))
    if(camemberts<100)
    {
      p<-1-(1-p)*(camemberts/100)
    }
  }
  return(1-p)
}

# Applying the bathtub in a surival function
setGeneric("Surv",function(Object){standardGeneric("Surv")})

setMethod("Surv","Leprechaun",function(Object){
  
  if(runif(1)>sizeSurvival(Object@age,Object@size,Object@camemberts)){
    Object@alive<-FALSE
    DEAD<<-c(DEAD,Object@ID)
  }
  
  return(Object)
  
})
# Simple function, simply adds 1 to the age
setGeneric("Age",function(Object){standardGeneric("Age")})

setMethod("Age","Leprechaun",function(Object){
  Object@age<-as.integer(Object@age+1)
  return(Object)
})

# Growth // Sizes change, this fact is known by many - if not all
setGeneric("Grow",function(Object){standardGeneric("Grow")})

setMethod("Grow","Leprechaun",function(Object){
  Object@size<-Object@size*runif(1,lowBoundGrowth,highBoundGrowth)
  return(Object)
})

# Retrieving the ID of an individual
setGeneric("IDretrieve",function(Object){standardGeneric("IDretrieve")})

setMethod("IDretrieve","Leprechaun",function(Object){
  return(Object@ID)
})

# Retrieving the size of an individual
setGeneric("Size",function(Object){standardGeneric("Size")})

setMethod("Size","Leprechaun",function(Object){
  return(Object@size)
})

# Retrieving the hunting quality of an individual
setGeneric("Hunting",function(Object){standardGeneric("Hunting")})

setMethod("Hunting","Leprechaun",function(Object){
  return(Object@hunting)
})

# Retrieving the sex of an individual
setGeneric("Sex",function(Object){standardGeneric("Sex")})

setMethod("Sex","Leprechaun",function(Object){
  return(Object@sex)
})

# Calculating the number of offspring for a females
setGeneric("Num_off",function(Object){standardGeneric("Num_off")})

setMethod("Num_off","Leprechaun",function(Object){
  sizeDeviation<-Object@size-MeanBirthSize*meanGrowth^Object@age    
  lambda<-exp(log(MeanRepro)+fertilitySelection1*sizeDeviation+fertilitySelection2*sizeDeviation^2+camembertSelection*((Object@camemberts)^(1/3)-survivalPenaltyForRepro))
  repro<-rpois(n=1,lambda=lambda)
  Object@ARS<-as.integer(repro)
  return(Object)
})

# Function for camembert attributions 
setGeneric("Food",function(Object,Camams){standardGeneric("Food")})

setMethod("Food","Leprechaun",function(Object,Camams){
  Object@camemberts<-as.integer(Camams)
  return(Object)
})
 

############### Function for printing values
print_info<-function(YR=YR,ALIVE=ALIVE,CID=CID,camembert=camembert){
  cat("\nAt the beginning of year:",YR,"\nThere are:",length(ALIVE),"Leprechauns\n-----------------\n")
  cat("ALIVE:",ALIVE,"\n")
  cat("Current ID:",CID,"\n")
  cat("\n-----------------\nThe camembert production is:",camembert,"\n")
  return()
}

