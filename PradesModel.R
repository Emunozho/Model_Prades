library(SoilR)
library(FME)
library(RColorBrewer)
library(expm)
library(parallel)


setwd("~/Documents/MSCA/DISEQ/Prades model") # set path.


# -----------------------------------------------------------------------------------------------------------------------------------------------------------
# Cores
# -----------------------------------------------------------------------------------------------------------------------------------------------------------
# Calculate the number of cores
no_cores=detectCores() #- 10

# Initiate cluster
cl <- makeCluster(no_cores)

# -----------------------------------------------------------------------------------------------------------------------------------------------------------
# Functions
# -----------------------------------------------------------------------------------------------------------------------------------------------------------
M=function(t,B,u){expm(t*B)%*%u} # Mass carbon.
Rt=function(t, B, u){diag(colSums(-B))%*%expm(t*B)%*%u} # Respiratory losses.

# create natrix B,
makeB=function(pars){
  B=diag(-1*pars[1:7])
  B[2,1]=pars["alpha21"]*pars["k1"]
  B[3,1]=pars["alpha31"]*pars["k1"]
  B[4,1]=pars["alpha41"]*pars["k1"]
  B[5,2]=pars["alpha52"]*pars["k2"]
  B[6,2]=pars["alpha62"]*pars["k2"]
  B[4,3]=pars["alpha43"]*pars["k3"]
  B[6,3]=pars["alpha63"]*pars["k3"]
  B[7,3]=pars["alpha73"]*pars["k3"]
  B[6,4]=pars["alpha64"]*pars["k4"]
  B[6,5]=pars["alpha65"]*pars["k5"]
  B[7,6]=pars["alpha76"]*pars["k6"]
  return(B)
}

ecoModel=function(pars,GPP){
  mod=Model(t=yr,A=makeB(pars),ivList=P0,inputFluxes=c(GPP,rep(0,6)),pass=TRUE)
  Ct=getC(mod)
  colnames(Ct)<-c("Fol","Wood","FL","CWD","SA","SB")
  return(data.frame(time=yr, Ct))
}

ecoCost=function(pars,GPP){
  out=ecoModel(pars,GPP)
  cost1=modCost(model=out,obs=obsFol,weight="none")
  cost2=modCost(model=out,obs=obsWood,cost=cost1,weight="none")
  cost3=modCost(model=out,obs=obsFL,cost=cost2,weight="std")
  cost4=modCost(model=out,obs=obsCWD,cost=cost3,weight="mean")
  cost5=modCost(model=out,obs=obsSC,cost=cost4,weight="none")
  return(cost5)
}

randomFit=function(GPP){
  ecoFit=modFit(f=ecoCost,p=inipars,GPP=GPP,lower=0,upper=c(rep(3,7),rep(1,9)), method="Marq")
  cp=coef(ecoFit)
  return(cp)
}

meanSATTs=function(u, B){
  x=-1*solve(B)%*%u
  SA=sum(-1*solve(B)%*%x)/(sum(x))
  TT=sum(x)/sum(u)
  return(c(x=sum(x),SA=SA,TT=TT))
}

makeMod=function(pars,GPP){
  B=makeB(pars)
  u=matrix(c(GPP,rep(0,6)),7,1)
  SATT=meanSATTs(u,B)
  return(SATT)
}

allObjects=ls()
clusterExport(cl, varlist=c(allObjects, "modFit", "ode","modCost","Model", "getC"))


# -----------------------------------------------------------------------------------------------------------------------------------------------------------
# Read and prepare data
# -----------------------------------------------------------------------------------------------------------------------------------------------------------
treat='Control' # Treatment: select 'Control' or 'Drought'
if (treat=='Control'){
  Cdata<-read.csv("Data/Prades_data_control.csv") 
} else if (treat=='Drought'){
  Cdata<-read.csv("Data/Prades_data_drought.csv")
} else {
  print("Select a treatment")
}

poolnames=c("Foliage", "Wood", "Roots", "Thin litter", "Thick litter", "Soil carbon A", "Soil carbon B")
poolnamesd=c("Foliage", "Wood", "Thin litter", "Thick litter", "Soil carbon A", "Soil carbon B")

pdf("Stocks.pdf",5,10)
par(mfrow=c(4,2))
for (i in 3:8){
  plot(x=Cdata[,1],y=Cdata[,i],
     col="#69b3a2",
     pch=20,
     cex=1,
     xlab="time [years]", ylab=expression("C stocks " ~(MgC~ha^{-1})),
     main=poolnamesd[i-2])
}
dev.off()

#frac_foliage=0.08; frac_wood=0.19+0.73

GPP=47.63 # [Mg ha-1] Interpolation for 2024 from data reported in Table 5 by Sabaté (2002).
GPP_sd=6.91 # [Mg ha-1] Interpolation for 2024 from data reported in Table 5 by Sabaté (2002).

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Cómo se sacaron estos estadísticos? datos diferentes a los que hay en el archivo PorceB2012.csv?

# Random variates for primary forests. Mean and sd values from Sierra et al. (2007, FEM)
nplots=33
set.seed(11); primaryplotAges=sample(100:150,nplots,replace=TRUE)
set.seed(22); folp=rnorm(nplots,mean=mean(Cdata$Leaves),sd=sd(Cdata$Leaves))
set.seed(33); woodp=rnorm(nplots,mean=mean(Cdata$Wood),sd=sd(Cdata$Wood))
#set.seed(44); finerootp=rnorm(nplots,mean=0.45*17.4, sd=0.45*2.6)
#set.seed(55); coarserootp=rnorm(nplots,mean=0.45*67.1, sd=0.45*16.4)
set.seed(66); finelitterp=rnorm(nplots,mean=mean(Cdata$Thin.litter,na.rm=TRUE),sd=sd(Cdata$Thin.litter,na.rm=TRUE))
set.seed(77); cwdp=rnorm(nplots,mean=mean(Cdata$Thick.litter,na.rm=TRUE),sd=sd(Cdata$Thick.litter,na.rm=TRUE))
set.seed(88); soilcAp=rnorm(nplots,mean=mean(Cdata$Soil.Hz.A,na.rm=TRUE),sd=sd(Cdata$Soil.Hz.A,na.rm=TRUE))
set.seed(99); soilcBp=rnorm(nplots,mean=mean(Cdata$Soil.Hz.B,na.rm=TRUE),sd=sd(Cdata$Soil.Hz.B,na.rm=TRUE))





# Observations from secondary forest chronosequence
obsFols=data.frame(time=Cdata$Year,Fol=Cdata$Leaves)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
obsWoods=data.frame(time=Cdata$Year,Wood=Cdata$Wood)
obsFLs=data.frame(time=Cdata$Year,FL=Cdata$Thin.litter)
obsCWDs=data.frame(time=Cdata$Year,CWD=Cdata$Thick.litter)
obsSAs=data.frame(time=Cdata$Year,SA=Cdata$Soil.Hz.A)
obsSBs=data.frame(time=Cdata$Year,SB=Cdata$Soil.Hz.B)
obsSTs=data.frame(time=Cdata$Year,ST=Cdata$Soil.Hz.A+Cdata$Soil.Hz.B)


# Merged observations
obsFol=rbind(obsFols,data.frame(time=primaryplotAges,Fol=folp))
obsWood=rbind(obsWoods,data.frame(time=primaryplotAges,Wood=woodp))
#obsFR=rbind(obsFRs,data.frame(time=primaryplotAges,FR=finerootp))
#obsCR=rbind(obsCRs,data.frame(time=primaryplotAges,CR=coarserootp))
obsFL=rbind(obsFLs,data.frame(time=primaryplotAges,FL=finelitterp))
obsCWD=rbind(obsCWDs,data.frame(time=primaryplotAges,CWD=cwdp))
obsSA=rbind(obsSAs,data.frame(time=primaryplotAges,SA=soilcAp))
obsSB=rbind(obsSBs,data.frame(time=primaryplotAges,SB=soilcBp))

yr=seq(0,150,0.5)
# Initial values
P0=c(Fol=0,Wood=0,RO=0,FL=0,CWD=0,SA=mean(obsSAs[,2],na.rm=TRUE),SB=mean(obsSBs[,2],na.rm=TRUE))

inipars=c(k1=0.5,k2=0.3,k3=0.1,k4=0.1,k5=0.1,k6=0.1,k7=0.1,
          alpha21=0.1,alpha31=0.1,alpha41=0.1,alpha51=0.1,alpha53=0.1,alpha62=0.1,alpha64=0.1,alpha75=0.1,alpha76=0.1)

inipars=c(k1=0.5,k2=0.3,k3=0.1,k4=0.1,k5=0.1,k6=0.1,k7=0.1,
          alpha21=0.1,alpha31=0.1,alpha41=0.1,alpha52=0.1,alpha62=0.1,alpha43=0.1,alpha63=0.1,alpha73=0.1,alpha64=0.1,alpha65=0.1,alpha76=0.1)

tau=seq(0,400,by=0.05)
pal=brewer.pal(7,"Dark2")


####
n=1000
set.seed(123); gpp=rnorm(n,mean=GPP,sd=GPP_sd)

pdf("Figures/GPPhist.pdf")
hist(gpp)
dev.off()

#----------------------------------------------------->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

