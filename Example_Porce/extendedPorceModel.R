library(SoilR)
library(FME)
library(RColorBrewer)
library(expm)
library(parallel)

# Calculate the number of cores
no_cores <- detectCores() #- 10

# Initiate cluster
cl <- makeCluster(no_cores)
setwd("~/Documents/MSCA/DISEQ/Prades_model/Example_Porce") # set path.

Cdata<-read.csv("PorcedB2012.csv")
poolnames=c("Foliage", "Wood", "Fine roots", "Coarse roots", "Fine litter", "Coarse woody debris", "Soil carbon")

# Proportion of components as reported by Zapata-Cuartas (thesis)
frac_foliage=0.08; frac_wood=0.19+0.73
GPP1=6.2711587*365/100 # Value of GPP from Ryu et al. and computed by Lina Estupiñan
GPP2=6.687468*365/100 # Value of GPP from Jung et al. and computed by Lina Estupiñan
GPP=mean(c(GPP1,GPP2))

# Random variates for primary forests. Mean and sd values from Sierra et al. (2007, FEM)
nplots=33
set.seed(11); primaryplotAges=sample(100:150, nplots, replace=TRUE)
set.seed(22); folp=rnorm(nplots, mean=frac_foliage*111.6, sd=frac_foliage*18.5)
set.seed(33); woodp=rnorm(nplots, mean=frac_wood*111.6, sd=frac_wood*18.5)
set.seed(44); finerootp=rnorm(nplots, mean=0.45*17.4, sd=0.45*2.6)
set.seed(55); coarserootp=rnorm(nplots, mean=0.45*67.1, sd=0.45*16.4)
set.seed(66); finelitterp=rnorm(nplots, mean=0.45*6.0, sd=0.45*0.4)
set.seed(77); cwdp=rnorm(nplots, mean=0.45*6.1, sd=0.45*1.3)
set.seed(88); soilcp=rnorm(nplots, mean=96.6, sd=2.5)

# Observations from secondary forest chronosequence
obsFols=data.frame(time=Cdata$Age, Fol=0.45*(frac_foliage*(Cdata$AGB+Cdata$Palm)+Cdata$HV))
obsWoods=data.frame(time=Cdata$Age, Wood=0.45*(frac_wood*(Cdata$AGB+Cdata$Palm)))
obsFRs=data.frame(time=Cdata$Age,FR=0.45*Cdata$FR)
obsCRs=data.frame(time=Cdata$Age,CR=0.45*Cdata$CR)
obsFLs=data.frame(time=Cdata$Age,FL=0.45*Cdata$FL)
obsCWDs=data.frame(time=Cdata$Age,CWD=0.45*Cdata$CWD)
obsS15s=data.frame(time=Cdata$Age,S15=Cdata$SC15)
obsS30s=data.frame(time=Cdata$Age,S30=Cdata$SC30)
obsSCs=data.frame(time=Cdata$Age,SC=Cdata$SC15+Cdata$SC30)

# Merged observations
obsFol=rbind(obsFols, data.frame(time=primaryplotAges, Fol=folp))
obsWood=rbind(obsWoods, data.frame(time=primaryplotAges, Wood=woodp))
obsFR=rbind(obsFRs, data.frame(time=primaryplotAges, FR=finerootp))
obsCR=rbind(obsCRs, data.frame(time=primaryplotAges, CR=coarserootp))
obsFL=rbind(obsFLs, data.frame(time=primaryplotAges, FL=finelitterp))
obsCWD=rbind(obsCWDs, data.frame(time=primaryplotAges, CWD=cwdp))
obsSC=rbind(obsSCs, data.frame(time=primaryplotAges, SC=soilcp))

yr=seq(0,150,0.5)
P0=c(Fol=0,Wood=0, FR=0, CR=0, FL=0, CWD=0,SC=mean(obsS15s[,2]+obsS30s[,2],na.rm=TRUE))

 inipars=c(k1=0.5,k2=0.3,k3=0.1,k4=0.1,k5=0.1,k6=0.1,k7=0.1,
           alpha21=0.1,alpha31=0.1,alpha41=0.1,alpha51=0.1,alpha53=0.1,alpha62=0.1,alpha64=0.1,alpha75=0.1,alpha76=0.1)
#inipars=c(k1=0.5,k2=0.3,k3=0.1,k4=0.1,k5=0.1,k6=0.1,k7=0.1,
#          alpha21=0.1,alpha32=0.1,alpha42=0.1,alpha51=0.1,alpha53=0.1,alpha62=0.1,alpha64=0.1,alpha75=0.1,alpha76=0.1,K=1.0)

tau=seq(0,400, by=0.05)
M=function(t,B,u){expm(t*B)%*%u}
Rt=function(t, B, u){diag(colSums(-B))%*%expm(t*B)%*%u}

pal=brewer.pal(7,"Dark2")

####
n=1000
set.seed(123); gppRyu=rnorm(n,mean=6.2711587, sd=sqrt(1.6605167))
set.seed(456); gppJung=rnorm(n, mean=6.687468, sd=sqrt(0.28476033))
gpp=((gppRyu+gppJung)/2)*365/100

pdf("Figures/GPPhist.pdf")
hist(gpp)
dev.off()

#ecoModel=function(pars, GPP){
#  derivs=function(time,pools,pars){
#    with(as.list(c(pars,pools)),{
##      dFol=(GPP*(exp(Fol)/(K+exp(Fol))))-k1*Fol
##      dFol=(GPP*(Fol/(K+Fol)))-k1*Fol
##      dFol=(GPP^K)-k1*Fol
#      dFol=GPP-k1*Fol
#      dWood=alpha21*k1*Fol-k2*Wood
#      dFR=alpha31*k1*Fol-k3*FR
##      dFR=alpha32*k2*Wood-k3*FR
#      dCR=alpha41*k1*Fol-k4*CR
##      dCR=alpha42*(1-alpha32)*k2*Wood-k4*CR
#      dFL=alpha51*k1*Fol+alpha53*k3*FR-k5*FL
#      dCWD=alpha62*k2*Wood+alpha64*k4*CR-k6*CWD
#      dSC=alpha75*k5*FL+alpha76*k6*CWD-k7*SC
#
#      return(list(c(dFol,dWood,dFR,dCR,dFL,dCWD,dSC)))
#    })
#  }
#  out=ode(y=P0,parms=pars,times=yr,func=derivs)
#  as.data.frame(out)
#}

makeB=function(pars){
  B=diag(-1*pars[1:7])
  B[2,1]=pars["alpha21"]*pars["k1"]
  B[3,1]=pars["alpha31"]*pars["k1"]
  B[4,1]=pars["alpha41"]*pars["k1"]
  B[5,1]=pars["alpha51"]*pars["k1"]
  B[5,3]=pars["alpha53"]*pars["k3"]
  B[6,2]=pars["alpha62"]*pars["k2"]
  B[6,4]=pars["alpha64"]*pars["k4"]
  B[7,5]=pars["alpha75"]*pars["k5"]
  B[7,6]=pars["alpha76"]*pars["k6"]
  return(B)
}

#===============================================================================


ecoModel=function(pars, GPP){
  mod=Model(t=yr,A=makeB(pars), ivList=P0, inputFluxes=c(GPP, rep(0,6)), pass=TRUE)
  Ct=getC(mod)
  colnames(Ct)<-c("Fol","Wood","FR","CR","FL","CWD","SC")
  return(data.frame(time=yr, Ct))
}

ecoCost=function(pars, GPP){
  out=ecoModel(pars, GPP)
  cost1=modCost(model=out,obs=obsFol, weight="none")
  cost2=modCost(model=out,obs=obsWood,cost=cost1, weight="none")
  cost3=modCost(model=out,obs=obsFR,cost=cost2, weight="none")
  cost4=modCost(model=out,obs=obsCR,cost=cost3, weight="none")
  cost5=modCost(model=out,obs=obsFL,cost=cost4, weight="std")
  cost6=modCost(model=out,obs=obsCWD,cost=cost5, weight="mean")
  cost7=modCost(model=out,obs=obsSC,cost=cost6, weight="none")
  return(cost7)
}

randomFit=function(GPP){
  ecoFit=modFit(f=ecoCost,p=inipars,GPP=GPP,lower=0,upper=c(rep(3,7),rep(1,9)), method="Marq")
  cp=coef(ecoFit)
  return(cp)
}

meanSATTs=function(u,B){
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

# modpars=parSapply(cl,gpp,FUN=randomFit)
modpars=t(modpars)
save(modpars, file="modpars.RData")
load("modpars.RData")

porceSATT=NULL
for(i in 1:length(gpp)){
  tmp=(makeMod(pars=modpars[i,], GPP=gpp[i]))
  porceSATT=rbind(porceSATT,tmp)
}


boxplot(porceSATT[,2:3], ylab="Age or transit time (years)", names=c("Mean system age", "Mean transit time"), outline=FALSE, ylim=c(0,60))
boxplot(modpars, outline=FALSE)
round(apply(porceSATT,2,mean),2);round(apply(porceSATT,2,sd),2) 

meanpars=apply(modpars, MARGIN = 2, FUN=mean)
sdpars=apply(modpars, MARGIN = 2, FUN=sd)
meanGPP=mean(gpp); meanU=matrix(c(meanGPP, rep(0,6)),7,1)
meanB=makeB(meanpars)
meanX=-1*solve(meanB)%*%meanU
meanFluxes=meanB%*%diag(as.numeric(meanX))
colSums(meanB%*%diag(as.numeric(1/diag(meanB)))) # First element is the proportion respired from foliage quoted in text

mSA=systemAge(A=meanB, u=meanU,a=tau)
mTT=transitTime(A=meanB, u=meanU, a=tau)
Mt=t(sapply(tau,FUN=M, B=meanB, u=meanU))
Resp=t(sapply(tau, FUN=Rt, B=meanB, u=meanU))

pdf("Figures/Fate.pdf")
par(mfrow=c(2,1),mar=c(4,4.5,0.1,0.1))
matplot(tau, Mt[,1:3], col=pal, lty=c(1,2,2), lwd=c(1,1,2), type="l", xlim=c(0,100),cex.lab=1.1,  xlab="", ylab=expression(paste("Carbon remaining (Mg C h",a^-1,")")), bty="n")
lines(tau, rowSums(Mt),lwd=2)
legend("topright", c("Total",poolnames[1:3]), col=c(1,pal[1:3]), lty=c(1,1,2,2),lwd=c(2,1,1,2), bty="n")
matplot(tau, Mt[,4:7], col=pal[4:7], lty=c(1,2,1,2),lwd=c(1,2,2,1), type="l", xlim=c(0,100),cex.lab=1.1,  xlab="Time since fixation (years)", ylab=expression(paste("Carbon remaining (Mg C h",a^-1,")")), bty="n")
legend("topright", poolnames[4:7], col=pal[4:7], lty=c(1,2,1,2),lwd=c(1,2,2,1), bty="n")
par(mfrow=c(1,1))
dev.off()

pdf("Figures/TransitTime.pdf")
par(mfrow=c(2,1), mar=c(4,4.5,1,1))
plot(tau,meanGPP*mTT$transitTimeDensity, xlim=c(0,10),type="l",lwd=2,cex.lab=1.1,xlab="", ylab=expression(paste("Carbon respired (Mg C h",a^-1, " y", r^-1, ")")),bty="n")
matlines(tau,Resp[,c(1,5)],lty=c(1,2),col=pal[c(1,5)])
abline(v=mTT$quantiles[2],col=2,lty=3)
legend("right", c("Total",poolnames[c(1,5)]), col=c(1,pal[c(1,5)]), lty=c(1,1,2),lwd=c(2,1,1), bty="n")
legend("topright", "a", bty="n")

par(mar=c(4,4.5,0,1))
matplot(tau,Resp[,c(-1,-5)], col=pal[c(-1,-5)], lty=c(2,2,1,2,1),lwd=c(1,2,1,2,1), type="l",xlim=c(0,20),ylim=c(0,0.1),xlab="Time since fixation (years)", ylab=expression(paste("Carbon respired (Mg C h",a^-1, " y", r^-1, ")")), bty="n")
legend("topright", c(poolnames[c(-1,-5)]), col=c(pal[c(-1,-5)]), lty=c(2,2,1,2,1),lwd=c(1,2,1,2,1), bty="n")
legend("topright", "b", bty="n")
par(mfrow=c(1,1))
dev.off()

ind50=which(tau <= mTT$quantiles[2])
ind95=which(tau <= mTT$quantiles[3])
bluepal=brewer.pal(3,"Blues")

pdf("Figures/areaTTdistribution.pdf")
par(mar=c(4,5,1,1))
plot(tau,meanGPP*mTT$transitTimeDensity, xlim=c(0,5),type="l",lwd=2,cex.lab=1.1,
     ylab=expression(paste("Carbon respired (Mg C h",a^-1, " y", r^-1, ")")),
     xlab="Time since C fixation (transit time in years)",bty="n")
polygon(x=c(tau[ind95], rev(tau[ind95])), y=c(meanGPP*mTT$transitTimeDensity[ind95], rep(0,length(ind95))),col=bluepal[2])
polygon(x=c(tau[ind50], rev(tau[ind50])), y=c(meanGPP*mTT$transitTimeDensity[ind50], rep(0,length(ind50))),col=bluepal[3])
legend("topright", c(paste("50% of GPP respired in less than ", round(mTT$quantiles[2],2), " years"),
                     paste("95% of GPP respired in less than ", round(mTT$quantiles[3],2), " years")),
       col=bluepal[c(3,2)], pch=15, bty="n")
dev.off()

plot(tau, mTT$transitTimeDensity, type="l", log="y",col=2, xlim=c(0,150), xlab="Age or transit time (years)", ylab="Density distribution", bty="n")
abline(v=mTT$meanTransitTime, col=2, lty=2)
lines(tau, mSA$systemAgeDensity, col=4)
abline(v=mSA$meanSystemAge, col=4, lty=2)
legend("topright", c("System age", "Transit time"), lty=1, col=c(4,2), bty="n")

##
avgModel=ecoModel(pars=meanpars, GPP=meanGPP)
allModels=array(0,c(length(yr), 8, n))
for(i in 1:n){
  allModels[,,i]<-as.matrix(ecoModel(pars=modpars[i,], GPP=gpp[i]))
}

meanpred=apply(allModels, c(1,2), mean)
sdpred=apply(allModels, c(1,2), sd)
minpred=apply(allModels, c(1,2), min)
maxpred=apply(allModels, c(1,2), max)

allX=NULL
for(i in 1:n){
  allX[i]=sum(-1*solve(makeB(modpars[i,]))%*%matrix(c(gpp[i],rep(0,6)),7,1))
}
round(mean(allX),1); round(sd(allX),1)

allR=NULL
allRa=NULL
allRh=NULL
for(i in 1:n){
  tmpB=makeB(modpars[i,])
  tmpX=-1*solve(makeB(modpars[i,]))%*%matrix(c(gpp[i],rep(0,6)),7,1)
  tmpr=colSums(tmpB)
  allR[i]=sum(tmpr*tmpX)
  allRa[i]=sum(tmpr[1:4]*tmpX[1:4])
  allRh[i]=sum(tmpr[5:7]*tmpX[5:7])
}
round(mean(allRa),1); round(sd(allRa),1)
round(mean(allRh),1); round(sd(allRh),1)
round(mean(allR),1); round(sd(allR),1)
mean(allRa)/mean(allR); mean(allRh)/mean(allR)

NPP=gpp+allRa
round(mean(NPP),1); round(sd(NPP),1)

meanPoolAgef=function(pars,GPP){
  B=makeB(pars)
  u=matrix(c(GPP,rep(0,6)),7,1)
  SA=systemAge(A=B,u=u,a=1,q=0.95)
  return(SA)
}

mPA=NULL
qSA=NULL
for(i in 1:n){
  tmp=meanPoolAgef(pars=modpars[i,],GPP=gpp[i])
  mPA=cbind(mPA,tmp$meanPoolAge)
  qSA=cbind(qSA,tmp$quantilesSystemAge)
}
round(apply(mPA,1,mean),2); round(apply(mPA,1,sd),2)
round(apply(qSA,1,mean),1); round(apply(qSA,1,sd),1)

TTf=function(pars,GPP){
  B=makeB(pars)
  u=matrix(c(GPP,rep(0,6)),7,1)
  TT=transitTime(A=B,u=u,a=1,q=c(0.5,0.95))
  return(TT)
}

q5TT=NULL
q95TT=NULL
for(i in 1:n){
  tmp=TTf(pars=modpars[i,],GPP=gpp[i])
  q5TT[i]=tmp$quantiles[1]
  q95TT[i]=tmp$quantiles[2]
}
round(mean(q5TT),2); round(sd(q5TT),2)
round(mean(q95TT),2); round(sd(q95TT),2)


pdf("Figures/avgModel.pdf")
matplot(avgModel[,1], avgModel[,-1], type="l",lty=1, col=pal, xlab="Time since recovery (years)", ylab="Carbon stocks (Mg C ha-1)", bty="n")
points(obsFol,col=pal[1], pch=19, cex=0.5)
points(obsWood, col=pal[2], pch=19, cex=0.5)
points(obsFR, col=pal[3], pch=19, cex=0.5)
points(obsCR, col=pal[4], pch=19, cex=0.5)
points(obsFL, col=pal[5], pch=19, cex=0.5)
points(obsCWD, col=pal[6], pch=19, cex=0.5)
points(obsSC, col=pal[7], pch=19, cex=0.5)
legend("topleft", legend=poolnames, lty=1, col=pal, bty="n")
dev.off()

pdf("Figures/dataModelFit.pdf")
par(mfrow=c(2,2), mar=c(4,4.5,0.1,0.1))
plot(meanpred[,1], meanpred[,2], type="l", col=pal[1], bty="n", ylim=c(0,15),cex.lab=1.1, ylab=expression(paste("Carbon stock (Mg C h",a^-1,")")), xlab="")
polygon(x=c(maxpred[,1],rev(minpred[,1])), y=c(maxpred[,2], rev(minpred[,2])), col="gray", border = pal[1])
points(obsFol,col=pal[1], pch=19, cex=0.5)
polygon(x=c(maxpred[,1],rev(minpred[,1])), y=c(maxpred[,6], rev(minpred[,6])), col="gray", border = pal[5])
points(obsFL, col=pal[5], pch=19, cex=0.5)
legend("topleft", "a", bty="n")

plot(meanpred[,1], meanpred[,3], type="l", col=pal[2], bty="n", ylim=c(0,150),cex.lab=1.1, ylab=expression(paste("Carbon stock (Mg C h",a^-1,")")), xlab="")
polygon(x=c(maxpred[,1],rev(minpred[,1])), y=c(maxpred[,3], rev(minpred[,3])), col="gray", border = pal[2])
points(obsWood, col=pal[2], pch=19, cex=0.5)
polygon(x=c(maxpred[,1],rev(minpred[,1])), y=c(maxpred[,5], rev(minpred[,5])), col="gray", border = pal[4])
points(obsCR, col=pal[4], pch=19, cex=0.5)
legend("topleft", "b", bty="n")

plot(meanpred[,1], meanpred[,7], col=pal[6], type="l", ylim=c(0,10),cex.lab=1.1, bty="n", ylab=expression(paste("Carbon stock (Mg C h",a^-1,")")), xlab="Time (yr)")
polygon(x=c(maxpred[,1],rev(minpred[,1])), y=c(maxpred[,7], rev(minpred[,7])), col="gray", border = pal[6])
points(obsCWD, col=pal[6], pch=19, cex=0.5)
polygon(x=c(maxpred[,1],rev(minpred[,1])), y=c(maxpred[,4], rev(minpred[,4])), col="gray", border = pal[3])
points(obsFR, col=pal[3], pch=19, cex=0.5)
legend("topleft", "c", bty="n")

plot(meanpred[,1], meanpred[,8], type="l", col=pal[7], bty="n", ylim=c(0,200),cex.lab=1.1, ylab=expression(paste("Carbon stock (Mg C h",a^-1,")")), xlab="Time (yr)")
polygon(x=c(maxpred[,1],rev(minpred[,1])), y=c(maxpred[,8], rev(minpred[,8])), col="gray", border = pal[7])
points(obsSC, col=pal[7], pch=19, cex=0.5)
legend("topright", legend=poolnames, lty=1, col=pal, bty="n")
legend("topleft", "d", bty="n")
par(mfrow=c(1,1))
dev.off()

pdf("Figures/modelDataFitboxplot.pdf")
par(mfrow=c(2,2), mar=c(4,4.5,0.1,0.1))
plot(meanpred[,1], meanpred[,2], type="l", col=pal[1], bty="n", ylim=c(0,15),xlim=c(0,100), ylab=expression(paste("Carbon stock (Mg C h",a^-1,")")), xlab="")
polygon(x=c(maxpred[,1],rev(minpred[,1])), y=c(maxpred[,2], rev(minpred[,2])), col="gray", border = pal[1])
points(obsFols,col=pal[1], pch=19, cex=0.5)
polygon(x=c(maxpred[,1],rev(minpred[,1])), y=c(maxpred[,6], rev(minpred[,6])), col="gray", border = pal[5])
points(obsFLs, col=pal[5], pch=19, cex=0.5)
boxplot(folp,add=TRUE, at=c(100),range=0,bty="n",col=pal[1])
boxplot(finelitterp,add=TRUE, at=c(100),range=0,bty="n",col=pal[5])
legend("topleft", "a", bty="n")

plot(meanpred[,1], meanpred[,3], type="l", col=pal[2], bty="n", ylim=c(0,150),xlim=c(0,100), ylab=expression(paste("Carbon stock (Mg C h",a^-1,")")), xlab="")
polygon(x=c(maxpred[,1],rev(minpred[,1])), y=c(maxpred[,3], rev(minpred[,3])), col="gray", border = pal[2])
points(obsWoods, col=pal[2], pch=19, cex=0.5)
polygon(x=c(maxpred[,1],rev(minpred[,1])), y=c(maxpred[,5], rev(minpred[,5])), col="gray", border = pal[4])
points(obsCRs, col=pal[4], pch=19, cex=0.5)
boxplot(woodp,add=TRUE, at=c(100),range=0,bty="n",col=pal[1])
boxplot(coarserootp,add=TRUE, at=c(100),range=0,bty="n",col=pal[5])
legend("topleft", "b", bty="n")

plot(meanpred[,1], meanpred[,7], col=pal[6], type="l", ylim=c(0,10),xlim=c(0,100), bty="n", ylab=expression(paste("Carbon stock (Mg C h",a^-1,")")), xlab="Time (yr)")
polygon(x=c(maxpred[,1],rev(minpred[,1])), y=c(maxpred[,7], rev(minpred[,7])), col="gray", border = pal[6])
points(obsCWDs, col=pal[6], pch=19, cex=0.5)
polygon(x=c(maxpred[,1],rev(minpred[,1])), y=c(maxpred[,4], rev(minpred[,4])), col="gray", border = pal[3])
points(obsFRs, col=pal[3], pch=19, cex=0.5)
boxplot(cwdp,add=TRUE, at=c(100),range=0,bty="n",col=pal[1])
boxplot(finerootp,add=TRUE, at=c(100),range=0,bty="n",col=pal[5])
legend("topleft", "c", bty="n")

plot(meanpred[,1], meanpred[,8], type="l", col=pal[7], bty="n", ylim=c(0,200),xlim=c(0,100), ylab=expression(paste("Carbon stock (Mg C h",a^-1,")")), xlab="Time (yr)")
polygon(x=c(maxpred[,1],rev(minpred[,1])), y=c(maxpred[,8], rev(minpred[,8])), col="gray", border = pal[7])
points(obsSCs, col=pal[7], pch=19, cex=0.5)
boxplot(soilcp,add=TRUE, at=c(100),range=0,bty="n",col=pal[1])
legend("topright", legend=poolnames, lty=1, col=pal, bty="n")
legend("topleft", "d", bty="n")
par(mfrow=c(1,1))
dev.off()

upxlim=c(10,100,100,100,10,100,200,200); agelab=c(rep(" ",6), "Age (yr)")
pdf("Figures/poolAges.pdf")
par(mfrow=c(4,2),mar=c(4,4,1,1))
for(i in 1:length(poolnames)){
  plot(tau,mSA$poolAgeDensity[,i],type="l",col=pal[i], xlim=c(0,upxlim[i]),cex.lab=1.1, ylab="Density", xlab=agelab[i], bty="n")
  legend("top",poolnames[i],bty="n", text.font = 3)
}
plot(tau,mSA$systemAgeDensity, type="l", xlim=c(0,200), xlab="Age (yr)", ylab="Density", bty="n")
legend("top", "Ecosystem", text.font=3, bty="n")
par(mfrow=c(1,1))
dev.off()

library(xtable)
parsymbols=c(paste("$k_",1:7,"$", sep=""),paste("$\\alpha_{",c("2,1","3,1","4,1","5,1","5,3","6,2","6,4","7,5","7,6"),"}$",sep=""))
partext=c("Cycling rate in foliage", #1
          "Cycling rate in wood", #2
          "Cycling rate in fine roots", #3
          "Cycling rate in coarse roots", #4
          "Cycling rate in fine litter", #5
          "Cycling rate in coarse woody debris", #6
          "Cycling rate in soil carbon", #7
          "Proportion transferred from foliage to wood",
          "Proportion transferred from foliage to fine roots", 
          "Proportion transferred from foliage to coarse roots",
          "Proportion transferred from foliage to fine litter",
          "Proportion transferred from fine roots to fine litter",
          "Proportion transferred from wood to coarse woody debris",
          "Proportion transferred from coarse roots to coarse woody debris",
          "Proportion transferred from fine litter to soil carbon",
          "Proportion transferred from coarse woody debris to soil carbon")
parstats=data.frame(Parameter=parsymbols, Description=partext,Mean=meanpars,SD=sdpars)
tabcap=paste("Mean and standard deviation (SD) of parameter values obtained from the", n, "iterations of the optimization procedure.
             Values of cycling rates are given in units of yr$^{-1}$, and values of transfer coefficients are unitless (proportion between
             0 and 1).")
print(xtable(parstats,align=rep("l",5), digits=3, caption=tabcap, label="tab:parstats"), file="../parstatsR1.tex",caption.placement="top",
      include.rownames=FALSE, include.colnames=TRUE, floating=TRUE, booktabs=TRUE,sanitize.text.function=function(x){x})

NPPstar=mean(c(12.59, 12.93))
sdNPPstar=sqrt(sum(c(0.9,0.96)))

stopCluster(cl)
