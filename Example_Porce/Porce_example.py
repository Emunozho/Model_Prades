# -*- coding: utf-8 -*-
"""
This code is the python version of the R script by Sierra et al. (2021, JE)

"""

# Libraries
import multiprocessing
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import random
import scipy
from scipy.linalg import expm

# Definition of functions

def M(t,B,u): # Verify after use.***********************************
    return  np.matmul(expm(t*B),u)

def Rt(t,B,u):  # Verify after use.***********************************
    return(np.matmul(np.diagonal(np.sum(-B))),expm(t*B),u)    

def makeB(pars):
    B=np.diagonal(-1*pars[1:7])
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


# Calculate number of cores
no_cores=multiprocessing.cpu_count()

# Read data
path='Documents/MSCA/DISEQ/Prades model/Example_Porce/'
Cdata=pd.read_csv('PorcedB2012.csv')


poolnames=["Foliage", "Wood", "Fine roots", "Coarse roots", "Fine litter", "Coarse woody debris", "Soil carbon"]

# Proportion of components as reported by Zapata-Cuartas (thesis)
frac_foliage=0.08; frac_wood=0.19+0.73
GPP1=6.2711587*365/100 # Value of GPP from Ryu et al. and computed by Lina Estupiñan
GPP2=6.687468*365/100 # Value of GPP from Jung et al. and computed by Lina Estupiñan
GPP=np.mean([GPP1,GPP2])


# Random variates for primary forests. Mean and sd values from Sierra et al. (2007, FEM)
nplots=33
random.seed(11); primaryplotAges=np.random.choice(np.arange(100,150), nplots, replace=True)
random.seed(22); folp=np.random.normal(frac_foliage*111.6,frac_wood*18.5,nplots)
random.seed(33); woodp=np.random.normal(frac_wood*111.6,frac_wood*18.5,nplots)
random.seed(44); finerootp=np.random.normal(0.45*17.4,0.45*2.6,nplots)
random.seed(55); coarserootp=np.random.normal(0.45*67.1,0.45*16.4,nplots)
random.seed(66); finelitterp=np.random.normal(0.45*6.0,0.45*0.4,nplots)
random.seed(77); cwdp=np.random.normal(0.45*6.1,0.45*1.3,nplots)
random.seed(88); soilcp=np.random.normal(96.6,2.5,nplots)

# Observations from secondary forest chronosequence
obsFols=pd.DataFrame({'time':Cdata.Age,'Fol':0.45*(frac_foliage*(Cdata.AGB+Cdata.Palm)+Cdata.HV)})
obsWoods=pd.DataFrame({'time':Cdata.Age,'Wood':0.45*(frac_wood*(Cdata.AGB+Cdata.Palm))})
obsFRs=pd.DataFrame({'time':Cdata.Age,'FR':0.45*Cdata.FR})
obsCRs=pd.DataFrame({'time':Cdata.Age,'CR':0.45*Cdata.CR})
obsFLs=pd.DataFrame({'time':Cdata.Age,'FL':0.45*Cdata.FL})
obsCWDs=pd.DataFrame({'time':Cdata.Age,'CWD':0.45*Cdata.CWD})
obsS15s=pd.DataFrame({'time':Cdata.Age,'S15':Cdata.SC15})
obsS30s=pd.DataFrame({'time':Cdata.Age,'S30':Cdata.SC30})
obsSCs=pd.DataFrame({'time':Cdata.Age,'SC':Cdata.SC15+Cdata.SC30})

pd.concat([obsFols,pd.DataFrame({'time':primaryplotAges,'Fol':folp})])

# Merged observations
obsFol=pd.concat([obsFols,pd.DataFrame({'time':primaryplotAges,'Fol':folp})])
obsWood=pd.concat([obsWoods,pd.DataFrame({'time':primaryplotAges,'Wood':woodp})])
obsFR=pd.concat([obsFRs,pd.DataFrame({'time':primaryplotAges,'FR':finerootp})])
obsCR=pd.concat([obsCRs,pd.DataFrame({'time':primaryplotAges,'CR':coarserootp})])
obsFL=pd.concat([obsFLs,pd.DataFrame({'time':primaryplotAges,'FL':finelitterp})])
obsCWD=pd.concat([obsCWDs,pd.DataFrame({'time':primaryplotAges,'CWD':cwdp})])
obsSC=pd.concat([obsSCs,pd.DataFrame({'time':primaryplotAges,'SC':soilcp})])

yr=np.arange(0,150.5,0.5)
#P0=[0,0,0,0,0,0,np.nanmean(obsS15s.S15+obsS30s.S30)]
P0=pd.DataFrame({'Fol':[0],'Wood':[0],'CR':[0],'FL':[0],'CWD':[0],'SC':[np.nanmean(obsS15s.S15+obsS30s.S30)]})


#inipars=c(k1=0.5,k2=0.3,k3=0.1,k4=0.1,k5=0.1,k6=0.1,k7=0.1,
#           alpha21=0.1,alpha31=0.1,alpha41=0.1,alpha51=0.1,alpha53=0.1,alpha62=0.1,alpha64=0.1,alpha75=0.1,alpha76=0.1)

v0=np.array([0.5,0.3,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]).reshape((1,16))
inipars=pd.DataFrame(v0,columns=['k1','k2','k3','k4','k5','k6','k7',
                                 'alpha21','alpha31','alpha41','alpha51','alpha53','alpha62','alpha64','alpha75','alpha76'])


tau=np.arange(0,400.05,0.05)


####
n=1000
random.seed(123); gppRyu=np.random.normal(6.2711587,np.sqrt(1.6605167),n)
random.seed(456); gppJung=np.random.normal(6.687468,np.sqrt(0.28476033),n)

gpp=((gppRyu+gppJung)/2)*365/100
#

fig=plt.hist(gpp)
plt.xlabel(r'GPP [MgC ha$^{−1}$]'); plt.ylabel(r'$f_{GPP}$ [ha MgC$^{−1}$]' ) 
plt.savefig('GPPhist.pdf')


# Discover what is Model.
#ecoModel=function(pars, GPP){
#  mod=Model(t=yr,A=makeB(pars), ivList=P0, inputFluxes=c(GPP, rep(0,6)), pass=TRUE)
#  Ct=getC(mod)
#  colnames(Ct)<-c("Fol","Wood","FR","CR","FL","CWD","SC")
#  return(data.frame(time=yr, Ct))
#}








