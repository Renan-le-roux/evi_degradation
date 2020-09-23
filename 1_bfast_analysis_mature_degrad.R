library(raster)
library(rgdal)
library(raster)
library(raster)
library(matrixStats)
library(rgdal)
library(reshape2)
library(ggplot2)
library(xts)
library(rts)
library(tidyverse)
library(ggpmisc)
library(foreach)
library(parallel)
library(doSNOW)
library(matrixStats)
library(bfast)
library(greenbrown)
library(rasterVis)
library(gridExtra)



bfast.saeason<-function(i,modis.evi.para.mx){
  # T0<-Sys.time()
  x<-modis.evi.para.mx[i,]
  ts.x<-ts(x,start = c(2002,1),end=c(2018,12),frequency = 12)
  na.nb<-sum(is.na(ts.x))
  pour<-(length(ts.x)-length(ts.x)*0.6)
  if(na.nb>pour){
    seasonal.degrad<-(rep(NA,204))
  }else{
    ts.x<-na.approx(ts.x)
    bf<-bfast(ts.x,h=0.15, season="harmonic", max.iter=10)
    niter<- length(bf$output)
    seasonal.deg<-bf$output[[niter]]$St
    bp<-bf$output[[niter]]$Vt.bp
    bp<-as.Date(time(ts.x))[bp]
    trend<-bf$output[[niter]]$Tt
    seasonal.degrad<-list(seasonal.deg,bp,trend)
  }
  # T1<-Sys.time()
  
  
  return(seasonal.degrad)
}

########################################Prepa data evi and impact#############
setwd("H:/Utilisateurs/leroux_r/post_doc_cnes/CNES/Travail/FEU_EVI/Para/")
para<-readOGR(dsn = "../../raw_data/Para.shp")
MF<-raster("../../raw_data/Map_biomas/mask_foret_para.tif")
para<-spTransform(para,CRSobj = proj4string(MF))
MF<-crop(MF,para)
MF<-mask(MF,para)
MF_mature<-raster("../../raw_data/Map_biomas/MF_mature1.tif")
MF_mature<-crop(MF_mature,para)
MF_mature<-mask(MF_mature,para)

###############Resample Data
LS<-list.files("EVI_max_monthly_aggregate_1km/",full.names = T,pattern = ".tif$")
modis.evi.para<-stack(LS)

degrad_map<-raster("G:/Utilisateurs/Leroux/data_papier_degrad/class_degrad_para_focal.tif")
RAPPORT<-round(res(modis.evi.para)/res(degrad_map))[1]

degrad_map<-aggregate(x = degrad_map,fact=RAPPORT,fun=modal,na.rm=T)
degrad_map<-projectRaster(degrad_map,crs = proj4string(modis.evi.para),method = "ngb")

degrad_map<-resample(degrad_map,modis.evi.para,"ngb")

MF<-projectRaster(MF,crs = proj4string(modis.evi.para),method = "ngb")
MF<-resample(MF,degrad_map,"ngb")

MF_mature<-projectRaster(MF_mature,crs = proj4string(modis.evi.para),method = "ngb")
MF_mature<-resample(MF_mature,degrad_map,"ngb")

degrad_map<-degrad_map*MF

NOM<-names(modis.evi.para)
NOM<-sapply(NOM,function(x) strsplit(x,"MODIS_EVI_PARA_sum_monthly_")[[1]][2])
NOM<-paste0(NOM,".01")

DAT<-as.Date(NOM,format="%Y.%m.%d")


names(modis.evi.para)<-DAT
para<-readOGR("../../raw_data/para.shp")
para<-spTransform(para,CRSobj = proj4string(modis.evi.para))
modis.evi.para.mx<-as.matrix(modis.evi.para)
sup_outlier<-quantile(modis.evi.para.mx,na.rm=T,c(0.01,0.99))
pos.sup.outlier<-which(modis.evi.para.mx<sup_outlier[1]|modis.evi.para.mx>sup_outlier[2])
modis.evi.para.mx[pos.sup.outlier]<-NA

###############################Classif profil mature#############

MF_mature.vc<-as.vector(MF_mature)
pos.mature<-which(MF_mature.vc==1)
modis.evi.para.mx.mature<-modis.evi.para.mx[pos.mature,]


strt = Sys.time()
n.cores <- detectCores()-5
cl<-makeCluster(n.cores)
registerDoSNOW(cl)

#####################################################Parallèle



###############################################On fait la boucle que sur les données où on des valeurs de précipitations 
# Stats.Para_EVI_test = 
res.bast_parallele.mature<-foreach(i = 1:nrow(modis.evi.para.mx.mature), 
                            # .export = c("stats_Cell"),
                            # .combine = cbind,
                            .packages=c("raster","parallel","doSNOW","bfast","zoo"),
                            .inorder=T) %dopar% {
                              bfast.saeason(i,modis.evi.para.mx.mature)  
                            }


stopCluster(cl)
res.bast_parallele.mature[[1]]
save(res.bast_parallele.mature,file = "res.bast_parallele.mature.RData")
##############################Exrtract seasonal profile##

###############################Classif profil degrad#############

degrad_map.vc<-as.vector(degrad_map)
# writeRaster(degrad_map,"degradation_map.rsp.tif",datatype='INT2S',overwrite=T)
pos.degrad<-which(!is.na(degrad_map.vc))
modis.evi.para.mx.degrad<-modis.evi.para.mx[pos.degrad,]



strt = Sys.time()
n.cores <- detectCores()-5
cl<-makeCluster(n.cores)
registerDoSNOW(cl)

#####################################################Parallèle



###############################################On fait la boucle que sur les données où on des valeurs de précipitations 
# Stats.Para_EVI_test = 
res.bast_parallele.degrad<-foreach(i = 1:nrow(modis.evi.para.mx.degrad), 
                                   # .export = c("stats_Cell"),
                                   # .combine = cbind,
                                   .packages=c("raster","parallel","doSNOW","bfast","zoo"),
                                   .inorder=T) %dopar% {
                                     bfast.saeason(i,modis.evi.para.mx.degrad)  
                                   }


stopCluster(cl)
res.bast_parallele.degrad[[1]]
save(res.bast_parallele.degrad,file ="res.bast_parallele.degrad.RData" )


#####################################If faut faire l'analyse pour le nombre de classe####"

