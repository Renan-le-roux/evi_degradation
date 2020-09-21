library(raster)
library(rgdal)
###########################Fire data papier degrad
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
########################################Prepa data evi and impact#############
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
plot(modis.evi.para)

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




plot(degrad_map,col=rainbow(9))

NOM<-names(modis.evi.para)
NOM<-sapply(NOM,function(x) strsplit(x,"MODIS_EVI_PARA_sum_monthly_")[[1]][2])
NOM<-paste0(NOM,".01")

DAT<-as.Date(NOM,format="%Y.%m.%d")


names(modis.evi.para)<-DAT
para<-readOGR("../../raw_data/para.shp")
para<-spTransform(para,CRSobj = proj4string(modis.evi.para))
modis.evi.para.mx<-as.matrix(modis.evi.para)
####################EVI FIRE ##################################
degrad_map.vc<-as.vector(degrad_map)
pos.fire<-which(degrad_map.vc==2)

modis.evi.para.mx.fire<-modis.evi.para.mx[pos.fire,]

bfast.saeason<-function(i,modis.evi.para.fire){
  # T0<-Sys.time()
  x<-modis.evi.para.fire[i,]
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
    
    
    seasonal.degrad<-list(seasonal.deg,bp)
  }
  # T1<-Sys.time()
  
  
  return(seasonal.degrad)
}



strt = Sys.time()
n.cores <- detectCores()-5
cl<-makeCluster(n.cores)
registerDoSNOW(cl)

#####################################################Parallèle



###############################################On fait la boucle que sur les données où on des valeurs de précipitations 
# Stats.Para_EVI_test = 
res.bast_parallele<-foreach(i = 1:nrow(modis.evi.para.mx.fire), 
                            # .export = c("stats_Cell"),
                            # .combine = cbind,
                            .packages=c("raster","parallel","doSNOW","bfast","zoo"),
                            .inorder=T) %dopar% {
                              bfast.saeason(i,modis.evi.para.mx.fire)  
                            }


stopCluster(cl)
res.bast_parallele[[1]]


res.bfast.xts<-list()

for (i in 1:length(res.bast_parallele)){
  print(i)
  X<-xts(res.bast_parallele[[i]][[1]], as.Date(time(res.bast_parallele[[i]][[1]])))
  
  
  # colnames(X)<-c("evi.seas","date")
  res.bfast.xts[[i]]<-X
  # res.date<-cbind(res.date,X$date)
}
res.bfast.xts<-do.call("cbind",res.bfast.xts)
res.date<-matrix(nrow = length(res.bast_parallele),ncol = 5)
for (i in 1:length(res.bast_parallele)){
  print(i)
  dat.X<-res.bast_parallele[[i]][[2]]
  dat.X<-as.character(dat.X)
  
  
  if(length(dat.X)==0){dat.X=c(NA,NA,NA,NA,NA)}
  if(length(dat.X)==1){dat.X=c(dat.X,NA,NA,NA,NA)}
  if(length(dat.X)==2){dat.X=c(dat.X,NA,NA,NA)}
  if(length(dat.X)==3){dat.X=c(dat.X,NA,NA)}
  if(length(dat.X)==4){dat.X=c(dat.X,NA)}
  if(length(dat.X)==5){dat.X=c(dat.X)}
  # colnames(X)<-c("evi.seas","date")
  res.date[i,]<-as.Date(dat.X)
  # res.date<-cbind(res.date,X$date)
}
hist(as.Date(as.vector(res.date)))
hist(as.Date(res.date[,2]),18)
date.vc<-NULL
for (i in 1 :ncol(res.date)){
  date.vc<-c(date.vc,res.date[,i])
}

date.vc<-as.data.frame(as.Date(date.vc))

date.vc$Year<-format(date.vc,"%Y")
plot(date.vc$`as.Date(date.vc)`,date.vc$Year)
pos.na<-which(is.na(date.vc[,1]))
date.vc<-date.vc[-pos.na,]
summ<-table(date.vc$Year)
plot(summ)
plot(as.Date(date.vc),)
str(res.date)
str(res.bfast.xts)
time(res.bfast.xts)
dim(res.bfast.xts)
subset_xts<-res.bfast.xts['2002/2018']
subset_xts.mx<-as.matrix(t(subset_xts))
dim(subset_xts.mx)
colnames(subset_xts.mx)
rownames(subset_xts.mx)<-1:nrow(subset_xts.mx)


kmeans.seasonal<-kmeans(na.omit(subset_xts.mx),5)
res <- kmeans(na.omit(subset_xts.mx), 2)$cluster
res
DF=cbind.data.frame(subset_xts.mx, 'clus'=NA)
DF[names(res),][,ncol(DF)] <- res
colnames(DF)[1:ncol(subset_xts.mx)]<-colnames(subset_xts.mx)
colnames(DF)[ncol(DF)]<-"class"
DF$class<-as.factor(DF$class)
DF.agg<-aggregate(.~class,DF,function(x) mean(x,na.rm=T))
table(DF$class)
DF.agg.melt<-melt(DF.agg)
head(DF.agg.melt)
colnames(DF.agg.melt)<-c("class","date","evi")
DF.agg.melt$month<-format(as.Date(DF.agg.melt$date),"%m")

DF.agg.melt.agg<-aggregate(evi~month+class,DF.agg.melt,mean)
gplot_seasggplot<-ggplot(DF.agg.melt.agg,aes(x = as.integer(month),y=evi,colour=class))+
  geom_line()
gplot_seasggplot
ggsave(gplot_seasggplot,filename = "")




vect<-as.vector(modis.evi.para.mx[,1])
vect[]<-NA
vect[pos.fire]<-DF$class





vect.clas.degrad<-setValues(modis.evi.para[[1]],vect)
plot(vect.clas.degrad,col=rainbow(2))

plot(vect.clas.degrad,col=rainbow(5))
writeRaster(vect.clas.degrad.agg,filename = "bfast_kmeans_5_agg_int.tif",datatype='INT2S')

#########################logging###############################

degrad_map.vc<-as.vector(degrad_map)
pos.logging<-which(degrad_map.vc==1)

modis.evi.para.mx.logging<-modis.evi.para.mx[pos.logging,]

bfast.saeason<-function(i,modis.evi.para.logging){
  # T0<-Sys.time()
  x<-modis.evi.para.logging[i,]
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
    
    
    seasonal.degrad<-list(seasonal.deg,bp)
  }
  # T1<-Sys.time()
  
  
  return(seasonal.degrad)
}



strt = Sys.time()
n.cores <- detectCores()-5
cl<-makeCluster(n.cores)
registerDoSNOW(cl)

#####################################################Parallèle



###############################################On fait la boucle que sur les données où on des valeurs de précipitations 
# Stats.Para_EVI_test = 
res.bast_parallele.logging<-foreach(i = 1:nrow(modis.evi.para.mx.logging), 
                            # .export = c("stats_Cell"),
                            # .combine = cbind,
                            .packages=c("raster","parallel","doSNOW","bfast","zoo"),
                            .inorder=T) %dopar% {
                              bfast.saeason(i,modis.evi.para.mx.logging)  
                            }


stopCluster(cl)
res.bast_parallele.logging[[1]]


res.bfast.xts<-list()

for (i in 1:length(res.bast_parallele)){
  print(i)
  X<-xts(res.bast_parallele[[i]][[1]], as.Date(time(res.bast_parallele[[i]][[1]])))
  
  
  # colnames(X)<-c("evi.seas","date")
  res.bfast.xts[[i]]<-X
  # res.date<-cbind(res.date,X$date)
}
res.bfast.xts<-do.call("cbind",res.bfast.xts)
res.date<-matrix(nrow = length(res.bast_parallele),ncol = 5)
for (i in 1:length(res.bast_parallele)){
  print(i)
  dat.X<-res.bast_parallele[[i]][[2]]
  dat.X<-as.character(dat.X)
  
  
  if(length(dat.X)==0){dat.X=c(NA,NA,NA,NA,NA)}
  if(length(dat.X)==1){dat.X=c(dat.X,NA,NA,NA,NA)}
  if(length(dat.X)==2){dat.X=c(dat.X,NA,NA,NA)}
  if(length(dat.X)==3){dat.X=c(dat.X,NA,NA)}
  if(length(dat.X)==4){dat.X=c(dat.X,NA)}
  if(length(dat.X)==5){dat.X=c(dat.X)}
  # colnames(X)<-c("evi.seas","date")
  res.date[i,]<-as.Date(dat.X)
  # res.date<-cbind(res.date,X$date)
}
hist(as.Date(as.vector(res.date)))
hist(as.Date(res.date[,2]),18)
date.vc<-NULL
for (i in 1 :ncol(res.date)){
  date.vc<-c(date.vc,res.date[,i])
}

date.vc<-as.data.frame(as.Date(date.vc))

date.vc$Year<-format(date.vc,"%Y")
plot(date.vc$`as.Date(date.vc)`,date.vc$Year)
pos.na<-which(is.na(date.vc[,1]))
date.vc<-date.vc[-pos.na,]
summ<-table(date.vc$Year)
plot(summ)
plot(as.Date(date.vc),)
str(res.date)
str(res.bfast.xts)
time(res.bfast.xts)
dim(res.bfast.xts)
subset_xts<-res.bfast.xts['2002/2018']
subset_xts.mx<-as.matrix(t(subset_xts))
dim(subset_xts.mx)
colnames(subset_xts.mx)
rownames(subset_xts.mx)<-1:nrow(subset_xts.mx)


kmeans.seasonal<-kmeans(na.omit(subset_xts.mx),5)
res <- kmeans(na.omit(subset_xts.mx), 2)$cluster
res
DF=cbind.data.frame(subset_xts.mx, 'clus'=NA)
DF[names(res),][,ncol(DF)] <- res
colnames(DF)[1:ncol(subset_xts.mx)]<-colnames(subset_xts.mx)
colnames(DF)[ncol(DF)]<-"class"
DF$class<-as.factor(DF$class)
DF.agg<-aggregate(.~class,DF,function(x) mean(x,na.rm=T))
table(DF$class)
DF.agg.melt<-melt(DF.agg)
head(DF.agg.melt)
colnames(DF.agg.melt)<-c("class","date","evi")
DF.agg.melt$month<-format(as.Date(DF.agg.melt$date),"%m")

DF.agg.melt.agg<-aggregate(evi~month+class,DF.agg.melt,mean)
gplot_seasggplot<-ggplot(DF.agg.melt.agg,aes(x = as.integer(month),y=evi,colour=class))+
  geom_line()
gplot_seasggplot
ggsave(gplot_seasggplot,filename = "")




vect<-as.vector(modis.evi.para.mx[,1])
vect[]<-NA
vect[pos.logging]<-DF$class





vect.clas.degrad<-setValues(modis.evi.para[[1]],vect)
plot(vect.clas.degrad,col=rainbow(2))

plot(vect.clas.degrad,col=rainbow(5))
writeRaster(vect.clas.degrad.agg,filename = "bfast_kmeans_5_agg_int.tif",datatype='INT2S')



