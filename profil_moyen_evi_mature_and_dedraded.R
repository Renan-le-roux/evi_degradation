library(raster)
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
library(rasterVis)
library(gridExtra)
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

########################Profil moyen forêt mature##########
degrad_map.vc<-as.vector(degrad_map)
MF_mature.vc<-as.vector(MF_mature)
pos.mature<-which(MF_mature.vc==1)


length(pos.mature)
modis.evi.para.mx.mature<-modis.evi.para.mx[pos.mature,]


date.modis<-as.Date(colnames(modis.evi.para.mx.mature),format="X%Y.%m.%d")

month<-sprintf("%02d",1:12)
month<-paste0("-",month,"-")
monthly.evi.mature<-NULL
for (i in 1:length(month)){
  print(i)
  pos.month<-grep(month[i],date.modis)
  modis.month<-modis.evi.para.mx.mature[,pos.month]
  modis.mean.month<-apply(modis.month,1,mean,na.rm=T)
  monthly.evi.mature<-cbind(monthly.evi.mature,modis.mean.month)
}

monthly.evi.mature<-as.data.frame(monthly.evi.mature)
colnames(monthly.evi.mature)<-sprintf('%02d',1:12)
res <- kmeans(na.omit(monthly.evi.mature), 5,iter.max = 100,nstart = 20)$cluster
DF=cbind.data.frame(monthly.evi.mature, 'clus'=NA)
DF[names(res),][,ncol(DF)] <- res
colnames(DF)[1:ncol(monthly.evi.mature)]<-colnames(monthly.evi.mature)
colnames(DF)[ncol(DF)]<-"class"
DF$class<-as.factor(DF$class)
DF.agg<-aggregate(.~class,DF,function(x) mean(x,na.rm=T))
table(DF$class)
DF.agg.melt<-melt(DF.agg)
head(DF.agg.melt)
colnames(DF.agg.melt)<-c("class","date","evi")
# DF.agg.melt$month<-format(as.Date(DF.agg.melt$date),"%m")

DF.agg.melt.agg<-aggregate(evi~date+class,DF.agg.melt,mean)
gplot_seasggplot.mature<-ggplot(DF.agg.melt.agg,aes(x = as.integer(date),y=evi,colour=class))+
  geom_line()+
  scale_color_identity(guide = "legend")
gplot_seasggplot.mature



vect<-as.vector(modis.evi.para.mx[,1])
vect[]<-NA
vect[pos.mature]<-DF$class


vect.clas.mature<-setValues(modis.evi.para[[1]],vect)
plot(vect.clas.mature,col=rainbow(10))

plot(vect.clas.mature,col=rainbow(5))
writeRaster(vect.clas.mature,filename = "results_raster/kmeans_foret_mature_moyen.tif",datatype='INT2S')


vect.clas.mature.gg<-rasterToPoints(vect.clas.mature)
vect.clas.mature.gg<-as.data.frame(vect.clas.mature.gg)
head(vect.clas.mature.gg)
colnames(vect.clas.mature.gg)<-c("long","lat","class")

carto_foret_mature<-ggplot()+
  geom_raster(data = vect.clas.mature.gg,aes(x=long,y=lat,fill=as.factor(class)))+
  geom_polygon(data=para, aes(x=long, y=lat, group=group), 
               fill=NA,color="grey50", size=1)+
  coord_equal()+
  scale_fill_identity(guide = "none")

carto_foret_mature

TT<-grid.arrange(gplot_seasggplot.mature,carto_foret_mature,ncol=2)
ggsave(TT,filename = "figure/bfast_season_result/mature_forest_5class_profil.jpeg") 



########################Profil moyen forêt degrad##########
degrad_map.vc<-as.vector(degrad_map)
pos.degrad<-which(degrad_map.vc>0)


length(pos.degrad)
modis.evi.para.mx.degrad<-modis.evi.para.mx[pos.degrad,]


date.modis<-as.Date(colnames(modis.evi.para.mx.degrad),format="X%Y.%m.%d")

month<-sprintf("%02d",1:12)
month<-paste0("-",month,"-")
monthly.evi.degrad<-NULL
for (i in 1:length(month)){
  print(i)
  pos.month<-grep(month[i],date.modis)
  modis.month<-modis.evi.para.mx.degrad[,pos.month]
  modis.mean.month<-apply(modis.month,1,mean,na.rm=T)
  monthly.evi.degrad<-cbind(monthly.evi.degrad,modis.mean.month)
}

monthly.evi.degrad<-as.data.frame(monthly.evi.degrad)
colnames(monthly.evi.degrad)<-sprintf('%02d',1:12)
res <- kmeans(na.omit(monthly.evi.degrad), 5,iter.max = 100,nstart = 20)$cluster
DF=cbind.data.frame(monthly.evi.degrad, 'clus'=NA)
DF[names(res),][,ncol(DF)] <- res
colnames(DF)[1:ncol(monthly.evi.degrad)]<-colnames(monthly.evi.degrad)
colnames(DF)[ncol(DF)]<-"class"
DF$class<-as.factor(DF$class)
DF.agg<-aggregate(.~class,DF,function(x) mean(x,na.rm=T))
table(DF$class)
DF.agg.melt<-melt(DF.agg)
head(DF.agg.melt)
colnames(DF.agg.melt)<-c("class","date","evi")
# DF.agg.melt$month<-format(as.Date(DF.agg.melt$date),"%m")

DF.agg.melt.agg<-aggregate(evi~date+class,DF.agg.melt,mean)
gplot_seasggplot.degrad<-ggplot(DF.agg.melt.agg,aes(x = as.integer(date),y=evi,colour=class))+
  geom_line()+
  scale_color_identity(guide = "legend")
gplot_seasggplot.degrad



vect<-as.vector(modis.evi.para.mx[,1])
vect[]<-NA
vect[pos.degrad]<-DF$class


vect.clas.degrad<-setValues(modis.evi.para[[1]],vect)
plot(vect.clas.degrad,col=rainbow(10))

plot(vect.clas.degrad,col=rainbow(5))
writeRaster(vect.clas.degrad,filename = "results_raster/kmeans_foret_degrad_moyen.tif",datatype='INT2S')


vect.clas.degrad.gg<-rasterToPoints(vect.clas.degrad)
vect.clas.degrad.gg<-as.data.frame(vect.clas.degrad.gg)
head(vect.clas.degrad.gg)
colnames(vect.clas.degrad.gg)<-c("long","lat","class")

carto_foret_degrad<-ggplot()+
  geom_raster(data = vect.clas.degrad.gg,aes(x=long,y=lat,fill=as.factor(class)))+
  geom_polygon(data=para, aes(x=long, y=lat, group=group), 
               fill=NA,color="grey50", size=1)+
  coord_equal()+
  scale_fill_identity(guide = "none")

carto_foret_degrad

TT<-grid.arrange(gplot_seasggplot.degrad,carto_foret_degrad,ncol=2)
ggsave(TT,filename = "figure/bfast_season_result/degrad_forest_5class_profil.jpeg") 


