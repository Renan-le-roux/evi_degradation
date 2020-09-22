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
library(rasterVis)
library(gridExtra)
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
    
    
    seasonal.degrad<-list(seasonal.deg,bp)
  }
  # T1<-Sys.time()
  
  
  return(seasonal.degrad)
}


###############################Classif profil mature#############

degrad_map.vc<-as.vector(degrad_map)
MF_mature.vc
pos.mature<-which(MF_mature.vc==1)


length(pos.mature)
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


res.bfast.xts<-list()

for (i in 1:length(res.bast_parallele.mature)){
  print(i)
  X<-xts(res.bast_parallele.mature[[i]][[1]], as.Date(time(res.bast_parallele.mature[[i]][[1]])))
  
  
  # colnames(X)<-c("evi.seas","date")
  res.bfast.xts[[i]]<-X
  # res.date<-cbind(res.date,X$date)
}
res.bfast.xts<-do.call("cbind",res.bfast.xts)
res.date<-matrix(nrow = length(res.bast_parallele.mature),ncol = 5)
for (i in 1:length(res.bast_parallele.mature)){
  print(i)
  dat.X<-res.bast_parallele.mature[[i]][[2]]
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



res <- kmeans(na.omit(subset_xts.mx), 5)$cluster
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
ggsave(gplot_seasggplot,filename = "figure/bfast_season_result/profil_mature_5C_forest.jpeg")




vect<-as.vector(modis.evi.para.mx[,1])
vect[]<-NA
vect[pos.mature]<-DF$class





vect.clas.degrad<-setValues(modis.evi.para[[1]],vect)
plot(vect.clas.degrad,col=rainbow(10))

plot(vect.clas.degrad,col=rainbow(5))
writeRaster(vect.clas.degrad,filename = "results_raster/kmeans_foret_mature5.tif",datatype='INT2S')


###############################Classif profil degrad#############

degrad_map.vc<-as.vector(degrad_map)
MF_degrad.vc
pos.degrad<-which(!is.na(degrad_map.vc))


length(pos.degrad)
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

res.bfast.xts<-list()

for (i in 1:length(res.bast_parallele.degrad)){
  print(i)
  X<-xts(res.bast_parallele.degrad[[i]][[1]], as.Date(time(res.bast_parallele.degrad[[i]][[1]])))
  
  
  # colnames(X)<-c("evi.seas","date")
  res.bfast.xts[[i]]<-X
  # res.date<-cbind(res.date,X$date)
}
res.bfast.xts<-do.call("cbind",res.bfast.xts)
res.date<-matrix(nrow = length(res.bast_parallele.degrad),ncol = 5)
for (i in 1:length(res.bast_parallele.degrad)){
  print(i)
  dat.X<-res.bast_parallele.degrad[[i]][[2]]
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



res <- kmeans(na.omit(subset_xts.mx), 5)$cluster
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
gplot_seasggplot.degrad<-ggplot(DF.agg.melt.agg,aes(x = as.integer(month),y=evi,colour=class))+
  geom_line()+
  scale_color_identity(guide = "legend")
gplot_seasggplot.degrad
ggsave(gplot_seasggplot,filename = "figure/bfast_season_result/profil_degrad_5C_forest.jpeg")




vect<-as.vector(modis.evi.para.mx[,1])
vect[]<-NA
vect[pos.degrad]<-DF$class


vect.clas.degrad<-setValues(modis.evi.para[[1]],vect)
plot(vect.clas.degrad,col=rainbow(10))

plot(vect.clas.degrad,col=rainbow(5))
writeRaster(vect.clas.degrad,filename = "results_raster/kmeans_foret_degrad5.tif",datatype='INT2S')


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
ggsave(TT,filename = "figure/bfast_season_result/degrad_forest_5class.jpeg") 
