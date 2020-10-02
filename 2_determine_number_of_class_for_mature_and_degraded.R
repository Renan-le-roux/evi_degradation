
library(factoextra)
library(NbClust)
library(xts)
library(raster)
library(gridExtra)
library(reshape2)

setwd("H:/Utilisateurs/leroux_r/post_doc_cnes/CNES/Travail/FEU_EVI/Para/")

###########################################Clustering degraded forest####################
load("res.bast_parallele.degrad.RData" )


res.bfast.xts<-list()

for (i in 1:length(res.bast_parallele.degrad)){
  print(i)
  X<-xts(res.bast_parallele.degrad[[i]][[1]], as.Date(time(res.bast_parallele.degrad[[i]][[1]])))
  
  
  # colnames(X)<-c("evi.seas","date")
  res.bfast.xts[[i]]<-X
  # res.date<-cbind(res.date,X$date)
}
res.bfast.xts<-do.call("cbind",res.bfast.xts)
subset_xts<-res.bfast.xts['2010']
subset_xts.mx<-as.matrix(t(subset_xts))
dim(subset_xts.mx)
colnames(subset_xts.mx)
rownames(subset_xts.mx)<-1:nrow(subset_xts.mx)
save(subset_xts.mx,file = "profil_bfast_degrad.RData")

kmeans.sample<-na.omit(subset_xts.mx)
kmeans.sample.df<-as.data.frame(kmeans.sample)
sasample(x = 1:nrow(kmeans.sample.df)))

res <- kmeans(na.omit(subset_xts.mx), 5)$cluster

elbow_method<-fviz_nbclust(kmeans.sample, kmeans, method = "wss")
elbow_method_lab<-elbow_method+ geom_vline(xintercept = 4, linetype = 2)+
  labs(subtitle = "Elbow method")


silhouette_methode<-fviz_nbclust(kmeans.sample, kmeans, method = "silhouette")
silhouette_methode_lab<-silhouette_methode+labs(subtitle = "Silhouette method")

gap_method<-fviz_nbclust(kmeans.sample, kmeans, nstart = 25,  method = "gap_stat", nboot = 50)
gap_method_lab<-gap_method+ labs(subtitle = "Gap statistic method")
TT<-grid.arrange(elbow_method_lab,silhouette_methode_lab,gap_method_lab,ncol=2)
ggsave(TT,filename = "figure/optimal_number_of_cluster_for_degraded_forest.jpeg")                 
  
###########################Figure degrad profile and map
load("profil_bfast_degrad.RData")
raster.degrad<-raster("results_raster/degradation_map.rsp.tif")
raster.degrad.vc<-as.vector(raster.degrad)
pos.nona.degrad<-which(!is.na(raster.degrad.vc))


res <- kmeans(na.omit(subset_xts.mx),4)$cluster
DF=cbind.data.frame(subset_xts.mx, 'clus'=NA)
DF[names(res),][,ncol(DF)] <- res
colnames(DF)[1:ncol(subset_xts.mx)]<-colnames(subset_xts.mx)
colnames(DF)[ncol(DF)]<-"class"
DF$class<-as.factor(DF$class)
DF.agg<-aggregate(.~class,DF,function(x) mean(x,na.rm=T))
save(DF.agg,file = "classif.profil.degrad.RData")
DF.agg.melt<-melt(DF.agg)
head(DF.agg.melt)
colnames(DF.agg.melt)<-c("class","date","evi")
DF.agg.melt$month<-format((DF.agg.melt$date))
DF.agg.melt.agg<-aggregate(evi~month+class,DF.agg.melt,mean)
DF.agg.melt.agg$month<-as.Date(DF.agg.melt.agg$month)
gplot_seasggplot<-ggplot(DF.agg.melt.agg,aes(x =month,y=evi,colour=class))+
  geom_line()+
  scale_color_identity(guide = "legend")+
  scale_x_date(date_labels = "%m",date_breaks = "1 month")+
  xlab("Month")+
  ylab("Seasonal profile")
gplot_seasggplot
ggsave(gplot_seasggplot,filename = "figure/bfast_season_result/kmeans_profile_degraded_forest_elbow4.jpeg")


vect<-raster.degrad.vc
vect[]<-NA
vect[pos.nona.degrad]<-DF$class


vect.clas.degrad<-setValues(raster.degrad[[1]],vect)
plot(vect.clas.degrad,col=rainbow(4))

plot(vect.clas.degrad,col=rainbow(5))
writeRaster(vect.clas.degrad,filename = "results_raster/kmeans_degrad_elbow4_2.tif",datatype='INT2S',overwrite=T)



###########################################Clustering mature forest####################
load("res.bast_parallele.mature.RData" )


res.bfast.xts<-list()

for (i in 1:length(res.bast_parallele.mature)){
  print(i)
  X<-xts(res.bast_parallele.mature[[i]][[1]], as.Date(time(res.bast_parallele.mature[[i]][[1]])))
  
  
  # colnames(X)<-c("evi.seas","date")
  res.bfast.xts[[i]]<-X
  # res.date<-cbind(res.date,X$date)
}
res.bfast.xts<-do.call("cbind",res.bfast.xts)
subset_xts<-res.bfast.xts['2010']
subset_xts.mx<-as.matrix(t(subset_xts))
dim(subset_xts.mx)
colnames(subset_xts.mx)
rownames(subset_xts.mx)<-1:nrow(subset_xts.mx)
save(subset_xts.mx,file = "profil_bfast_mature.RData")

kmeans.sample<-na.omit(subset_xts.mx)
kmeans.sample.df<-as.data.frame(kmeans.sample)
NR<-sample(x = 1:nrow(kmeans.sample.df),size = 10000,replace = F)
sample.kmeans<-kmeans.sample.df[NR,]

elbow_method<-fviz_nbclust(sample.kmeans, kmeans, method = "wss")
elbow_method_lab<-elbow_method+ geom_vline(xintercept = 2, linetype = 2)+
  labs(subtitle = "Elbow method")


silhouette_methode<-fviz_nbclust(sample.kmeans, kmeans, method = "silhouette")
silhouette_methode_lab<-silhouette_methode+labs(subtitle = "Silhouette method")

gap_method<-fviz_nbclust(sample.kmeans, kmeans, nstart = 25,  method = "gap_stat", nboot = 50)
gap_method_lab<-gap_method+ labs(subtitle = "Gap statistic method")
TT<-grid.arrange(elbow_method_lab,silhouette_methode_lab,ncol=2)
ggsave(TT,filename = "figure/optimal_number_of_cluster_for_matureed_forest.jpeg")                 

###########################Figure mature profile and map
load("profil_bfast_mature.RData")
raster.mature<-raster("results_raster/MF_mature.rsp.tif")
raster.mature.vc<-as.vector(raster.mature)
pos.nona.mature<-which(!is.na(raster.mature.vc))

res <- kmeans(na.omit(subset_xts.mx),2)$cluster
DF=cbind.data.frame(subset_xts.mx, 'clus'=NA)
DF[names(res),][,ncol(DF)] <- res
colnames(DF)[1:ncol(subset_xts.mx)]<-colnames(subset_xts.mx)
colnames(DF)[ncol(DF)]<-"class"
DF$class<-as.factor(DF$class)
DF.agg<-aggregate(.~class,DF,function(x) mean(x,na.rm=T))
save(DF.agg,file = "classif.profil.mature.RData")
DF.agg.melt<-melt(DF.agg)
head(DF.agg.melt)
colnames(DF.agg.melt)<-c("class","date","evi")
DF.agg.melt$month<-as.Date(DF.agg.melt$date)
DF.agg.melt.agg<-aggregate(evi~month+class,DF.agg.melt,mean)

gplot_seasggplot<-ggplot(DF.agg.melt.agg,aes(x = month,y=evi,colour=class))+
  geom_line()+
  scale_color_identity(guide = "legend")+
  scale_x_date(date_labels = "%m",date_breaks = "1 month")+
  xlab("Month")+
  ylab("Seasonal profile")
gplot_seasggplot
ggsave(gplot_seasggplot,filename = "figure/bfast_season_result/kmeans_profile_mature_forest_elbow2.jpeg")


vect<-raster.mature.vc
vect[]<-NA
vect[pos.nona.mature]<-DF$class
length(which(!is.na(vect)))

vect.clas.mature<-setValues(raster.mature[[1]],vect)
plot(vect.clas.mature,col=rainbow(4))

plot(vect.clas.mature,col=rainbow(5))
writeRaster(vect.clas.mature,filename = "results_raster/kmeans_mature_elbow2.tif",datatype='INT2S')

