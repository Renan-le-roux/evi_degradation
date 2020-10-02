
library(factoextra)
library(NbClust)
library(xts)
library(raster)
library(gridExtra)
library(reshape2)
library(rgdal)
setwd("H:/Utilisateurs/leroux_r/post_doc_cnes/CNES/Travail/FEU_EVI/Para/")

#################################################Site 1##############################
###########################################Clustering degraded forest####################
load("profil_bfast_degrad.RData" )
classif.degrad<-raster("results_raster/kmeans_degrad_elbow4.tif")


classif.degrad.vc<-as.vector(classif.degrad)
pos.nona<-which(!is.na(classif.degrad.vc))
study_site<-readOGR("../../raw_data/study_site.shp")
st<-study_site[1,]
st<-spTransform(st,CRSobj = proj4string(classif.degrad))
kmeans.sample<-na.omit(subset_xts.mx)
dim(kmeans.sample)
n<-ncol(kmeans.sample)
rast.ref <- stack(replicate(n, classif.degrad))
r<-matrix(nrow = length(classif.degrad.vc),ncol=ncol(kmeans.sample))
r[pos.nona,]<-kmeans.sample
# save(subset_xts.mx,file = "profil_bfast_degrad.RData")
profile_raster.degrad.site1<-setValues(rast.ref,r)
profile_raster.degrad.site1<-crop(profile_raster.degrad.site1,st)
kmeans.sample.mx.site1<-as.matrix(profile_raster.degrad.site1)
kmeans.sample.df<-as.data.frame(kmeans.sample.mx.site1)
kmeans.sample.df<-na.omit(kmeans.sample.df)
colnames(kmeans.sample.df)<-colnames(subset_xts.mx)

elbow_method<-fviz_nbclust(kmeans.sample.df, kmeans, method = "wss")
elbow_method_lab<-elbow_method+ geom_vline(xintercept = 4, linetype = 2)+
  labs(subtitle = "Elbow method")


silhouette_methode<-fviz_nbclust(kmeans.sample, kmeans, method = "silhouette")
silhouette_methode_lab<-silhouette_methode+labs(subtitle = "Silhouette method")

gap_method<-fviz_nbclust(kmeans.sample, kmeans, nstart = 25,  method = "gap_stat", nboot = 50)
gap_method_lab<-gap_method+ labs(subtitle = "Gap statistic method")
TT<-grid.arrange(elbow_method_lab,silhouette_methode_lab,gap_method_lab,ncol=2)
ggsave(TT,filename = "figure/optimal_number_of_cluster_for_degraded_forest.jpeg")                 

###########################Figure degrad profile and map
raster.degrad<-raster("results_raster/degradation_map.rsp.tif")
raster.degrad.vc<-as.vector(raster.degrad)
pos.nona.degrad<-which(!is.na(raster.degrad.vc))


res <- kmeans(kmeans.sample.df,4)$cluster
DF=cbind.data.frame(kmeans.sample.df, 'clus'=NA)
DF[names(res),][,ncol(DF)] <- res
colnames(DF)[1:ncol(kmeans.sample.df)]<-colnames(kmeans.sample.df)
colnames(DF)[ncol(DF)]<-"class"
DF$class<-as.factor(DF$class)
DF.agg<-aggregate(.~class,DF,function(x) mean(x,na.rm=T))
DF.agg.melt<-melt(DF.agg)
head(DF.agg.melt)
colnames(DF.agg.melt)<-c("class","date","evi")
DF.agg.melt$month<-format(as.Date(DF.agg.melt$date),"%m")
DF.agg.melt.agg<-aggregate(evi~month+class,DF.agg.melt,mean)
gplot_seasggplot<-ggplot(DF.agg.melt.agg,aes(x = as.integer(month),y=evi,colour=class))+
  geom_line()+
  scale_color_identity(guide = "legend")
gplot_seasggplot
ggsave(gplot_seasggplot,filename = "figure/bfast_season_result/kmeans_profile_degraded_forest_elbow4_site1.jpeg")


vect<-as.vector(profile_raster.degrad.site1[[1]])
vect[which(!is.na(vect))]<-DF$class


vect.clas.degrad<-setValues(profile_raster.degrad.site1[[1]],vect)
plot(vect.clas.degrad,col=rainbow(4))

plot(vect.clas.degrad,col=rainbow(5))
writeRaster(vect.clas.degrad,filename = "results_raster/kmeans_degrad_elbow4_site1.tif",datatype='INT2S')




###########################################Clustering mature forest####################
load("profil_bfast_mature.RData" )
classif.mature<-raster("results_raster/kmeans_mature_elbow2.tif")


classif.mature.vc<-as.vector(classif.mature)
pos.nona<-which(!is.na(classif.mature.vc))
study_site<-readOGR("../../raw_data/study_site.shp")
st<-study_site[1,]
st<-spTransform(st,CRSobj = proj4string(classif.mature))
kmeans.sample<-na.omit(subset_xts.mx)
dim(kmeans.sample)
n<-ncol(kmeans.sample)
rast.ref <- stack(replicate(n, classif.mature))
r<-matrix(nrow = length(classif.mature.vc),ncol=ncol(kmeans.sample))
r[pos.nona,]<-kmeans.sample
# save(subset_xts.mx,file = "profil_bfast_mature.RData")
profile_raster.mature.site1<-setValues(rast.ref,r)
profile_raster.mature.site1<-crop(profile_raster.mature.site1,st)
kmeans.sample.mx.site1<-as.matrix(profile_raster.mature.site1)
kmeans.sample.df<-as.data.frame(kmeans.sample.mx.site1)
kmeans.sample.df<-na.omit(kmeans.sample.df)
colnames(kmeans.sample.df)<-colnames(subset_xts.mx)

# elbow_method<-fviz_nbclust(kmeans.sample.df, kmeans, method = "wss")
# elbow_method_lab<-elbow_method+ geom_vline(xintercept = 4, linetype = 2)+
#   labs(subtitle = "Elbow method")
# 
# 
# silhouette_methode<-fviz_nbclust(kmeans.sample, kmeans, method = "silhouette")
# silhouette_methode_lab<-silhouette_methode+labs(subtitle = "Silhouette method")
# 
# gap_method<-fviz_nbclust(kmeans.sample, kmeans, nstart = 25,  method = "gap_stat", nboot = 50)
# gap_method_lab<-gap_method+ labs(subtitle = "Gap statistic method")
# TT<-grid.arrange(elbow_method_lab,silhouette_methode_lab,gap_method_lab,ncol=2)
# ggsave(TT,filename = "figure/optimal_number_of_cluster_for_matureed_forest.jpeg")                 

###########################Figure mature profile and map
raster.mature<-raster("results_raster/degradation_map.rsp.tif")
raster.mature.vc<-as.vector(raster.mature)
pos.nona.mature<-which(!is.na(raster.mature.vc))


res <- kmeans(kmeans.sample.df,2)$cluster
DF=cbind.data.frame(kmeans.sample.df, 'clus'=NA)
DF[names(res),][,ncol(DF)] <- res
colnames(DF)[1:ncol(kmeans.sample.df)]<-colnames(kmeans.sample.df)
colnames(DF)[ncol(DF)]<-"class"
DF$class<-as.factor(DF$class)
DF.agg<-aggregate(.~class,DF,function(x) mean(x,na.rm=T))
DF.agg.melt<-melt(DF.agg)
head(DF.agg.melt)
colnames(DF.agg.melt)<-c("class","date","evi")
DF.agg.melt$month<-format(as.Date(DF.agg.melt$date),"%m")
DF.agg.melt.agg<-aggregate(evi~month+class,DF.agg.melt,mean)
gplot_seasggplot<-ggplot(DF.agg.melt.agg,aes(x = as.integer(month),y=evi,colour=class))+
  geom_line()+
  scale_color_identity(guide = "legend")
gplot_seasggplot
ggsave(gplot_seasggplot,filename = "figure/bfast_season_result/kmeans_profile_mature_forest_elbow2_site1.jpeg")


vect<-as.vector(profile_raster.mature.site1[[1]])
vect[which(!is.na(vect))]<-DF$class


vect.clas.mature<-setValues(profile_raster.mature.site1[[1]],vect)
plot(vect.clas.mature,col=rainbow(4))

plot(vect.clas.mature,col=rainbow(5))
writeRaster(vect.clas.mature,filename = "results_raster/kmeans_mature_elbow4_site1.tif",datatype='INT2S')





#################################################Site 2##############################
###########################################Clustering degraded forest####################
load("profil_bfast_degrad.RData" )
classif.degrad<-raster("results_raster/kmeans_degrad_elbow4.tif")


classif.degrad.vc<-as.vector(classif.degrad)
pos.nona<-which(!is.na(classif.degrad.vc))
study_site<-readOGR("../../raw_data/study_site.shp")
st<-study_site[2,]
st<-spTransform(st,CRSobj = proj4string(classif.degrad))
kmeans.sample<-na.omit(subset_xts.mx)
dim(kmeans.sample)
n<-ncol(kmeans.sample)
rast.ref <- stack(replicate(n, classif.degrad))
r<-matrix(nrow = length(classif.degrad.vc),ncol=ncol(kmeans.sample))
r[pos.nona,]<-kmeans.sample
# save(subset_xts.mx,file = "profil_bfast_degrad.RData")
profile_raster.degrad.site2<-setValues(rast.ref,r)
profile_raster.degrad.site2<-crop(profile_raster.degrad.site2,st)
kmeans.sample.mx.site2<-as.matrix(profile_raster.degrad.site2)
kmeans.sample.df<-as.data.frame(kmeans.sample.mx.site2)
kmeans.sample.df<-na.omit(kmeans.sample.df)
colnames(kmeans.sample.df)<-colnames(subset_xts.mx)

# elbow_method<-fviz_nbclust(kmeans.sample.df, kmeans, method = "wss")
# elbow_method_lab<-elbow_method+ geom_vline(xintercept = 4, linetype = 2)+
#   labs(subtitle = "Elbow method")
# 
# 
# silhouette_methode<-fviz_nbclust(kmeans.sample, kmeans, method = "silhouette")
# silhouette_methode_lab<-silhouette_methode+labs(subtitle = "Silhouette method")
# 
# gap_method<-fviz_nbclust(kmeans.sample, kmeans, nstart = 25,  method = "gap_stat", nboot = 50)
# gap_method_lab<-gap_method+ labs(subtitle = "Gap statistic method")
# TT<-grid.arrange(elbow_method_lab,silhouette_methode_lab,gap_method_lab,ncol=2)
# ggsave(TT,filename = "figure/optimal_number_of_cluster_for_degraded_forest.jpeg")                 
# 
# ###########################Figure degrad profile and map
raster.degrad<-raster("results_raster/degradation_map.rsp.tif")
raster.degrad.vc<-as.vector(raster.degrad)
pos.nona.degrad<-which(!is.na(raster.degrad.vc))


res <- kmeans(kmeans.sample.df,4)$cluster
DF=cbind.data.frame(kmeans.sample.df, 'clus'=NA)
DF[names(res),][,ncol(DF)] <- res
colnames(DF)[1:ncol(kmeans.sample.df)]<-colnames(kmeans.sample.df)
colnames(DF)[ncol(DF)]<-"class"
DF$class<-as.factor(DF$class)
DF.agg<-aggregate(.~class,DF,function(x) mean(x,na.rm=T))
DF.agg.melt<-melt(DF.agg)
head(DF.agg.melt)
colnames(DF.agg.melt)<-c("class","date","evi")
DF.agg.melt$month<-format(as.Date(DF.agg.melt$date),"%m")
DF.agg.melt.agg<-aggregate(evi~month+class,DF.agg.melt,mean)
gplot_seasggplot<-ggplot(DF.agg.melt.agg,aes(x = as.integer(month),y=evi,colour=class))+
  geom_line()+
  scale_color_identity(guide = "legend")
gplot_seasggplot
ggsave(gplot_seasggplot,filename = "figure/bfast_season_result/kmeans_profile_degraded_forest_elbow4_site2.jpeg")


vect<-as.vector(profile_raster.degrad.site2[[1]])
vect[which(!is.na(vect))]<-DF$class


vect.clas.degrad<-setValues(profile_raster.degrad.site2[[1]],vect)
plot(vect.clas.degrad,col=rainbow(4))

plot(vect.clas.degrad,col=rainbow(5))
writeRaster(vect.clas.degrad,filename = "results_raster/kmeans_degrad_elbow4_site2.tif",datatype='INT2S',overwrite=T)




###########################################Clustering mature forest####################
load("profil_bfast_mature.RData" )
classif.mature<-raster("results_raster/kmeans_mature_elbow2.tif")


classif.mature.vc<-as.vector(classif.mature)
pos.nona<-which(!is.na(classif.mature.vc))
study_site<-readOGR("../../raw_data/study_site.shp")
st<-study_site[2,]
st<-spTransform(st,CRSobj = proj4string(classif.mature))
kmeans.sample<-na.omit(subset_xts.mx)
dim(kmeans.sample)
n<-ncol(kmeans.sample)
rast.ref <- stack(replicate(n, classif.mature))
r<-matrix(nrow = length(classif.mature.vc),ncol=ncol(kmeans.sample))
r[pos.nona,]<-kmeans.sample
# save(subset_xts.mx,file = "profil_bfast_mature.RData")
profile_raster.mature.site2<-setValues(rast.ref,r)
profile_raster.mature.site2<-crop(profile_raster.mature.site2,st)
kmeans.sample.mx.site2<-as.matrix(profile_raster.mature.site2)
kmeans.sample.df<-as.data.frame(kmeans.sample.mx.site2)
kmeans.sample.df<-na.omit(kmeans.sample.df)
colnames(kmeans.sample.df)<-colnames(subset_xts.mx)

# elbow_method<-fviz_nbclust(kmeans.sample.df, kmeans, method = "wss")
# elbow_method_lab<-elbow_method+ geom_vline(xintercept = 4, linetype = 2)+
#   labs(subtitle = "Elbow method")
# 
# 
# silhouette_methode<-fviz_nbclust(kmeans.sample, kmeans, method = "silhouette")
# silhouette_methode_lab<-silhouette_methode+labs(subtitle = "Silhouette method")
# 
# gap_method<-fviz_nbclust(kmeans.sample, kmeans, nstart = 25,  method = "gap_stat", nboot = 50)
# gap_method_lab<-gap_method+ labs(subtitle = "Gap statistic method")
# TT<-grid.arrange(elbow_method_lab,silhouette_methode_lab,gap_method_lab,ncol=2)
# ggsave(TT,filename = "figure/optimal_number_of_cluster_for_matureed_forest.jpeg")                 

###########################Figure mature profile and map
raster.mature<-raster("results_raster/degradation_map.rsp.tif")
raster.mature.vc<-as.vector(raster.mature)
pos.nona.mature<-which(!is.na(raster.mature.vc))


res <- kmeans(kmeans.sample.df,2)$cluster
DF=cbind.data.frame(kmeans.sample.df, 'clus'=NA)
DF[names(res),][,ncol(DF)] <- res
colnames(DF)[1:ncol(kmeans.sample.df)]<-colnames(kmeans.sample.df)
colnames(DF)[ncol(DF)]<-"class"
DF$class<-as.factor(DF$class)
DF.agg<-aggregate(.~class,DF,function(x) mean(x,na.rm=T))
DF.agg.melt<-melt(DF.agg)
head(DF.agg.melt)
colnames(DF.agg.melt)<-c("class","date","evi")
DF.agg.melt$month<-format(as.Date(DF.agg.melt$date),"%m")
DF.agg.melt.agg<-aggregate(evi~month+class,DF.agg.melt,mean)
gplot_seasggplot<-ggplot(DF.agg.melt.agg,aes(x = as.integer(month),y=evi,colour=class))+
  geom_line()+
  scale_color_identity(guide = "legend")
gplot_seasggplot
ggsave(gplot_seasggplot,filename = "figure/bfast_season_result/kmeans_profile_mature_forest_elbow2_site2.jpeg")


vect<-as.vector(profile_raster.mature.site2[[1]])
vect[which(!is.na(vect))]<-DF$class


vect.clas.mature<-setValues(profile_raster.mature.site2[[1]],vect)
plot(vect.clas.mature,col=rainbow(4))

plot(vect.clas.mature,col=rainbow(5))
writeRaster(vect.clas.mature,filename = "results_raster/kmeans_mature_elbow4_site2.tif",datatype='INT2S',overwrite=T)



