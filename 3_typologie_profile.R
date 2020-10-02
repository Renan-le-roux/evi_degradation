library(raster)
setwd("H:/Utilisateurs/leroux_r/post_doc_cnes/CNES/Travail/FEU_EVI/Para/")
######################Classif forêts matures

classif.mature<-raster("results_raster/kmeans_mature_elbow2.tif")
plot(classif.mature,col=rainbow(10))
classif.mature.vc<-as.vector(classif.mature)
classif.degrad<-raster("results_raster/kmeans_degrad_elbow4.tif")
classif.degrad.vc<-as.vector(classif.degrad)
raster.degrad<-raster("results_raster/degradation_map.rsp.tif")
raster.degrad.vc<-as.vector(raster.degrad)
length(raster.degrad.vc)
pos.non.na.degrad<-which(!is.na(raster.degrad.vc))
length(pos.non.na.degrad)

plot(classif.degrad,col=rainbow(4))
pos.nona.mature<-which(!is.na(classif.mature.vc))
pos.nona.degrad<-which(!is.na(classif.degrad.vc))
# plot the raster

# create sample points and add them to the plot
xy = rasterToPoints(classif.degrad)
nrow(xy)

compare_classif<-NULL

for (i in 1:nrow(xy)){
  
print(i)
dis.fr.point<-distanceFromPoints(classif.mature, xy[i,1:2])
dis.fr.point<-dis.fr.point[pos.nona.mature]
pos.min<-which.min(dis.fr.point)
val.mat<-classif.mature.vc[pos.nona.mature][pos.min]
val.deg<-xy[i,3]
res<-cbind(val.mat,val.deg)
compare_classif<-rbind(compare_classif,res)
}
# show output of both procedures

head(compare_classif)
compare_classif.df<-as.data.frame(compare_classif)
# compare_classif.df[,1]<-1
# compare_classif.df$classif_degrad<-classif.degrad.vc[which(!is.na(raster.degrad.vc))]
head(compare_classif.df)
unique_combination<-unique(compare_classif.df[,c(1,2)])
nrow(unique_combination)

compare_classif.df$combi<-paste0(compare_classif.df[,1],compare_classif.df[,2])
unique(compare_classif.df$combi)
#####################################Associate profile with class#################
load("classif.profil.degrad.RData")
profil.degrad<-DF.agg
load("classif.profil.mature.RData")
profil.mature<-DF.agg
######################################################

################################Associate profile with typo#####
typologie_profile_forets_degradees<-NULL
for (i in 1:length(unique(compare_classif.df$combi))){
  pos.typo<-unique(compare_classif.df$combi)[i]
  pos.typo.class<-which(compare_classif.df$combi==pos.typo)[1]
  
  
  class.degrad<-compare_classif.df$val.deg[pos.typo.class]
  class.mature<-compare_classif.df$val.mat[pos.typo.class]
  
  pos.profil.mature<-which(profil.mature$class==class.mature)
  prof.mature.typo<-profil.mature[pos.profil.mature,]
  pos.profil.degrad<-which(profil.degrad$class==class.degrad)
  prof.degrad.typo<-profil.degrad[pos.profil.degrad,]
  prof.mature.typo<-as.data.frame(prof.mature.typo)
  compare_classif.df$combi[pos.typo.class]
  prof.mature.typo$typo<-pos.typo
  prof.degrad.typo$typo<-pos.typo
  prof.typo<-rbind(prof.mature.typo,prof.degrad.typo)
  prof.typo$class<-c("mature","degraded")
  typologie_profile_forets_degradees<-rbind(typologie_profile_forets_degradees,prof.typo)
}




typologie_profile_forets_degradees.melt<-melt(typologie_profile_forets_degradees,id.vars = c("class","typo"))
colnames(typologie_profile_forets_degradees.melt)<-c("class","typo","month","seasonal")
typologie_profile_forets_degradees.melt$month<-as.Date(typologie_profile_forets_degradees.melt$month)
graph.typo.gg<-ggplot()+
  geom_line(data = typologie_profile_forets_degradees.melt,aes(x=month,y=seasonal,color=class))+
  facet_wrap(~typo,nrow=2)+
  # scale_color_identity(guide = "legend")+
  scale_x_date(date_labels = "%m",date_breaks = "1 month")+
  xlab("Month")+
  ylab("Seasonal profile")
graph.typo.gg
ggsave(graph.typo.gg,filename = "figure/bfast_season_result/profil_typo4_2mature.jpeg")

###################☻Cartographie typologie##############
typo.vc<-raster.degrad.vc
typo.vc[]<-NA
typo.vc[pos.nona.degrad]<-as.integer(compare_classif.df$combi)
# typo.vc<-as.integer(typo.vc)
rast.typo<-setValues(raster.degrad,typo.vc)
plot(rast.typo,col=rainbow(19))
writeRaster(rast.typo,filename = "results_raster/rastertypo4.tif",overwrite=T)
