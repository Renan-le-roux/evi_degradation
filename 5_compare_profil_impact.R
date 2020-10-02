library(factoextra)
library(NbClust)
library(xts)
library(raster)
library(gridExtra)
library(reshape2)
library(rgdal)
library(FactoMineR)
##########################Description of classe by degradation and climate#############
cumul_degrad_map<-stack("results_raster/tot_degrad_1km.grd")
profile.degraded<-raster("results_raster/kmeans_degrad_elbow4_2.tif")
cumul_degrad.mx<-as.matrix(cumul_degrad_map)

profile.degraded.vc<-as.vector(profile.degraded)
class_impact<-cbind.data.frame(profile.degraded.vc,cumul_degrad.mx)
class_impact$profile.degraded.vc<-as.factor(class_impact$profile.degraded.vc)
colnames(class_impact)<-c("class",names(cumul_degrad_map))



class_impact.agg<-aggregate(.~class,class_impact,mean)
class_impact.agg.melt<-melt(class_impact.agg)
colnames(class_impact.agg.melt)<-c("class","impact","Number")
 class_impact.agg.melt$impact<-factor(class_impact.agg.melt$impact,level=c("logging","burned","logged_burned"))
impact_per_profil<-ggplot()+
    geom_bar(data = class_impact.agg.melt,aes(x=class,y=Number,fill=impact),stat = "identity")+
    scale_fill_manual(values = c("bisque3","red","firebrick4"))+
    xlab("profil class")+
    ylab("Average number of impact per 1km²")
impact_per_profil
ggsave(impact_per_profil,filename =   "figure/barplot_impact_profil_non_dodge.jpeg")
