

# Required packages
library(ggplot2)
library(FSA) 
library(rcompanion) 
library(dplyr)
library(ggfortify)
library(tibble)
library(reshape2)

## Load initial dataset 'metadata'. This contains sample information, phenotypic observation such as biomass and strigolactone
meta<-read.csv("Rice_metadata.csv", row.names = 1)


#### 1. Biomass (Figure 1a) ####

### Subsetting biomass data
colnames(meta)
bm<-subset(meta[,c(8,1:4)],Compartment=="Root")
bm$Genotype2<-factor(bm$Genotype, 
                     c("IAC165", "IAC1246", "GWD", "Dullo","Bina","Sonk","Kinko", "TN1","Bhas","SC","Shiokari","d14","d3","d10","d17","d27"))


### Kruskal wallis test, compare biomass depending on soil type
bm_kw<-kruskal.test(bm$total_fresh_biomass_g_per_plant~bm$Soil, data=bm)
chi_squared<-bm_kw$statistic
Df<-bm_kw$parameter
p<-bm_kw$p.value
Adj_p<-p.adjust(bm_kw$p.value,method="BH" )
bm_kw_res<-cbind(chi_squared,Df,p,Adj_p)
rownames(bm_kw_res)<-"Total biomass_by_soil" 


### Duun test, compare biomass among genotypes in each soil

## in forest soil
fobm_dunn<-subset(bm, Soil=="Forest")
PT<-dunnTest(fobm_dunn$total_fresh_biomass_g_per_plant~Genotype, data=fobm_dunn, method = "bh")
fobm_dunn_res<-PT$res
fobm_dunn_res_letter<-cldList(comparison = fobm_dunn_res$Comparison,p.value = fobm_dunn_res$P.adj,threshold  = 0.05)


## in field soil
fibm_dunn<-subset(bm, Soil=="Forest")
PT<-dunnTest(fibm_dunn$total_fresh_biomass_g_per_plant~Genotype, data=fibm_dunn, method = "bh")
fibm_dunn_res<-PT$res
fibm_dunn_res_letter<-cldList(comparison = fibm_dunn_res$Comparison,p.value = fibm_dunn_res$P.adj,threshold  = 0.05)



### plot 

ggplot(bm) + 
  geom_boxplot(aes(x=Genotype2, y=total_fresh_biomass_g_per_plant, fill=Soil), outlier.colour = NA)+ 
  labs(x="", y = "Total fresh biomass g/plant") + scale_fill_manual(values=c("#a6611a","#018571"))+
  geom_point(aes(x=Genotype2, y=total_fresh_biomass_g_per_plant, fill=Soil), alpha = 0.3, shape = 21, 
             position = position_jitterdodge())+
  theme(axis.text.x = element_text(size = 13, colour = "black",angle=90, hjust=1), 
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(face = "bold", size = 13, vjust = 3), 
        legend.text = element_text(size = 13), legend.title = element_text(size = 13,face = "bold"))




#### 2. Strigolactone (Figure 1b and 1c) ####
#As SLs were not detected in field soil, further SL-related studies are only samples obtained from forest soil

### PCA
colnames(meta)

foSLs<-subset(meta, Soil_compartment=="Fo_RS"&SL_analysis=="yes") 
foSLs_pca<-prcomp(foSLs[11:13], scale. = T)
autoplot(foSLs_pca, data=foSLs_pca)
autoplot(foSLs_pca, data = foSLs_pca, loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3, label.colour = 'blue')


### Dunn test, compare level of SLs among genotypes. 
foSLs_dunn<-foSLs %>% filter(!Genotype%in%c('Kinko','TN1','Bhas','SC','d10','d17','d27')) # these genotypes are removed due to zeros (non-detected) which is not appropriate for testing 
foSLs_dunn<-foSLs_dunn[,c(11:13,1)]

indices=2 #number of indices that I am testing, As level of MeO5DS was not significantly different among genotype, it was excluded. 
Z<-as.data.frame(matrix(NA, 36, indices)) #results list =36
P.unadj<-as.data.frame(matrix(NA, 36, indices)) #results list =36
P.adj<-as.data.frame(matrix(NA, 36, indices)) #results list =36
foSLs_dunn_letter<-as.data.frame(matrix(NA, 9, indices)) #results list =9

for(i in 1:indices) {
  PT<-dunnTest(foSLs_dunn[,i]~Genotype, data=foSLs_dunn, method = "bh")
  Z[,i]<-PT$res$Z
  P.unadj[,i]<-PT$res$P.unadj
  P.adj[,i]<-PT$res$P.adj
  PT2<-PT$res
  cl<-cldList(comparison = PT2$Comparison, p.value = PT2$P.adj,threshold  = 0.05)
  foSLs_dunn_letter[,i]<-cl$Letter
}

foSLs_dunn_res<-cbind(Z,P.unadj, P.adj)
rownames(foSLs_dunn_res)<- PT$res$Comparison
colnames(foSLs_dunn_res)<-c("Z_orobanchol","Z_4DO",
                             "unadjusted.P_orobanchol","unadjusted.P_4DO",
                             "adjusted.P_orobanchol","adjusted.P_4DO")
rownames(foSLs_dunn_letter) <- cl$Group
colnames(foSLs_dunn_letter) <- colnames(foSLs_dunn[1:indices])


### SL plot
foSLs$Genotype<-droplevels(as.factor(foSLs$Genotype))
foSLs$Genotype2<-factor(foSLs$Genotype, 
                         c("IAC165", "IAC1246", "GWD", "Dullo","Bina","Sonk","Kinko", "TN1","Bhas","SC","Shiokari","d3","d14","d10","d17","d27"))

ggplot(foSLs)+
  geom_boxplot(aes(x=Genotype2, y=X4DO_pmol_g), outlier.colour = NA)+labs(x="Rice genotype", y = "4DO pmol/gFW") + 
  geom_point(aes(x=Genotype2, y=X4DO_pmol_g), alpha = 0.3, shape = 21) +
  theme(axis.text.x = element_text(size = 13, colour = "black",angle=90, hjust=1), 
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(face = "bold", size = 13, vjust = 3), 
        legend.text = element_text(size = 13), legend.title = element_text(size = 13,face = "bold"))

