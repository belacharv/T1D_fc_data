###########################################
## VALIDATION OF MANUALLY GATED CLUSTERS ##
###########################################

getwd()
library(tidyverse)
library(ggpubr)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
clusters <- list()
selected <- c(10,7,12,3,15,11)
all_selected_clusters<- tibble()
setwd("exported/exported_cd3/260411/260411_cluster_results")
for (i in 1:length(selected)){
  print(i)
  print(selected[i])
  cluster_name <- paste("Pheno_Cluster_",selected[i],".csv",sep="")
  clust <- read_csv(cluster_name)
  colnames(clust)[5] <- paste("Cluster_",selected[i],sep="")
  clusters[[i]] <- clust
  #colnames(clust)[5] <- paste("Cluster_",selected[i],sep="")
  if (i == 1){
    all_selected_clusters <- clust
  }
  else{
    all_selected_clusters <- full_join(all_selected_clusters,clust[,-c(1,2,4)],by="sample_var")
  }
  #full_join(all_selected_clusters,)
}
#cluster <- read_csv("exported/exported_cd3/260411/260411_cluster_results/Pheno_Cluster_7.csv")
setwd("..")
setwd("..")
setwd("..")
setwd("exported_gdTcells/260406/260406_cluster_results")
gd_clusters <- tibble()
sel_gd <- c(6,10)
for (i in 1:2){
  cluster_name <- paste("Pheno_Cluster_",sel_gd[i],".csv",sep="")
  clust <- read_csv(cluster_name)
  colnames(clust)[5] <- paste("Cluster_gd_",sel_gd[i],sep="")
  clusters[[i+5]] <- clust
  if (i == 1){
    gd_clusters <- clust
  }
  else{gd_clusters <- full_join(gd_clusters,clust[,-c(1,2,4)],by="sample_var")}
}


#### old 260405 data ####
# 12, 15 ## CD25+ CD25hi

# 45 tcrgd-
# 47 tcrgd-cd25+
# 73, 75, 77
# 85, 87
### 1,12,15,21,42,59,45,73,75,77,85,84


#### new 260415 after data ###
# 14, 17, 20 cluster10
# 26, 29 cluster 7
# 51 cluster 12
# 54, 57, 48 cluster 3
############
#yd in T cells clust
# 100 cluster 11
# 92,94, 96 cluster 15
# 98 cluster 15
############
#yd in yd
# 108, 110 in cluster10!
# 118, 120, 104 in cluster6
# 14,17,20,26,29,51,54,57,48,99,92,94,96,98,108,110,118,120,104
# 10,10,10,7,7,12,3,3,3,11,15,15,15,15,


setwd("..")
setwd("..")
setwd("..")
setwd("..")
manual_gating <- read_csv("results/260415_all_gated_extra_dia_motol_after.csv")
for (i in 1:ncol(manual_gating)) {
  print(paste(i,substr(colnames(manual_gating)[i],66,200)))
}
manual_gating_sel <- tibble(manual_gating[c(1,14,17,20,26,29,51,54,57,48,100,92,94,96,98,107,109,117,119,103)])
colnames(manual_gating_sel) <- c("filename","cluster10_1","cluster10_2","cluster10_3","cluster7_1","cluster7_2",
                                 "cluster12","cluster3_1","cluster3_2","cluster3_3","cluster11","cluster15_1",
                                 "cluster15_2","cluster15_3","cluster15_4","cluster10gd_1","cluster10gd_2",
                                 "cluster6gd_1","cluster6gd_2","cluster6gd_3")
#manual_gating <- read_csv("results/260405_all_gated_extra_dia_motol1.csv")
#colnames(manual_gating)[21]
#manual_gating_cl7 <- tibble(manual_gating[,c(1,21)])
#colnames(manual_gating_cl7)[1] <- "filename"
cluster <- all_selected_clusters %>%
  mutate(filename = substr(all_selected_clusters$sample_var,18,22)) 
cluster$filename <- cluster$filename %>%
  str_remove("_") %>%
  str_remove(" ") %>%
  str_remove("C")
cluster$filename
data <- full_join(cluster,manual_gating_sel,"filename")
head(data)
#colnames(data)[7] = "freq_gate"
col0 <- c("#7CD6FC","#fa7065","#be2115")

ggtheme <- function() {
  theme(
    axis.text.x = element_text(size = 20,color = "#000"),
    axis.text.y = element_text(size = 20,color = "#000"),
    #axis.title = element_text(size = 17),
    #plot.title = element_text(size = 24),
    text = element_text(size = 20, colour = "black", family = "Arial"),
    axis.title = element_text(size = 20))
  #legend.text = element_text(size = 20),
  #legend.key.size = unit(10, units = "points"))
  
} 
# 10,10,10,7,7,12,3,3,3,11,15,15,15,15,
cols1 <- c(5,5,5,6,6,7,8,8,8,10,9,9,9,9)
cols2 <- c(12:30)
setwd("results/gating_plots/260415/correlations")
plotting <- function(freq_clust,freq_gate,data,clust, gate,k){
  ggplot(data = data, aes(x = freq_clust, y = freq_gate, color = group_var, label = filename)) +
  geom_point(size = 2) +
 # geom_text(nudge_y = 1, size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "grey40", linetype = "dashed") +
  stat_cor(aes(group = 1),method = "pearson", label.x.npc = "left", label.y.npc = "top") +  # requires ggpubr
  labs(
    x = paste("pheno ",colnames(data)[clust],sep=""),
    y = paste("gated ",colnames(data)[gate],sep=""),
   # title = "Phenograph cluster vs. manual gate frequency per sample",
    color = "Group"
  ) +
  theme_bw()+
  scale_color_manual(values=col0)+
  scale_x_continuous(limits = c(0,NA))+
  scale_y_continuous(limits = c(0,NA))+
  ggtheme()#+
  #geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed")
  ggsave(paste(k,"_correlation_manual_gating_phenograph.svg",sep=""),width = 16,height = 8,units = "cm")
}
for (j in 1:length(cols1)){
  plotting(data[[cols1[j]]],data[[cols2[j]]],data,cols1[j],cols2[j],j)
}

getwd()


cluster_gd <- gd_clusters %>%
  mutate(filename = substr(gd_clusters$sample_var,18,22)) 
cluster_gd$filename <- cluster_gd$filename %>%
  str_remove("_") %>%
  str_remove(" ") %>%
  str_remove("C")
cluster_gd$filename
data_gd <- full_join(cluster_gd,manual_gating_sel[c(1,16:20)],"filename")
colsgd <- c(6,6,5,5,5)
for (gd in 1:5){
  plotting(data_gd[[colsgd[gd]]],data_gd[[gd+7]],data_gd,colsgd[gd],gd+7,paste("gd",gd,sep = ""))
}

setwd("..")
setwd("..")
setwd("..")
setwd("..")
setwd("..")
getwd()
setwd("intracellular_pannels/dia_intra01/exported/exported_tregs/260411/260411_cluster_results")
# clust 6 =   2_27
# clust 5 = 2_35
# clust 9 = 2_38
# clust 10 = 2_39

#selected <- c(6,5,9,10)
selected <- 4
tregs_clusters <- tibble()
for (i in 1:length(selected)){
  print(i)
  print(selected[i])
  cluster_name <- paste("Pheno_merged6_Cluster_",selected[i],".csv",sep="")
  clust <- read_csv(cluster_name)
  colnames(clust)[5] <- paste("Cluster_",selected[i],sep="")
  clusters[[i]] <- clust
  #colnames(clust)[5] <- paste("Cluster_",selected[i],sep="")
  if (i == 1){
    tregs_clusters <- clust
  }
  else{
    tregs_clusters <- full_join(tregs_clusters,clust[,-c(1,2,4)],by="sample_var")
  }
  #full_join(all_selected_clusters,)
}
setwd("..")
setwd("..")
setwd("..")
setwd("..")
#setwd("results/gating_plots/260405/")
getwd()
manual_tregs <- read_csv("260404_intra_tregs_dia_motol.csv")
manual_tregs_sel <- manual_tregs[c(1,27)]
colnames(manual_tregs_sel)
colnames(manual_tregs_sel) <- c("filename","cluster4")
tregs_clusters <- tregs_clusters %>%
  mutate(filename = substr(gd_clusters$sample_var,18,22)) 
tregs_clusters$filename <- tregs_clusters$filename %>%
  str_remove("_") %>%
  str_remove(" ") %>%
  str_remove("C")
tregs_clusters$filename
colnames(tregs_clusters)
data_tregs <- full_join(tregs_clusters,manual_tregs_sel,"filename")
plotting(data_tregs[[5]],data_tregs[[10]],data_tregs,5,10,paste("tregs",gd,sep = ""))


#   geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed")
plotting(data[[cols1[7]]],data[[cols2[7]]],data,cols1[7],cols2[7],77)
ggplot(data = data, aes(x =data[[cols1[7]]], y = data[[cols2[7]]], color = group_var, label = filename)) +
  geom_point(size = 2) +
  # geom_text(nudge_y = 1, size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "grey40", linetype = "dashed") +
  stat_cor(aes(group = 1),method = "pearson", label.x.npc = "left", label.y.npc = "top") +  # requires ggpubr
  labs(
    x = paste("pheno ",colnames(data)[cols1[7]],sep=""),
    y = paste("gated ",colnames(data)[cols2[7]],sep=""),
    # title = "Phenograph cluster vs. manual gate frequency per sample",
    color = "Group"
  ) +
  theme_bw()+
  scale_color_manual(values=col0)+
  scale_x_continuous(limits = c(0,NA))+
  scale_y_continuous(limits = c(0,NA))+
  ggtheme()#+
  #geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed")
ggsave(paste(77,"_correlation_manual_gating_phenograph.svg",sep=""),width = 16,height = 8,units = "cm")
getwd()


plotting(data_tregs[[5]],data_tregs[[10]],data_tregs,5,10,paste("tregs",gd,sep = ""))
ggplot(data = data_tregs, aes(x = data_tregs[[5]], y = data_tregs[[7]], color = group_var, label = filename)) +
  geom_point(size = 2) +
  # geom_text(nudge_y = 1, size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "grey40", linetype = "dashed") +
  stat_cor(aes(group = 1),method = "pearson", label.x.npc = "left", label.y.npc = "top") +  # requires ggpubr
  labs(
    x = paste("pheno ",colnames(data_tregs)[5],sep=""),
    y = paste("gated ",colnames(data_tregs)[7],sep=""),
    # title = "Phenograph cluster vs. manual gate frequency per sample",
    color = "Group"
  ) +
  theme_bw()+
  scale_color_manual(values=col0)+
  scale_x_continuous(limits = c(0,NA))+
  scale_y_continuous(limits = c(0,NA))+
  ggtheme()#+
#geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed")
ggsave(paste(66,"_correlation_manual_gating_phenograph.svg",sep=""),width = 16,height = 8,units = "cm")
