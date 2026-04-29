library(cyCONDOR)
library(ggplot2)
library(ggrastr)
library(cowplot)
library(readxl)
library(tidyverse)
library(flowCore)
library(RColorBrewer)
library(pheatmap)

##### loading the data #########################################################
getwd()
directory_default = "/home/charvatb/PRIMUS/data/48_lab/Project Diabetes Motol/FACS validation/240311 Analysis BCH"
setwd(directory_default)
panel = readline(prompt = "Choose panel [extra/intra01/intra02]: ")
extra
if (panel == "extra" ) {
  directory = paste(directory_default,"/extracellular_pannel/exported/",sep = "")
} else if (panel == "intra01") {
  directory = paste(directory_default,"/intracellular_pannels/dia_intra01/exported/",sep = "")
} else if (panel == "intra02") {
  directory = paste(directory_default,"/intracellular_pannels/dia_intra02/exported/",sep = "")
}
#directory = directory_default
folder_name = readline(prompt = "Name of the folder: ")
exported_cd3
#exported_gdTcells
#exported_cd3 #exported_tregs
directory = paste(directory,folder_name,sep="")
setwd(directory)
getwd()
if (!dir.exists(format(Sys.Date(), "%y%m%d"))) {
  dir.create(format(Sys.Date(), "%y%m%d"), recursive = TRUE)
}
dir_clust = paste(directory,format(Sys.Date(), "%y%m%d"),sep="/")


use_default = readline(prompt = "Use default patients data? [y/n]: ");
y
if (use_default == "y"){
  data_path = "/home/charvatb/PRIMUS/data/48_lab/Project Diabetes Motol/FACS validation/240311 Analysis BCH/patients_data/241114_patients.xlsx"
} else {data_path =  readline(prompt = "Insert path to your data: ");}


#my data "/home/charvatb/PRIMUS/data/48_lab/Project Diabetes Motol/FACS validation/240311 Analysis BCH/patients_data/241114_patients.xlsx"
# change the adress to a place, where are your files
#directory <- paste(getwd(),"/extracellular_pannel/exported/exported_cd3/241122/", sep = "")

#data1 <- readxl::read_excel(paste(getwd(),"patients_data/241114_patients.xlsx",sep = "/"))
data1 = readxl::read_excel(data_path)
# editing the table for further use
# adding columns disease, timepoint, group etc.
data1 <- data1 %>%
  mutate(Disease = ifelse(substr(SampleID,1,1)==2,"0","1")) %>%
  mutate(Disease = as.numeric(Disease)) %>%
  mutate(Timepoint = ifelse(substr(SampleID,4,4) == "b","1","0")) %>%
  mutate(Timepoint = as.numeric(Timepoint)) %>%
  mutate(Group = ifelse(Disease == 1,ifelse(Timepoint == 1,"DiaT1","DiaT0"),"Control")) %>%
  mutate(sample_number = substr(SampleID,1,3)) %>%
  mutate(sample_number = as.numeric(sample_number)) %>%
  mutate(filename = `Sample:`)

fcs_files <- list.files(directory, pattern = "\\.fcs$", full.names = FALSE)
samplename_length = nchar(data1$`Sample:`[1])-4

##### EDIT !!!!!! ####
file_name_start = substr(fcs_files[1],1,14)
file_name_end = substr(fcs_files[1],15+samplename_length,nchar(fcs_files[1]))
#the column filename has to fit perfectly with the name, of the .fcs file you want
# to load, otherwise the cyCONDOR object will have 0 columns
data1$filename <- data1$filename %>%
  str_replace(".fcs",file_name_end)
data1 <- data1 %>%
  mutate(filename = paste(file_name_start,filename,sep=""))

# Assuming 'my_data.fcs' is your FCS file
setwd(directory)
flow_object <- read.FCS(fcs_files[1],truncate_max_range = FALSE) # Example keyword for parameter name

# Create a data frame to see both together
channel_info <- data.frame(
  channel = flow_object@parameters@data$name,
  marker = flow_object@parameters@data$desc
)
print("These are channels present in your panel")
print(channel_info)
print("These channels will be removed from the clustering automatically")
remove_channels = c("FSC-H", "FSC-A", "SSC-H", "SSC-A", "SSC-B-H","SSC-B-A", "Time","Via","FlowAI")
print(remove_channels)

################################################################################
library(flowCore)

# Read one of your FCS files directly
fcs <- flow_object

# CD25 raw distribution
hist(exprs(flow_object))
for (chan in 1:10){
  png(paste("hist",chan,".png",sep=""))
  h <- hist(exprs(flow_object)[,channel_info$channel[chan]],breaks=100,main = channel_info$marker[chan])
  dev.off()
}



################################################################################

remove_final = c(remove_channels,c("CD3",#,
                                   #"CD4",
                                   #"CD8", 
                                   "B220",
                                   # extra
                                   "iNKT",
                                   #"CD25",
                                   "CD16",
                                   #"TCRgd",
                                   "HLA"
                                   # intra01
                                   #"CD8","B220","CD4",
                                   #"Foxp3","CD127"
))

# not necessary to do it all the time, just for the first time with new 
# data it is crucial

write.csv(data1, file.path(directory,"patients_path.csv"))


# debug
# Verify files are where you think they are
# print(directory)
# print(list.files(directory, pattern = "\\.fcs$"))
# 
# Check if the filenames in your annotation match actual files
files_exist <- data1$filename %in% list.files(directory, pattern = "\\.fcs$")
print(data1$filename[!files_exist])  # Show missing files
# Print actual vs expected filenames
print("Expected filenames:")
print(data1$filename)
print("Actual files:")
print(list.files(directory, pattern = "\\.fcs$"))
he = list.files(directory, pattern = "\\.fcs$")
if (data1$filename[1] == he[6]){
  print("huuuu")
}
print(data1$filename[1])
print(he[6])
################################################################################
####### Creating condor object  ################################################
condor <- prep_fcd(data_path = directory, 
                   max_cell = 1000, 
                   useCSV = FALSE, 
                   transformation = "auto_logi", 
                   remove_param = c(remove_final), 
                   #anno_table = paste(directory,"patients_path.csv",sep=""), 
                   anno_table = file.path(directory, "patients_path.csv"),
                   filename_col = "filename",
                   truncate_max_range = FALSE,
                   ignore.text.offset = TRUE
)


# library(flowCore)
# lgcl <- logicleTransform(w = 0.5, t = 262144, m = 4.5, a = 0)
#condor$expr$orig$CD25 <- lgcl(condor$expr$orig$CD25)

hist
for (chan in 1:10){
  png(paste("condor_hist",chan,".png",sep=""))
  print(colnames(condor$expr$orig[chan]))
  h <- hist(condor$expr$orig[chan],breaks=100,main =colnames(condor$expr$orig)[chan])
  dev.off()
}

#?prep_fcd

class(condor)

condor <- runPCA(fcd = condor, 
                 data_slot = "orig")


condor <- runUMAP(fcd = condor, 
                  input_type = "pca", 
                  data_slot = "orig"
)

plot_dim_red(fcd= condor,  
             reduction_method = "umap", 
             reduction_slot = "pca_orig", 
             param = "Group", 
             title = "UMAP colored by group"
)
ggsave(filename = "plot.svg",path = dir_clust,width=12,height = 12,units = "cm")


################################################################################
#### BATCH EFFECT ##############################################################

# Colored by run, displaying the batch effect
plot_dim_red(fcd= condor,  
             reduction_method = "umap", 
             reduction_slot = "pca_orig", 
             param = "Run_no", 
             title = "UMAP before edit"
)
ggsave(filename = "plot1.svg",path = dir_clust,width=12,height = 12,units = "cm")

# Deleting batch effect
condor <- harmonize_PCA(fcd = condor, 
                        batch_var = c("Run_no"), 
                        data_slot = "orig")

#condor$pca$norm[1:10, 1:5]

condor <- runUMAP(fcd = condor, 
                  input_type = "pca", 
                  data_slot = "norm",
                  prefix= NULL)


plot_dim_red(fcd= condor,  
             reduction_method = "umap", 
             reduction_slot = "pca_norm", 
             param = "Run_no", 
             title = "Harmonized UMAP")
ggsave(filename = "plot2.svg",path = dir_clust,width=12,height = 12,units = "cm")

plot_dim_red(fcd= condor,  
             reduction_method = "umap", 
             reduction_slot = "pca_norm", 
             param = "Group", 
             title = "Harmonized UMAP")
ggsave(filename = "plot3.svg",path = dir_clust,width=12,height = 12,units = "cm")


# split UMAP plots (showing groups from metadata, under param)
plot_dim_red(fcd = condor,
             expr_slot = NULL,
             reduction_method = "umap",
             reduction_slot = "pca_norm",
             cluster_slot = NULL,
             param = "Group",
             facet_by_variable = T,
             title = "UMAP")
ggsave(filename = "plot4.svg",path = dir_clust,width=16,height = 8,units = "cm")



################################################################################
### CLUSTERING #################################################################
setwd(dir_clust)
getwd()
condor <- runPhenograph(fcd = condor, 
                        input_type = "pca", 
                        data_slot = "norm", 
                        k = 60, 
                        seed = 91)
n_clusters_p <- nlevels(condor[["clustering"]][["phenograph_pca_norm_k_60"]][["Phenograph"]])
FlowSOM_n_clusters = readline(prompt = paste("Phenograph found",n_clusters_p,"clusters. How many clusters do you want when running FlowSOM? "));
6
FlowSOM_n_clusters = as.integer(FlowSOM_n_clusters)

condor <- runFlowSOM(fcd = condor, 
                     input_type = "pca", 
                     data_slot = "norm", 
                     nClusters = FlowSOM_n_clusters, 
                     seed = 91, 
                     ret_model = TRUE)

FlowSOM_name = paste("FlowSOM_pca_norm_k_",FlowSOM_n_clusters,sep="")

plot_dim_red(fcd = condor,
             expr_slot = NULL,
             reduction_method = "pca",
             reduction_slot = "orig",
             cluster_slot = NULL,
             param = "Group",
             title = "PCA")

plot_dim_red(fcd = condor,
             expr_slot = NULL,
             reduction_method = "umap",
             reduction_slot = "pca_norm",
             cluster_slot = NULL,
             param = "Group",
             title = "UMAP")




palette <- c(
  "#67160F",
  "#C0392B",
  "#D68342",
  "#EDDF51",
  "#DEC4B1",
  "#7DBB3A",
  "#1E8449",
  "#98EBD6",
  "#67B0F4",
  "#1A5276",
  "#69215B",
  "#DB6EC6",
  "#DCC3F3",
  "#7D4A11",
  "#696664",
  "#0D0802"
)



plot_dim_red(fcd = condor,
             expr_slot = NULL,
             reduction_method = "umap",
             reduction_slot = "pca_norm",
             cluster_slot = "phenograph_pca_norm_k_60",
             param = "Phenograph",
             title = "UMAP Phenograph")+
  #scale_color_manual(values = rainbow(15))
  scale_color_manual(values = palette)
ggsave(filename = "plot5.svg",width=12,height = 12,units = "cm")

plot_dim_red(fcd = condor,
             expr_slot = NULL,
             reduction_method = "umap",
             reduction_slot = "pca_norm",
             cluster_slot = FlowSOM_name,
             param = "FlowSOM",
             title = "UMAP FlowSOM")
ggsave(filename = "plot6.svg",path = dir_clust,width=12,height = 12,units = "cm")



condor <- metaclustering(fcd = condor,
                        cluster_slot = "phenograph_pca_norm_k_60",
                        cluster_var = "Phenograph",
                        cluster_var_new = "metaclusters",
                        metaclusters = c("1" = "CD4+ AgExp 01",
                                         "2" = "CD4+ AgExp 02",
                                         "3" = "CD8+ Naive 01",
                                         "4" = "CD8+ Naive 02",
                                         "5" = "CD8+ AgExp 01",
                                         "6" = "CD4+ Naive 01",
                                         "7" = "CD4+ Naive 02",
                                         "8" = "γδ 01",
                                         "9" = "CD8+ Naive 03",
                                         "10" = "CD4+ Naive 03",
                                         "11" = "γδ 02",
                                         "12" = "CD8+ Naive 04",
                                         "13" = "CD8+ AgExp 02",
                                         "14" = "γδ low",
                                         "15"="γδ 03"
                        ))

plot_dim_red(fcd = condor,
             expr_slot = NULL,
             reduction_method = "umap",
             reduction_slot = "pca_norm",
             cluster_slot = "phenograph_pca_norm_k_60",
             #cluster_var = "metaclusters",
             param = "metaclusters",
             title = "UMAP Phenograph")+
  #scale_color_manual(values = rainbow(15))
  scale_color_manual(values = palette)
ggsave(filename = "plot7.svg",width=12,height = 12,units = "cm")


condor$clustering$phenograph_pca_norm_k_60$metaclusters <- factor(
  condor$clustering$phenograph_pca_norm_k_60$metaclusters,
  levels = sort(unique(as.character(
    condor$clustering$phenograph_pca_norm_k_60$metaclusters
  )))
)
levels(condor$clustering$phenograph_pca_norm_k_60$metaclusters)


plot_dim_red(fcd = condor,
             expr_slot = NULL,
             reduction_method = "umap",
             reduction_slot = "pca_norm",
             cluster_slot = "phenograph_pca_norm_k_60",
             #cluster_var = "metaclusters",
             param = "metaclusters",
             title = "UMAP Phenograph")+
  #scale_color_manual(values = rainbow(15))
  scale_color_manual(values = palette)
ggsave(filename = "plot8.svg",width=12,height = 12,units = "cm")


palette1 <- c(
  "#BFBABA",
  "#BFBABA",
  "#BFBABA",
  "#BFBABA",
  "#BFBABA",
  "#BFBABA",
  "#BFBABA",
  "#BFBABA",
  "#BFBABA",
  "#BFBABA",
  "#BFBABA",
  "#DB6EC6",
  "#DCC3F3",
  "#7D4A11",
  "#BFBABA",
  "#0D0802"
)


plot_dim_red(fcd = condor,
             expr_slot = NULL,
             reduction_method = "umap",
             reduction_slot = "pca_norm",
             cluster_slot = "phenograph_pca_norm_k_60",
             #cluster_var = "metaclusters",
             param = "metaclusters",
             title = "UMAP Phenograph")+
  #scale_color_manual(values = rainbow(15))
  scale_color_manual(values = palette1)
ggsave(filename = "plot9.svg",width=12,height = 12,units = "cm")

# 3 + 7 + 10
palette2 <- c(
  "#BFBABA",
  "#BFBABA",
  "#BFBABA",
  "#EDDF51",
  "#DEC4B1",
  "#BFBABA",
  "#BFBABA",
  "#98EBD6",
  "#BFBABA","#BFBABA","#BFBABA","#BFBABA","#BFBABA","#BFBABA","#BFBABA"
)


plot_dim_red(fcd = condor,
             expr_slot = NULL,
             reduction_method = "umap",
             reduction_slot = "pca_norm",
             cluster_slot = "phenograph_pca_norm_k_60",
             #cluster_var = "metaclusters",
             param = "metaclusters",
             title = "UMAP Phenograph")+
  #scale_color_manual(values = rainbow(15))
  scale_color_manual(values = palette2)
ggsave(filename = "plot10.svg",width=12,height = 12,units = "cm")

palette3 <- c(
  "#BFBABA",
  "#BFBABA",
  "#BFBABA",
  "#BFBABA",
  "#BFBABA",
  "#BFBABA",
  "#BFBABA",
  "#BFBABA",
  "#BFBABA",
  "#BFBABA",
  "#BFBABA",
  "#BFBABA",
  #"#DB6EC6",
  "#DCC3F3",
  "#7D4A11",
  "#BFBABA",
  "#0D0802"
)

plot_dim_red(fcd = condor,
             expr_slot = NULL,
             reduction_method = "umap",
             reduction_slot = "pca_norm",
             cluster_slot = "phenograph_pca_norm_k_60",
             #cluster_var = "metaclusters",
             param = "metaclusters",
             title = "UMAP Phenograph")+
  #scale_color_manual(values = rainbow(15))
  scale_color_manual(values = palette3)
ggsave(filename = "plot11.svg",width=12,height = 12,units = "cm")



palette4 <- c(
  "#BFBABA",
  "#BFBABA",
  "#BFBABA",
  "#BFBABA",
  "#BFBABA",
  "#BFBABA",
  "#BFBABA",
  "#BFBABA",
  "#BFBABA",
  "#BFBABA",
  "#69215B",
  "#BFBABA",
  "#BFBABA",
  "#BFBABA",
  "#BFBABA"
)


plot_dim_red(fcd = condor,
             expr_slot = NULL,
             reduction_method = "umap",
             reduction_slot = "pca_norm",
             cluster_slot = "phenograph_pca_norm_k_60",
             #cluster_var = "metaclusters",
             param = "metaclusters",
             title = "UMAP Phenograph")+
  #scale_color_manual(values = rainbow(15))
  scale_color_manual(values = palette4)
ggsave(filename = "plot12.svg",width=12,height = 12,units = "cm")


palette5 <- c(
  "#BFBABA",
  "#BFBABA",
  "#BFBABA",
  "#EDDF51",
  "#DEC4B1",
  "#BFBABA",
  "#BFBABA",
  "#BFBABA",
  #"#98EBD6",
  "#BFBABA",
  "#BFBABA",
  "#69215B",
  "#BFBABA",
  "#DCC3F3",
  "#7D4A11",
  "#BFBABA"
)

plot_dim_red(fcd = condor,
             expr_slot = NULL,
             reduction_method = "umap",
             reduction_slot = "pca_norm",
             cluster_slot = "phenograph_pca_norm_k_60",
             #cluster_var = "metaclusters",
             param = "metaclusters",
             title = "UMAP Phenograph")+
  #scale_color_manual(values = rainbow(15))
  scale_color_manual(values = palette5)
ggsave(filename = "plot13.svg",width=12,height = 12,units = "cm")


palette6 <- c(
  "#BFBABA",
  "#BFBABA",
  "#BFBABA",
  "#EDDF51",
  "#DEC4B1",
  "#BFBABA",
  "#BFBABA",
  "#BFBABA",
  #"#98EBD6",
  "#BFBABA",
  "#BFBABA",
  "#69215B",
  "#BFBABA",
  "#BFBABA",
  "#BFBABA",
  "#BFBABA",
  "#DCC3F3",
  "#7D4A11",
  "#BFBABA"
)

plot_dim_red(fcd = condor,
             expr_slot = NULL,
             reduction_method = "umap",
             reduction_slot = "pca_norm",
             cluster_slot = "phenograph_pca_norm_k_60",
             #cluster_var = "metaclusters",
             param = "metaclusters",
             title = "UMAP Phenograph")+
  #scale_color_manual(values = rainbow(15))
  scale_color_manual(values = palette6)
ggsave(filename = "plot15.svg",width=12,height = 12,units = "cm")

#condor <- metaclustering(fcd = condor, 
#                         cluster_slot = FlowSOM_name, 
#                        cluster_var = "FlowSOM")
#                          cluster_var_new = "metaclusters", 
#                          metaclusters = c("1" = "i",
#                                           "2" = "ii",
#                                           "3" = "iii", 
#                                           "4" = "iv",
#                                           "5" = "v",
#                                           "6" = "vi",
#                                           "7" = "vii", 
#                                           "8" = "viii",
#                                           "9" = "ix",
#                                           "10" = "x",
#                                           "11" = "xi", 
#                                           "12" = "xii", 
#                                           "13" = "xiii",
#                                           "14" = "xiv",
#                                           "15" = "xv"))

clustersExpression <- function(marker,clusterMethod,cluster_var) {
  plot_marker_boxplot(fcd = condor,
                      expr_slot = "orig", 
                      cluster_slot = clusterMethod, 
                      cluster_var = cluster_var,
                      group_var = "Group",
                      sample_var = "SampleID", 
                      marker = marker,
                      fun = "mean",
                      facet_ncol = 3,
                      dot_size = 0.5)
  ggsave(filename = paste(j,"_clusters_",marker,"_expr.svg",sep=""),path = dir_clust,width=16,height = 8,units = "cm")
}

clusteringMethods <- c("phenograph_pca_norm_k_60",FlowSOM_name)
clusterVar <- c("Phenograph","FlowSOM")
for (j in 1:2) {
  for (i in 1:ncol(condor$expr$orig)) {
    name <- colnames(condor$expr$orig)[i]
    #clusteringMethod <- clusteringMethods[j]
    #if (name == "FlowAI") { next }
    clustersExpression(name,clusteringMethods[j],clusterVar[j])
  } 
}

################################################################################
### PLOT EXPRESSION OF EACH MARKER IN UMAP #####################################

plot_expression <- function(marker) {
  plot_dim_red(fcd = condor,
               expr_slot = "orig",
               reduction_method = "umap",
               reduction_slot = "pca_norm",
               cluster_slot = NULL,
               param = marker, 
               order = T,
               title = paste(marker,"expression"))
  ggsave(filename = paste("plot_expr_",i,".svg",sep=""),path = dir_clust,width=12,height = 12,units = "cm")
}

## for cycle across all the markers present in the population ##################
for (i in 1:ncol(condor$expr$orig)) {
  name <- colnames(condor$expr$orig)[i]
  plot_expression(name)
} 

plot_marker_HM(fcd = condor,
               expr_slot = "orig",
               marker_to_exclude = c("FSC-A","SSC-A"),
               cluster_slot = "phenograph_pca_norm_k_60",
               cluster_var = "Phenograph")

#ranran = data.frame(condor[["clustering"]][["phenograph_pca_norm_k_60"]])
#heatmap(ranran)
### heatmaps ###
#library(pheatmap)
library(dplyr)

# Extract expression data and cluster assignments
expr_data <- condor$expr$orig
cluster_data_pheno <- condor$clustering$phenograph_pca_norm_k_60$Phenograph
cluster_data_flowSOM <- condor$clustering[[2]]$FlowSOM
# 
# create_heatmap <- function(expr_data,clust_data,name_clust){
#   # Calculate median expression per cluster for each marker
#   cluster_medians <- expr_data %>%
#     as.data.frame() %>%
#     mutate(cluster = clust_data) %>%
#     group_by(cluster) %>%
#     summarise(across(everything(), median)) %>%
#     as.data.frame()
#   # Set cluster as rownames
#   rownames(cluster_medians) <- cluster_medians$cluster
#   cluster_medians$cluster <- NULL
#   # creating heatmap
#   svg(paste("04_heatmap_",name_clust,".svg",sep=""), width = 800, height = 600, res = 100)
#   coul <- colorRampPalette(brewer.pal(9, "YlOrRd"))(5)
#   heatmap_data <- t(sapply(cluster_medians, as.numeric ))
#   h = heatmap(heatmap_data, Colv = NA, Rowv = NA, scale="row", col = coul)
#   dev.off()
# }


expr_data <- condor$expr$orig  # or condor_object$raw.data

# Get cluster assignments from phenograph results
# If stored in condor object: condor_object$clustering$phenograph
# generate_heatmaps <- function(pheno, name){
#   #clusters <- tibble(condor$clustering$phenograph)  # Adjust path as needed
#   clusters <- tibble(condor$clustering[[pheno]])
#   print(clusters)
#   clusters$Phenograph <- as.numeric(clusters$Phenograph)
#   # Calculate mean expression for each cluster
#   cluster_ids <- unique(clusters$Phenograph)
#   n_clusters <- length(cluster_ids)
#   n_markers <- ncol(expr_data)
#   
#   # Create matrix of cluster centroids (mean expression per cluster)
#   cluster_centroids <- matrix(0, nrow = n_clusters, ncol = n_markers)
#   rownames(cluster_centroids) <- paste0("Cluster_", cluster_ids)
#   colnames(cluster_centroids) <- colnames(expr_data)
#   
#   for (i in seq_along(cluster_ids)) {
#     cluster_cells <- clusters == cluster_ids[i]
#     cluster_centroids[i, ] <- colMeans(expr_data[cluster_cells, , drop = FALSE])
#   }
#   
#   # Calculate distance matrix between clusters
#   # Using Euclidean distance by default
#   
#   cluster_dist <- as.matrix(dist(cluster_centroids, method = "euclidean"))
#   rownames(cluster_dist) <- paste0("Cluster_", cluster_ids)
#   colnames(cluster_dist) <- paste0("Cluster_", cluster_ids)
#   
#   # Optional: Try different distance metrics
#   # Correlation distance
#   cluster_dist_cor <- as.matrix(as.dist(1 - cor(t(cluster_centroids))))
#   rownames(cluster_dist_cor) <- paste0("Cluster_", cluster_ids)
#   colnames(cluster_dist_cor) <- paste0("Cluster_", cluster_ids)
#   
#   # Create heatmap with clusters ordered by similarity
#   # pheatmap will automatically cluster columns (clusters) by distance
#   p2 <- pheatmap(
#     cluster_centroids,
#     color = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100),
#     scale = "column",           # Scale markers (rows) to z-scores
#     cluster_rows = TRUE,     # Cluster markers
#     cluster_cols = TRUE,     # Cluster columns (this orders clusters by distance)
#     clustering_distance_cols = "euclidean",  # Distance metric for clusters
#     clustering_method = "complete",          # Hierarchical clustering method
#     show_colnames = TRUE,
#     show_rownames = TRUE,
#     angle_col = 45, 
#     main = "Marker Expression Across Phenograph Clusters",
#     fontsize_row = 20,
#     fontsize_col = 20,
#     border_color = NA
#   )
#   name_plot <- paste("heatmap",name,"2.png",sep="")
#   ggsave(name_plot,plot = p2,width=24,height = 24, units = "cm")
# }
# 
# generate_heatmaps("phenograph_pca_norm_k_60","_og_")



generate_heatmaps <- function(pheno, name){
  clusters <- tibble(condor$clustering[[pheno]])
  
  # use metaclusters instead of Phenograph
  cluster_ids <- unique(clusters$metaclusters)
  n_clusters <- length(cluster_ids)
  n_markers <- ncol(expr_data)
  
  cluster_centroids <- matrix(0, nrow = n_clusters, ncol = n_markers)
  rownames(cluster_centroids) <- as.character(cluster_ids)  # metacluster names as row labels
  colnames(cluster_centroids) <- colnames(expr_data)
  
  for (i in seq_along(cluster_ids)) {
    cluster_cells <- clusters$metaclusters == cluster_ids[i]  # match on metaclusters column
    cluster_centroids[i, ] <- colMeans(expr_data[cluster_cells, , drop = FALSE])
  }
  
  p2 <- pheatmap(
    cluster_centroids,
    color = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100),
    scale = "column",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    show_colnames = TRUE,
    show_rownames = TRUE,
    angle_col = 45,
    main = "Marker Expression Across Phenograph Clusters",
    fontsize_row = 20,
    fontsize_col = 20,
    border_color = NA
  )
  name_plot <- paste("heatmap", name, "2.svg", sep="")
  ggsave(name_plot, plot = p2, width=24, height=24, units="cm")
}

generate_heatmaps("phenograph_pca_norm_k_60","_og_")




plot_dim_red(fcd = condor,
             expr_slot = NULL,
             reduction_method = "umap",
             reduction_slot = "pca_norm",
             cluster_slot = NULL,
             param = "Group",
             facet_by_variable = T,
             title = "UMAP")




plot_frequency_barplot(fcd = condor,
                       cluster_slot = FlowSOM_name,
                       cluster_var = "FlowSOM",
                       group_var = "Group",
                       title = "Stacked barplot of cluster frequencies")

plots <- plot_frequency_boxplot(fcd = condor,
                                cluster_slot = FlowSOM_name, 
                                cluster_var = "FlowSOM",
                                sample_var = "expfcs_filename", 
                                group_var = "Group", 
                                numeric = T)

plot_grid(plotlist = plots, ncol = 3)


################################################################################
##### Extracting table for each cluster ########################################
##### Creating plots for each cluster ##########################################


##### FLOWSOM CLUSTERING PLOTS ##############################################
clusters_table <- plots[["Cluster_1"]][["data"]]
save_clusters <- paste(dir_clust,"/",format(Sys.Date(), "%y%m%d"),"_cluster_results", sep="")
if (!dir.exists(save_clusters)) {
  dir.create(save_clusters, recursive = TRUE)
}
setwd(save_clusters)

clusters_list_of_tables1 <- list()
for (i in 1:FlowSOM_n_clusters) {
  name <- paste("Cluster_",i,sep="")
  newCluster <- plots[[name]][["data"]]
  clusters_list_of_tables1[[i]] <- newCluster
  filename <- paste("FlowSOM_",name,".csv",sep="")
  write.csv(newCluster,file.path(save_clusters,filename))
}


comparisons <- list(c("Control", "DiaT0"),c("DiaT0","DiaT1"), c("Control", "DiaT1"))
expandingAxis <- function(maxValue) {
  expand <- 0
  ifelse(maxValue > 50,
         ifelse(maxValue >80,expand <- 15,expand <- 10),
         ifelse(maxValue >10,expand <- 5,expand <- 3))
  return(expand)
}
col0 <- c("#7CD6FC","#fa7065","#be2115")


ggtheme <- function() {
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 30),
    #axis.title = element_text(size = 17),
    #plot.title = element_text(size = 24),
    text = element_text(size = 24, colour = "black", family = "Arial"),
    axis.title = element_blank())
  #legend.text = element_text(size = 20),
  #legend.key.size = unit(10, units = "points"))
  
} 

plotting <- function(population,data,name,i,l) 
{
  
  ggplot(aes(x = group_var, y = population), data=data) +
    geom_boxplot(alpha = 0.2,aes(fill = group_var)) +
    geom_dotplot(binaxis='y',
                 dotsize = 1.5,
                 stackdir='center',
                 aes(fill = group_var)) +
    
    labs(title = name) +
    theme_classic() +
    ggpubr::stat_compare_means(label.x= 1.5,size = 10, vjust = -0.05, 
                               label = "p.format",comparisons = comparisons)  +
    #ylab("Frequency") +
    #xlab("Condition") +
    ggtheme() +
    scale_y_continuous(limits = c(0,NA), expand = expansion(mult = c(0.05, 0.25))) +
    #scale_y_continuous(limits = c(0,NA), expand = c(0.05,0,0,expandingAxis(max(population)))) +
    guides(fill = FALSE)+
    scale_fill_manual(values = col0)
  Filename <- paste(c(l,"_plot",i,".svg"),collapse = "")
  ggsave(filename=Filename,width=9,height = 12,units = "cm")
}
col10 <- c("#c52718","#74c7f3")


for (j in 1:FlowSOM_n_clusters) {
  cluster <- clusters_list_of_tables1[[j]] #value
  name <- paste("FlowSOM_Cluster_",j,sep="")
  plotting(cluster$value,cluster,name,j,1)
}

##### PHENOGRAPH CLUSTERING PLOTS ##############################################
plots <- plot_frequency_boxplot(fcd = condor,
                                cluster_slot = "phenograph_pca_norm_k_60", 
                                cluster_var = "Phenograph",
                                sample_var = "expfcs_filename", 
                                group_var = "Group", 
                                numeric = T)
# Get number of Phenograph clusters
#n_clusters_p <- nlevels(condor[["clustering"]][["phenograph_pca_norm_k_60"]][["Phenograph"]])
clusters_list_of_tables2 <- list()
for (i in 1:n_clusters_p) {
  name <- paste("Cluster_",i,sep="")
  newCluster <- plots[[name]][["data"]]
  clusters_list_of_tables2[[i]] <- newCluster
  filename <- paste("Pheno_",name,".csv",sep="")
  write.csv(newCluster,file.path(save_clusters,filename))
}

for (j in 1:n_clusters_p) {
  cluster <- clusters_list_of_tables2[[j]] #value
  name <- paste("Cluster_",j,sep="")
  plotting(cluster$value,cluster,name,j,2)
}

# library(dplyr)
# 
# # Compute centroids, cluster them, and merge
# centroids <- condor$expr %>% 
#   as.data.frame() %>%
#   mutate(cluster = condor$clustering$phenograph_pca_norm_k_60) %>%
#   group_by(cluster) %>%
#   summarise(across(everything(), mean))
# 
# hc <- hclust(dist(centroids[,-1]))
# merged <- cutree(hc, k = 10)  # Or use h = for height cutoff
# 
# # Map back to original clusters
# condor$clustering$phenograph_pca_norm_k_60 <- merged[match(condor$clustering$phenograph_pca_norm_k_60, centroids$cluster)]





################################################################################
# Extract cluster assignments as numeric
original_clusters <- as.numeric(as.character(condor$clustering$phenograph_pca_norm_k_60$Phenograph))

# Compute centroids, cluster them, and merge
centroids <- condor$expr %>% 
  as.data.frame() %>%
  mutate(cluster = original_clusters) %>%
  group_by(cluster) %>%
  summarise(across(everything(), mean))

hc <- hclust(dist(centroids[,-1]))

# Visualize to find good height
plot(hc, main = "Cluster Dendrogram", xlab = "Cluster ID", sub = "")
abline(h = 2.2, col = "red")

# Check current clusters
cat("Current clusters:", length(unique(original_clusters)), "\n")

# Cut by number of clusters
merged <- cutree(hc, k = 3)

# Map back to original clusters
# Create a mapping from old cluster IDs to new merged cluster IDs
cluster_mapping <- data.frame(
  old_cluster = centroids$cluster,
  new_cluster = merged
)

# Apply mapping to all cells
condor$clustering$phenograph_pca_norm_k_60_merged$Phenograph <- cluster_mapping$new_cluster[
  match(original_clusters, cluster_mapping$old_cluster)
]

cat("Merged clusters:", length(unique(condor$clustering$phenograph_pca_norm_k_60_merged)), "\n")


# After creating the merged clusters, format it as a data frame
condor$clustering$phenograph_pca_norm_k_60_merged <- data.frame(
  Phenograph = cluster_mapping$new_cluster[
    match(original_clusters, cluster_mapping$old_cluster)
  ],
  Description = "merged_k5"
)

# Set rownames to match original
rownames(condor$clustering$phenograph_pca_norm_k_60_merged) <- 
  rownames(condor$clustering$phenograph_pca_norm_k_60)

# Now plot
p <- plot_dim_red(fcd = condor,
                  expr_slot = NULL,
                  reduction_method = "umap",
                  reduction_slot = "pca_norm",
                  cluster_slot = "phenograph_pca_norm_k_60_merged",
                  param = "Phenograph",
                  title = "UMAP Phenograph Merged (k=5)")
p$data[["poi"]] <- factor(p$data[["poi"]])
p + 
  #scale_color_manual(values = c("#D32248", "#D35422","#D3AD22", "#22D3AD", "#2248D3"))+
  scale_color_manual(values = c("#2BC7F3","#2B64F3","#2BF3BA"))+
  guides(color = guide_legend(override.aes = list(size = 5)))
ggsave(filename = "plot444.svg",width=12,height = 12,units = "cm")

expr_data <- condor$expr$orig
cluster_data_pheno_merged <- condor$clustering$phenograph_pca_norm_k_60_merged$Phenograph
#create_heatmap(expr_data,cluster_data_pheno_merged,"pheno_merged")


plots <- plot_frequency_boxplot(fcd = condor,
                                cluster_slot = "phenograph_pca_norm_k_60_merged", 
                                cluster_var = "Phenograph",
                                sample_var = "expfcs_filename", 
                                group_var = "Group", 
                                numeric = T)

clusters_list_of_tables2 <- list()
for (i in 1:n_clusters_p) {
  name <- paste("Cluster_",i,sep="")
  newCluster <- plots[[name]][["data"]]
  clusters_list_of_tables2[[i]] <- newCluster
  filename <- paste("Pheno_merged_",name,".csv",sep="")
  write.csv(newCluster,file.path(save_clusters,filename))
}

for (j in 1:5) {
  cluster <- clusters_list_of_tables2[[j]] #value
  name <- paste("Cluster_merged_",j,sep="")
  plotting(cluster$value,cluster,name,j,3)
}


library(pheatmap)

library(RColorBrewer)

# Assuming you have:
# - condor_object: your Cycondor object
# - cluster_assignments: vector of cluster labels from phenograph

# Extract marker expression data from condor object
# Adjust this based on your condor object structure
# Common structures: condor_object$expr.data or condor_object$raw.data
expr_data <- condor$expr$orig  # or condor_object$raw.data

# Get cluster assignments from phenograph results
# If stored in condor object: condor_object$clustering$phenograph
generate_heatmaps <- function(pheno, name){
  #clusters <- tibble(condor$clustering$phenograph)  # Adjust path as needed
  clusters <- tibble(condor$clustering[[pheno]])
  print(clusters)
  clusters$Phenograph <- as.numeric(clusters$Phenograph)
  # Calculate mean expression for each cluster
  cluster_ids <- unique(clusters$Phenograph)
  n_clusters <- length(cluster_ids)
  n_markers <- ncol(expr_data)
  
  # Create matrix of cluster centroids (mean expression per cluster)
  cluster_centroids <- matrix(0, nrow = n_clusters, ncol = n_markers)
  rownames(cluster_centroids) <- paste0("Cluster_", cluster_ids)
  colnames(cluster_centroids) <- colnames(expr_data)
  
  for (i in seq_along(cluster_ids)) {
    cluster_cells <- clusters == cluster_ids[i]
    cluster_centroids[i, ] <- colMeans(expr_data[cluster_cells, , drop = FALSE])
  }
  
  # Calculate distance matrix between clusters
  # Using Euclidean distance by default
  
  cluster_dist <- as.matrix(dist(cluster_centroids, method = "euclidean"))
  rownames(cluster_dist) <- paste0("Cluster_", cluster_ids)
  colnames(cluster_dist) <- paste0("Cluster_", cluster_ids)
  
  # Create heatmap with pheatmap
  p1 <- pheatmap(
    cluster_dist,
    color = colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(100),
    display_numbers = TRUE,  # Show distance values
    number_format = "%.2f",  # Format numbers to 2 decimal places
    fontsize_number = 10,
    cluster_rows = TRUE,     # Hierarchical clustering of rows
    cluster_cols = TRUE,     # Hierarchical clustering of columns
    main = "Distance Matrix Between Phenograph Clusters",
    cellwidth = 30,
    cellheight = 30,
    border_color = "grey60"
  )
  name_plot <- paste("heatmap",name,"1.png",sep="")
  ggsave(name_plot,plot = p1,width=24,height = 24, units = "cm")
  
  
  # Optional: Try different distance metrics
  # Correlation distance
  cluster_dist_cor <- as.matrix(as.dist(1 - cor(t(cluster_centroids))))
  rownames(cluster_dist_cor) <- paste0("Cluster_", cluster_ids)
  colnames(cluster_dist_cor) <- paste0("Cluster_", cluster_ids)
  p3 <- pheatmap(
    cluster_dist_cor,
    color = colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(100),
    display_numbers = TRUE,
    number_format = "%.2f",
    fontsize_number = 10,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    main = "Correlation Distance Between Phenograph Clusters",
    cellwidth = 30,
    cellheight = 30,
    border_color = "grey60"
  )
  name_plot <- paste("heatmap",name,"3.png",sep="")
  ggsave(name_plot,plot = p3,width=24,height = 24, units = "cm")
  
  
  
  # Create heatmap with clusters ordered by similarity
  # pheatmap will automatically cluster columns (clusters) by distance
  p2 <- pheatmap(
    cluster_centroids,
    color = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100),
    scale = "column",           # Scale markers (rows) to z-scores
    cluster_rows = TRUE,     # Cluster markers
    cluster_cols = TRUE,     # Cluster columns (this orders clusters by distance)
    clustering_distance_cols = "euclidean",  # Distance metric for clusters
    clustering_method = "complete",          # Hierarchical clustering method
    show_colnames = TRUE,
    show_rownames = TRUE,
    angle_col = 45, 
    main = "Marker Expression Across Phenograph Clusters",
    fontsize_row = 20,
    fontsize_col = 20,
    border_color = NA
  )
  name_plot <- paste("heatmap",name,"2.svg",sep="")
  ggsave(name_plot,plot = p2,width=24,height = 24, units = "cm")
}

generate_heatmaps("phenograph_pca_norm_k_60","_og_")
generate_heatmaps("phenograph_pca_norm_k_60_merged","_merged_")

################################################################################
condor$clustering$phenograph_pca_norm_k_60$Phenograph





################################################################################
# 
# #### EXCLUDING CONTAMINATION OF CD4+ from yd T cells
# if (!dir.exists(format(Sys.Date(), "clean"))) {
#   dir.create(format(Sys.Date(), "clean"), recursive = TRUE)
# }
# setwd("clean")
# getwd()
# # ── 1. Check which column holds the cluster labels ──────────────────────────
# head(condor$clustering)          # inspect the clustering slot
# # typically the column is named something like "Phenograph" or "phenograph_cluster"
# 
# # ── 2. Identify cells to KEEP (exclude cluster 15) ──────────────────────────
# keep_cells <- rownames(condor$clustering$phenograph_pca_norm_k_60)[
#   condor$clustering$phenograph_pca_norm_k_60$Phenograph != 15
# ]
# 
# cat("Cells before removal:", nrow(condor$expr$orig), "\n")
# cat("Cells after removal: ", length(keep_cells), "\n")
# 
# # ── 3. Subset all slots of the condor object ────────────────────────────────
# condor_clean <- condor  # keep original intact
# 
# condor_clean$expr$orig        <- condor$expr$orig[keep_cells, ]
# condor_clean$pca$orig         <- condor$pca$orig[keep_cells, ]
# condor_clean$clustering$Phenograph_orig <- 
#   condor$clustering$Phenograph_orig[keep_cells, , drop = FALSE]
# 
# # If you have an annotation/sample table slot, subset that too
# # if (!is.null(condor$anno)) {
# #   print(nrow(condor_clean$anno)
# #   condor_clean$anno <- condor$anno[keep_cells, ]
# # }
# #anno_clean <- tibble(condor_clean$anno$cell_anno)
# condor_clean$anno$cell_anno <- condor$anno$cell_anno[keep_cells, ]
# cat("Sanity check — expr rows:", nrow(condor_clean$expr$orig), "\n")
# 
# # ── 4. Re-run PCA on clean data ──────────────────────────────────────────────
# condor_clean <- runPCA(fcd       = condor_clean,
#                        data_slot = "orig")
# 
# # ── 5. Re-run UMAP on clean data ─────────────────────────────────────────────
# condor_clean <- runUMAP(fcd        = condor_clean,
#                         input_type = "pca",
#                         data_slot  = "orig")
# 
# # ── 6. Re-run Phenograph clustering on clean data ───────────────────────────
# # condor_clean <- runPhenograph(fcd       = condor_clean,
# #                               input_type = "pca",
# #                               data_slot  = "orig",
# #                               k          = 60)   # use same k as your original run
# # 
# # 
# # 
# # 
# # 
# # # ── 8. Save the clean object ─────────────────────────────────────────────────
# # saveRDS(condor_clean, file = "condor_gdT_clean.rds")
# 
# ################################################################################
# 
# plot_dim_red(fcd= condor_clean,  
#              reduction_method = "umap", 
#              reduction_slot = "pca_orig", 
#              param = "Group", 
#              title = "UMAP colored by group"
# )
# ggsave(filename = "clean_plot.svg",path = dir_clust,width=12,height = 12,units = "cm")
# 
# 
# ################################################################################
# #### BATCH EFFECT ##############################################################
# 
# # Colored by run, displaying the batch effect
# plot_dim_red(fcd= condor_clean,  
#              reduction_method = "umap", 
#              reduction_slot = "pca_orig", 
#              param = "Run_no", 
#              title = "UMAP before edit"
# )
# ggsave(filename = "clean_plot1.svg",path = dir_clust,width=12,height = 12,units = "cm")
# 
# # Deleting batch effect
# condor_clean <- harmonize_PCA(fcd = condor_clean, 
#                               batch_var = c("Run_no"), 
#                               data_slot = "orig")
# 
# #condor_clean$pca$norm[1:10, 1:5]
# 
# condor_clean <- runUMAP(fcd = condor_clean, 
#                         input_type = "pca", 
#                         data_slot = "norm",
#                         prefix= NULL)
# 
# 
# plot_dim_red(fcd= condor_clean,  
#              reduction_method = "umap", 
#              reduction_slot = "pca_norm", 
#              param = "Run_no", 
#              title = "Harmonized UMAP")
# ggsave(filename = "clean_plot2.svg",path = dir_clust,width=12,height = 12,units = "cm")
# 
# plot_dim_red(fcd= condor_clean,  
#              reduction_method = "umap", 
#              reduction_slot = "pca_norm", 
#              param = "Group", 
#              title = "Harmonized UMAP")
# ggsave(filename = "clean_plot3.svg",path = dir_clust,width=12,height = 12,units = "cm")
# 
# 
# # split UMAP plots (showing groups from metadata, under param)
# plot_dim_red(fcd = condor_clean,
#              expr_slot = NULL,
#              reduction_method = "umap",
#              reduction_slot = "pca_norm",
#              cluster_slot = NULL,
#              param = "Group",
#              facet_by_variable = T,
#              title = "UMAP")
# ggsave(filename = "clean_plot4.svg",path = dir_clust,width=16,height = 8,units = "cm")
# 
# condor_clean <- runPhenograph(fcd       = condor_clean,
#                               input_type = "pca",
#                               data_slot  = "norm",
#                               k          = 60)   # use same k as your original run
# 
# 
# clustersExpression <- function(marker,clusterMethod,cluster_var) {
#   plot_marker_boxplot(fcd = condor_clean,
#                       expr_slot = "orig", 
#                       cluster_slot = clusterMethod, 
#                       cluster_var = cluster_var,
#                       group_var = "Group",
#                       sample_var = "SampleID", 
#                       marker = marker,
#                       fun = "mean",
#                       facet_ncol = 3,
#                       dot_size = 0.5)
#   ggsave(filename = paste("clean_clusters_",marker,"_expr.svg",sep=""),path = dir_clust,width=16,height = 8,units = "cm")
# }
# 
# clusteringMethods <- c("phenograph_pca_norm_k_60")
# clusterVar <- c("Phenograph")
# for (i in 1:ncol(condor_clean$expr$orig)) {
#   name <- colnames(condor_clean$expr$orig)[i]
#   print(colnames(condor_clean$expr$orig)[i])
#   #clusteringMethod <- clusteringMethods[j]
#   #if (name == "FlowAI") { next }
#   clustersExpression(name,"phenograph_pca_norm_k_60","Phenograph")
# } 
# 
# 
# plot_expression <- function(marker) {
#   plot_dim_red(fcd = condor_clean,
#                expr_slot = "orig",
#                reduction_method = "umap",
#                reduction_slot = "pca_norm",
#                cluster_slot = NULL,
#                param = marker, 
#                order = T,
#                title = paste(marker,"expression"))
#   ggsave(filename = paste("clean_plot_expr_",i,".svg",sep=""),path = dir_clust,width=12,height = 12,units = "cm")
# }
# 
# ## for cycle across all the markers present in the population ##################
# for (i in 1:ncol(condor_clean$expr$orig)) {
#   name <- colnames(condor_clean$expr$orig)[i]
#   plot_expression(name)
# } 
# 
# plot_marker_HM(fcd = condor_clean,
#                expr_slot = "orig",
#                marker_to_exclude = c("FSC-A","SSC-A"),
#                cluster_slot = "phenograph_pca_norm_k_60",
#                cluster_var = "Phenograph")
# 
# #ranran = data.frame(condor_clean[["clustering"]][["phenograph_pca_norm_k_60"]])
# #heatmap(ranran)
# ### heatmaps ###
# #library(pheatmap)
# library(dplyr)
# 
# # Extract expression data and cluster assignments
# expr_data <- condor_clean$expr$orig
# cluster_data_pheno <- condor_clean$clustering$phenograph_pca_norm_k_60$Phenograph
# cluster_data_flowSOM <- condor_clean$clustering[[2]]$FlowSOM
# 
# 
# expr_data <- condor_clean$expr$orig  # or condor_clean_object$raw.data
# 
# # Get cluster assignments from phenograph results
# # If stored in condor_clean object: condor_clean_object$clustering$phenograph
# generate_heatmaps <- function(pheno, name){
#   #clusters <- tibble(condor_clean$clustering$phenograph)  # Adjust path as needed
#   clusters <- tibble(condor_clean$clustering[[pheno]])
#   print(clusters)
#   clusters$Phenograph <- as.numeric(clusters$Phenograph)
#   # Calculate mean expression for each cluster
#   cluster_ids <- unique(clusters$Phenograph)
#   n_clusters <- length(cluster_ids)
#   n_markers <- ncol(expr_data)
#   
#   # Create matrix of cluster centroids (mean expression per cluster)
#   cluster_centroids <- matrix(0, nrow = n_clusters, ncol = n_markers)
#   rownames(cluster_centroids) <- paste0("Cluster_", cluster_ids)
#   colnames(cluster_centroids) <- colnames(expr_data)
#   
#   for (i in seq_along(cluster_ids)) {
#     cluster_cells <- clusters == cluster_ids[i]
#     cluster_centroids[i, ] <- colMeans(expr_data[cluster_cells, , drop = FALSE])
#   }
#   
#   # Calculate distance matrix between clusters
#   # Using Euclidean distance by default
#   
#   cluster_dist <- as.matrix(dist(cluster_centroids, method = "euclidean"))
#   rownames(cluster_dist) <- paste0("Cluster_", cluster_ids)
#   colnames(cluster_dist) <- paste0("Cluster_", cluster_ids)
#   
#   # Optional: Try different distance metrics
#   # Correlation distance
#   cluster_dist_cor <- as.matrix(as.dist(1 - cor(t(cluster_centroids))))
#   rownames(cluster_dist_cor) <- paste0("Cluster_", cluster_ids)
#   colnames(cluster_dist_cor) <- paste0("Cluster_", cluster_ids)
#   
#   # Create heatmap with clusters ordered by similarity
#   # pheatmap will automatically cluster columns (clusters) by distance
#   p2 <- pheatmap(
#     cluster_centroids,
#     color = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100),
#     scale = "column",           # Scale markers (rows) to z-scores
#     cluster_rows = TRUE,     # Cluster markers
#     cluster_cols = TRUE,     # Cluster columns (this orders clusters by distance)
#     clustering_distance_cols = "euclidean",  # Distance metric for clusters
#     clustering_method = "complete",          # Hierarchical clustering method
#     show_colnames = TRUE,
#     show_rownames = TRUE,
#     angle_col = 45, 
#     main = "Marker Expression Across Phenograph Clusters",
#     fontsize_row = 20,
#     fontsize_col = 20,
#     border_color = NA
#   )
#   name_plot <- paste("heatmap",name,"2.png",sep="")
#   ggsave(name_plot,plot = p2,width=24,height = 24, units = "cm")
# }
# 
# generate_heatmaps("phenograph_pca_norm_k_60","_og_")
# 
# 
# 
# plots <- plot_frequency_boxplot(fcd = condor_clean,
#                                 cluster_slot = "phenograph_pca_norm_k_60", 
#                                 cluster_var = "Phenograph",
#                                 sample_var = "expfcs_filename", 
#                                 group_var = "Group", 
#                                 numeric = T)
# # Get number of Phenograph clusters
# #n_clusters_p <- nlevels(condor_clean_clean[["clustering"]][["phenograph_pca_norm_k_60"]][["Phenograph"]])
# clusters_list_of_tables2 <- list()
# for (i in 1:n_clusters_p) {
#   name <- paste("Cluster_",i,sep="")
#   newCluster <- plots[[name]][["data"]]
#   clusters_list_of_tables2[[i]] <- newCluster
#   filename <- paste("Pheno_",name,".csv",sep="")
#   write.csv(newCluster,file.path(save_clusters,filename))
# }
# unique(condor_clean$clustering$phenograph_pca_norm_k_60$Phenograph)
# 
# for (j in 1:length(unique(condor_clean$clustering$phenograph_pca_norm_k_60$Phenograph))) {
#   cluster <- clusters_list_of_tables2[[j]] #value
#   name <- paste("Cluster_",j,sep="")
#   plotting(cluster$value,cluster,name,j,2)
# }
# 
# plot_dim_red(fcd = condor_clean,
#              expr_slot = NULL,
#              reduction_method = "umap",
#              reduction_slot = "pca_norm",
#              cluster_slot = "phenograph_pca_norm_k_60",
#              param = "Phenograph",
#              title = "UMAP Phenograph")
# ggsave(filename = "clean_plot6.svg",path = dir_clust,width=12,height = 12,units = "cm")
# ggsave(filename = "clean_plot6.png",path = directory,width=12,height = 12,units = "cm")

# plot_dim_red(fcd = condor_clean,
#              expr_slot = NULL,
#              reduction_method = "umap",
#              reduction_slot = "pca_norm",
#              cluster_slot = FlowSOM_name,
#              param = "FlowSOM",
#              title = "UMAP FlowSOM clusters")
# ggsave(filename = "clean_plot7.svg",path = dir_clust,width=12,height = 12,units = "cm")
# #ggsave(filename = "plot7.svg",path = directory,width=12,height = 12,units = "cm")
