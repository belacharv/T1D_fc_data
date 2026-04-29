library(readxl)
library(tidyverse)
library(ggpubr)
library("ggbeeswarm")
getwd()
directory <- "/home/charvatb/PRIMUS/data/48_lab/Project Diabetes Motol/FACS validation/240311 Analysis BCH/extracellular_pannel/"
setwd(directory)
gdTcells <- read_xlsx(paste(directory,"241120_gdTcells_extra_dia_panel.xlsx",sep=""))
naive_ef <- read_xlsx(paste(directory,"241120_naive_vs_eff.xlsx",sep=""))
tem_tcm <- read_xlsx(paste(directory,"241120_gating_tem_tcm_extra_dia_panel.xlsx",sep=""))
cd8low <- read_xlsx(paste(directory,"241120_cd8low_extra_dia_panel.xlsx",sep=""))
other <- read_xlsx(paste(directory,"241120_other_gating_extra_dia_panel.xlsx",sep=""))
patients <- read_xlsx("/home/charvatb/PRIMUS/data/48_lab/Project Diabetes Motol/FACS validation/240311 Analysis BCH/patients_data/metadata_v07.xlsx")
VN_gating <- read_xlsx(paste(directory,"241116_gating_VN_extra_dia_motol.xlsx",sep=""))
new_yd <- read_csv(paste(directory,"results/260403_yd_T_cells_extra_dia_motol.csv",sep=""))
cxcr3_cd56 <- read_csv(paste(directory,"results/260403_cd56_cxcr3_extra_dia_motol.csv",sep=""))
extra_clust <- read_csv(paste(directory,"results/260403_extra_clust_after_gating_dia_motol.csv",sep=""))
all_together <- read_csv(paste(directory,"results/260404_all_gated_extra_dia_motol1.csv",sep=""))
all_together1 <- read_csv(paste(directory,"results/260414_all_gated_extra_dia_motol1.csv",sep=""))
all_together1 <- read_csv(paste(directory,"results/260415_all_gated_extra_dia_motol_before.csv",sep=""))
all_together2 <- read_csv(paste(directory,"results/260428_all_gated_extra_dia_motol_after.csv",sep=""))
#naive_ef$CD4_eff_naive_ratio <- naive_ef$`FlowAIGoodEvents/Lymphocytes/Single Cells/Single Cells/Live cells/CD3+ CD16-/CD4+/CD45RA+CCR7+ neg | Freq. of Parent`/naive_ef$`FlowAIGoodEvents/Lymphocytes/Single Cells/Single Cells/Live cells/CD3+ CD16-/CD4+/CD45RA+CCR7+ | Freq. of Parent`
edit_table <- function(data) {
  colnames(data) <- colnames(data) %>%
    str_remove(" Freq. of") %>%
    #str_remove("CD8_naive_") %>%
    str_remove("TEMRA_") %>%
    str_remove("FlowAIGoodEvents/Lymphocytes/Single Cells/Single Cells/Live cells/") %>%
    str_remove(" CD16-/")
  #colnames(data) <- make.unique(colnames(data), sep = "_dup")
  colnames(data)[1] <- "sample"
  data <- data %>%
    mutate_at(vars(2:ncol(data)),as.numeric) %>%
    mutate(SampleID = substr(sample,4,length(sample)))
  data$SampleID <- data$SampleID %>%
    str_remove(" ") %>%
    #str_remove(".") %>%
    str_remove(".fcs")
  data <- data %>%
    mutate(Disease = ifelse(substr(SampleID,1,1) == 2,"0","1")) %>%
    mutate(Disease = as.numeric(Disease)) %>%
    mutate(Timepoint = ifelse(substr(SampleID,4,4) == "b","1","0")) %>%
    mutate(Timepoint = as.numeric(Timepoint)) %>%
    mutate(Group = ifelse(Disease == 1,ifelse(Timepoint == 1,"DiaT1","DiaT0"),"Control")) %>%
    #mutate(SampleID = substr(SampleID,1,3)) %>%
    mutate(id_pat = substr(SampleID,1,3))
  #mutate(Disease = substr(SampleID,1,1))
  return(data)
}
data1 <-edit_table(gdTcells)
data2 <- edit_table(naive_ef)
data3 <- edit_table(tem_tcm)
data4 <- edit_table(cd8low)
data5 <- edit_table(other)
data6 <- edit_table(VN_gating)

edit_table2 <- function(data) {
  colnames(data) <- colnames(data) %>%
    str_remove(" Freq. of") %>%
    #str_remove("CD8_naive_") %>%
    str_remove("TEMRA_") %>%
    str_remove("FlowAIGoodEvents/Lymphocytes/Single Cells/Single Cells/Live cells/") %>%
    str_remove(" CD16-/")
  colnames(data)[1] <- "sample"
  data <- data %>%
    mutate_at(vars(2:ncol(data)),as.numeric) %>%
    mutate(SampleID = sample)
  data <- data %>%
    mutate(Disease = ifelse(substr(SampleID,1,1) == 2,"0","1")) %>%
    mutate(Disease = as.numeric(Disease)) %>%
    mutate(Timepoint = ifelse(substr(SampleID,4,4) == "b","1","0")) %>%
    mutate(Timepoint = as.numeric(Timepoint)) %>%
    mutate(Group = ifelse(Disease == 1,ifelse(Timepoint == 1,"DiaT1","DiaT0"),"Control")) %>%
    mutate(id_pat = as.numeric(substr(SampleID,1,3)))
  return(data)
}

data7 <- edit_table2(new_yd)
data8 <- edit_table2(extra_clust)
data9 <- edit_table2(cxcr3_cd56)
data10 <- edit_table2(all_together1)
data11 <- edit_table(all_together2)
col0 <- c("#7CD6FC","#fa7065","#be2115")
col1 <- c("#7CD6FC","#fa7065")
col1 <- c("#00ADD4","#DC001E")

comparisons <- list(c("Control", "DiaT0"),c("DiaT0","DiaT1"), c("Control", "DiaT1"))  

#data1a <- filter(data1,!(data1$Disease == 1 & data1$Timepoint == 1))

expandingAxis <- function(maxValue) {
  expand <- 0
  ifelse(maxValue > 50,
         ifelse(maxValue >80,expand <- 15,expand <- 10),
         ifelse(maxValue >10,expand <- 5,expand <- 3))
  return(expand)
}

setwd(paste(directory,"results/gating_plots/",sep=""))
if (!dir.exists(format(Sys.Date(), "%y%m%d"))) {
  dir.create(format(Sys.Date(), "%y%m%d"), recursive = TRUE)
}
setwd(format(Sys.Date(), "%y%m%d"))
dir_save = getwd()

ggtheme <- function() {
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 8),
    text = element_text(size = 20, colour = "black", family = "Arial"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank())
  #legend.text = element_text(size = 20),
  #legend.key.size = unit(10, units = "points"))
  
} 

ggtheme <- function() {
  theme(
    axis.text = element_text(size = 26),
    axis.title = element_text(size = 26),
    plot.title = element_text(size = 8),
    text = element_text(size = 26, colour = "black", family = "Arial"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank())
  #legend.text = element_text(size = 20),
  #legend.key.size = unit(10, units = "points"))
  
} 




plotting <- function(population,format,i,data,j) 
{
  p <- ggplot(aes(x = Group, y = population), data=data) +
    geom_boxplot(alpha = 0.2,aes(fill = Group)) +
    geom_dotplot(binaxis='y',
                 dotsize = 1.5,
                 stackdir='center',
                 aes(fill = Group)) +
    labs(title = colnames(data)[i], y=population) +
    theme_classic() +
    ggpubr::stat_compare_means(label.x= 1.5,size = 9, vjust = -0.05, 
                               label = "p.format",comparisons = comparisons)  +
    xlab("Condition") +
    ggtheme() +
    scale_y_continuous(limits = c(0,NA), expand = expansion(mult = c(0.05, 0.25))) +
    guides(fill = FALSE)+
    scale_fill_manual(values = col0)
  Filename <- paste(c(j,"_plot_",i,format),collapse = "")
  ggsave(plot = p,filename=Filename,path=dir_save,width=9,height = 12,units = "cm")
}


col10 <- c("#74c7f3","#c52718")


# for (i in 2:(ncol(data1)-5)){
#   data1[[i]] <- as.numeric(data1[[i]])
# }

# for (i in 2:(ncol(data1)-5)){
#   plotting(data1[[i]],".png",i,data1)
# }

#dir_save <- paste(dir_save,"/axis_limit",sep="")
data10$CD4_Naive_log <- log((data10$`CD3+CD4+/CD45RA+CCR7+ neg | Parent`)/(data10$`CD3+CD4+/CD45RA+CCR7+ | Parent`))
data10$CD8_Naive_log <- log(data10$`CD3+CD8+/CD45RA+CCR7+ neg | Parent`/data10$`CD3+CD8+/CD45RA+CCR7+ | Parent`)
data11$CD4_Naive_log <- log(data11$`CD3+CD4+/CD45RA+CCR7+ neg | Parent`/data11$`CD3+CD4+/CD45RA+CCR7+ | Parent`)
data11$CD8_Naive_log <- log(data11$`CD3+CD8+/CD45RA+CCR7+ neg | Parent`/data11$`CD3+CD8+/CD45RA+CCR7+ | Parent`)
#all_data <- list(data1,data2,data3,data4,data5,data6,data7,data8,data9)
all_data <- list(data11)
for (j in 1:length(all_data)){
  for (i in 2:(ncol(all_data[[j]])-7)){
    plotting(all_data[[j]][[i]],".png",i,all_data[[j]],j)
  }
  plotting(all_data[[j]][[ncol(all_data[[j]])-2]],".png",i+1,all_data[[j]],j)
  plotting(all_data[[j]][[ncol(all_data[[j]])-1]],".png",i+2,all_data[[j]],j)
}
setwd("..")
setwd("..")
#selected = c(2:ncol(all_data[1]))

for (j in 1:length(all_data)){
  for (i in 2:(ncol(all_data[[j]])-7)){
  plotting(all_data[[j]][[i]],".svg",i,all_data[[j]],j)
  }
  plotting(all_data[[j]][[ncol(all_data[[j]])-2]],".svg",i+1,all_data[[j]],j)
  plotting(all_data[[j]][[ncol(all_data[[j]])-1]],".svg",i+2,all_data[[j]],j)
}


dir_save <- paste(dir_save,"/axis_limit",sep="")
plotting <- function(population,format,i,data,j) 
{
  p <- ggplot(aes(x = Group, y = population), data=data) +
    geom_boxplot(alpha = 0.2,aes(fill = Group)) +
    geom_dotplot(binaxis='y',
                 dotsize = 1.5,
                 stackdir='center',
                 aes(fill = Group)) +
    labs(title = colnames(data)[i], y=population) +
    theme_classic() +
    ggpubr::stat_compare_means(label.x= 1.5,size = 9, vjust = -0.05, 
                               label = "p.format",comparisons = comparisons)  +
    xlab("Condition") +
    ggtheme() +
    scale_y_continuous(limits = c(0,ifelse(max(population)>85,100,NA))) +
    guides(fill = FALSE)+
    scale_fill_manual(values = col0)
  Filename <- paste(c(j,"_plot_",i,format),collapse = "")
  ggsave(plot = p,filename=Filename,path=dir_save,width=9,height = 12,units = "cm")
}

for (j in 1:length(all_data)){
  for (i in 2:(ncol(all_data[[j]])-7)){
    plotting(all_data[[j]][[i]],".png",i,all_data[[j]],j)
  }
  plotting(all_data[[j]][[ncol(all_data[[j]])-2]],".png",i+1,all_data[[j]],j)
  plotting(all_data[[j]][[ncol(all_data[[j]])-1]],".png",i+2,all_data[[j]],j)
}
#setwd("..")
#setwd("..")
#selected = c(2:ncol(all_data[1]))

for (j in 1:length(all_data)){
  for (i in 2:(ncol(all_data[[j]])-7)){
    plotting(all_data[[j]][[i]],".svg",i,all_data[[j]],j)
  }
  plotting(all_data[[j]][[ncol(all_data[[j]])-2]],".svg",i+1,all_data[[j]],j)
  plotting(all_data[[j]][[ncol(all_data[[j]])-1]],".svg",i+2,all_data[[j]],j)
}

plotting <- function(population,format,i,data,j) 
{
  p <- ggplot(aes(x = Group, y = population), data=data) +
    geom_boxplot(alpha = 0.2,aes(fill = Group)) +
    geom_dotplot(binaxis='y',
                 dotsize = 1.5,
                 stackdir='center',
                 aes(fill = Group)) +
    labs(title = colnames(data)[i], y=population) +
    theme_classic() +
    ggpubr::stat_compare_means(label.x= 1.5,size = 9, vjust = -0.05, 
                               label = "p.format",comparisons = comparisons)  +
    xlab("Condition") +
    ggtheme() +
    scale_y_continuous(limits = c(NA,ifelse(max(population)>85,100,NA))) +
    guides(fill = FALSE)+
    scale_fill_manual(values = col0)
  Filename <- paste(c(j,"_plot_",i,format),collapse = "")
  ggsave(plot = p,filename=Filename,path=dir_save,width=9,height = 12,units = "cm")
}


plotting(data10[[121]],".svg",400,data10,3)
plotting(data10[[129]],".svg",401,data10,3)
plotting(data10[[130]],".svg",402,data10,3)
plotting(data11[[129]],".svg",403,data11,3)
plotting(data11[[130]],".svg",404,data11,3)
plotting(data11[[131]],".svg",405,data11,3)


data2$CD4_eff_naive_ratio <- data2$`CD3+CD4+/CD45RA+CCR7+ neg | Parent`/data2$`CD3+CD4+/CD45RA+CCR7+ | Parent`
data2$CD8_eff_naive_ratio <- data2$`CD3+CD8+/CD45RA+CCR7+ neg | Parent`/data2$`CD3+CD8+/CD45RA+CCR7+ | Parent`
sel2 = c(3,4,6,7,19,20)
for (sele in sel2){
  plotting(data2[sele],".svg",sele,data2,2)
}


plotting(data3[[12]],".svg",12,data3,3)
plotting(data4[[5]],".svg",5,data4,4)



p <- ggplot(aes(x = Group, y = data4[[5]]), data=data4) +
  geom_boxplot(alpha = 0.2,aes(fill = Group)) +
  geom_dotplot(binaxis='y',
               dotsize = 1.5,
               stackdir='center',
               aes(fill = Group)) +
  labs(title = colnames(data4)[5], y=population) +
  theme_classic() +
  ggpubr::stat_compare_means(label.x= 1.5,size = 7, vjust = -0.05, 
                             label = "p.format",comparisons = comparisons)  +
  xlab("Condition") +
  ggtheme() +
  scale_y_continuous(limits = c(0,90), expand = expansion(mult = c(0.05, 0.25))) +
  guides(fill = FALSE)+
  scale_fill_manual(values = col0)
Filename <- paste(c(04,"_plot_",500,".svg"),collapse = "")
ggsave(plot = p,filename=Filename,path=dir_save,width=9,height = 12,units = "cm")

