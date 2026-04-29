library(readxl)
library(tidyverse)
library(ggpubr)
getwd()
directory <- "/home/charvatb/PRIMUS/data/48_lab/Project Diabetes Motol/FACS validation/240311 Analysis BCH/intracellular_pannels/dia_intra01/"
setwd(directory)
cd4 <- read_xlsx(paste(directory,"251005_intra01_gating.xlsx",sep=""))
cd4 <- cd4[,-c(2,3,4,5,6,7)]
intra_new <- read_csv(paste(directory,"260428_intra_tregs_edit3.csv",sep=""))
patients <- read_xlsx("/home/charvatb/PRIMUS/data/48_lab/Project Diabetes Motol/FACS validation/240311 Analysis BCH/patients_data/metadata_v07.xlsx")

 
colnames(patients)[1] = "id_pat"
setwd(paste(directory,"results/gating_plots/",sep=""))
if (!dir.exists(format(Sys.Date(), "%y%m%d"))) {
  dir.create(format(Sys.Date(), "%y%m%d"), recursive = TRUE)
}
dir_save = paste(directory,"results/gating_plots/",format(Sys.Date(), "%y%m%d"),sep="")

colnames(cd4)
#str_remove(colnames(cd4),"FlowAIGoodEvents/Lymphocytes/Single Cells/Single Cells/SSC-B-H, SSC-H subset/Via, FSC-A subset/")
edit_table <- function(data) {
  colnames(data) <- colnames(data) %>%
    str_remove(" Freq. of") %>%
    #str_remove("CD8_naive_") %>%
    #str_remove("TEMRA_") %>%
    str_remove("FlowAIGoodEvents/Lymphocytes/Single Cells/Single Cells/SSC-B-H, SSC-H subset/Via, FSC-A subset/") %>%
    str_remove(" CD8-")
  colnames(data)[1] <- "sample"
  data <- data %>%
    mutate_at(vars(2:ncol(data)),as.numeric) %>%
    mutate(SampleID = substr(sample,4,8))
  data$SampleID <- data$SampleID %>%
    str_remove(" ") %>%
    str_remove("f")
  data <- data %>%
    mutate(Disease = ifelse(substr(SampleID,1,1) == 2,"0","1")) %>%
    mutate(Disease = as.numeric(Disease)) %>%
    mutate(Timepoint = ifelse(substr(SampleID,4,4) == "b","1","0")) %>%
    mutate(Timepoint = as.numeric(Timepoint)) %>%
    mutate(Group = ifelse(Disease == 1,ifelse(Timepoint == 1,"Dia T1","Dia T0"),"Control")) %>%
    mutate(SampleID = substr(SampleID,1,3)) %>%
    mutate(id_pat = as.numeric(SampleID))
  #mutate(Disease = substr(SampleID,1,1))
  return(data)
}

data1 <- edit_table(cd4)


edit_table2 <- function(data) {
  colnames(data) <- colnames(data) %>%
    str_remove(" Freq. of") %>%
    #str_remove("TEMRA_") %>%
    #str_remove("FlowAIGoodEvents/Lymphocytes/Single Cells/Single Cells/Live cells/") %>%
    str_remove("FlowAIGoodEvents/Lymphocytes/Single Cells/Single Cells/SSC-B-H, SSC-H subset/Via, FSC-A subset/") %>%
    str_remove(" CD8-/")
  colnames(data)[1] <- "sample"
  data <- data %>%
    mutate_at(vars(2:ncol(data)),as.numeric) %>%
    mutate(SampleID = sample)
  data <- data %>%
    mutate(Disease = ifelse(substr(SampleID,1,1) == 2,"0","1")) %>%
    mutate(Disease = as.numeric(Disease)) %>%
    mutate(Timepoint = ifelse(substr(SampleID,4,4) == "b","1","0")) %>%
    mutate(Timepoint = as.numeric(Timepoint)) %>%
    mutate(Group = ifelse(Disease == 1,ifelse(Timepoint == 1,"Dia T1","Dia T0"),"Control")) %>%
    mutate(id_pat = as.numeric(substr(SampleID,1,3)))
  return(data)
}

data2 <- edit_table(intra_new)
col0 <- c("#7CD6FC","#7CA5FC","#927CFC")
col0 <- c("#7CD6FC","#fa7065","#be2115")
col1 <- c("#7CD6FC","#fa7065")
col1 <- c("#00ADD4","#DC001E")

comparisons <- list(c("Control", "Dia T0"),c("Dia T0","Dia T1"), c("Control", "Dia T1"))  

data1a <- filter(data1,!(data1$Disease == 1 & data1$Timepoint == 1))
data1aa <- data.frame(data1a)
expandingAxis <- function(maxValue) {
  expand <- 0
  ifelse(maxValue > 50,
         ifelse(maxValue >80,expand <- 15,expand <- 10),
         ifelse(maxValue >10,expand <- 5,expand <- 3))
  return(expand)
}


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


plotting <- function(data,population,format,i,j) 
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
  Filename <- paste(c(j,"_plot",i,format),collapse = "")
  ggsave(plot = p,filename=Filename,path=dir_save,width=9,height = 12,units = "cm")
}


col10 <- c("#74c7f3","#c52718")

data_all <- list(data1,data2)
for (j in 1:length(data_all)){
  for (i in 2:(ncol(data_all[[j]])-5)){
    data_all[[j]][[i]] <- as.numeric(data_all[[j]][[i]])
  }
}
  
for (j in 1:length(data_all)){
  for (i in 2:(ncol(data_all[[j]])-5)){
    plotting(data_all[[j]],data_all[[j]][[i]],".png",i,j)
  }
}

selected = c(2:(ncol(data2)-5)
             #4,18,19,22,27,29,30,31,35,38,39,44,45
             )
for (sel in selected){
  plotting(data2,data2[[sel]],".svg",sel,2)
}



plotting <- function(data,population,format,i,j) 
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
    scale_y_continuous(limits = c(0,100)) +
    guides(fill = FALSE)+
    scale_fill_manual(values = col0)
  Filename <- paste(c(j,"_plotaxis",i,format),collapse = "")
  ggsave(plot = p,filename=Filename,path=dir_save,width=9,height = 12,units = "cm")
}

for (j in 1:length(data_all)){
  for (i in 2:(ncol(data_all[[j]])-5)){
    plotting(data_all[[j]],data_all[[j]][[i]],".png",i,j)
  }
}

selected = c(2:(ncol(data2)-5)
             #4,18,19,22,27,29,30,31,35,38,39,44,45
)
for (sel in selected){
  plotting(data2,data2[[sel]],".svg",sel,2)
}




################################################################################

p <- ggplot(aes(x = Group, y = data1[[3]]), data=data1) +
  geom_boxplot(alpha = 0.2,aes(fill = Group)) +
  geom_dotplot(binaxis='y',
               dotsize = 1.5,
               stackdir='center',
               aes(fill = Group)) +
  labs(title = colnames(data1)[3], y=population) +
  theme_classic() +
  ggpubr::stat_compare_means(label.x= 1.5,size = 9, vjust = -0.05, 
                             label = "p.format",comparisons = comparisons)  +
  xlab("Condition") +
  ggtheme() +
  scale_y_continuous(limits = c(0,95), expand = expansion(mult = c(0.05, 0.25))) +
  guides(fill = FALSE)+
  scale_fill_manual(values = col0)
Filename <- paste(c("plot",333,".svg"),collapse = "")
ggsave(plot = p,filename=Filename,path=dir_save,width=9,height = 12,units = "cm")


### divided by group ###########################################################

ggplot(aes(x = Group, y = data2$`CD3+ CD16-/CD8+/CD45RA+CCR7+ neg | Parent`), data=data2) +
  #ylim(0,30) +
  geom_boxplot(alpha = 0.2,aes(fill = Group)) +
  
  # Add points that will be colored by group
  geom_dotplot(binaxis='y',
               dotsize = 0.7,
               stackdir='center',
               aes(fill = Group)) +
  #geom_jitter(aes(colour = Group)) +
  
  labs(title = "Effector CD8+ T cells") +
  # Select a nice theme
  theme_classic() +
  ggpubr::stat_compare_means(label.x= 1.5,size = 6, vjust = -0.05, 
                             label = "p.format",comparisons = comparisons)  +
  ylab("Frequency") +
  xlab("Condition") +
  facet_wrap(~Group, scales = "free", ncol = 9) +
  scale_y_continuous(limits = c(0,NA), expand = c(0.05,0,0,10)) +
  guides(fill = FALSE)+
  scale_fill_manual(values = col0)
Filename <- paste(c("plot",4,".png"),collapse = "")
ggsave(filename=Filename,path=directory,width=9,height = 12,units = "cm")



corplotFunction <- function(pop1,population,data,name,i) {
  data %>%
    ggplot(aes(x=pop1,y = population,color = Group)) +
    ylim(0,NA) +
    geom_point()+
    theme_classic() +
    labs(title = name,x=name, y = 'frequency population') +
    geom_smooth(method = "lm", se = FALSE)+
    scale_fill_manual(values = col10)
}
