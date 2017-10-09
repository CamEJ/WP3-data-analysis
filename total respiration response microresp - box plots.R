
# plotting total respiration response microresp

library(ggplot2)


setwd("C:/Users/Camilla/Dropbox/Data & analysis/WP3 Slurry disturbance/Plotting Data/PhenotypicData-R")

pod <- read.table("MResp-BioReps_UnTransformed.txt", sep = "\t", header=T, row.names=1)
meta <- read.table("Mresp-metadata.txt", header = TRUE, row.names = 1, sep='\t')



pod[,"tot"] <- rowSums(pod)

g = pod
g$treatment = meta$Treatment
g$time = meta$timepoint

g[["time"]] <- setFactorOrder(g[["time"]], c("T0", "T1", "T2", "T3", "T4", "T5","T6","T7", "T8", "T9","T10","T11","T12", "T13"))

# order of treatments. 
g[["treatment"]] <- setFactorOrder(g[["treatment"]], c("Control", "Slurry", "Flood", "Flood+Slurry"))
# 

write.csv(g, "totalSubsInducedRespirationResponse.csv")

c <- ggplot(g, aes(factor(treatment), tot, fill = factor(treatment))) +
  
  ## + geom_boxplot so it knows what type of plot
  # and put colour = black to make lines of box black. 
  
  geom_boxplot(colour="black") +
  scale_fill_manual(values=c("black", "chocolate4", "slateblue", "olivedrab"))
c4 = c + facet_wrap(~time, ncol = 3)


rep <- c4 + labs(fill="    Treatment ") +
  ylab(expression(paste('Total substrate induced respiration (', mu, 'g ', C-CO[2],' ', g^-1, 'soil ', hr^-1, ')', '\n' )))+
  
  ## specify labels for axes and plot title if required
  
  theme_bw() +
  
  
  ## change text size and placing of axes and titles
  ## ensure , at end of each line 
  ## legend.key.size changes size of legend and required 'grid' package to be loaded for unit() to work
  
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_text(size=16, vjust=0.5, colour = "black"),
        axis.title.y=element_text(size=18, vjust=1, colour="black"),
        legend.text=element_text(size=15, vjust=0.5),
        legend.title=element_blank(),
        legend.key.size=unit(1.2, "cm"),
        axis.title.x=element_blank(),
        legend.direction = "horizontal",
        strip.text.x = element_text(size = 18, colour = "black"),# change font of facet label
        strip.background =  element_rect(fill = "white"),
        legend.position = "top"
  ) # remove grey backgroup of facet label

rep
# === different tings with legend 
rep + guides(fill = guide_legend(ncol=2))

rep + guides(fill=FALSE)

# --------------------- scatter plot of the same


  
  c <- ggplot(g, aes(factor(treatment), tot, fill = factor(treatment))) +
  
  ## + geom_boxplot so it knows what type of plot
  # and put colour = black to make lines of box black. 
  
    geom_point(shape=21, size = 5) +
  scale_fill_manual(values=c("black", "chocolate4", "slateblue", "olivedrab"))
c4 = c + facet_wrap(~time, ncol = 3)


c4
rep <- c4 + labs(fill="    Treatment ") +
  ylab(expression(paste('Total substrate induced respiration (', mu, 'g ', C-CO[2],' ', g^-1, 'soil ', hr^-1, ')', '\n' )))+
  
  ## specify labels for axes and plot title if required
  
  theme_bw() +
  
  
  ## change text size and placing of axes and titles
  ## ensure , at end of each line 
  ## legend.key.size changes size of legend and required 'grid' package to be loaded for unit() to work
  
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_text(size=16, vjust=0.5, colour = "black"),
        axis.title.y=element_text(size=18, vjust=1, colour="black"),
        legend.text=element_text(size=15, vjust=0.5),
        legend.title=element_blank(),
        legend.key.size=unit(1.2, "cm"),
        axis.title.x=element_blank(),
        #legend.direction = "horizontal",
        strip.text.x = element_text(size = 18, colour = "black"),# change font of facet label
        strip.background =  element_rect(fill = "white"),
        legend.position = c(1, 0), legend.justification = c(1, 0)
  ) # remove grey backgroup of facet label

rep

rep + guides(fill = guide_legend(ncol=2))
  



# - --------------------------------- trying some other stuff that didnt really work. 




# Basic dot plot
p<-ggplot(ToothGrowth, aes(x=dose, y=len)) + 
  geom_dotplot(binaxis='y', stackdir='center')
p
# Change dotsize and stack ratio
ggplot(ToothGrowth, aes(x=dose, y=len)) + 
  geom_dotplot(binaxis='y', stackdir='center',
               stackratio=1.5, dotsize=1.2)
# Rotate the dot plot
p + coord_flip()

# Basic dot plot
p<-ggplot(g, aes(x=treatment, y=tot)) + 
  geom_dotplot(binaxis='y', stackdir='center') +
  facet_wrap(~time, ncol = 3)
p
# Change dotsize and stack ratio
p = ggplot(g, aes(x=treatment, y=tot, colour = treatment)) + 
  geom_dotplot(binaxis='y', stackdir='center',
               stackratio=1.5, dotsize=1.2, fill = "grey69")+
  facet_wrap(~time, ncol = 3)



# Function to produce summary statistics (mean and +/- sd)
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}
#Use a custom summary function :
  
  p + stat_summary(fun.data=data_summary, colour = "black", size=1)
  