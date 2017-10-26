
# dot plots of alpha diversity


library(ggplot2)
library(grid)
library(ggthemr)
library(dplyr)


# read in table

# # what my data looks like
head(do)
# percov     sobs     chao simpsoneven invsimpson
# T0_Tmt1_a 92.6538 3752.244 7376.146    0.006636   24.89716
# T0_Tmt1_b 92.8413 3675.903 7056.404    0.005307   19.50796
# T0_Tmt1_c 93.2834 3580.545 6430.525    0.005305   18.99483
# T0_Tmt2_a 92.3209 3995.823 7688.366    0.013165   52.60369
# T0_Tmt2_b 92.1468 4143.379 7679.125    0.009722   40.27926
# T0_Tmt2_c 92.5942 3963.788 7351.110    0.010625   42.11221


# =============== plotting alpha diveristy by day 


head(do)
do$id = row.names(do)

# 1. Do first for inv simpson diversity


# make averages of dataframe (called 'do' here ) using dplyr

AvsDiv <- do %>% 
  group_by(timepoint, treatment) %>% # define how your biological replicates are grouped. If only one treatment then just group_by(MyTreatment)
  summarise(av = mean(invsimpson), # mean of
            sdv = sd(invsimpson),
            phs = first(phase),
            total = n())


# order you want timepoints to occur run setFactorOrder.R function in R basics folder (github)

AvsDiv[["timepoint"]] <- setFactorOrder(AvsDiv[["timepoint"]], c("T0", "T2", "T3",  "T5","T7", "T8", "T11", "T13"))
AvsDiv[["treatment"]] <- setFactorOrder(AvsDiv[["treatment"]], c("Control", "Slurry", "Flood", "Flood+Slurry"))



Div = ggplot(AvsDiv, aes(x=timepoint, y=av, ymin=av-sdv, ymax=av+sdv))+
  geom_pointrange(size=0.75) + # lines (error lines)
  geom_point(colour="black", shape=21, size = 7) + # shapes = c(21, 11)
  aes(fill = factor(treatment)) # colour of fill be treatment 


divPlot = Div + theme_bw()   + labs(y = "Average diversity\n(Simpson's index)\n", fill = "Treatment") +
  theme(axis.title.x=element_blank(), 
        axis.text.y=element_text(size=13, colour="black"),
        axis.text.x=element_text(size=13, colour="black"),
        axis.title.y=element_text(size=13, vjust=1, colour="black"),
        legend.title=element_text(size=14, vjust=1),
        legend.text=element_text(size=13, vjust=0.5),
        legend.key.size=unit(1.1, "cm"),
        legend.position="top"
  )

# Second plot simpsons evenness

head(do)
AvsEven <- do %>% 
  group_by(timepoint, treatment) %>%
  summarise(av = mean(simpsoneven), 
            sdv = sd(simpsoneven),
            phs = first(phase),
            total = n())
# order you want timepoints to occur
AvsEven[["timepoint"]] <- setFactorOrder(AvsEven[["timepoint"]], c("T0", "T2", "T3",  "T5","T7", "T8", "T11", "T13"))

# order of treatments. 
AvsEven[["treatment"]] <- setFactorOrder(AvsEven[["treatment"]], c("Control", "Slurry", "Flood", "Flood+Slurry"))



even = ggplot(AvsEven, aes(x=timepoint, y=av, ymin=av-sdv, ymax=av+sdv))+
  geom_pointrange(size=0.75) +
  geom_point(colour="black", shape=21, size = 7) +
  aes(fill = factor(treatment)) 


evenPlot = even + theme_bw()   + labs(y = "Average evenness\n", fill = "Treatment") +
  theme(axis.title.x=element_blank(), 
        axis.text.y=element_text(size=13, colour="black"),
        axis.text.x=element_text(size=13, colour="black"),
        axis.title.y=element_text(size=13, vjust=1, colour="black"),
        legend.title=element_text(size=14, vjust=1),
        legend.text=element_text(size=13, vjust=0.5),
        legend.key.size=unit(1.1, "cm"),
        legend.position="top"
  )

# third plot - no of observed OTUS (richness)


AvsSob <- do %>% 
  group_by(timepoint, treatment) %>%
  summarise(av = mean(sobs), 
            sdv = sd(sobs),
            phs = first(phase),
            total = n())
# order you want timepoints to occur
AvsSob[["timepoint"]] <- setFactorOrder(AvsSob[["timepoint"]], c("T0", "T2", "T3",  "T5","T7", "T8", "T11", "T13"))

# order of treatments. 
AvsSob[["treatment"]] <- setFactorOrder(AvsSob[["treatment"]], c("Control", "Slurry", "Flood", "Flood+Slurry"))



sob = ggplot(AvsSob, aes(x=timepoint, y=av, ymin=av-sdv, ymax=av+sdv))+
  geom_pointrange(size=0.75) +
  geom_point(colour="black", shape=21, size = 7) +
  aes(fill = factor(treatment)) 


sobPlot = sob + theme_bw()   + labs(y = "Average richness\n(OTUs observed)\n", fill = "Treatment") +
  theme(axis.title.x=element_blank(), 
        axis.text.y=element_text(size=13, colour="black"),
        axis.text.x=element_text(size=13, colour="black"),
        axis.title.y=element_text(size=13, vjust=1, colour="black"),
        legend.title=element_text(size=14, vjust=1),
        legend.text=element_text(size=13, vjust=0.5),
        legend.key.size=unit(1.1, "cm"),
        legend.position="top"
  )


# =================== combine the three

library(ggpubr)


ggarrange(sobPlot, divPlot, evenPlot, common.legend = TRUE, nrow=3, align="hv")


