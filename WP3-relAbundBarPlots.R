## Total abundance plot
library(RColorBrewer)
library(ggplot2)

# read in data
AllT <- read.csv("NoSingPhyloseqOTU_RelAbund_1T0.csv", header = TRUE)

# run setFactorOrder function - save in R basics folder
# this keeps subplots in correct order when plotted

AllT[["timepoint"]] <- setFactorOrder(AllT[["timepoint"]], c("Slurry", "T0", "T2", "T3", "T5", "T7", "T8", "T11", "T13"))

AllT[["Treatment"]] <- setFactorOrder(AllT[["Treatment"]], c("Control", "Slurry", "Flood", "Flood+Slurry"))

# ================= CHOOSING COLOURS ========================================================
# using colour brewer make a palette
mypal <- colorRampPalette( brewer.pal( 11 , "Spectral" ) )
intercalate <- function(n){ #my crude attempt to shuffle the colors
  c(rbind(1:(n/2), n:(n/2+1))) #it will only work for even numbers
}


## if you want to define colours yourself, use ggthemr
library(ggthemr)

## 30 colours
abund_cols <- c("darkred", "red", "orange1","firebrick", "tomato","goldenrod1","darkorange1",
                "yellow", "chartreuse2", "darkgreen","slategray3","lightslateblue", "blueviolet", 
                "blue",  "aquamarine3","olivedrab3",   "green", "darkslateblue",
                "magenta4","darkorchid4", "orchid","mediumturquoise","navyblue","lightcoral","orangered4",
                'lightgoldenrod1', "deeppink4", "cornflowerblue","black", "grey18", "blue4", "grey",  "violet" )



# you have to add a colour at the start of your palette for outlining boxes, we'll use a grey:
abundCols <- c("#555555", abund_cols)
# remove previous effects:
ggthemr_reset()
# Define colours for your figures with define_palette

DCols <- define_palette(
  swatch = abundCols, # colours for plotting points and bars
  gradient = c(lower = abundCols[1L], upper = abundCols[2L]), #upper and lower colours for continuous colours
  background = "white" #defining a grey-ish background 
)

# set the theme for your figures:
ggthemr(DCols)
#====================================================================================
## now plot
# change as appropritae:
# x=
# y=
# fill = 
# facet_wrap(SamplingTime) # this is if you want subplots of tmt or day eg.

ggplot(AllT, aes(x=smallID, y=Abundance, fill=Phylum)) + 
  geom_bar(stat="identity") + 
  facet_wrap(~timepoint, scales="free") + 
theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust =1, colour= "black", size = 10),
        axis.text.y=element_text(colour= "black", size = 10),
        axis.title.y=element_text(colour= "black"),
        axis.title.x=element_blank(),
        legend.text=element_text(colour="black", size = 10),
        legend.title=element_text(colour="black", size = 12),
        strip.text.x = element_text(size = 14, colour = "black"),# change font of facet label
        strip.background =  element_rect(fill = "white") ) + # remove grey backgroup of facet label) +
  labs(fill = "                             Phylum\n", y = "Relative Abundance") 


# trying to make bars same width



#= ================some wrestling to make bars same width =========================


T0y2 <- subset(AllT, (timepoint %in% c("T0","T2")))
T0y2[["Treatment"]] <- setFactorOrder(T0y2[["Treatment"]], c("Control", "Slurry", "Flood", "Flood+Slurry"))

slurry = subset(AllT, (timepoint %in% c("Slurry")))

Trest <- subset(AllT, !(timepoint %in% c("Slurry", "T0","T2")))
Trest[["Treatment"]] <- setFactorOrder(Trest[["Treatment"]], c("Control", "Slurry", "Flood", "Flood+Slurry"))

#Trest$Sample = factor(Trest$Sample, level =  unique(Trest$Sample[order(Trest$Treatment)]))
#Trest$Sample = factor(Trest$Sample, level =  Trest$Sample[order(Trest$Treatment)])

# orders = c("Control", "Slurry", "Flood","Flood+Slurry")
# T3 to T13 PLOT: 
full = ggplot(Trest, aes(x=smallID, y=Abundance, fill=Phylum)) + 
  geom_bar(stat="identity") + 
  facet_wrap(~timepoint, scales="free") + 
  theme_bw() +
  #scale_x_discrete(breaks=arrange(Trest, group)$Treatment, labels=arrange(Trest, group)$Treatment) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, colour= "black", size = 11),
        axis.text.y=element_text(colour= "black", size = 10),
        axis.title.y=element_text(colour= "black"),
        axis.title.x=element_blank(),
        legend.text=element_text(colour="black", size = 12),
        legend.title=element_blank(),
        legend.position="bottom",
        strip.text.x = element_text(size = 13, colour = "black"),# change font of facet label
        strip.background =  element_rect(fill = "white") ) + # remove grey backgroup of facet label) +
  labs(y = "Relative Abundance") +  guides(fill = guide_legend(ncol = 4))

# T0 & 2 PLOT: 
T02 = ggplot(T0y2, aes(x=smallID, y=Abundance, fill=Phylum)) + 
  geom_bar(stat="identity") + 
  facet_wrap(~timepoint, scales="free") + 
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust =1,colour= "black", size = 11),
        axis.text.y=element_text(colour= "black", size = 10),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        legend.text=element_text(colour="black", size = 12),
        legend.title=element_blank(),
        strip.text.x = element_text(size = 13, colour = "black"),# change font of facet label
        strip.background =  element_rect(fill = "white") ) + # remove grey backgroup of facet label) +
  labs(fill = " ", y = "Relative Abundance") + guides(fill=FALSE) 

T02 #+ guides(fill=FALSE)

# plotting slurry bar
slu = ggplot(slurry, aes(x=Sample, y=Abundance, fill=Phylum)) + 
  geom_bar(stat="identity") + 
  #facet_wrap(~timepoint, scales="free") + 
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, colour= "black", size = 13),
        axis.text.y=element_text(colour= "black", size = 13),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        legend.text=element_text(colour="black", size = 12),
        legend.title=element_blank(),
        strip.text.x = element_text(size = 13, colour = "black"),# change font of facet label
        strip.background =  element_rect(fill = "white") ) + # remove grey backgroup of facet label) +
  labs(fill = "                             Phylum\n", y = "Relative Abundance") 

slu + guides(fill=FALSE)


## trying to put in one figure:

install.packages("ggpubr")
library("ggpubr")

ggarrange(slu, T02, full, 
          ncol = 2, nrow=2,
          common.legend=TRUE)


# == colour messings 

library(ggthemr)

## 30 colours
abund_cols <- c("black", "darkred", "red", "orange1","firebrick", "tomato","goldenrod1","darkorange1",
                "yellow", "chartreuse2", "darkgreen","lightslateblue", "blueviolet", 
                "blue", "slategray3", "olivedrab3",  "aquamarine3", "green", "darkslateblue",
                "magenta4","darkorchid4", "orchid","navyblue","mediumturquoise","lightcoral","orangered4",
                'lightgoldenrod1', "deeppink4","cornflowerblue","blue4", "grey",  "violet" )



# you have to add a colour at the start of your palette for outlining boxes, we'll use a grey:
abundCols <- c("#555555", abund_cols)
# remove previous effects:
ggthemr_reset()
# Define colours for your figures with define_palette

DCols <- define_palette(
  swatch = abundCols, # colours for plotting points and bars
  gradient = c(lower = abundCols[1L], upper = abundCols[2L]), #upper and lower colours for continuous colours
  background = "white" #defining a grey-ish background 
)

# set the theme for your figures:
ggthemr(DCols)
# now plot away!! 

                    
