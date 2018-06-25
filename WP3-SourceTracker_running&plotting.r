# Running SourceTracker version 1.0.1
# https://github.com/danknights/sourcetracker
# Download tar.gz and then extract this and set folder as wkdir
# Open example.r and run script the first time to check everything works 
# correctly using their example data. 

# ================== getting data ready ====================================# 

# you will need 1) a metadatafile and 2) an OTU table

# 1. metadata file formate
# this must have:
  # a) a column called 'SourceSink' saying if your sample is a source or a sink.
  # b) a column called 'Env' saying which environment your sample is from
  # c) a column called Description which you can use to describe your samples
  # d) any samples you don't now want to be analysed as a source or sink should 
  # have NA in these 3 columns. 

# head(metadata)
# SlurryAdd   Phase Treatment SourceSink    Env Description
# Slurry1        <NA>    <NA>      <NA>     source slurry      slurry
# Slurry2        <NA>    <NA>      <NA>     source slurry      slurry
# T0_Tmt1_a     Minus NoFlood   Control     source   soil        soil
# T0_Tmt1_b     Minus NoFlood   Control     source   soil        soil
# T0_Tmt1_c     Minus NoFlood   Control     source   soil        soil
# T0_Tmt2_a      Plus NoFlood    Slurry       sink   tmt2    0_slurry

# 2. OTU table.

# This needs to have same sample names as in your metadata file
# They suggest using trimmed OTU to remove low abundance OTUs
# I ran this with OTUs trimmed to > 0.005% as per Bokulich et al 2013. 
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3531572/

# The script below expects that initially your samples will be columns and 
# OTUs will be rows (& they will therefore be transformed). 
# your OTU table should not include any taxo information. 


# ====================== loading data & running scrips ================================

# load sample metadata
metadata <- read.table('metadata2slurry.txt',sep='\t',h=T,row.names=1)

# load OTU table
# This 'read.table' command is designed for a 
# QIIME-formatted OTU table.
# namely, the first line begins with a '#' sign
# and actually _is_ a comment; the second line
# begins with a '#' sign but is actually the header
otus <- read.table('OTUtable_2slurry.txt',sep='\t', header=T,row.names=1,check=F,skip=1,comment='')
head(otus)
dim(otus)
otus2 = otus[,c(1:8,15:90)] # removing extra T0 samples & also taxo info as not in their eg data
head(otus2) # check
otus <- t(as.matrix(otus2))

# extract only those samples in common between the two tables
common.sample.ids <- intersect(rownames(metadata), rownames(otus))
otus <- otus[common.sample.ids,]
metadata <- metadata[common.sample.ids,]
# double-check that the mapping file and otu table
# had overlapping samples
if(length(common.sample.ids) <= 1) {
    message <- paste(sprintf('Error: there are %d sample ids in common '),
                    'between the metadata file and data table')
    stop(message)
}

# extract the source environments and source/sink indices
train.ix <- which(metadata$SourceSink=='source')
test.ix <- which(metadata$SourceSink=='sink')
envs <- metadata$Env
if(is.element('Description',colnames(metadata))) desc <- metadata$Description


# load SourceTracker package
source('src/SourceTracker.r')

# tune the alpha values using cross-validation (this is slow!)
# tune.results <- tune.st(otus[train.ix,], envs[train.ix])
# alpha1 <- tune.results$best.alpha1
# alpha2 <- tune.results$best.alpha2
# note: to skip tuning, run this instead:
alpha1 <- alpha2 <- 0.001

# train SourceTracker object on training data
st <- sourcetracker(otus[train.ix,], envs[train.ix])

# Error in sum(x[i, ]) : invalid 'type' (character) of argument
# I initially got this error when I had only one source sample.
# added a second and it worked a charm

# Estimate source proportions in test data (read note below before doing this, re 'full results'
results <- predict(st,otus[test.ix,], alpha1=alpha1, alpha2=alpha2)

# if you want to get OTU ids need to add full.results=TRUE argument to this command. 
#https://groups.google.com/forum/#!msg/qiime-forum/2Xn1niCHWaI/ob37bsE5UeUJ;context-place=topic/qiime-forum/EstvrrRmDfY
results <- predict(st,otus[test.ix,], alpha1=alpha1, alpha2=alpha2, full.results=TRUE)

# Estimate leave-one-out source proportions in training data 
results.train <- predict(st, alpha1=alpha1, alpha2=alpha2)

# plot results
labels <- sprintf('%s %s', envs,desc)
plot(results, labels[test.ix], type='pie')

# other plotting functions
# plot(results, labels[test.ix], type='bar')
plot(results, labels[test.ix], type='dist')
#plot(results.train, labels[train.ix], type='pie')
# plot(results.train, labels[train.ix], type='bar')
# plot(results.train, labels[train.ix], type='dist')

# plot results with legend
plot(results, labels[test.ix], type='pie', include.legend=TRUE, env.colors=c('#47697E','#5B7444','#CC6666','#79BEDB','#885588'))



# ===================reading data out to plot more clearly======================

PieNos = results$proportions

PieNos = as.data.frame(PieNos)

# in order to group and summarise samples (ie average) we will add metadata back
# and this will allow grouping of samples
  metadata$id = row.names(metadata) 
  PieNos$id = row.names(PieNos)
  zz <- merge(PieNos, metadata, by = "id")

  require(dplyr)
  # grouping and summarising variables within a df called 'zz'
  # to get mean and sd of each
  
              PieAvs <- zz %>% # choose name of new df (here it's PieAvs) to which new data will be put
                group_by(Description) %>% # group by sample
                summarise(Slurry = mean(slurry),  # then use summarise to work out mean etc for each diff group
                          Soil = mean(soil), # where each of these will be a new column in new df 
                          unknown = mean(Unknown),
                          time = first(timepoint),
                          treatment = first(Treatment),
                          varSlurry = sd(slurry), # and cols of sdev
                          varSoil = sd(soil), 
                          varUnknown = sd(Unknown))
 
              write.csv(PieAvs, "averagesOutPutPiesSourceTracker.csv")
              
# ============================= and now to plotting =============================== # 
setwd("~/Desktop/sourcetracker-1.0.1")
pies <- read.csv('averagesOutPutPiesSourceTracker.csv',h=T,row.names=1)

library(ggplot2) 
library(reshape2)
              
# df = pies[,1:6]  
# mdf = as.data.frame(df)
# df[["time"]] <- setFactorOrder(df[["time"]], c("T0", "T2", "T3", "T5", "T7", "T8", "T11", "T13"))
# 
# row.names(mdf)=mdf$Description
# 
# mdf[["time"]] <- setFactorOrder(mdf[["time"]], c("T0", "T2", "T3", "T5", "T7", "T8", "T11", "T13"))
# 
# mmdf = melt(mdf)
# mmdf[["time"]] <- setFactorOrder(mmdf[["time"]], c("T0", "T2", "T3", "T5", "T7", "T8", "T11", "T13"))


pie <- read.table('PieDataordered.txt',sep='\t',h=T)

pie[["time"]] <- setFactorOrder(pie[["time"]], c("T0", "T2", "T3", "T5", "T7", "T8", "T11", "T13"))
# set factor order of id based on order of time
pie$id = factor(pie$id, level =  unique(pie$id[order(pie$time)]))
# =============== subset by treatment to make one plot for each tmt ============== #
library(dplyr)
PS = subset(pie, (treatment %in% c("Slurry")))
PS$Day = factor(PS$Day)

PSF =  subset(pie, !(treatment %in% c("Slurry")))
PSF$Day = factor(PSF$Day)

# set labels
myLabs = c("0 ", "3 ", "6 ", "15 ", "29 ",
           "92 ", "114", "140")

# plot slurry (tmt2)

Splot = ggplot(PS, aes(x=Day, y=prop*100, fill=source)) + # time prop by 100 to give %
  geom_bar(stat="identity") + # add bars
  geom_text( # adding labels to end of bars
    aes(label = c("26.9%", "", "", "24.3%","", "", " 2.8%", "", ""," 1.0%", "", ""," 0.9%", "", ""," 0.3%","", "", " 0.1%", "", "",  " 0.2%", "", "")),    
    size = 5, hjust = -0.05,  position = "stack", inherit.aes = TRUE, colour="black") +
  coord_flip()+ # flip axes
  theme_classic() + #change theme
  theme(axis.text.x=element_text(size=15, colour="black"),
        axis.text.y=element_text(size=15, colour="black"),
        axis.title.y=element_text(colour="black", size = 15),
        legend.text=element_text(colour="black", size = 15),
        legend.title=element_text(colour="black", size = 16),
        axis.title.x=element_text(colour="black", size = 15),
        legend.key.size = unit(1.2, "cm"),
        legend.position = "top",
        plot.margin = unit(c(0.5,1.5,0.5,0.5), "cm"))+
  labs(y = "\nPercentage abundance attributed to source", x="Days post slurry amendment", fill="                    Source\n ")+
  scale_x_discrete(labels= myLabs) +
  guides(fill = guide_legend(title.position = "top"))

Splot
# percentage labels are hidden by clipping so now
# use grid to turn off clipping
library(grid)
gt <- ggplot_gtable(ggplot_build(Splot))
gt$layout$clip[gt$layout$name == "panel"] <- "off"
grid.draw(gt)


# ============= plot flood + slurry data now =================

Flood = ggplot(PSF, aes(x=Day, y=prop*100, fill=source)) + 
  geom_bar(stat="identity") + 
 geom_text(
  aes(label = c(rep("",6), "12.5%", "", "", " 5.4%","", "", " 2.4%", "", ""," 0.6%", "", ""," 0.2%", "", ""," 0.4%","", "")),
  size = 5, hjust = -0.05,  position = "stack", inherit.aes = TRUE, colour="black") +
  coord_flip()+ # flip axes
  theme_classic() + #change theme
  theme(axis.text.x=element_text(size=15, colour="black"),
        axis.text.y=element_text(size=15, colour="black"),
        axis.title.y=element_text(colour="black", size = 15),
        legend.text=element_text(colour="black", size = 15),
        legend.title=element_text(colour="black", size = 16),
        axis.title.x=element_text(colour="black", size = 15),
        legend.key.size = unit(1.2, "cm"),
        legend.position = "top",
        plot.margin = unit(c(0.5,1.5,0.5,0.5), "cm"))+
  labs(y = "\nPercentage abundance attributed to source", x="Days post slurry amendment", fill="                    Source\n ")+
  scale_x_discrete(labels= myLabs) +
  guides(fill = guide_legend(title.position = "top"))

Flood

library(grid)
gt <- ggplot_gtable(ggplot_build(Flood))
gt$layout$clip[gt$layout$name == "panel"] <- "off"
grid.draw(gt)




# =================== specifying colours using ggthemr =================================#
library(ggthemr)

WP3_colsA <- c("yellow4", "chocolate4", "slateblue", "olivedrab")
# add background white (#55555)
WP3_cols <- c("#555555", WP3_colsA)

# define palette
WP3Cols <- define_palette(
  swatch = WP3_cols, # colours for plotting points and bars
  gradient = c(lower = WP3_cols[1L], upper = WP3_cols[2L]), #upper and lower colours for continuous colours
  background = "white" #defining a grey-ish background 
)

# set new col theme as current

ggthemr(WP3Cols)   
