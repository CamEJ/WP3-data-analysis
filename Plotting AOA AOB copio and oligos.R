
# =========== getting out and plotting specific phylotypes ===========
# NOTE this is using 0.005% trimmed OTU table. 
# as recommended by Bokulvich et al. 

yy <- read.csv("Trimmed.0.05.phyloSeq.csv", header = TRUE)
head(yy)
# run setFactorFunction.R
# function for this saved in github R basics folder
yy[["timepoint"]] <- setFactorOrder(yy[["timepoint"]], c("T0", "T2", "T3", "T5", "T7", "T8", "T11", "T13"))

yy[["Treatment"]] <- setFactorOrder(yy[["Treatment"]], c("Control", "Slurry", "Flood", "Flood+Slurry"))

# =========================== AOB ==============================
#genera: Nitrosomonas, Nitrosococcus, Nitrobacter and Nitrococcus

## look Nitrosomonas
Nitrosomonas = subset_taxa(SubRel, Genus=="Nitrosomonas(100)")
# none

Nitrobacter = subset_taxa(SubRel, Genus=="Nitrobacter(100)")
# none 
# ====================== Nitrospira genus ========================== 
SubRelNitroSpira = subset_taxa(SubRel, Genus=="Nitrospira(100)")

dat1 <- psmelt(SubRelNitroSpira) #write out so can use ggplot
write.csv(dat1, file='NitrospiraGenus.csv')
# read back in
NP <- read.csv("NitrospiraGenus.csv", header = TRUE)
# function for this saved in github R basics folder
NP[["timepoint"]] <- setFactorOrder(NP[["timepoint"]], c("Slurry", "T0", "T2", "T3", "T5", "T7", "T8", "T11", "T13"))
NP[["Treatment"]] <- setFactorOrder(NP[["Treatment"]], c("Control", "Slurry", "Flood", "Flood+Slurry"))
NP[["Phase"]] <- setFactorOrder(NP[["Phase"]], c("PreFlood", "Flood", "Recovery"))

abund2 <- NP %>% 
  group_by(smallID) %>%
  summarise(abund = sum(Abundance), 
            Tp = first(timepoint),
            Tmt = first(Treatment),
            Phase = first(Phase),
            total = n())


nsp <- ggplot(abund2, aes(x=Tp,y=abund, fill = Tmt, colour="black")) +
  geom_boxplot(colour="black")
nsp

# make it pretty(ier)

nsp2 <- nsp + labs(fill="Treatment") + # update legend titles.
  facet_wrap(~Phase, scales = "free_x") +
  ylab(expression(paste(italic(" Nitrospira "), "spp." ))) +
  theme_bw()  +
  theme(axis.text.x=element_text(size=13, colour="black")) +
  theme(axis.text.y=element_text(size=12, colour="black")) +
  theme(axis.title.y=element_text(size=14, colour="black")) +
  theme(legend.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14)) +
  guides(color = FALSE) +
  #guides(fill = guide_legend(override.aes = list(size=1))) +
  theme(legend.key.size = unit(1.5, "cm"),
        #strip.text.x = element_text(size = 14, colour = "black"),# change font of facet label
        #strip.background =  element_rect(fill = "white"),
        strip.text.x = element_blank(),
        strip.background =  element_blank(),
        axis.title.x = element_blank()) 
nitros = nsp2

## ============================ look Nitrosospira genus =================================
microSubRelNitrosopira = subset_taxa(SubRel, Genus=="Nitrosospira(100)")

#write out so can use ggplot
dat <- psmelt(microSubRelNitrosopira)
write.csv(dat, file='NitrosospiraGenus.csv')
# read back in
NSP <- read.csv("NitrosospiraGenus.csv", header = TRUE)

# set factor orders
# function for this saved in github R basics folder
NSP[["timepoint"]] <- setFactorOrder(NSP[["timepoint"]], c("Slurry", "T0", "T2", "T3", "T5", "T7", "T8", "T11", "T13"))

NSP[["Treatment"]] <- setFactorOrder(NSP[["Treatment"]], c("Control", "Slurry", "Flood", "Flood+Slurry"))
NSP[["Phase"]] <- setFactorOrder(NSP[["Phase"]], c("PreFlood", "Flood", "Recovery"))


# now make new df called abundtot
# grouping samples by name, adding up Abundance data
# and putting first occurence of SamplingTime
# and Treatment.A as will be same for each group in my case
library(dplyr)
abundTot <- NSP %>% 
  group_by(smallID) %>%
  summarise(abund = sum(Abundance), 
            Tp = first(timepoint),
            Tmt = first(Treatment),
            Phase = first(Phase),
            total = n())


nsp <- ggplot(abundTot, aes(x=Tp,y=abund, fill = Tmt, colour="black")) +
  geom_boxplot(colour="black")
nsp

# make it pretty(ier)

nsp2 <- nsp + labs(fill="Treatment") + # update legend titles.
  facet_wrap(~Phase, scales = "free_x") +
  ylab(expression(paste(italic(" Nitrosospira "), "spp." ))) +
  theme_bw()  +
  theme(axis.text.x=element_text(size=13, colour="black")) +
  theme(axis.text.y=element_text(size=12, colour="black")) +
  theme(axis.title.y=element_text(size=14, colour="black")) +
  theme(legend.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14)) +
  guides(color = FALSE) +
  #guides(fill = guide_legend(override.aes = list(size=1))) +
  theme(legend.key.size = unit(1.5, "cm"),
        #strip.text.x = element_text(size = 14, colour = "black"),# change font of facet label
        #strip.background =  element_rect(fill = "white"),
        strip.text.x = element_blank(),
        strip.background =  element_blank(),
        axis.title.x = element_blank()) 
nitroso = nsp2

# ================== AOA ===================================
# as many of these were unclassified and i found them searching 
# rep seqs on silva db, I have to subset on OTU id

OTU_thaum <- subset(yy, (OTU %in% c("Otu000007","Otu000016","Otu000023", "Otu000026")))


# melt to readout
write.csv(OTU_thaum, file='ThaumarchaetoaSpp.csv') # write out
# read back in
OTU_thaum <- read.csv("ThaumarchaetoaSpp.csv", header = TRUE)

bet = OTU_thaum

bet[["timepoint"]] <- setFactorOrder(bet[["timepoint"]], c("Slurry", "T0", "T2", "T3", "T5", "T7", "T8", "T11", "T13"))
bet[["Treatment"]] <- setFactorOrder(bet[["Treatment"]], c("Control", "Slurry", "Flood", "Flood+Slurry"))
bet[["Phase"]] <- setFactorOrder(bet[["Phase"]], c("PreFlood", "Flood", "Recovery"))

# now make new df called abundtot
# grouping samples by name, adding up Abundance data
# and putting first occurence of SamplingTime

abundTot <- bet %>% 
  group_by(Sample) %>%
  summarise(abund = sum(Abundance), 
            Tp = first(timepoint),
            Tmt = first(Treatment),
            Phase = first(Phase),
            total = n())

Ac <- ggplot(abundTot, aes(x=Tp,y=abund, fill = Tmt, colour="black")) +
  geom_boxplot(colour="black")
Ac

# make it pretty(ier)

Ac2 <- Ac + labs(fill="Treatment") + # update legend titles.
  facet_wrap(~Phase, scales = "free_x") +
  ylab(expression(paste(italic(" Thaumarchaeota"), "spp."))) +
  theme_bw()  +
  theme(axis.text.x=element_text(size=13, colour="black")) +
  theme(axis.text.y=element_text(size=12, colour="black")) +
  theme(axis.title.y=element_text(size=14, colour="black")) +
  theme(legend.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14)) +
  guides(color = FALSE) +
  #guides(fill = guide_legend(override.aes = list(size=1))) +
  theme(legend.key.size = unit(1.5, "cm"),
        strip.text.x = element_text(size = 14, colour = "black"),# change font of facet label
        strip.background =  element_rect(fill = "white"),
        axis.title.x = element_blank()) 
Ac2


# =========== put the three in one plot ===========
library("ggpubr")

ggarrange(Ac2, nitros, nitroso, common.legend = TRUE, nrow=3,align="hv")



#====================copio and oligotrophs ======================

library(dplyr)
# Acidobacteria(100)
## look phylum Acidobacteria(100)
microSubRelAcido = subset_taxa(SubRel, Phylum=="Acidobacteria(100)")
BarAcido <- plot_bar(microSubRelAcido, fill="Class", title="Acidobacteria")
## changes scale/day so cant see difference really so just use all together
BarAcido + facet_wrap(~timepoint, scales="free") + theme_bw() +
  theme(axis.text.x=element_text(angle=90, hjust=1, size=9, vjust=0.4))

dat <- psmelt(microSubRelAcido)
write.csv(dat, file='AcidobacteriaPhylum.csv')
# read back in
acid <- read.csv("AcidobacteriaPhylum.csv", header = TRUE)

acid[["timepoint"]] <- setFactorOrder(acid[["timepoint"]], c("Slurry", "T0", "T2", "T3", "T5", "T7", "T8", "T11", "T13"))
acid[["Treatment"]] <- setFactorOrder(acid[["Treatment"]], c("Control", "Slurry", "Flood", "Flood+Slurry"))
acid[["Phase"]] <- setFactorOrder(acid[["Phase"]], c("PreFlood", "Flood", "Recovery"))

# now make new df called abundtot
# grouping samples by name, adding up Abundance data
# and putting first occurence of SamplingTime

abundTot <- acid %>% 
  group_by(Sample) %>%
  summarise(abund = sum(Abundance), 
            Tp = first(timepoint),
            Tmt = first(Treatment),
            Phase = first(Phase),
            total = n())

Ac <- ggplot(abundTot, aes(x=Tp,y=abund, fill = Tmt, colour="black")) +
  geom_boxplot(colour="black")
Ac

# make it pretty(ier)

Ac2 <- Ac + labs(fill="Treatment") + # update legend titles.
  facet_wrap(~Phase, scales = "free_x") +
  ylab(expression(paste(italic(" Acidobacteria")))) +
  theme_bw()  +
  theme(axis.text.x=element_text(size=13, colour="black")) +
  theme(axis.text.y=element_text(size=12, colour="black")) +
  theme(axis.title.y=element_text(size=13, colour="black")) +
  theme(legend.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14)) +
  guides(color = FALSE) +
  #guides(fill = guide_legend(override.aes = list(size=1))) +
  theme(legend.key.size = unit(1.5, "cm"),
        strip.text.x = element_text(size = 14, colour = "black"),# change font of facet label
        strip.background =  element_rect(fill = "white"),
        axis.title.x = element_blank()) 

Ac2

acido = Ac2

# =============================================== Bacteroidetes(100) ========================
## look phylum Bacteroidetes(100)
microSubRelBacto = subset_taxa(SubRel, Phylum=="Bacteroidetes(100)")

dat <- psmelt(microSubRelBacto)
write.csv(dat, file='BacteriodetesPhylum.csv')
# read back in
bac <- read.csv("BacteriodetesPhylum.csv", header = TRUE)

#Bacteroidia = subset(bact, Class=="Bacteroidia(100)" )
#rest = subset(bact, (Class %in% c("Flavobacteria(100) ","Sphingobacteria(100)")))

#bac = Bacteroidia
#bac= rest
bac[["timepoint"]] <- setFactorOrder(bac[["timepoint"]], c("Slurry", "T0", "T2", "T3", "T5", "T7", "T8", "T11", "T13"))
bac[["Treatment"]] <- setFactorOrder(bac[["Treatment"]], c("Control", "Slurry", "Flood", "Flood+Slurry"))
bac[["Phase"]] <- setFactorOrder(bac[["Phase"]], c("PreFlood", "Flood", "Recovery"))

# now make new df called abundtot
# grouping samples by name, adding up Abundance data
# and putting first occurence of SamplingTime

abundTot <- bac %>% 
  group_by(Sample) %>%
  summarise(abund = sum(Abundance), 
            Tp = first(timepoint),
            Tmt = first(Treatment),
            Phase = first(Phase),
            total = n())

Ac <- ggplot(abundTot, aes(x=Tp,y=abund, fill = Tmt, colour="black")) +
  geom_boxplot(colour="black")
Ac

# make it pretty(ier)

Ac2 <- Ac + labs(fill="Treatment") + # update legend titles.
  facet_wrap(~Phase, scales = "free_x") +
  ylab(expression(paste(italic(" Bacteroidetes ") ))) +
  theme_bw()  +
  theme(axis.text.x=element_text(size=13, colour="black")) +
  theme(axis.text.y=element_text(size=12, colour="black")) +
  theme(axis.title.y=element_text(size=13, colour="black")) +
  theme(legend.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14)) +
  guides(color = FALSE) +
  #guides(fill = guide_legend(override.aes = list(size=1))) +
  theme(legend.key.size = unit(1.5, "cm"),
        #strip.text.x = element_text(size = 14, colour = "black"),# change font of facet label
        #strip.background =  element_rect(fill = "white"),
        strip.text.x = element_blank(),
        strip.background =  element_blank(),
        axis.title.x = element_blank()) 
Ac2

bacto = Ac2

FyS = Ac2



# ================ Betaproteobacteria ===================================
## look class Betaproteobacteria(100)
beta = subset(yy, Class=="Betaproteobacteria(100)")
bet = beta
microSubRelBetaP = subset_taxa(SubRel, Class=="Betaproteobacteria(100)")
dat <- psmelt(microSubRelBetaP) # melt to readout
write.csv(dat, file='BetaproteobacteriaClass.csv') # write out
# read back in
bet <- read.csv("BetaproteobacteriaClass.csv", header = TRUE)

bet[["timepoint"]] <- setFactorOrder(bet[["timepoint"]], c("Slurry", "T0", "T2", "T3", "T5", "T7", "T8", "T11", "T13"))
bet[["Treatment"]] <- setFactorOrder(bet[["Treatment"]], c("Control", "Slurry", "Flood", "Flood+Slurry"))
bet[["Phase"]] <- setFactorOrder(bet[["Phase"]], c("PreFlood", "Flood", "Recovery"))

# now make new df called abundtot
# grouping samples by name, adding up Abundance data
# and putting first occurence of SamplingTime

abundTot <- bet %>% 
  group_by(Sample) %>%
  summarise(abund = sum(Abundance), 
            Tp = first(timepoint),
            Tmt = first(Treatment),
            Phase = first(Phase),
            total = n())

Ac <- ggplot(abundTot, aes(x=Tp,y=abund, fill = Tmt, colour="black")) +
  geom_boxplot(colour="black")
Ac

# make it pretty(ier)

Ac2 <- Ac + labs(fill="Treatment") + # update legend titles.
  facet_wrap(~Phase, scales = "free_x") +
  ylab(expression(paste(italic(" Betaproteobacteria")))) +
  theme_bw()  +
  theme(axis.text.x=element_text(size=13, colour="black")) +
  theme(axis.text.y=element_text(size=11, colour="black")) +
  theme(axis.title.y=element_text(size=13, colour="black")) +
  theme(legend.text = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14)) +
  guides(color = FALSE) +
  #guides(fill = guide_legend(override.aes = list(size=1))) +
  theme(legend.key.size = unit(1.5, "cm"),
        #strip.text.x = element_text(size = 14, colour = "black"),# change font of facet label
        #strip.background =  element_rect(fill = "white"),
        strip.text.x = element_blank(),
        strip.background =  element_blank(),
        axis.title.x = element_blank()) 
Ac2

betaPro = Ac2

# =========== in one 
library("ggpubr")

ggarrange(acido, bacto, betaPro, common.legend = TRUE, nrow=3,align="hv")
