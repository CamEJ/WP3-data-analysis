# Comparing means - alpha diversity WP3 data

## using ouput from 
#mothur > summary.single(shared=stability.opti_mcc.shared, calc=nseqs-coverage-sobs-simpsoneven-invsimpson-simpsoneven-chao, subsample=n)
# using shared file with singletons removed by remove.rare nseqs=1

library(dunn.test)

setwd("C:/Users/Camilla/Dropbox/Data & analysis/WP3 Slurry disturbance/Microbial Community Analyses/DNA_16S/combining all tps/No singletons")

pd = read.table("stability.opti_mcc.0.03.pick.groups.ave-std.summary", sep = "\t", header=T)
meta = read.table("metadataCombined-S.txt", sep = "\t", header=T)
rownames(meta) = meta$group

# my simpson even was is diff file to rest of alpha div metrics hence
# the second table called pd. 

d = pd[2:89,] # choose just rows with ave and cut out std rows 
              # also cut out slurry (row 1)
rownames(d) = d$group # specify row names

d$percov <- with(d, coverage *100) # make a new column with coverage expressed as %

d1 = d[c(1:6,13:88),] # cut out second T0 samples
# now make df with just the metrics you want to plot

do = subset(d1, select=c(percov,sobs,chao,simpsoneven,invsimpson))

head(do)

# =============== cut two dropped samples from metadata file


rejects = c("T7_Tmt3_a","T3_Tmt3_c")

meta1 = subset(meta, !rownames(meta) %in% rejects)
meta1 = meta1[c(1:6,13:88),] # cut out second T0 samples

# now metadata file and df with alpha div in shuld have same no of obs. 

# == add metadata to data frame of alpha metrics. 

do$phase = meta1$Phase
do$treatment = meta1$Treatment


# ======================================= coverage:

site = do$treatment # define your variables for test
# how samples are grouped into tmts

H <- do$percov # define variable
tapply(H, INDEX=site, FUN=shapiro.test)
# Some ps < 0.05 so will do kruskal

kruskal.test(H ~ site) # did kruskal too as some borderline p vals with normality test
# p-value = 0.03014

dunn.test(H, site,
          method="bh") 

# Kruskal-Wallis rank sum test
# 
# data: H and site
# Kruskal-Wallis chi-squared = 8.9373, df = 3, p-value = 0.03
# 
# 
# Comparison of H by site   (Benjamini-Hochberg)                
# Col Mean-|
# Row Mean |    Control      Flood   Flood+Sl
# ---------+---------------------------------
#   Flood |   1.509712
#         |     0.0983
#         |
# Flood+Sl|   1.754414   0.173977
#         |     0.0794     0.4309
#         |
# Slurry  |  -0.848500  -2.268634  -2.539973
#         |     0.2377     0.0349     0.0333


# ============================== sobs:

site = do$treatment

H <- do$sobs
tapply(H, INDEX=site, FUN=shapiro.test)
# if no p > 0.05 then normally distributed and can continue to aov
kruskal.test(H ~ site) # did kruskal too as some borderline p vals with normality test
# p-value = 0.002
Haov = aov(H ~ site)
summary(Haov)
dunn.test(H, site,
          method="bh") 

# Kruskal-Wallis chi-squared = 14.6849, df = 3, p-value = 0
# Comparison of H by site                            
# (Benjamini-Hochberg)                              
# Col Mean- |
# Row Mean  |    Control      Flood   Flood+Sl
# ---------+---------------------------------
# Flood   |  -2.420418
#         |    0.0155*
#         |
# Flood+Sl|  -3.445236  -0.852913
#         |    0.0017*     0.2362
#         |
# Slurry  |  -0.830318   1.677759   2.676510
#         |     0.2032     0.0700    0.0112*


# ============================== chao1:

site = do$treatment

H <- do$chao
tapply(H, INDEX=site, FUN=shapiro.test)
# one p < 0.05 so do kruskal
kruskal.test(H ~ site)  # p = 0.04349
# also check with aov out of interest as more strict
Haov = aov(H ~ site)
summary(Haov) # p is 0.09

# ============================== simp even :

site = do$treatment

H <- do$simpsoneven
tapply(H, INDEX=site, FUN=shapiro.test)
# few p < 0.05 so do kruskal
kruskal.test(H ~ site) # p = 9.608e-06
# do post hoc dunn test
dunn.test(H, site,
          method="bh") 

# Kruskal-Wallis chi-squared = 25.9846, df = 3, p-value = 0 
# Comparison of H by site   (Benjamini-Hochberg)                          
# Col Mean-|
# Row Mean |    Control      Flood   Flood+Sl
# ---------+---------------------------------
#     Flood |  -1.627625
#           |     0.0622
#           |
# Flood+Sl  |  -4.905092  -2.922411
#           |    0.0000*    0.0035*
#           |
# Slurry    |  -3.251592  -1.280687   1.894703
#           |    0.0017*     0.1002     0.0436

# ============================== simp div:

site = do$treatment

H <- do$invsimpson
tapply(H, INDEX=site, FUN=shapiro.test)
# few p < 0.05 so do kruskal
kruskal.test(H ~ site) # p = 3.841e-06
# do post hoc dunn test
dunn.test(H, site,
          method="bh") 
#ouput:
# Col Mean-|
# Row Mean |    Control      Flood   Flood+Sl
# ---------+---------------------------------
# Flood     |  -2.122270
#           |    0.0203*
#           |
# Flood+Sl  |  -5.233317  -2.755632
#           |    0.0000*    0.0059*
#           |
#   Slurry  |  -3.509793  -0.561060   2.455806
#           |    0.0040*     0.2874    0.0105*

