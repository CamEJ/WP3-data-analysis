## PAM regression analysis
## from: 

# http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0101238#s9
# Hypothetical example

#Procrustes analysis associated with ANOVA
# Hypothetical datasets
#******************************

# Soil Microbial structure raw data [X] PLFA data set

# Soil microbial functioning raw data [Y] ( denitrifiers (nirS, nirK, nosZ), nitrifiers (amoA), methanotrophic gene (pmoA) )

# Land use types (LU) factor with four levels: Original forest fragment, Silvipastoral system, improved pasture, unimproved pasture


###################################################################
###################################################################
######################################################################

# Running the principal components analyses and making PC matrices from [X] and [Y]
#*************************************************************

#  load vegan R package:

library(vegan)

#  PCA of [X] (PLFA data set):
#*****************************
sharedsubFile2 = read.table('stability.opti_mcc.0.03.subsample.shared')
sharedsubFile = t(sharedsubFile2)
rownames(sharedsubFile) = sharedsubFile[,1]
colnames(sharedsubFile) = row.names(metaFile)
dim(sharedsubFile) # [1] 87682    90
#sharedsubFile = sharedsubFile[,2:90] # all samples
sharedsubFile = sharedsubFile[4:87682,]
sharedsubFile = sharedsubFile[,c(2:8,15:90)] # removing extra T0 samples
class(sharedsubFile) <- "numeric"
head(sharedsubFile)
sharedsubFile = sharedsubFile[,2:83] # cut out slurry sample

# for now, I Kknow its a little random, but will cut to 5000 OTUs
dim(sharedsubFile)
sharedsubFileSMall = sharedsubFile[1:5000,]

#now transform as I think they want samples as rows
sharedsubFileSMall = t(sharedsubFileSMall)
class(sharedsubFileSMall) <- "numeric"
X = sharedsubFileSMall

head(sharedsubFileSMall)
X.log<- log(X+1) # transformation
X.pca<- rda(X.log)

# Extracting and obtaining principal components of X

X.2axes<- scores(X.pca, display = c("sites"), choice = c(1,2)) # 2 PCA axes 
X.2axes<-as.matrix(X.2axes) # 2 axes matrix


# PCA of [Y] ( denitrifiers (nirS, nirK, nosZ), nitrifiers (amoA), methanotrophic gene (pmoA) ))
# In my case this is microresp data - undtransformed and with correct samples cut out. 
# to match miseq samples
#********************************
Tmt = read.table('PAM-treatment.txt', header=T, sep='\t')
rownames(Tmt) = Tmt[,1]
Tmt = Tmt[,2]
Tmt = as.data.frame(Tmt)
class(LU) <- "numeric"
LU = Tmt


#*****************************************
Mresp = read.table('MicroRespForPAM.txt', header=T, sep='\t')
rownames(Mresp) = Mresp[,1]
Mresp = Mresp[,2:12]
Y = Mresp

Y.log<- log(Y+1)
Y.pca <-rda(Y.log)

# Extracting and obtaining principal components matrices of Y


Y.2axes<- scores(Y.pca, display = c("sites"), choice = c(1,2)) # 2 PCA axes 
Y.2axes<-as.matrix(Y.2axes) # 2 axes matrix

#####################################
####################################
####################################
########################################
##########################################


# Run the Procrustes relationships between X and Y axes matrices
#***************************************************************

# Between 2 axes matrices (Figure 6a. Main text)

bet2axes<-procrustes(X.2axes,Y.2axes) # Notice both matrices must be the same number of columns


#######################################################
########################################################
########################################################

# Obtaining the PAMs (Procrustes association metric)

PAM2axes<-residuals(bet2axes) # Relationship between SMC-SMF based on  2 axes matrices


######################################
############################################
################################################

# Run an ANOVA with Land use (LU) as fixed factor and Soil microbial structure-soil microbial functioning, i.e PAM2axes as response
#********************************************************************************************************************************
datafr<-data.frame(cbind(PAM2axes,LU)) 
                   
                   ANOVA.proc<-aov(PAM2axes~Tmt,datafr)
                   
                   ANOVA.proc # In our hypothetical example, the one-way ANOVA output shows that at least two types of land 
                   #use differ statistically (F = 7.047, P = 0.003)
                   
                   #########################################################3
                   ################################
                   ###################################
                   ##################################
                   # Call:
                   #   aov(formula = PAM2axes ~ Tmt, data = datafr)
                   # 
                   # Terms:
                   #   Tmt Residuals
                   # Sum of Squares   35.44655 110.72643
                   # Deg. of Freedom         3        78
                   # 
                   # Residual standard error: 1.191457
                   # Estimated effects may be unbalanced
                   summary(ANOVA.proc)
                   # Df Sum Sq Mean Sq F value   Pr(>F)    
                   # Tmt          3  35.45   11.82   8.323 7.19e-05 ***
                   #   Residuals   78 110.73    1.42                     
                   # ---
                   #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
                   # Running a multiple comparisons mean test to assess different land use types in regard to their effect 
                   # on the relationship between soil microbial structure and functioning.
                   #****************************************************************************
                   
                   
                   #load agricolae R package
                   library(agricolae)
                   HSD.test(ANOVA.proc,"Tmt",group=TRUE) # Tukey 95%
                   # didnt get their HSD to work, just used normal one I use
                   TukeyHSD(ANOVA.proc)
                   #Output:
                   # $Tmt
                   # diff        lwr       upr     p adj
                   # Flood-Control         0.19382448 -0.8157081 1.2033571 0.9578891
                   # Flood+Slurry-Control  1.39505862  0.4197575 2.3703597 0.0018517
                   # Slurry-Control        1.37941237  0.4764590 2.2823658 0.0007813
                   # Flood+Slurry-Flood    1.20123415  0.1265047 2.2759636 0.0223225
                   # Slurry-Flood          1.18558790  0.1760553 2.1951205 0.0147303
                   # Slurry-Flood+Slurry  -0.01564625 -0.9909474 0.9596549 0.9999729
                   
                   