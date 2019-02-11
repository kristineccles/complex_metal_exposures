############################################################################
# Comparison of transformation and range normalization
# Written R version 3.4.3
###########################################################################
# Load Libraries
library(Hmisc)
library(ade4)
library (adegraphics)
library(lattice)
library(sp)
library(adegenet)
library(spdep)
library(adespatial)
library(rgdal)
library(factoextra)
library(ggplot2)

# Set up colour pallet 
coul = c("#8dd3c7", "#ffd92f", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69")
#########################################################################
## RAW DATA ##
#Import Data
complete_data = read.csv("grid_complete_data.csv")
hist(complete_data[,10:31])

amphib = subset(complete_data, Species == "Amphibian")
gull =  subset(complete_data, Species == "Gull"| Species == "Tern")
mammal =  subset(complete_data, Species == "Fisher" | Species == "Marten"| Species == "Mink" | Species == "Muskrat"| Species == "River Otter")

# EDA of original data
hist(amphib[,10:31], nclass = 10)
hist(gull[,10:31], nclass = 10)
hist(mammal[,10:31], nclass = 10)

##########################################################################
#### Transform Data ####
#Sqrt Tranform
sqrt_df = sqrt(complete_data[,10:31])
summary(sqrt_df)
hist(sqrt_df)
sqrt_data = cbind(complete_data[,1:9], sqrt_df)

# Histogram by speices
amphib2 = subset(sqrt_data, Species == "Amphibian")
gull2 =  subset(sqrt_data, Species == "Gull"| Species == "Tern")
mammal2 =  subset(sqrt_data, Species == "Fisher" | Species == "Marten"| Species == "Mink" | Species == "Muskrat"| Species == "River Otter")

hist(amphib2[,10:31], nclass = 10)
hist(gull2[,10:31], nclass = 10)
hist(mammal2[,10:31], nclass = 10)

# Log10 Transform
log10_df = log10(complete_data[,10:31]+1)
summary(log10_df)
hist(log10_df)
log10_data = cbind(complete_data[,1:9], log10_df)

# Histogram by speices
amphib3 = subset(log10_data, Species == "Amphibian")
gull3 =  subset(log10_data, Species == "Gull"| Species == "Tern")
mammal3 =  subset(log10_data, Species == "Fisher" | Species == "Marten"| Species == "Mink" | Species == "Muskrat"| Species == "River Otter")

hist(amphib3[,10:31], nclass = 10)
hist(gull3[,10:31], nclass = 10)
hist(mammal3[,10:31], nclass = 10)

##################################################################################
#### Normalize between 0 and 1 data by speices and lifestage ####
# Untransformed 
fac = interaction(factor(complete_data$Species),factor(complete_data$Lifestage), drop = T)
fnorm.var = function(y){ return((y - min(y)) / (max(y) - min(y)))}
fnorm.tab = function(tab) {apply(tab, 2, fnorm.var)}
range_norm_df=unsplit(lapply(split(as.data.frame(complete_data[,10:31]), fac), 
                             FUN = function(x) as.data.frame(fnorm.tab(x))), fac)
summary(range_norm_df)
# Replace NA's with zeros
range_norm_df=rapply( range_norm_df, f=function(x) ifelse(is.nan(x),0,x), how="replace")
complete=as.data.frame(cbind(complete_data[,1:9],range_norm_df))

# Histogram by speices
amphib1 = subset(sqrt, Species == "Amphibian")
gull1 =  subset(sqrt, Species == "Gull"| Species == "Tern")
mammal1 =  subset(sqrt, Species == "Fisher" | Species == "Marten"| Species == "Mink" | Species == "Muskrat"| Species == "River Otter")

hist(amphib1[,10:31], nclass = 10)
hist(gull1[,10:31], nclass = 10)
hist(mammal1[,10:31], nclass = 10)

# SQRT 
fac = interaction(factor(sqrt_data$Species),factor(sqrt_data$Lifestage), drop = T)
fnorm.var = function(y){ return((y - min(y)) / (max(y) - min(y)))}
fnorm.tab = function(tab) {apply(tab, 2, fnorm.var)}
range_norm_df=unsplit(lapply(split(as.data.frame(sqrt_data[,10:31]), fac), 
                             FUN = function(x) as.data.frame(fnorm.tab(x))), fac)
summary(range_norm_df)
# Replace NA's with zeros
range_norm_df=rapply( range_norm_df, f=function(x) ifelse(is.nan(x),0,x), how="replace")
sqrt=as.data.frame(cbind(sqrt_data[,1:9],range_norm_df))

# Histogram by speices
amphib4 = subset(sqrt, Species == "Amphibian")
gull4 =  subset(sqrt, Species == "Gull"| Species == "Tern")
mammal4 =  subset(sqrt, Species == "Fisher" | Species == "Marten"| Species == "Mink" | Species == "Muskrat"| Species == "River Otter")

hist(amphib4[,10:31], nclass = 10)
hist(gull4[,10:31], nclass = 10)
hist(mammal4[,10:31], nclass = 10)

# LOG10
fac = interaction(factor(log10_data$Species),factor(log10_data$Lifestage), drop = T)
fnorm.var = function(y){ return((y - min(y)) / (max(y) - min(y)))}
fnorm.tab = function(tab) {apply(tab, 2, fnorm.var)}
range_norm_df=unsplit(lapply(split(as.data.frame(log10_data[,10:31]), fac), 
                             FUN = function(x) as.data.frame(fnorm.tab(x))), fac)
summary(range_norm_df)
# Replace NA's with zeros
range_norm_df=rapply( range_norm_df, f=function(x) ifelse(is.nan(x),0,x), how="replace")
log10=as.data.frame(cbind(log10_data[,1:9],range_norm_df))

#Plot log10 Transformed data
amphib4 = subset(log10, Species == "Amphibian")
gull4 =  subset(log10, Species == "Gull"| Species == "Tern")
mammal4 =  subset(log10, Species == "Fisher" | Species == "Marten"| Species == "Mink" | Species == "Muskrat"| Species == "River Otter")

hist(amphib4[,10:31], nclass = 10)
hist(gull4[,10:31], nclass = 10)
hist(mammal4[,10:31], nclass = 10)

#########################################################################################
#### PCA and BCA on Untransformed Unnormalized data ####
# Run the PCA
pca1 = dudi.pca(df = complete_data[,10:31], scannf = FALSE, nf = 2)
g1_1=s.arrow(pca1$co, plot=FALSE, plabels=list(cex=0.75,optim=TRUE),psub.text = "Untransformed Unnormalized")
g2_1=s.class(pca1$li,complete_data$Species, col=coul, plabels=list(optim=TRUE),psub.text = "Untransformed Unnormalized")
ADEgS(list(g1_1,g2_1))

# BCA: Site
bca_site = bca(x = pca1, fac = complete_data$Site, scannf = FALSE, nf = 4)
s.arrow(bca_site$co, plabels=list(cex=0.75,optim=TRUE))
rt_between_species=randtest(bca_site)
rt_between_species
# 53% of total inertia comes from the differences between sites

# BCA: Speices
bca_species = bca(x = pca1, fac = complete_data$Species, scannf = FALSE, nf = 2)
bca_species
s.arrow(bca_species$co)
rt_between_species=randtest(bca_species)
rt_between_species
# 34% of total inertia comes from the differences between species

# BCA: Year
bca_species = bca(x = pca1, fac = as.factor(complete_data$Year), scannf = FALSE, nf = 2)
bca_species
s.arrow(bca_species$co)
rt_between_species=randtest(bca_species)
rt_between_species
# 7.3% of total inertia comes from the differences between species

#########################################################################################
#### PCA and BCA on Untransformed, range normalized data ####
# Run the PCA
pca2 = dudi.pca(df = complete[,10:31], scannf = FALSE, nf = 2)
g1_2=s.arrow(pca2$co, plot=FALSE, plabels=list(cex=0.75,optim=TRUE),psub.text = "Untransformed Range Normalized")
g2_2=s.class(pca2$li,complete$Species, col=coul, plabels=list(optim=TRUE),psub.text = "Untransformed Range Normalized")
ADEgS(list(g1_2,g2_2))

#BCA: Site
bca_site = bca(x = pca2, fac = complete$Site, scannf = FALSE, nf = 4)
#s.arrow(bca_site$co, plabels=list(cex=0.75,optim=TRUE))
rt_between_species=randtest(bca_site)
rt_between_species
# 39% of total inertia comes from the differences between sites

#BCA: Speices
bca_species = bca(x = pca2, fac = complete$Species, scannf = FALSE, nf = 2)
bca_species
s.arrow(bca_species$co)
rt_between_species=randtest(bca_species)
rt_between_species
# 18.7% of total inertia comes from the differences between species

#BCA: Year
bca_Year = bca(x = pca2, fac = as.factor(complete$Year), scannf = FALSE, nf = 2)
bca_Year
s.arrow(bca_Year$co)
rt_between_Year=randtest(bca_Year)
rt_between_Year
# 7.3% of total inertia comes from the differences between Year

############################################################################################################
#### PCA and BCA on sqrt Transform, range normalized data ####
pca3 = dudi.pca(df = sqrt[,10:31], scannf = FALSE, nf = 2)
g1_3=s.arrow(pca3$co, plot=FALSE, plabels=list(cex=0.75,optim=TRUE),psub.text = "SQRT Tranformation Range Normalized")
g2_3=s.class(pca3$li,sqrt$Species, col=coul, plabels=list(optim=TRUE),psub.text = "SQRT Tranformation Range Normalized")
ADEgS(list(g1_3,g2_3))

#BCA: Site
bca_site = bca(x = pca3, fac = sqrt$Site, scannf = FALSE, nf = 4)
#s.arrow(bca_site$co, plabels=list(cex=0.75,optim=TRUE))
rt_between_species=randtest(bca_site)
rt_between_species
# 41.6% of total inertia comes from the differences between sites

#BCA: Speices
bca_species = bca(x = pca3, fac = sqrt$Species, scannf = FALSE, nf = 2)
bca_species
s.arrow(bca_species$co)
rt_between_species=randtest(bca_species)
rt_between_species
# 21% of total inertia comes from the differences between species

#BCA: Year
bca_year = bca(x = pca3, fac = as.factor(sqrt$Year), scannf = FALSE, nf = 2)
bca_year
s.arrow(bca_year$co)
rt_bca_year=randtest(bca_year)
rt_bca_year
# 21% of total inertia comes from the differences between species

######################################################################################
### PCA and BCA on log10 transformed range normalized data ###
pca4 = dudi.pca(df = log10[,10:31], scannf = FALSE, nf = 2)
g1_4=s.arrow(pca4$co, plot=FALSE, plabels=list(cex=0.75,optim=TRUE),psub.text = "Log10 Tranformation Range Normalized")
g2_4=s.class(pca4$li,log10$Species, col=coul, plabels=list(optim=TRUE),psub.text ="Log10 Tranformation Range Normalized")
ADEgS(list(g1_4,g2))

#BCA: Site
bca_site = bca(x = pca4, fac = log10$Site, scannf = FALSE, nf = 4)
#s.arrow(bca_site$co, plabels=list(cex=0.75,optim=TRUE))
rt_between_species=randtest(bca_site)
rt_between_species
# 41% of total inertia comes from the differences between sites

#BCA: Speices
bca_species = bca(x = pca4, fac = log10$Species, scannf = FALSE, nf = 2)
bca_species
s.arrow(bca_species$co)
rt_between_species=randtest(bca_species)
rt_between_species
# 20% of total inertia comes from the differences between species

#BCA: Year
bca_Year = bca(x = pca4, fac = as.factor(log10$Year), scannf = FALSE, nf = 2)
bca_Year
s.arrow(bca_Year$co)
rt_between_Year=randtest(bca_Year)
rt_between_Year
# 7.5% of total inertia comes from the differences between Year
###########################################################################################
#Combine Plots
ADEgS(list(g1_1,g2_1, g1_2,g2_2, g1_3,g2_3, g1_4,g2_4), layout=c(4, 2))

##########################################################################################
# Lifestage
#########################################################################################
# Set up colour pallet 
coul = c("#8dd3c7", "#ffd92f", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69")

#### PCA and BCA on Untransformed Unnormalized data ####
# Run the PCA
pca1 = dudi.pca(df = complete_data[,10:31], scannf = FALSE, nf = 2)
g1_1=s.arrow(pca1$co, plot=FALSE, plabels=list(cex=0.75,optim=TRUE),psub.text = "Untransformed Unnormalized")
g2_1=s.class(pca1$li,complete_data$Lifestage, col=coul, plabels=list(optim=TRUE),psub.text = "Untransformed Unnormalized")
ADEgS(list(g1_1,g2_1))

# BCA: Site
bca_site = bca(x = pca1, fac = complete_data$Site, scannf = FALSE, nf = 4)
s.arrow(bca_site$co, plabels=list(cex=0.75,optim=TRUE))
rt_between_species=randtest(bca_site)
rt_between_species
# 53% of total inertia comes from the differences between sites

# BCA: Speices
bca_species = bca(x = pca1, fac = complete_data$Lifestage, scannf = FALSE, nf = 2)
bca_species
s.arrow(bca_species$co)
rt_between_species=randtest(bca_species)
rt_between_species
# 34% of total inertia comes from the differences between species

# BCA: Year
bca_species = bca(x = pca1, fac = as.factor(complete_data$Year), scannf = FALSE, nf = 2)
bca_species
s.arrow(bca_species$co)
rt_between_species=randtest(bca_species)
rt_between_species
# 7.3% of total inertia comes from the differences between species

#########################################################################################
#### PCA and BCA on Untransformed, range normalized data ####
# Run the PCA
pca2 = dudi.pca(df = complete[,10:31], scannf = FALSE, nf = 2)
g1_2=s.arrow(pca2$co, plot=FALSE, plabels=list(cex=0.75,optim=TRUE),psub.text = "Untransformed Range Normalized")
g2_2=s.class(pca2$li,complete$Lifestage, col=coul, plabels=list(optim=TRUE),psub.text = "Untransformed Range Normalized")
ADEgS(list(g1_2,g2_2))

#BCA: Site
bca_site = bca(x = pca2, fac = complete$Site, scannf = FALSE, nf = 4)
#s.arrow(bca_site$co, plabels=list(cex=0.75,optim=TRUE))
rt_between_species=randtest(bca_site)
rt_between_species
# 39% of total inertia comes from the differences between sites

#BCA: Speices
bca_species = bca(x = pca2, fac = complete$Lifestage, scannf = FALSE, nf = 2)
bca_species
s.arrow(bca_species$co)
rt_between_species=randtest(bca_species)
rt_between_species
# 18.7% of total inertia comes from the differences between species

#BCA: Year
bca_Year = bca(x = pca2, fac = as.factor(complete$Year), scannf = FALSE, nf = 2)
bca_Year
s.arrow(bca_Year$co)
rt_between_Year=randtest(bca_Year)
rt_between_Year
# 7.3% of total inertia comes from the differences between Year

############################################################################################################
#### PCA and BCA on sqrt Transform, range normalized data ####
pca3 = dudi.pca(df = sqrt[,10:31], scannf = FALSE, nf = 2)
g1_3=s.arrow(pca3$co, plot=FALSE, plabels=list(cex=0.75,optim=TRUE),psub.text = "SQRT Tranformation Range Normalized")
g2_3=s.class(pca3$li,sqrt$Lifestage, col=coul, plabels=list(optim=TRUE),psub.text = "SQRT Tranformation Range Normalized")
ADEgS(list(g1_3,g2_3))

#BCA: Site
bca_site = bca(x = pca3, fac = sqrt$Site, scannf = FALSE, nf = 4)
#s.arrow(bca_site$co, plabels=list(cex=0.75,optim=TRUE))
rt_between_species=randtest(bca_site)
rt_between_species
# 41.6% of total inertia comes from the differences between sites

#BCA: Speices
bca_species = bca(x = pca3, fac = sqrt$Lifestage, scannf = FALSE, nf = 2)
bca_species
s.arrow(bca_species$co)
rt_between_species=randtest(bca_species)
rt_between_species
# 21% of total inertia comes from the differences between species

#BCA: Year
bca_year = bca(x = pca3, fac = as.factor(sqrt$Year), scannf = FALSE, nf = 2)
bca_year
s.arrow(bca_year$co)
rt_bca_year=randtest(bca_year)
rt_bca_year
# 21% of total inertia comes from the differences between species

######################################################################################
### PCA and BCA on log10 transformed range normalized data ###
pca4 = dudi.pca(df = log10[,10:31], scannf = FALSE, nf = 2)
g1_4=s.arrow(pca4$co, plot=FALSE, plabels=list(cex=0.75,optim=TRUE),psub.text = "Log10 Tranformation Range Normalized")
g2_4=s.class(pca4$li,log10$Lifestage, col=coul, plabels=list(optim=TRUE),psub.text ="Log10 Tranformation Range Normalized")
ADEgS(list(g1_4,g2))

#BCA: Site
bca_site = bca(x = pca4, fac = log10$Site, scannf = FALSE, nf = 4)
#s.arrow(bca_site$co, plabels=list(cex=0.75,optim=TRUE))
rt_between_species=randtest(bca_site)
rt_between_species
# 41% of total inertia comes from the differences between sites

#BCA: Speices
bca_species = bca(x = pca4, fac = log10$Lifestage, scannf = FALSE, nf = 2)
bca_species
s.arrow(bca_species$co)
rt_between_species=randtest(bca_species)
rt_between_species
# 20% of total inertia comes from the differences between species

#BCA: Year
bca_Year = bca(x = pca4, fac = as.factor(log10$Year), scannf = FALSE, nf = 2)
bca_Year
s.arrow(bca_Year$co)
rt_between_Year=randtest(bca_Year)
rt_between_Year
# 7.5% of total inertia comes from the differences between Year
###########################################################################################
#Combine Plots
ADEgS(list(g1_1,g2_1, g1_2,g2_2, g1_3,g2_3, g1_4,g2_4), layout=c(4, 2))





