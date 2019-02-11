###############################################################

###############################################################
##### Import Data ####
#Take the aboslute values of the imputed data
data=read.csv("grid_complete_data.csv")

#### Load Libraries ####
library(ade4)
library(adespatial)
library(adegenet)
library(sp)
library(spdep)
library(maptools)
library(adegraphics)
library(RColorBrewer)
library(readOGR)
library(rgdal)

###############################################################
#### Set up data ####
#Normalize between 0 and 1 data by speices and lifestage
fac = interaction(factor(data$Species),factor(data$Lifestage), drop = T)
fnorm.var = function(y){ return((y - min(y)) / (max(y) - min(y)))}
fnorm.tab = function(tab) {apply(tab, 2, fnorm.var)}
range_norm_df=unsplit(lapply(split(as.data.frame(data[,9:31]), fac), 
                             FUN = function(x) as.data.frame(fnorm.tab(x))), fac)
summary(range_norm_df)
# Replace NA's with zeros
range_norm_df=rapply( range_norm_df, f=function(x) ifelse(is.nan(x),0,x), how="replace")
norm_data=as.data.frame(cbind(data[,2:8],range_norm_df))

# Subets data
amphib = subset(norm_data, Species == "Amphibian")
gull =  subset(norm_data, Species == "Gull"| Species == "Tern")
mammal =  subset(norm_data, Species == "Fisher" | Species == "Marten"| Species == "Mink" | Species == "Muskrat"| Species == "River Otter")

#################################################################
#### Amphibians ####
# PCA
pca1= dudi.pca(df = amphib[,9:30], center = TRUE, scannf = FALSE, nf = 4)
g1=s.corcircle(pca1$co, plot=FALSE)
g2=s.label(pca1$li, plot=FALSE)
ADEgS(list(g1,g2))

# Between group analysis: SITE
amphib_site=bca(pca1,as.factor(amphib$Lat),scannf=FALSE)
#plot(amphib_site, row.pellipse.col=adegpar()$ppalette$quali(6))
rt_between_site=randtest(amphib_site)
rt_between_site

# Spatial Analysis 
# Aggregate the x y coordinates for sites
xy=coordinates(amphib[,3:2])
sitexy=aggregate(amphib[, 3:2], list(amphib$Site), mean)

# Extract information by site from BCA
dim(amphib_site$tab)
scores_by_site=amphib_site$li

# Attach Coordinates to sites
scores_xy=cbind(sitexy, scores_by_site)
write.csv(scores_xy, "amphib_site_scores_xy.csv")

# Make the Maps
rivers = readOGR(dsn = "C:/Users/keccl081/Dropbox/Chapter4_PCA/All_Species/PCA_All","rivers_clip")
r2 = plot(rivers, col="#CBF6FF", border="#CBF6FF")

amphib_map_bca=s.value(scores_xy[,2:3],amphib_site$li, symbol="circle", pgrid.draw = FALSE, 
                       Sp=rivers, pSp.col="#CBF6FF", pSp.border="#CBF6FF",
                       ylim = c(55, 62), xlim = c(-115, -109))
amphib_bi = s.arrow(amphib_site$co)

#################################################################
#### Gulls ####
# PCA
pca1= dudi.pca(df = gull[,8:30], center = TRUE, scannf = FALSE, nf = 4)
g1=s.corcircle(pca1$co, plot=FALSE)
g2=s.label(pca1$li, plot=FALSE)
ADEgS(list(g1,g2))

# Between group analysis: SITE
gull_site=bca(pca1,as.factor(gull$Lat),scannf=FALSE)
#plot(gull_site, row.pellipse.col=adegpar()$ppalette$quali(6))
rt_between_site=randtest(gull_site)
rt_between_site

# Spatial Analysis 
# Aggregate the x y coordinates for sites
xy=coordinates(gull[,3:2])
sitexy=aggregate(gull[, 3:2], list(gull$Site), mean)

# Extract information by site from BCA
dim(gull_site$tab)
scores_by_site=gull_site$li

# Attach Coordinates to sites
scores_xy=cbind(sitexy, scores_by_site)
write.csv(scores_xy, "gull_site_scores_xy.csv")

# Make the Maps
rivers = readOGR(dsn = "C:/Users/keccl081/Dropbox/Chapter4_PCA/All_Species/PCA_All","rivers_clip")
#plot(rivers, col="#CBF6FF", border="#CBF6FF")

gull_map_bca=s.value(scores_xy[,2:3],gull_site$li, symbol="circle", pgrid.draw = FALSE, 
                       Sp=rivers, pSp.col="#CBF6FF", pSp.border="#CBF6FF",
                       ylim = c(46, 62), xlim = c(-120, -101))
gulls_bi = s.arrow(gull_site$co)

#################################################################
#### Mammals ####
# PCA
pca1= dudi.pca(df = mammal[,9:30], center = TRUE, scannf = FALSE, nf = 4)
g1=s.corcircle(pca1$co, plot=FALSE)
g2=s.label(pca1$li, plot=FALSE)
ADEgS(list(g1,g2))

# Between group analysis: SITE
mammal_site=bca(pca1,as.factor(mammal$Lat),scannf=FALSE)
#plot(mammal_site, row.pellipse.col=adegpar()$ppalette$quali(6))
rt_between_site=randtest(mammal_site)
rt_between_site

# Spatial Analysis 
# Aggregate the x y coordinates for sites
xy=coordinates(mammal[,3:2])
sitexy=aggregate(mammal[, 3:2], list(mammal$Site), mean)

# Extract information by site from BCA
dim(mammal_site$tab)
scores_by_site=mammal_site$li

# Attach Coordinates to sites
scores_xy=cbind(sitexy, scores_by_site)
write.csv(scores_xy, "mammal_site_scores_xy.csv")

# Make the Maps
rivers = readOGR(dsn = "C:/Users/keccl081/Dropbox/Chapter4_PCA/All_Species/PCA_All","rivers_clip")
#plot(rivers, col="#CBF6FF", border="#CBF6FF")

mammal_map_bca=s.value(scores_xy[,2:3],mammal_site$li, symbol="circle", pgrid.draw = FALSE, 
                     Sp=rivers, pSp.col="#CBF6FF", pSp.border="#CBF6FF",
                     ylim = c(53, 59), xlim = c(-120, -109))
mammals_bi = s.arrow(mammal_site$co)

###################################################
# Plot all figures together

ADEgS(list(amphib_map_bca, amphib_bi, 
           gull_map_bca, gulls_bi, 
           mammal_map_bca, mammals_bi), layout = c(3, 2))
