########################################################
# Spatial PCA
# written in R version 3.4.3
########################################################
# Library
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
coul = c("#8dd3c7", "#ffd92f", "#bebada", "#fb8072", "#80b1d3", "#fdb461", "#b3de69")
#########################################################
## RAW DATA ##
#Import Data
complete_data = read.csv("grid_complete_data.csv")

#### Transform Data ####
#Sqrt Tranform
sqrt_df = sqrt(complete_data[,10:31])
sqrt_data = cbind(complete_data[,1:9], sqrt_df)

# SQRT 
fac = interaction(factor(sqrt_data$Species),factor(sqrt_data$Lifestage), drop = T)
fnorm.var = function(y){ return((y - min(y)) / (max(y) - min(y)))}
fnorm.tab = function(tab) {apply(tab, 2, fnorm.var)}
range_norm_df=unsplit(lapply(split(as.data.frame(sqrt_data[,10:31]), fac), 
                             FUN = function(x) as.data.frame(fnorm.tab(x))), fac)
summary(range_norm_df)
# Replace NA's with zeros
range_norm_df=rapply(range_norm_df, f=function(x) ifelse(is.nan(x),0,x), how="replace")
sqrt=as.data.frame(cbind(sqrt_data[,1:9],range_norm_df))

#### PCA and BCA on sqrt Transform, range normalized data ####
pca3 = dudi.pca(df = sqrt[,10:31], scannf = FALSE, nf = 4)
g1_3=s.arrow(pca3$co, plot=FALSE, plabels=list(cex=0.75,optim=TRUE),psub.text = "SQRT Tranformation Range Normalized")
g2_3=s.class(pca3$li,sqrt$Species, col=coul, plabels=list(optim=TRUE),psub.text = "SQRT Tranformation Range Normalized")
ADEgS(list(g1_3,g2_3))

#BCA: Site
bca_site = bca(x = pca3, fac = sqrt$Site, scannf = FALSE, nf = 4)
#s.arrow(bca_site$co, plabels=list(cex=0.75,optim=TRUE))
rt_between_species=randtest(bca_site)
rt_between_species
# 41.9% of total inertia comes from the differences between sites

s.arrow(bca_site$co)
write.csv(bca_site$C1, "sqrt_bca_loadings.csv")

###################################################################
#### Map the BCA components ####
# Extract data for mapping
# Aggregate the x y coordinates for sites
sitexy=aggregate(sqrt[,3:4], list(sqrt$Site), mean)

# Extract information by site from BCA
dim(bca_site$tab)
scores_by_site=bca_site$li

# Attach Coordinates to sites
scores_xy=cbind(sitexy, scores_by_site)
write.csv(scores_xy, "site_scores_sqrt_norm.csv")

# Set up Spatial DataFrame
rivers = readOGR(dsn = "C:/Users/Kristin/Dropbox/chapter4_pca/All_Species/PCA_All", layer = "rivers_clip")
rivers_convert <- fortify(rivers)
g1.map.bca=s.value(scores_xy[,3:2],bca_site$li, symbol="circle", 
                   pgrid.draw = FALSE,ppoints.cex=0.75,
                   Sp=rivers, pSp.col="#CBF6FF", pSp.border="#CBF6FF", ylim = c(52, 61), xlim = c(-122, -106))

# Map Species Locations
ggplot(complete_data, aes(x=Long, y=Lat)) +
  xlab("Longitude")+
  ylab("Latitude")+
  geom_polygon(data=rivers_convert, aes(long, lat, group = group, fill = hole), colour = alpha("lightgrey", 1/2), size = 0.7) + 
  scale_fill_manual(values = c("lightgrey", "white"), guide=FALSE)+
  geom_point(aes(color = factor(complete_data$Species)),size=4)+
  scale_color_manual(values=c("#8dd3c7", "#ffd92f", "#bebada", "#fb8072", "#80b1d3", "#fdb461", "#b3de69"), name="Species")+
  theme_minimal()

####################################################################
#### Spatia PCA ####
# Defining spatial weights
# Create spatial coordinates
xy=coordinates(scores_xy[,3:2])
# To explore other spatial weights matrices
#listw.explore()

#No weights on the edges
#nb <- chooseCN(coordinates(xy), type = 2, plot.nb = FALSE) 
#lw <- nb2listw(nb, style = 'W', zero.policy = TRUE)

# Weighted edges
nb <- chooseCN(coordinates(xy), type = 2, plot.nb = FALSE) 
distnb <- nbdists(nb, xy) 
fdist <- lapply(distnb, function(x) 1 - x/max(dist(xy))) 
lw <- nb2listw(nb, style = 'W', glist = fdist, zero.policy = TRUE)

# Plot neighbourhoods
g1=s.label(scores_xy[,3:2],nb=nb, pub.edge.col="red", ppoints.cex=0.25)

# Moran's I IDW
moran.randtest(scores_xy[,"Axis1"], listw=lw,nrepet=999)
moran.plot(scores_xy[,"Axis1"], listw=lw) # clustered

moran.randtest(scores_xy[,"Axis2"], listw=lw,nrepet=999)
moran.plot(scores_xy[,"Axis2"], listw=lw) #Clustered

moran.randtest(scores_xy[,"Axis3"], listw=lw,nrepet=999)
moran.plot(scores_xy[,"Axis3"], listw=lw) #clustered 

moran.randtest(scores_xy[,"Axis4"], listw=lw,nrepet=999)
moran.plot(scores_xy[,"Axis4"], listw=lw) # Clustered 

#Moran's eigenvector maps using gabriel neighbourhood
(me=mem(lw))
s.value(scores_xy[,3:2], me[,c(1:4)], ppoints.cex=0.75, ylim = c(52, 61), xlim = c(-120, -108))
scalo1=scalogram(scores_xy[,"Axis1"], me, nblocks=10)
plot(scalo1)
scalo2=scalogram(scores_xy[,"Axis2"], me, nblocks=10)
plot(scalo2)
scalo3=scalogram(scores_xy[,"Axis3"], me, nblocks=10)
plot(scalo3)
scalo4=scalogram(scores_xy[,"Axis4"], me, nblocks=10)
plot(scalo4)


### multispqti
ms1 = multispati(bca_site, lw)
ms1_weight_c1=as.data.frame(ms1$c1)
plot(ms1)
s.value(sitexy[,3:2], ms1$li, porigin.include = FALSE)
s.value(sitexy[,3:2], ms1$li[,1], porigin.include = FALSE, ppoints.cex=0.5, ylim = c(52, 61), xlim = c(-120, -108),
        Sp=rivers, pSp.col="#CBF6FF", pSp.border="#CBF6FF")
s.value(sitexy[,3:2], ms1$li[,2], porigin.include = FALSE, ppoints.cex=0.5, ylim = c(52, 61), xlim = c(-120, -108),
        Sp=rivers, pSp.col="#CBF6FF", pSp.border="#CBF6FF")
s.value(sitexy[,3:2], ms1$li[,3], porigin.include = FALSE, ppoints.cex=0.5, ylim = c(52, 61), xlim = c(-120, -108),
        Sp=rivers, pSp.col="#CBF6FF", pSp.border="#CBF6FF")
s.value(sitexy[,3:2], ms1$li[,4], porigin.include = FALSE, ppoints.cex=0.5, ylim = c(52, 61), xlim = c(-120, -108),
        Sp=rivers, pSp.col="#CBF6FF", pSp.border="#CBF6FF")
s.arrow(ms1$c1)

write.csv(ms1$c1, "multispati_loadings_site_all.csv")

