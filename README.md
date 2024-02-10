# Microarray-Data-Analysis-and-Subgroup-Identification-of-Medulloblastoma-Patients
## Install CRAN libraries
install.packages("RColorBrewer")
install.packages("rgl")
## Install CRAN libraries
install.packages("RColorBrewer")
install.packages("rgl")
install.packages("NMF")
BiocManager::install()
## Load libraries
library(RColorBrewer)
library(NMF)
library(rgl)
## Set a working directory – I suggest a folder below, but you’ll have to
## create it first
setwd("Microarray - Subgrouping Tumor cells")
## Load the data
load("filt.exp (1).RData")
## Check the dimensions of this filtered data using dim()
dim(filt.exp) # Gives the answer in rows first then columns
## Explore patterns of variation using principal component analysis
## Do the analysis and store it in an object called 'filt.exp.pca'
filt.exp.pca <- prcomp(t(filt.exp), scale=T)
## Let's look at the amount of variation in the first few PCs
plot(filt.exp.pca)
summary(filt.exp.pca)
## Let's plot the first two principal components
plot(filt.exp.pca$x, pch=19, main="62 medulloblastoma samples - PCA plot")
## Heatmap of the data
heatmap(as.matrix(filt.exp),
col=colorRampPalette(c("darkblue", "white", "darkred"))(255),
labRow=NA, labCol=NA, scale="none")
estim.r <- nmf(filt.exp, rank = 2:6, nrun=50, seed=1234)
## Look to maximise cophenetic coefficient, dispersion and silhouette
##(consensus), whilst also maximising cluster number
plot(estim.r)
## Let's look at the consensus map
## This shows the reproducibility with which samples are assigned to the same cluster
consensusmap(estim.r)
## Lets tell the NMF algorigthm to identify four groups
filt.nmf <- nmf(filt.exp, rank = 4, nrun=20, seed=1)
## Coefficient maps show the metagene values for each sample
coefmap(filt.nmf)
## Assign groups by predicting group membership on the metagene patterns
groups <- predict(filt.nmf)
names(groups) <- colnames(filt.exp)
## Save to disk in a format that can be read by Excel
write.table(groups, file="groupAssignments.csv", sep=",")
## Plot silhouette plot of assigned groups
x <- silhouette(dist(t(filt.exp))^2, x=as.numeric(groups))
plot(x)
## Fill in the numbers of each group and the colours as appropriate
plot(x, col=c(rep("darkgreen",26),rep("blue",9),
rep("yellow2",12),rep("red",15)))
## Compare with a 5 group solution
filt.nmf.5 <- nmf(filt.exp, rank = 5, nrun=20, seed=1)
groups.5 <- predict(filt.nmf.5)
## Plot silhouette plot of assigned groups
x.5 <- silhouette(dist(t(filt.exp))^2, x=as.numeric(groups.5))
## Assignment of group number is random, so rearrange numbers to match the plot
plot(x.5)
plot(x.5, col=c(rep("darkgreen",9),rep("blue",14),
rep("yellow2",15),rep("red",11), rep("purple",13)))
cols <- ifelse(groups==1, "darkgreen",
ifelse(groups == 2, "blue",
ifelse(groups == 3, "yellow2","red")))
plot(filt.exp.pca$x, pch=19, col=cols)
## And plot in 3D – play around with the orientation of the plot by
## clicking on the plot and moving the mouse to rotate
plot3d(filt.exp.pca$x[,1:3], col=cols, type="s", size=1.2)
#########################################
## Do the 4 groups that we’ve identified fit with what you see in the
## PCA plot?
#########################################
q()
