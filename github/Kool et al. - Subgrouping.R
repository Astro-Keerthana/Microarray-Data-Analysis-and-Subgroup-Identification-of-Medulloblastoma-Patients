## Install CRAN libraries
install.packages("RColorBrewer")
install.packages("rgl")
install.packages("NMF")
BiocManager::install()
BiocManager::install("limma")


## Load libraries
library(RColorBrewer)
library(NMF)
library(rgl)


## Set a working directory – I suggest a folder below, but you’ll have to 
## create it first
setwd("E:\\dESKTOP\\Keerthana\\Bsc (Hons) in Biotechnology\\3Y1S\\Bioinformatics\\Assignment\\R\\Microarray - Subgrouping Tumor cells")

#Part 2: Microarray exploratory analysis and subgroup identification
#When investigating high-dimensional dataset for the first time, it’s important to 
## explore the data to see what sort of patterns are present. In this section, you 
## will import the data, explore structure in the data and identify subgroups using 
# Non-negative matrix factorisation (NMF). The quality of the assigned subgroups will be 
# tested using silhouette plots and principal component analysis (PCA).
#Since the normalisation process is time-consuming and needs a lot of processing power, 
# I have pre-computed this step for you. I have downloaded the data, 
# normalised using gcrma (see lecture notes) and selected the 1500 most variable probes 
# for us to discover groups of related tumour samples.

## Load the data
load("E:\\dESKTOP\\Keerthana\\Bsc (Hons) in Biotechnology\\3Y1S\\Bioinformatics\\Assignment\\R\\filt.exp (1).RData")
## Check the dimensions of this filtered data using dim()
dim(filt.exp) # Gives the answer in rows first then columns


### Exploratory Investigation

## Explore patterns of variation using principal component analysis
## Do the analysis and store it in an object called 'filt.exp.pca'
filt.exp.pca <- prcomp(t(filt.exp), scale=T)

## Let's look at the amount of variation in the first few PCs
plot(filt.exp.pca)

summary(filt.exp.pca)

## Look at row 2 and 3 - proportion of variance explained and cumulative proportion
## So what this is saying is that if we look at the first three components (PC1, PC2 and PC3),
## we are capturing almost 50% of the total variation in the dataset

## Let's plot the first two principal components
plot(filt.exp.pca$x, pch=19, main="62 medulloblastoma samples - PCA plot")


#################
## There is no wrong answer here - how many groups do you think 
## there are by eye - give a range if you think it is unclear?
##################

## Heatmaps can be useful to display patterns of expression from multiple
## probes across multiple samples. Here we show genes with low expression ## blue, high expression red. 
# Each column is a sample and each row is a probe

## Heatmap of the data
heatmap(as.matrix(filt.exp), 
        col=colorRampPalette(c("darkblue", "white", "darkred"))(255), 
        labRow=NA, labCol=NA, scale="none")

#######################
## The dendrogram can give clues as to the number of clusters - how many do 
## you see (again, no wrong answer)
########################

## The above plots have shown us how to guess the number of clusters by
## eye, but we need a more formal, unbiased way to assess optimal cluster
## number.
## The optimal choice of cluster number is a balance, whereby
## cluster numbers would ideally be maximised while maintaining cluster
## quality.

### Cluster using NMF

## Identifying cluster numbers by examining dendrograms is 
## imprecise and not ideal, since this sort of high dimensional 
## data is noisy - using NMF can help us to pick out the major
## biological effects in the data (metagenes) which can be used to 
## identify optimal cluster numbers

###########################################
## How does NMF work and what is a metagene?
###########################################

### Cluster using NMF – try a range of metagene numbers
### This takes up to 10 minutes!
### So take a break if you wish

estim.r <- nmf(filt.exp, rank = 2:6, nrun=50, seed=1234)
## Look to maximise cophenetic coefficient, dispersion and silhouette
##(consensus), whilst also maximising cluster number
plot(estim.r)

## Let's look at the consensus map
## This shows the reproducibility with which samples are assigned to the same cluster
consensusmap(estim.r)

#####################
## How many clusters do you think is appropriate, bearing in mind 
## we want to maximise numbers of groups, whilst also maximising quality 
## of clusters – also bear in mind that we know there are at least 3!
#####################


## I’m selecting 4 in subsequent analyses, but am comparing with 5
## groups, since that's what Kool and colleagues selected.
## If you think differently, 
## that’s fine, speak to me to discuss 
## and we can move forward with your group number

## Lets tell the NMF algorigthm to identify four groups
filt.nmf <- nmf(filt.exp, rank = 4, nrun=20, seed=1)

## Coefficient maps show the metagene values for each sample
coefmap(filt.nmf)

## Assign groups by predicting group membership on the metagene patterns
groups <- predict(filt.nmf)

names(groups) <- colnames(filt.exp)
## Save to disk in a format that can be read by Excel
write.table(groups, file="groupAssignments.csv", sep=",")

## Silhouette plots give a silhouette score for each cluster and each sample
## A perfectly clustered sample would have a silhouette score of 1
## A sample that lay exactly between two groups would have a silhouette
## score of 0
## A sample most likely assigned to the wrong group would have a 
## negative silhouette score

## Plot silhouette plot of assigned groups
x <- silhouette(dist(t(filt.exp))^2, x=as.numeric(groups))
plot(x)

#####################
## What is a silhouette plot showing?
## How many samples are in each cluster group?
## Are the clusters good quality?
## Are there any samples that are poorly clustered?
#####################

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

#########################################
## Which do you think gives the better quality clusters, 4 or 5 groups?
## Do you agree or disagree with Kool and colleagues that there are 5
## groups?
#########################################


## Replot PCA, colouring the samples by the assigned group

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


## When you’ve finished, save your work by using the q() command. Click yes when prompted if you want to save, making sure that you save your script and your workspace.

#########################################
## The first part of your report should focus on the data exploration that 
## we’ve done today. Remember I am not interested in the sequence of 
## commands you’ve used, rather that you’ve understood the results we’ve 
## generated
#########################################

