# Downloaded openjdk (8 version?) msi installer, install; Custom setup: Set or override JAVA_HOME variable

# if packages are not installed, we have to install them
#install.packages("readxl")
#install package for chemoinformatics
#install.packages("rcdk")
#install.packages("rcdklibs")


# import library for excel reading
library(readxl)
library(rcdk)
library(rcdklibs)
library(fingerprint)
library(proxy)

# check the version of the library
cdk.version()

# Intro
# read molecule in a form of SMILES
smile <- 'c1ccccc1CC(=O)C(N)CC1CCCCOC1'
mol <- parse.smiles(smile)[[1]]

# custom function to visualize the molecule
toshow <- function(mol,pathsd){
  result = tryCatch({
    factory <- .jnew("org.openscience.cdk.depict.DepictionGenerator")$withAtomColors()
    factory$withSize(1000,1000)#$getStyle("cow")
    temp1 <- paste0(pathsd)
    result<-factory$depict(mol)$writeTo(temp1)
  }, warning = function(w) {
    result=NULL
  }, error = function(e) {
    result=NULL
  })
  return(result)
}
# visualize the mol and save it in png file
toshow(mol,paste0(getwd(),"/mol_check.png"))

# toy example: predict molecules' boiling points from mol descriptors
# take a look at the available descriptor categories.
dc <- get.desc.categories()
dc
# get the names of the descriptors for a topological category
dn <- get.desc.names(dc[4])
dn

#The CDK package also provides an example data set, called bpdata which contains 277 molecules, 
#in SMILES format and their associated boiling points (BP) in Kelvin. The data.frame has two columns, 
#viz., the SMILES and the BP. Molecules names are used as row names
data(bpdata)
head(bpdata)
dim(bpdata)

mols <- parse.smiles(bpdata[,1])
descNames <- c(
  'org.openscience.cdk.qsar.descriptors.molecular.KierHallSmartsDescriptor',
  'org.openscience.cdk.qsar.descriptors.molecular.APolDescriptor',
  'org.openscience.cdk.qsar.descriptors.molecular.HBondDonorCountDescriptor')
descs <- eval.desc(mols, descNames)
class(descs)
head(descs)

#build a linear regression model. 
# First, remove NA’s, correlated and constant columns. The code is shown below, but 
# since it involves a stochastic element, we will not run it for this example. 
# If we were to perform feature selection, then this type of reduction would have to be performed.
descs <- descs[, !apply(descs, 2, function(x) any(is.na(x)) )]
descs <- descs[, !apply( descs, 2, function(x) length(unique(x)) == 1 )]
r2 <- which(cor(descs)^2 > .6, arr.ind=TRUE)
r2 <- r2[ r2[,1] > r2[,2] , ]
descs <- descs[, -unique(r2[,2])]

# modeling
model <- lm(BP ~ khs.sCH3 + khs.sF + apol + nHBDon, data.frame(bpdata, descs))
summary(model)
# plot the result
plot(bpdata$BP, predict(model, descs),
     xlab="Observed BP", ylab="Predicted BP",
     pch=19, xlim=c(100, 700), ylim=c(100, 700))
abline(0,1, col='red')

# work with real dataset of molecules with and without antibacterial properties
# read xlsx file, skip first row, create dataframe
df <- read_excel("C:/Users/akuznetsova/PycharmProjects/TA_bioinf/data/mmc1.xlsx", sheet='S1B', skip=1)
head(df)
# check no missing values
sum(is.na(df))
# see unique target values
unique(df$Activity)
# convert "Active" to 1 and "Inactive" to 0
binary <- c(1,0)
names(binary) = unique(df$Activity)
df$Act_bin = binary[df$Activity]
head(df)

# get ECFP fingerprints for the molecule
fp_mol <- get.fingerprint(mol, type ='extended', depth=6, size=1024)

# repeat the same for the list of molecules
smiles <- df$SMILES
mols <- sapply(smiles, parse.smiles)

# get fingerprints for all the molecules
fps <- lapply(mols, get.fingerprint, type ='extended', depth=6, size=1024)
# write fingerprints from fingerprint object to a matrix; 
# Each element is 1 or 0 (1's being specified by the positions in each fingerprint vector)
fps_matrix <- fp.to.matrix(fps)
print(dim(fps_matrix))

# calculate Tanimoto similarity between first two mols
simil(list(fps_matrix[1,], fps_matrix[2,]), method="tanimoto") #, a = "numeric", b = "numeric")
# calculate a pairwise similarity matrix using the Tanimoto metric
fp.sim <- fingerprint::fp.sim.matrix(fps, method='tanimoto')
# Since R’s hclust method requires a distance matrix, 
# we convert the similarity matrix to a distance matrix
fp.dist <- 1 - fp.sim
# cluster mols
cls <- hclust(as.dist(fp.dist))
plot(cls, main='A Clustering of the a/b dataset', labels=FALSE)