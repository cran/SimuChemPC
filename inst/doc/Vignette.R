#######################################################
# Data Initialization
#######################################################
library(rcdk)
library(SimuChemPC)
mols <- load.molecules( c('sample.sdf') )
dn <- get.desc.names()
dn[47]
features <- eval.desc(mols, dn[47])
get.properties(mols[[1]])
get.property(mols[[1]], "Potency[nM]")
get.property(mols[[1]], "cr_index")
p = data.frame()

for (i in 1:length(mols))
{
	Potency = get.property(mols[[i]], "Potency[nM]");
	p = rbind(p, data.frame(Potency))
}

datafile <- cbind(features, p)
dim(datafile)
datafile[1:5,]
#######################################################
# Utility Functions
#######################################################
len <- dim(datafile)[1]
len
half <- dim(datafile)[1] / 2
half
col <- dim(datafile)[2]
col

xTrain <- datafile[1:half, 1:(col-1)]
yTrain <- as.matrix(as.numeric(array(datafile[1:half, col])))
xTest <- datafile[(half+1):len, (1:col-1)]
yTest <- as.matrix(as.numeric(array(datafile[(half+1):len, col])))

loghyper = trainChemPC(xTrain, yTrain)
loghyper

predictChemPC(xTrain, yTrain, xTest, loghyper, method="RA")
predictChemPC(xTrain, yTrain, xTest, loghyper, method="NN")
predictChemPC(xTrain, yTrain, xTest, loghyper, method="GP")
predictChemPC(xTrain, yTrain, xTest, loghyper, method="EI")
#######################################################
# Simulation set-up
#######################################################
dataX <- datafile[, 1:(col-1)]
dataY <- as.matrix(as.numeric(array(datafile[, col])))
dim(dataX)
dim(dataY)
rank <- SimuChemPC(dataX, dataY, "RA", 5)
rank <- SimuChemPC(dataX, dataY, "NN", 5)
rank <- SimuChemPC(dataX, dataY, "GP", 5)
rank <- SimuChemPC(dataX, dataY, "EI", 5)
#######################################################
#R Session Information
#######################################################
sessionInfo()