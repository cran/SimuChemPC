\name{trainChemPC}
\alias{trainChemPC}
\alias{trainChemPC, character list, character list}
\title{trainChemPC}
\description{This function performs learning step of GP and EI.}
\usage{
trainChemPC( trainData, targetVector)
}

\arguments{
  \item{trainData}{\code{trainData} specifies a data frame including an array of a training data with dimension of n*m. This data is used in to find a hyperlog vector.}
  \item{targetVector}{\code{targetVector} is a one dimensional array with n rows which is equal to number of rows in trainData.}
}

\details{
This function performs training step of GP or EI by finding a loghyper using gpr package. A loghyper can be used to predict potent compounds.
}

\value{
It returns a vector that holds a calculated loghyper.
}

\references{
	\code{1.}Predicting Potent Compounds via Model-Based Global Optimization, Journal of Chemical Information and Modeling, 2013, 53 (3), pp 553-559, M Ahmadi, M Vogt, P Iyer, J Bajorath, H Froehlich.
	\code{2.}Software MOE is used to calculate the numerical descriptors in data sets. Ref:  http://www.chemcomp.com/MOE-Molecular_Operating_Environment.htm
	\code{3.}ChEMBL was the source of the compound data and potency annotations in data sets. Ref: https://www.ebi.ac.uk/chembl/
}

\author{Mohsen Ahmadi}
\examples{
library(gpr)
library(SimuChemPC)
a = as.data.frame(array(1:10, dim=c(2,5)))
b = array(1:2)
myloghyper = trainChemPC( a, b)
}
\keyword{train, chemical, potent compounds, constraint global optimization, expected potency improvement, gpr, gaussian process}
