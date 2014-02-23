\name{trainChemPC}
\alias{trainChemPC}
\alias{trainChemPC, character list, character list}
\title{trainChemPC}
\description{Apply Gaussian Process Regression to learn the model.}
\usage{
trainChemPC(xTrain, yTrain)
}

\arguments{
  \item{xTrain}{ m * n martrix of train data.}
  \item{yTrain}{ m * 1 matrix of target values consist of potencies, pIC50 or other measurements of compound affinities that are desired to be maximized.}
}

\details{
This function performs training step of GP or EI by finding loghyper parameters.
}

\value{
It returns a vector that holds calculated loghyper parameters.
}

\references{
	\code{1.}Predicting Potent Compounds via Model-Based Global Optimization, Journal of Chemical Information and Modeling, 2013, 53 (3), pp 553-559, M Ahmadi, M Vogt, P Iyer, J Bajorath, H Froehlich.
	\code{2.}Software MOE is used to calculate the numerical descriptors in data sets. Ref:  http://www.chemcomp.com/MOE-Molecular_Operating_Environment.htm
	\code{3.}ChEMBL was the source of the compound data and potency annotations in data sets. Ref: https://www.ebi.ac.uk/chembl/
}

\author{Mohsen Ahmadi}
\examples{
x = as.data.frame(array(1:100, dim=c(20,5)))
y = as.matrix(as.numeric(array(1:20, dim=c(20,1))))
loghyper = trainChemPC(x, y)
}
\keyword{train, chemical, potent compounds, constraint global optimization, expected potency improvement, gaussian process}
