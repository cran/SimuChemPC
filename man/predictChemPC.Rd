\name{predictChemPC}
\alias{predictChemPC}
\alias{predictChemPC, character list, character list, vector, character list, vector}
\title{predictChemPC}
\description{This function performs a prediction method from 4 predicting potent compounds methods of this package. These methods are RA, EI, NN and GP.}
\usage{
predictChemPC(xTrain, yTrain, xTest, loghyper, method="RA")
}

\arguments{
  \item{xTrain}{ m * n martrix of train data.}
  \item{yTrain}{ m * 1 matrix of target values consist of potencies, pIC50 or other measurements of compound affinities that are desired to be maximized.}
  \item{xTest}{ j * n matrix of test data.}
  \item{loghyper}{ 3 * 1 matrix of loghyper parameters which is the output of trainChemPC function. } 
  \item{method}{ One of "EI", "GP", "NN" or "RA". }
}

\details{
This function withholds 4 methods to predict potent compounds.


\code{method} is one of:
	\code{EI} A compound for which maximum expected potency improvement is reached.
	\code{GP} A compound holding maximum predicted potency in test data is selected.
	\code{NN} A compound that is nearest (Tonimito Coefficient as distance measure) to the most potent compound in training data is selected.
	\code{RA} As it's name suggests, a compound is selected randomly.

\code{Feature selection}
Feature selection employed in this package is based on Spearman Rank Correlation such that
before each training step those attributes in which revealed a significant Spearman rank correlation
with the logarithmic potency values (q-value < 5%) of the training data are selected. Q-values
are computed from original p-values via the multiple testing correction method by Benjamini and Hochberg.

}

\value{
It returns index of most potent compound in original test set w.r.t. selected method.
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
xstar = as.data.frame(array(5:105, dim=c(20,5)))
loghyper = trainChemPC(x, y)
index = predictChemPC(x, y, xstar, loghyper, method="RA")
}
\keyword{predict, prediction, chemical, potent compounds, constraint global optimization, expected potency improvement, gaussian process}
