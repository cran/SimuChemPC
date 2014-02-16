\name{SimuChemPC}
\alias{SimuChemPC}
\alias{SimuChemPC, character list, character list, character list, character list,integer}
\title{SimuChemPC}
\description{ This function executes a simulation to compare 4 methods for predicting potent compounds. These methods are EI, GP, NN and RA.}
\usage{
SimuChemPC(dataX, dataY, method="RA", experiment=1)
}

\arguments{
  \item{dataX}{ m * n martrix of data (features/descriptors).}
  \item{dataY}{ m * 1 matrix of target values (potencies).} 
  \item{method}{ One of "EI", "GP", "NN" or "RA". }
  \item{experiment}{ An integer value that indicates a number by which experiment repeats. In our published experiment it was set to 25. }
}

\details{
This function withholds 4 simulation methods to predict potent compounds .
\code{method} can be RA, NN , EI or GP. The explanation of the abbreviations is listed below.

	\code{RA selection:} One compound will be selected randomly and added to train data each time.
	
	\code{NN selection:} The compound which is nearest (based on Tonimito Coefficient) to the most potent compound in training data is selected and added to train data.
	
	\code{EI selection} A compound for which maximum expected potency improvement is reached, is selected and then it is added to train data.
	
	\code{GP selection} A compound holding maximum potency in test data is selected.

\code{Feature selection}
Feature selection employed in this package is based on Spearman Rank Correlation such that
before each training step those attributes in which revealed a significant Spearman rank correlation
with the logarithmic potency values (q-value < 5%) of the training data are selected. Q-values
are computed from original p-values via the multiple testing correction method by Benjamini and Hochberg.

\code{The purpose of simulation step}
Simulation step is employed to select the compound(in the case where input files are chemical compounds)
in which maximal expected potency improvement is met. Subsequently, this compound is added to train data 
and simulation continues until all test data are consumed. Finally, the number of simulation steps is determined which
the algorithm used to select the most potent compound in the "original" test set.

\code{In this code, given our data sets (chemical compounds), we do the followings:}

	1. We split our data into two distinguish parts namely Train and Test data

	2. We do normalization on both parts

	3. We employ a specific feature selection algorithm (i.e. Multiple Testing Correction) to overcome high dimensionality

	4. Then we benefit Gaussian Process Regression in order to learn our model iteratively such that in each iteration training data are trained, the model is learnt and prediction is done for test data. One compound holding specific property will be added to train data and the progress will repeat until no test data is left.

	Result of this work is accepted in the Journal of Chemical Information and Modeling within the subject  "Predicting Potent Compounds via Model-Based Global Optimization".

}

\value{returns a matrix (m * experiment) of original potencies in test set.
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
SimuChemPC(x, y, "RA", 5)
}
\keyword{chemical, potent compounds, constraint global optimization, expected potency improvement, gaussian process}
