% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BayesSpikeSolver.R
\name{BayesSpikeSolver}
\alias{BayesSpikeSolver}
\title{Create a Solver class object using the Bayes Spike Solver}
\usage{
BayesSpikeSolver(mtx.assay = matrix(), targetGene, candidateRegulators,
  nOrderings = 10, quiet = TRUE)
}
\arguments{
\item{mtx.assay}{An assay matrix of gene expression data}

\item{candidateRegulators}{The designated set of transcription factors that could be associated
with the target gene}

\item{nOrderings}{An integer denoting the number of random starts to use for the Bayes Spike
method (default = 10)}

\item{quiet}{A logical denoting whether or not the solver should print output}

\item{target.gene}{A designated target gene that should be part of the mtx.assay data}
}
\value{
A Solver class object with Bayes Spike as the solver
}
\description{
Create a Solver class object using the Bayes Spike Solver
}
\examples{
solver <- BayesSpikeSolver()
}
\seealso{
\code{\link{solve.BayesSpike}}, \code{\link{getAssayData}}

Other Solver class objects: \code{\link{EnsembleSolver}},
  \code{\link{LassoPVSolver}}, \code{\link{LassoSolver}},
  \code{\link{PearsonSolver}},
  \code{\link{RandomForestSolver}},
  \code{\link{RidgeSolver}}, \code{\link{Solver-class}},
  \code{\link{SpearmanSolver}},
  \code{\link{SqrtLassoSolver}}
}
