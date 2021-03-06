% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TReNA.R
\name{getSolverObject}
\alias{getSolverObject}
\title{Get the Solver Object from a TReNA Object}
\usage{
getSolverObject(obj)
}
\arguments{
\item{obj}{An object of class TReNA}
}
\value{
The Solver object contained by the given TReNA object
}
\description{
Get the Solver Object from a TReNA Object
}
\examples{

# Create a LassoSolver object using the included Alzheimer's data and retrieve the solver object
load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
solver <- TReNA(mtx.sub, solver = "lasso")
mtx <- getSolverObject(solver)

}
