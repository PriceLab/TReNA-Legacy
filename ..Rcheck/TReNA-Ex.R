pkgname <- "TReNA"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('TReNA')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("LassoPVSolver")
### * LassoPVSolver

flush(stderr()); flush(stdout())

### Name: LassoPVSolver
### Title: Create a Solver class object using the LASSO P-Value solver
### Aliases: LassoPVSolver

### ** Examples

solver <- LassoPVSolver()



cleanEx()
nameEx("LassoSolver")
### * LassoSolver

flush(stderr()); flush(stdout())

### Name: LassoSolver
### Title: Create a Solver class object using the LASSO solver
### Aliases: LassoSolver

### ** Examples

solver <- LassoSolver()



cleanEx()
nameEx("PearsonSolver")
### * PearsonSolver

flush(stderr()); flush(stdout())

### Name: PearsonSolver
### Title: Create a Solver class object using Pearson correlation
###   coefficients as the solver
### Aliases: PearsonSolver

### ** Examples

solver <- PearsonSolver()



cleanEx()
nameEx("RidgeSolver")
### * RidgeSolver

flush(stderr()); flush(stdout())

### Name: RidgeSolver
### Title: Create a Solver class object using the Ridge Regression solver
### Aliases: RidgeSolver

### ** Examples

solver <- RidgeSolver()



cleanEx()
nameEx("SpearmanSolver")
### * SpearmanSolver

flush(stderr()); flush(stdout())

### Name: SpearmanSolver
### Title: Create a Solver class object using Spearman correlation
###   coefficients as the solver
### Aliases: SpearmanSolver

### ** Examples

solver <- SpearmanSolver()



cleanEx()
nameEx("SqrtLassoSolver")
### * SqrtLassoSolver

flush(stderr()); flush(stdout())

### Name: SqrtLassoSolver
### Title: Create a Solver class object using the Square Root LASSO solver
### Aliases: SqrtLassoSolver

### ** Examples

solver <- SqrtLassoSolver()



cleanEx()
nameEx("getSolverName-LassoPVSolver-method")
### * getSolverName-LassoPVSolver-method

flush(stderr()); flush(stdout())

### Name: getSolverName,LassoPVSolver-method
### Title: Get LassoPV Solver name
### Aliases: getSolverName,LassoPVSolver-method

### ** Examples

solver <- LassoPVSolver()
getSolverName(solver)



cleanEx()
nameEx("getSolverName-LassoSolver-method")
### * getSolverName-LassoSolver-method

flush(stderr()); flush(stdout())

### Name: getSolverName,LassoSolver-method
### Title: Get Lasso Solver name
### Aliases: getSolverName,LassoSolver-method

### ** Examples

solver <- LassoSolver()
getSolverName(solver)



cleanEx()
nameEx("getSolverName-PearsonSolver-method")
### * getSolverName-PearsonSolver-method

flush(stderr()); flush(stdout())

### Name: getSolverName,PearsonSolver-method
### Title: Get Pearson Solver name
### Aliases: getSolverName,PearsonSolver-method

### ** Examples

solver <- PearsonSolver()
getSolverName(solver)



cleanEx()
nameEx("getSolverName-RidgeSolver-method")
### * getSolverName-RidgeSolver-method

flush(stderr()); flush(stdout())

### Name: getSolverName,RidgeSolver-method
### Title: Get Ridge Solver name
### Aliases: getSolverName,RidgeSolver-method

### ** Examples

solver <- RidgeSolver()
getSolverName(solver)



cleanEx()
nameEx("getSolverName-SpearmanSolver-method")
### * getSolverName-SpearmanSolver-method

flush(stderr()); flush(stdout())

### Name: getSolverName,SpearmanSolver-method
### Title: Get Spearman Solver Name
### Aliases: getSolverName,SpearmanSolver-method

### ** Examples

solver <- SpearmanSolver()
getSolverName(solver)



cleanEx()
nameEx("getSolverName-SqrtLassoSolver-method")
### * getSolverName-SqrtLassoSolver-method

flush(stderr()); flush(stdout())

### Name: getSolverName,SqrtLassoSolver-method
### Title: Get SqrtLasso Solver name
### Aliases: getSolverName,SqrtLassoSolver-method

### ** Examples

solver <- SqrtLassoSolver()
getSolverName(solver)



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
