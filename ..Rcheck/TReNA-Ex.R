pkgname <- "TReNA"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('TReNA')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
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
nameEx("getSolverName-LassoSolver-method")
### * getSolverName-LassoSolver-method

flush(stderr()); flush(stdout())

### Name: getSolverName,LassoSolver-method
### Title: Get Lasso Solver name
### Aliases: getSolverName,LassoSolver-method

### ** Examples

solver <- LassoSolver()
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
