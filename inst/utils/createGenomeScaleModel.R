# Grab ensemble score table for all genes in a gene list
#----------------------------------------------------------------------------------------------------
createGenomeScaleModel <- function(mtx.assay, gene.list, genome.db.uri, project.db.uri,
                                   size.upstream=1000, size.downstream=1000,
                                   solver.list = "all.solvers", gene.cutoff = 0.1, num.cores = NULL){

    footprint.filter <- FootprintFilter(mtx.assay = mtx.assay)
    trena <- TReNA(mtx.assay, solver = "ensemble")

    # Setup the parallel structure with a default of half the cores
    if(is.null(num.cores)){
        num.cores <- detectCores()/2}  
    cl <- makeForkCluster(cores = num.cores)
    registerDoParallel(cl)
    #    browser()

#    dumb.test <- foreach(i = 1:10) %dopar% {
#        i+1}
    
    full.result.list <- foreach(i = 1:length(gene.list), .packages='TReNA') %dopar% {

        # Designate the target gene and grab the tfs
        target.gene <- gene.list[[i]]        
        tfs <- getCandidates(footprint.filter,
                             target.gene = target.gene,
                             genome.db.uri = genome.db.uri,
                             project.db.uri = project.db.uri,
                             size.upstream = size.upstream,
                             size.downstream = size.downstream)

        # Solve the trena problem using the supplied values and the ensemble solver
        solve(trena, target.gene, tfs,
                              extraArgs = list("solver.list" = solver.list,
                                               "gene.cutoff" = gene.cutoff))
    }

    # Name the list after the genes supplied
    names(full.result.list) <- gene.list
    return(full.result.list)

} # createGenomeScaleModel
#----------------------------------------------------------------------------------------------------