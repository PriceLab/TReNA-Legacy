.Cook <- setClass('Cook',
                   slots=c(name="character",
                           recipes="list",
                           quiet="logical")
                  )
#------------------------------------------------------------------------------------------------------------------------
setGeneric("getCooksName",             signature="obj", function(obj) standardGeneric ("getCooksName"))
setGeneric("getRecipes",               signature="obj", function(obj) standardGeneric ("getRecipes"))
setGeneric("createModelsFromRecipes",  signature="obj", function(obj) standardGeneric ("createModelsFromRecipes"))
#------------------------------------------------------------------------------------------------------------------------
Cook <- function(name=NA_character_, recipes=list(), quiet=TRUE)
{
   .Cook(name=name, recipes=recipes, quiet=quiet)

} # Cook constructor
#------------------------------------------------------------------------------------------------------------------------
setMethod("getCooksName", "Cook",
          function(obj){
             return(obj@name)
          })
#------------------------------------------------------------------------------------------------------------------------
setMethod("getRecipes", "Cook",
          function(obj){
             return(obj@recipes)
          })
#------------------------------------------------------------------------------------------------------------------------
setMethod("createModelsFromRecipes", "Cook",

      function(obj){
          result <- list()
          for(recipe in obj@recipes){
             model <- .createGeneModel(recipe)
             name <- recipe@name
             result[[name]] <- model
             } # for recipe
          return(result)
          })

#------------------------------------------------------------------------------------------------------------------------
.createGeneModel <- function(recipe)
{
   candidateFilterSpec <- getCandidateFilterSpec(recipe)
   solverSpec <-getSolverSpec(recipe)

   filterType <- candidateFilterSpec$filterType

   if(filterType == "DNaseFootprints"){
      genome.uri <- candidateFilterSpec$genomeDB
      footprint.uri <- candidateFilterSpec$linkingDataURI
      cFilter <- FootprintFilter(genome.uri, footprint.uri,
                                 candidateFilterSpec$geneCenteredSpec,
                                 candidateFilterSpec$regionsSpec)
      }

   # candidate regulators may have been provided explicitly. far more commonly they
   # are to be calculated by the specified CandidateFilter

   if(is.na(solverSpec$candidateRegulators)){
      candidateInfo <- getCandidates(cFilter)
      tfs <- sort(unique(unlist(candidateInfo$tfs)))
      solverSpec$candidateRegulators <- tfs
      }

   mtx <- .loadAndPrepareAssayMatrix(solverSpec$matrixName)
   solverName <- solverSpec$solver
   targetGene <-

   solver <- NULL

   switch(solverName,
          randomForest={
             printf("randomForest")
             solver <- RandomForestSolver(mtx=mtx,
                                          targetGene=solverSpec$targetGene,
                                          candidateRegulators=solverSpec$candidateRegulators);
             },
          stop(sprintf(printf("Cook.R::.createGeneModel, unrecognized solver: %s", solverName)))
          )

   if(!is.null(solver)){
      solver.out <- run(solver)
      return(solver.out)
      }
   else{
      return(NA)
      }

} # .createGeneModel
#------------------------------------------------------------------------------------------------------------------------
.loadAndPrepareAssayMatrix <- function(matrixName)
{
   if(matrixName == "rosmap") {
      load("~/github/dora/datasets/AD/rosmap_counts_matrix_normalized_geneSymbols_25031x638.RData")
      mtx.asinh <- asinh(mtx)
      medians <- apply(mtx.asinh, 1, median)
      deleters <- as.integer(which(medians <= 0.1))
      mtx <- mtx.asinh[-deleters,]
      return(mtx)
      }

   stop(sprintf("unrecognized assayMatrix name: '%s'", matrixName));

} # .loadAndPrepareAssayMatrix
#------------------------------------------------------------------------------------------------------------------------
.loadCandidateFilter <- function(filterName)
{
   if(tolower(filterName) == "brainhintfootprints"){
      genome.db.uri    <- "postgres://whovian/hg38"              # has gtf and motifsgenes tables
      footprint.db.uri <- "postgres://whovian/brain_hint"        # has hits and regions tables
      fpFilter <- FootprintFilter()
      return(fpFilter)
      }
   else if(tolower(filterName) == "geneontology"){
      if(!exists("org.Hs.eg.db"))
         require("org.Hs.eg.db")
      goFilter <- GeneOntologyFilter(org.Hs.eg.db)
      return(goFilter)
      }

   stop(sprintf("unrecognized candidateFilter name: '%s'", filterName));

} # .loadCandidateFilter
#------------------------------------------------------------------------------------------------------------------------
.createSingleGeneModel <- function(assay.mtx, target.gene, regions, candidateFilter)
{
   regions.string <- paste(regions, collapse="; ")
   printf("--- %s createGeneModel for %s, %s", date(), target.gene, regions.string)

   print("cgm 1")
   if(!target.gene %in% rownames(assay.mtx)){
      msg <- sprintf("no expression data for %s", target.gene);  # todo: pass this back as payload
      print(msg)
      return(list(model=data.frame(), regulatoryRegions=data.frame(), msg=msg))
      }

   print("cgm 2")
   tbl.fp <- data.frame()
   print("cgm 3")

   for(region in regions){
      region.parsed <- extractChromStartEndFromChromLocString(region)
      chrom <- region.parsed$chrom
      start <- region.parsed$start
      end <-   region.parsed$end
      printf("region parsed: %s:%d-%d", chrom, start, end);
      tbl.fp.new <- getFootprintsInRegion(fpf, chrom, start, end)
      tbl.fp <- rbind(tbl.fp, tbl.fp.new)
      filter.args <- list(chrom=chrom, start=start, end=end,
                          region.score.threshold=700,
                          motif.min.match.percentage=90,
                          tableName="wgEncodeRegDnaseClustered")
      tbl.regions <- getCandidates(dhsFilter, filter.args)
      }

   print("cgm 4")
   printf("region: %s", region)
   if(nrow(tbl.fp) == 0){
      msg <- printf("no footprints found within regions: %s", regions.string)
      print(msg)
      return(list(model=data.frame(), regulatoryRegions=data.frame(), msg=msg))
      }
   print("cgm 5")

   printf("range in which fps are requested: %d", end - start)
   printf("range in which fps are reported:  %d", max(tbl.fp$end) - min(tbl.fp$start))

   tbl.fptf <- mapMotifsToTFsMergeIntoTable(fpf, tbl.fp)
   candidate.tfs <- sort(unique(tbl.fptf$tf))
   candidate.tfs <- intersect(rownames(assay.mtx), candidate.tfs)
   goi <- sort(unique(c(target.gene, candidate.tfs)))

   mtx.matched <- assay.mtx[goi,]

   trena.all <- TReNA(mtx.matched, solver="ensemble")

   solvers <- c("lasso", "ridge", "randomForest",  "lassopv", "pearson", "spearman") # "sqrtlasso",
   tbl.model <- solve(trena.all, target.gene, candidate.tfs, extraArgs = list(solver.list=solvers, gene.cutoff=0.1))

   tbl.regRegions <- mergeTFmodelWithRegulatoryRegions(tbl.model, tbl.fp, tbl.fptf, target.gene)

   gene.info <- subset(tbl.tss, gene_name==target.gene)[1,]

   if(gene.info$strand  == "+"){
      gene.start <- gene.info$start
      tbl.regRegions$distance <- gene.start - tbl.regRegions$start
   }else{
      gene.start <- gene.info$endpos
      tbl.regRegions$distance <-  tbl.regRegions$start - gene.start
      }

   printf("after distances calculated, distances to %s tss:  %d - %d", target.gene,
          min(tbl.regRegions$distance), max(tbl.regRegions$distance))
   print(gene.info)
   msg <- sprintf("%d putative TFs found", nrow(tbl.regRegions))
   print(msg)
   tbl.model$target.gene <- target.gene

      # make sure that every tf in the model appears in the regulatoryRegions tf column
      # some regulatory tfs were lost when two motifs were mapped to the same chrom loc

   tfs.model <- sort(tbl.model$gene)
   tfs.reg   <- sort(unique(unlist(strsplit(tbl.regRegions$tf, ";"))))

   stopifnot(all(tfs.model == tfs.reg))

   invisible(list(model=tbl.model, regulatoryRegions=tbl.regRegions, msg=msg))

} # .createSingleGeneModel
#------------------------------------------------------------------------------------------------------------------------
