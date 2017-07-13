library(TReNA)
library(RSQLite)
library(RPostgreSQL)
library(RUnit)
library(MotifDb); mdb <- MotifDb # convenience variable
library(splitstackshape)
#------------------------------------------------------------------------------------------------------------------------
fp.file <- "/Users/paul/s/work/priceLab/ohsu-aquaporin/aqp4.sqlite"
footprint.db.uri <- sprintf("sqlite://%s", fp.file);
fp.uri <- sprintf("sqlite://%s", fp.file)
genome.db.uri <- "postgres://bddsrds.globusgenomics.org/hg38"
#------------------------------------------------------------------------------------------------------------------------
tbl.mg <- read.table("~/github/TReNA/inst/extdata/motifGenes.tsv", sep="\t", as.is=TRUE, header=TRUE)
#------------------------------------------------------------------------------------------------------------------------
# if we need to do direct inspeciont
#   db.fp <- dbConnect(dbDriver("SQLite"), fp.file)
#------------------------------------------------------------------------------------------------------------------------
tss <- 26865884   # minus strand
tssUpstream <- 500
tssDownstream <- 500
target.gene <- "AQP4"
#------------------------------------------------------------------------------------------------------------------------
if(!exists("mtx")){
  load("~/github/projects/examples/microservices/trenaGeneModel/datasets/coryAD/rosmap_counts_matrix_normalized_geneSymbols_25031x638.RData")
  #load(system.file(package="TReNA", "extdata/ampAD.154genes.mef2cTFs.278samples.RData"))
  mtx <- asinh(mtx)
  mtx.var <- apply(mtx, 1, var)
  deleters <- which(mtx.var < 0.01)
  if(length(deleters) > 0)   # 15838 x 638
     mtx <- mtx[-deleters,]
  mtx <<- mtx
  }
#------------------------------------------------------------------------------------------------------------------------
genome.name <- "hg38"
tbl.snp <- data.frame(target.gene=rep("AQP4", 5),
                      chromosome=rep("chr18", 5),
                      loc=c(26864410, 26865469, 26855623, 26855854, 26850565),
                      snp=c("rs3763040", "rs3875089", "rs335929", "rs3763043", "rs9951307"),
                      shoulder=rep(2000, 5),
                      genome=rep(genome.name, 5),
                      stringsAsFactors=FALSE)
#------------------------------------------------------------------------------------------------------------------------
find.footprints <- function(chrom, start, end)
{
   fp <- FootprintFinder(genome.db.uri, fp.uri, quiet=TRUE)
   tbl.fp <- getFootprintsInRegion(fp, chrom, start, end)


} # find.footprints
#------------------------------------------------------------------------------------------------------------------------
find.snps.in.footprints <- function()
{
   for(r in 1:nrow(tbl.snp)){
      chrom <- tbl.snp$chromosome[r]
      loc   <- tbl.snp$loc[r]
      shoulder <- 20
      query <- sprintf("select * from regions where chrom='%s' and start < %d and endpos > %d",
                       chrom, loc+shoulder, loc-shoulder)
      printf("---- %s: %d", tbl.snp$snp[r], tbl.snp$loc[r])
      print(dbGetQuery(db.fp, query))
      } # for r

} # find.snps.in.footprints
#------------------------------------------------------------------------------------------------------------------------
getFootprintCandidates <- function(target.gene, tssUpstream, tssDownstream)
{
   regionsSpec <- "chr18:26865459-26865479"
   #regionsSpec <- list()
   #geneCenteredSpec <- list(targetGene=target.gene, tssUpstream=tssUpstream, tssDownstream=tssDownstream)
   geneCenteredSpec <- list();


   filterSpec <- list(genomeDB=genome.db.uri,
                      footprintDB=footprint.db.uri,
                      geneCenteredSpec=geneCenteredSpec,
                      regionsSpec=regionsSpec)

   filter <- FootprintFilter(filterSpec$genomeDB,
                             filterSpec$footprintDB,
                             filterSpec$geneCenteredSpec,
                             filterSpec$regionsSpec,
                             quiet=TRUE)

   getCandidates(filter)

} # getFootprintCandidates
#------------------------------------------------------------------------------------------------------------------------
getDHSCandidates <- function(target.gene, tssUpstream, tssDownstream)
{
   genome <- "hg38"
   regionsSpec <- "chr18:26865459-26865479"
   #regionsSpec <- NA_character_;
   #geneCenteredSpec <- list(targetGene=target.gene, tssUpstream=tssUpstream, tssDownstream=tssDownstream)
   geneCenteredSpec <- list()
   #variants <- NA_character_
   variants <- "rs3875089"

   filterSpec <- list(filterType="EncodeDNaseClusters",
                      genomeName=genome,
                      encodeTableName="wgEncodeRegDnaseClustered",
                      pwmMatchPercentageThreshold=75L,
                      geneInfoDB="postgres://whovian/gtf",
                      geneCenteredSpec=geneCenteredSpec,
                      regionsSpec=regionsSpec,
                      variants=variants)

   filter <- with(filterSpec, HumanDHSFilter(genomeName,
                                             encodeTableName=encodeTableName,
                                             pwmMatchPercentageThreshold=85L,
                                             geneInfoDatabase.uri=geneInfoDB,
                                             regionsSpec=regionsSpec,
                                             geneCenteredSpec=geneCenteredSpec,
                                             quiet=FALSE))

   getCandidates(filter)

} # getDHSCandidates
#------------------------------------------------------------------------------------------------------------------------
#createModel <- function(target.gene, tssUpstream, tssDownstream)
createModel <- function(target.gene, chrom, start, end, variants=NA_character_)
{
   #geneCenteredSpec <- list(targetGene=target.gene, tssUpstream=tssUpstream, tssDownstream=tssDownstream)
   geneCenteredSpec <- list()

   #regionsSpec <- list()
   regionsSpec <- sprintf("%s:%d-%d", chrom, start, end)

   filterSpec <- list(genomeDB=genome.db.uri,
                      footprintDB=footprint.db.uri,
                      geneCenteredSpec=geneCenteredSpec,
                      regionsSpec=regionsSpec)

   filter <- FootprintFilter(filterSpec$genomeDB,
                             filterSpec$footprintDB,
                             filterSpec$geneCenteredSpec,
                             filterSpec$regionsSpec,
                             quiet=TRUE)

   x.fp <- getCandidates(filter)


   genome <- "hg38"
   #regionsSpec <- "chr18:26865459-26865479"
   #regionsSpec <- NA_character_;
   #geneCenteredSpec <- list(targetGene=target.gene, tssUpstream=tssUpstream, tssDownstream=tssDownstream)
   #geneCenteredSpec <- list()


   filterSpec <- list(filterType="EncodeDNaseClusters",
                      genomeName=genome,
                      encodeTableName="wgEncodeRegDnaseClustered",
                      pwmMatchPercentageThreshold=75L,
                      geneInfoDB="postgres://whovian/gtf",
                      geneCenteredSpec=geneCenteredSpec,
                      regionsSpec=regionsSpec,
                      variants=variants)

   filter <- with(filterSpec, HumanDHSFilter(genomeName,
                                             encodeTableName=encodeTableName,
                                             pwmMatchPercentageThreshold=90L,
                                             geneInfoDatabase.uri=geneInfoDB,
                                             regionsSpec=regionsSpec,
                                             geneCenteredSpec=geneCenteredSpec,
                                             variants=variants,
                                             quiet=FALSE))

   x.dhs <- getCandidates(filter)

   tfs.fp <- intersect(rownames(mtx), x.fp$tfs)

   if(length(tfs.fp) > 0){
      solver.wt <- RandomForestSolver(mtx, targetGene=target.gene, candidateRegulators=tfs.fp)
      model.fp <- run(solver.wt)
      }
   else{
     printf("no footprint-derived tfs")
     x.fp <- list(tfs=NA_character_, tbl=data.frame())
     model.fp <- list()
     }

   tfs.dhs <- intersect(rownames(mtx), x.dhs$tfs)

   solver.wt <- RandomForestSolver(mtx, targetGene=target.gene, candidateRegulators=tfs.dhs)
   model.dhs <- run(solver.wt)


   return(list(candidates.fp=x.fp, model.fp=model.fp, candidates.dhs=x.dhs, model.dhs=model.dhs))

} # createModel
#------------------------------------------------------------------------------------------------------------------------
createDHSModel <- function(target.gene, chrom, start, end, pwmMatchThreshold=80L, variants=NA_character_)
{
   regionsSpec <- sprintf("%s:%d-%d", chrom, start, end)
   geneCenteredSpec <- list()
   genome <- "hg38"

   filterSpec <- list(filterType="EncodeDNaseClusters",
                      genomeName=genome,
                      encodeTableName="wgEncodeRegDnaseClustered",
                      pwmMatchPercentageThreshold=pwmMatchThreshold,
                      geneInfoDB="postgres://whovian/gtf",
                      geneCenteredSpec=geneCenteredSpec,
                      regionsSpec=regionsSpec,
                      variants=variants)

   filter <- with(filterSpec, HumanDHSFilter(genomeName,
                                             encodeTableName=encodeTableName,
                                             pwmMatchPercentageThreshold=90L,
                                             geneInfoDatabase.uri=geneInfoDB,
                                             regionsSpec=regionsSpec,
                                             geneCenteredSpec=geneCenteredSpec,
                                             variants=variants,
                                             quiet=FALSE))

   x.dhs <- getCandidates(filter)

   tfs.dhs <- intersect(rownames(mtx), x.dhs$tfs)

   solver.wt <- RandomForestSolver(mtx, targetGene=target.gene, candidateRegulators=tfs.dhs)
   model.dhs <- run(solver.wt)

   return(list(candidates=x.dhs, model=model.dhs))

} # createDHSModel
#------------------------------------------------------------------------------------------------------------------------
# 18:26865469 A/G
test.createModel <- function()
{
    # chr18:26,860,742-26,882,685: includes all footprints and dhs clusters
    # chr18:26,860,992-26,870,492: a 10kb region enriched with brain hint footprints
    # chr18:26,864,303-26,866,143: a 1.8kb region including tss
   start <- 26864303 #26860992  #26860742
   end   <- 26866143 #26870492 # 26882685


   m1 <- createModel("AQP4", "chr18", start, end)
   m1.mut <- createModel("AQP4", "chr18", start, end, variants="rs3875089")

    # chr18:26865450-26865480: 30 bases around rs3875089
   start <- 26865450
   end   <- 26865480
   m2 <- createModel("AQP4", "chr18", start, end)
   m2.mut <- createModel("AQP4", "chr18", start, end, variants="rs3875089")

  browser()
  if(length(m$model.fp) > 0){
     head(m$model.dhs$edges)
      #        IncNodePurity  gene.cor
      # TEAD1      39.052237 0.7605383
      # ATF7       19.679876 0.7163899
      # SMAD9      15.450478 0.6694617
      # SP3        13.706970 0.6837160
      # NFE2L2      9.666957 0.6886920
      # GLI2        9.514838 0.5750392

     m$candidates.dhs$tbl[grep("TEAD", m$candidates.fp$tbl$tf),]
     #      chrom    start      end motifName length strand score1  score2   score3                      tf
     # 5851 chr18 26865890 26865897  MA0808.1      8      -      9 12.3265 6.51e-05 TEAD3;TEAD1;TEAD2;TEAD4
     } # if fp results

  head(m$model.dhs$edges)
  m$candidates.dhs$tbl[grep("TEAD", m$candidates.dhs$tbl$tf),]
  tbl.around.snp2 <- subset(m$candidates.dhs$tbl, motifStart<tbl.snp$loc[2] & motifEnd > tbl.snp$loc[2])
  tfs.big <- rownames(subset(m$model.dhs$edges, IncNodePurity > 5))
     # which of these high value tfs have a motif that includes snp2?
  tbl.around.snp2[unlist(lapply(tfs.big, function(tf) grep(tf, tbl.around.snp2$tfs))),]

    #      motifName chrom motifStart motifEnd strand motifScore motifRelativeScore      match regulatoryRegionStart regualtoryRegionEnd   regulatorySequence variant                                                                                                                                                                                 tfs
    # 17    MA0090.2 chr18   26865465 26865474      +   5.653203          0.8174939 AGCATCCCTT              26865450            26865480 AGGATTTGGCTAAAAAG...      wt                                                                                                                                                             TEAD1;TEAD2;TEAD3;TEAD4
    # 111   MA0808.1 chr18   26865466 26865473      +   5.564402          0.8368466   GCATCCCT              26865450            26865480 AGGATTTGGCTAAAAAG...      wt                                                                                                                                                             TEAD3;TEAD1;TEAD2;TEAD4
    # 112   MA0809.1 chr18   26865465 26865474      +   5.470297          0.8060669 AGCATCCCTT              26865450            26865480 AGGATTTGGCTAAAAAG...      wt                                                                                                                                                             TEAD4;TEAD1;TEAD2;TEAD3
    # 16    MA0084.1 chr18   26865468 26865476      -   4.535714          0.7016575  TTTGGCTAA              26865450            26865480 AGGATTTGGCTAAAAAG...      wt                                                                        SRY;SOX30;SOX15;SOX7;SOX17;SOX18;SOX8;SOX9;SOX10;SOX5;SOX6;SOX13;SOX4;SOX11;SOX12;SOX1;SOX2;SOX3;SOX14;SOX21
    # 16.1  MA0084.1 chr18   26865468 26865476      -   4.535714          0.7016575  TTTGGCTAA              26865450            26865480 AGGATTTGGCTAAAAAG...      wt                                                                        SRY;SOX30;SOX15;SOX7;SOX17;SOX18;SOX8;SOX9;SOX10;SOX5;SOX6;SOX13;SOX4;SOX11;SOX12;SOX1;SOX2;SOX3;SOX14;SOX21
    # 56    MA0597.1 chr18   26865467 26865475      +   4.070352          0.7155477  CATCCCTTT              26865450            26865480 AGGATTTGGCTAAAAAG...      wt                                                                                                          THAP1;THAP10;THAP11;PRKRIR;THAP2;THAP3;THAP4;THAP5;THAP6;THAP7;THAP8;THAP9
    # 130   MA0886.1 chr18   26865463 26865472      -   5.230989          0.7422321 GCTAAAAAGC              26865450            26865480 AGGATTTGGCTAAAAAG...      wt EMX2;NANOG;NOTO;VENTX;BSX;HHEX;HLX;EN1;EN2;EMX1;DLX1;DLX2;DLX3;DLX4;DLX5;DLX6;DBX1;DBX2;VAX1;VAX2;TLX1;TLX2;TLX3;BARX1;BARX2;HMX1;HMX2;HMX3;MSX1;MSX2;LBX1;LBX2;BARHL1;BARHL2;NR2E1
    # 87    MA0710.1 chr18   26865463 26865472      -   5.096429          0.7505068 GCTAAAAAGC              26865450            26865480 AGGATTTGGCTAAAAAG...      wt NOTO;NANOG;VENTX;BSX;HHEX;HLX;EN1;EN2;EMX1;EMX2;DLX1;DLX2;DLX3;DLX4;DLX5;DLX6;DBX1;DBX2;VAX1;VAX2;TLX1;TLX2;TLX3;BARX1;BARX2;HMX1;HMX2;HMX3;MSX1;MSX2;LBX1;LBX2;BARHL1;BARHL2;NR2E1
    # 81    MA0699.1 chr18   26865463 26865472      -   4.596502          0.7295830 GCTAAAAAGC              26865450            26865480 AGGATTTGGCTAAAAAG...      wt LBX2;NANOG;NOTO;VENTX;BSX;HHEX;HLX;EN1;EN2;EMX1;EMX2;DLX1;DLX2;DLX3;DLX4;DLX5;DLX6;DBX1;DBX2;VAX1;VAX2;TLX1;TLX2;TLX3;BARX1;BARX2;HMX1;HMX2;HMX3;MSX1;MSX2;LBX1;BARHL1;BARHL2;NR2E1
    # 127   MA0879.1 chr18   26865463 26865472      -   4.369845          0.7033855 GCTAAAAAGC              26865450            26865480 AGGATTTGGCTAAAAAG...      wt DLX1;NANOG;NOTO;VENTX;BSX;HHEX;HLX;EN1;EN2;EMX1;EMX2;DLX2;DLX3;DLX4;DLX5;DLX6;DBX1;DBX2;VAX1;VAX2;TLX1;TLX2;TLX3;BARX1;BARX2;HMX1;HMX2;HMX3;MSX1;MSX2;LBX1;LBX2;BARHL1;BARHL2;NR2E1
    # 122   MA0876.1 chr18   26865464 26865471      -   4.064815          0.7189870   CTAAAAAG              26865450            26865480 AGGATTTGGCTAAAAAG...      wt BSX;NANOG;NOTO;VENTX;HHEX;HLX;EN1;EN2;EMX1;EMX2;DLX1;DLX2;DLX3;DLX4;DLX5;DLX6;DBX1;DBX2;VAX1;VAX2;TLX1;TLX2;TLX3;BARX1;BARX2;HMX1;HMX2;HMX3;MSX1;MSX2;LBX1;LBX2;BARHL1;BARHL2;NR2E1

   # next up: do with the tbl.snp[2,] variant.   does these tfs.big drop out?

} # test.createModel
#------------------------------------------------------------------------------------------------------------------------
createModels <- function(spec, mtx)
{
   stopifnot(all(c("target.gene", "chromosome", "loc", "snp", "shoulder", "genome") %in% names(spec)))

   regionsSpec <- with(spec, sprintf("%s:%d-%d", chromosome, loc-shoulder, loc+shoulder))

   dhsFilterSpec <- list(filterType="EncodeDNaseClusters",
                  genomeName=spec$genome,
                  encodeTableName="wgEncodeRegDnaseClustered",
                  pwmMatchPercentageThreshold=80L,
                  geneInfoDB="postgres://whovian/gtf",
                  geneCenteredSpec=list(),
                  regionsSpec=c(regionsSpec),
                  variants=spec$snp)

   hdf.wt <- with(dhsFilterSpec,
            HumanDHSFilter(genomeName=genomeName,
                           encodeTableName=encodeTableName,
                           pwmMatchPercentageThreshold=pwmMatchPercentageThreshold,
                           geneCenteredSpec=geneCenteredSpec,
                           geneInfoDatabase.uri=geneInfoDB,
                           regionsSpec=regionsSpec,
                           quiet=TRUE))

   x.wt <- getCandidates(hdf.wt)

      # find any motif matches, in regulatory regions, which contain the variant site
   subset(x.wt$tbl, motifStart < spec$loc & motifEnd > spec$loc)

   hdf.mut <- with(dhsFilterSpec,
            HumanDHSFilter(genomeName=genomeName,
                           encodeTableName=encodeTableName,
                           pwmMatchPercentageThreshold=pwmMatchPercentageThreshold,
                           geneCenteredSpec=geneCenteredSpec,
                           geneInfoDatabase.uri=geneInfoDB,
                           regionsSpec=regionsSpec,
                           variants=variants,
                           quiet=TRUE))
   x.mut <- getCandidates(hdf.mut)

     # find any motif matches, in regulatory regions, which contain the variant site
   #subset(x.wt$tbl, motifStart <= spec$loc & motifEnd >= spec$loc)
   #subset(x.mut$tbl, motifStart <= spec$loc & motifEnd >= spec$loc)

   tfs.wt <- intersect(rownames(mtx), x.wt$tfs)
   tfs.mut <- intersect(rownames(mtx), x.mut$tfs)


   identical.tfs <- (length(tfs.wt) == length(tfs.mut)) & all(tfs.wt %in% tfs.mut)
   printf("%d wt tfs, %d mut tfs, identical? %s", length(tfs.wt), length(tfs.mut), identical.tfs)

   model.wt <- data.frame()
   model.mut <- data.frame()


      if(length(tfs.wt) > 0){
        solver.wt <- RandomForestSolver(mtx, targetGene=spec$target.gene, candidateRegulators=tfs.wt)
        model.wt <- run(solver.wt)
        }
      else{
         printf("no candidate tfs for wt")
         }

      if(length(tfs.mut) > 0){
         solver.mut <- RandomForestSolver(mtx, targetGene=spec$target.gene, candidateRegulators=tfs.mut)
         model.mut <- run(solver.mut)
         }
      else{
         printf("no candidate tfs for wt")
         }

   list(wt.candidates=x.wt, mut.candidates=x.mut, wt.model=model.wt, mut.model=model.mut)

} # createModels
#------------------------------------------------------------------------------------------------------------------------
explore.aqp4 <- function()
{
   start <- 26864303 #26860992  #26860742
   end   <- 26866143 #26870492 # 26882685

   start.300bp <- 26865352
   end.300bp   <- 26865619

   start <- start.300bp
   end   <- end.300bp

   #m1 <- createDHSModel("AQP4", "chr18", start, end, 95L)
   tv <- TReNA.Viz(portRange=11011:11051)

   # addBedTrackFromHostedFile(tv,
   #                           trackName="brain HINT",
   #                           uri="http://pshannon.systemsbiology.net/annotations/brain_hint.bed.gz",
   #                           index.uri="http://pshannon.systemsbiology.net/annotations/brain_hint.bed.gz.tbi",
   #                           displayMode="SQUISHED",
   #                           color="blue")

   # addBedTrackFromHostedFile(tv,
   #                           trackName="skin HINT",
   #                           uri="http://pshannon.systemsbiology.net/annotations/skin_hint.bed.gz",
   #                           index.uri="http://pshannon.systemsbiology.net/annotations/skin_hint.bed.gz.tbi",
   #                           displayMode="SQUISHED",
   #                           color="blue")

    addBedTrackFromHostedFile(tv,
                             trackName="EncodeDHSclustered",
                             uri="http://pshannon.systemsbiology.net/annotations/dhsClusters_hg38.bed.gz",
                             index.uri="http://pshannon.systemsbiology.net/annotations/dhsClusters_hg38.bed.gz.tbi",
                             displayMode="SQUISHED",
                             color="blue")


   addBedTrackFromHostedFile(tv,
                            trackName="AQP4 snps",
                            uri="http://pshannon.systemsbiology.net/annotations/aqp4-chr18-snps.bed",
                            index.uri=NA,
                            displayMode="SQUISHED",
                            color="red")

   #tbl.bed <- m1$candidates$tbl[, c("chrom", "motifStart", "motifEnd", "motifName", "motifScore")]
   #addBedTrackFromDataFrame(tv, "DHS motifs", tbl.bed, color="green")

   tv

} # explore.aqp4
#----------------------------------------------------------------------------------------------------
displayAllTracks <- function(chrom, start, end)
{
   if(!exists("tv"))
      tv <<- explore.aqp4()

   if(!exists("db"))
      db <<- dbConnect(PostgreSQL(), user= "trena", password="trena", host="whovian")

   dbNames <- c("brain_hint", "brain_wellington",
                "lymphoblast_hint", "lymphoblast_wellington",
                "skin_hint", "skin_wellington")

   actual.dbNames <- grep("_", dbGetQuery(db, "select datname from pg_database")[, 1], value=TRUE)
   dbNames <- intersect(dbNames, actual.dbNames)

   genome.db.uri <- "postgres://whovian/hg38"
   # chrom <- "chr18"
   #  # based on all 5 snps
   # start.17kb.all.snps <- 26849000
   # end.17kb.all.snps   <- 26866000
   #  # based on tss
   # start <- tss - 2000
   # end   <- tss + 2000
   # start <- 26847227
   # end   <- 26869162


   # start.81kb <- 26824211
   # end.81kb <- 26905511


   # start <- start.81kb
   # end   <- end.81kb

   showRegion(tv, sprintf("%s:%d-%d", chrom, start, end))

   tbl.regions <- data.frame();
   load("aqp4_fp.RData")
   tbl.cory <- fp$AQP4$tbl[, c(3, 4, 5, 1, 15)]
   tbl.cory$source <- "cory"
   tbl.regions <- rbind(tbl.regions, tbl.cory)

   addBedTrackFromDataFrame(tv, "coryFP", tbl.cory, displayMode="COLLAPSED", color="cyan")

   for(dbName in dbNames){
      printf("--- querying %s", dbName)
      fp.db.uri <- sprintf("postgres://whovian/%s", dbName)
      fpf <- FootprintFinder(genome.db.uri, fp.db.uri)
      tbl.fp <- getFootprintsInRegion(fpf, chrom, start, end)
      printf("%40s: %5d footprints", dbName, nrow(tbl.fp))
      if(nrow(tbl.fp) > 0){
         tbl.bed <- tbl.fp[, c("chrom", "start", "endpos", "name", "score2")]
         tbl.bed$source <- dbName
         tbl.regions <- rbind(tbl.regions, tbl.bed)
         addBedTrackFromDataFrame(tv, dbName, tbl.bed, displayMode="COLLAPSED", color="green")
         }
      closeDatabaseConnections(fpf)
      } # for dbName


   addBedTrackFromDataFrame(tv, "fp.all", tbl.regions, color="orange")
   save(tbl.regions, file="all.aqp4.footprints.RData")
   invisible(tbl.regions)

} # displayAllTracks
#----------------------------------------------------------------------------------------------------
getCandidates <- function(tbl.regions)
{
   i <- grep("endpos", colnames(tbl.regions))
   if(length(i) > 0)
      colnames(tbl.regions)[i] <- "end"
   motifMatcher <- MotifMatcher(name="mm", genomeName="hg38", quiet=FALSE)
   #x <- findMatchesByChromosomalRegion(motifMatcher, tbl.regions, pwmMatchMinimumAsPercentage=92)
   tbl.regions.single <- data.frame(chrom="chr18", start=26863884, end=26867884, stringsAsFactors=FALSE)
   x92 =  findMatchesByChromosomalRegion(motifMatcher, tbl.regions.single, pwmMatchMinimumAsPercentage=92)
   x70 =  findMatchesByChromosomalRegion(motifMatcher, tbl.regions.single, pwmMatchMinimumAsPercentage=70)
   print(system.time(x2 <- findMatchesByChromosomalRegion(motifMatcher, tbl.regions.single, pwmMatchMinimumAsPercentage=72)))

   tfs <- intersect(rownames(mtx), x92$tfs)
   solver.wt <- RandomForestSolver(mtx, targetGene="AQP4", candidateRegulators=tfs)
   model.wt <- run(solver.wt)

   snps <- tbl.snp$snp
   x92.snps =  findMatchesByChromosomalRegion(motifMatcher, tbl.regions.single, pwmMatchMinimumAsPercentage=92, variants=snps)



} # getCandidates
#----------------------------------------------------------------------------------------------------
newFriday.rs3875089 <- function()
{
   mm <- MotifMatcher(name="mm", genomeName="hg38", quiet=FALSE)
   rsid <- "rs3875089"
   loc <- subset(tbl.snp, snp==rsid)$loc
   tbl.regions.noSeq <- data.frame(chrom="chr18",
                                   start=loc - 10,
                                   end=loc + 10,
                                   stringsAsFactors=FALSE)
   tbl.wt <- getSequence(mm, tbl.regions.noSeq)
   tbl.mut <- getSequence(mm, tbl.regions.noSeq, rsid)
   checkTrue(tbl.wt$seq != tbl.mut$seq)

   x.wt  <- findMatchesByChromosomalRegion(mm, tbl.regions.noSeq, pwmMatchMinimumAsPercentage=80)
   x.mut <- findMatchesByChromosomalRegion(mm, tbl.regions.noSeq, pwmMatchMinimumAsPercentage=80, rsid)

} # newFriday.rs3875089
#----------------------------------------------------------------------------------------------------
identifyMotifs <- function(chrom="chr18", start=26865324, end=26866747)
{
   mm <- MotifMatcher(name="mm", genomeName="hg38", quiet=FALSE)
   tbl.regions.noSeq <- data.frame(chrom=chrom,
                                   start=start, # min(tbl.snp$loc) - 100,
                                   end=end, #max(tbl.snp$loc) + 5000,
                                   stringsAsFactors=FALSE)
   tbl.wt <- getSequence(mm, tbl.regions.noSeq)
   tbl.mut <- getSequence(mm, tbl.regions.noSeq, tbl.snp$snp)
   checkEquals(nchar(tbl.wt$seq), nchar(tbl.mut$seq))
   for(i in 1:nchar(tbl.mut$seq)){
     wt.base  <- substr(tbl.wt$seq, i, i)
     mut.base <- substr(tbl.mut$seq, i, i)
     if(wt.base != mut.base)
        printf("variant at position %d: %s -> %s", i, wt.base, mut.base)
     } # for i

   x.wt <- findMatchesByChromosomalRegion(mm, tbl.regions.noSeq, pwmMatchMinimumAsPercentage=80)
   x.mut <- findMatchesByChromosomalRegion(mm, tbl.regions.noSeq, pwmMatchMinimumAsPercentage=80, variants=tbl.snp$snp)

   tbl.wt.freq <- as.data.frame(table(x.wt$tbl$motifName))
   tbl.wt.freq <- tbl.wt.freq[order(tbl.wt.freq$Freq, decreasing=TRUE),]
   colnames(tbl.wt.freq) <- c("motif", "wtCount")

   tbl.mut.freq <- as.data.frame(table(x.mut$tbl$motifName))
   tbl.mut.freq <- tbl.mut.freq[order(tbl.mut.freq$Freq, decreasing=TRUE),]
   colnames(tbl.mut.freq) <- c("motif", "mutCount")
     # "MA0751.1"  lost in
   tbl.counts <- merge(tbl.wt.freq, tbl.mut.freq, by="motif", all=TRUE);
   tbl.counts$mutCount[is.na(tbl.counts$mutCount)] <- 0
   tbl.counts$wtCount[is.na(tbl.counts$wtCount)] <- 0

   tbl.counts$diff <- tbl.counts$wtCount - tbl.counts$mutCount
   tbl.counts <- tbl.counts[order(abs(tbl.counts$diff), decreasing=TRUE),]
   tbl.mg <- read.table("~/github/TReNA/inst/extdata/motifGenes.tsv", sep="\t", as.is=TRUE, header=TRUE)
   tbl.counts$motif <- as.character(tbl.counts$motif)
   tfs <- unlist(lapply(tbl.counts$motif, function(m) paste(subset(tbl.mg, motif==m)$tf.gene, collapse=";")))
   tbl.counts$tfs <- tfs
   mapMotifToTF <- function(motif){
      motif.stem <- strsplit(motif, ".", fixed=TRUE)[[1]][1]
      motif.stem.regex <- sprintf("%s\\.", motif.stem)
      tfs <- unique(mcols(query(mdb, motif.stem.regex))$geneSymbol)
      paste(tfs, collapse=";")
      }
   mdb.tfs <- unlist(lapply((tbl.counts$motif), mapMotifToTF))
   tbl.counts$motifDB.tfs <- mdb.tfs
   save(tbl.counts, file="aqp4.tbl.counts.5kbUpstream.RData")

   wt.ma0090.locs <- unique(with(subset(x.wt$tbl, motifName=="MA0090.2"), sprintf("%s:%d-%d", chrom, motifStart, motifEnd)))
   mut.ma0090.locs <- unique(with(subset(x.mut$tbl, motifName=="MA0090.2"), sprintf("%s:%d-%d", chrom, motifStart, motifEnd)))
   setdiff(wt.ma0090.locs, mut.ma0090.locs) # [1] "chr18:26855847-26855856" "chr18:26865465-26865474"

   subset(x.wt$tbl, (motifStart==26855847 | motifStart==26865465) & motifName=="MA0090.2")
    #      motifName chrom motifStart motifEnd strand motifScore motifRelativeScore      match chromStart chromEnd                                      seq status                      tf
    #                                                                                      |
    # 3462  MA0090.2 chr18   26855847 26855856      +   5.734231          0.8292111 AATATTCCAG   26850465 26865569 CCCATATATATGCTCACAATTGATAATTATTCTAATG...     wt TEAD1;TEAD2;TEAD3;TEAD4
    #                                                                                   |
    # 3481  MA0090.2 chr18   26865465 26865474      +   5.653203          0.8174939 AGCATCCCTT   26850465 26865569 CCCATATATATGCTCACAATTGATAATTATTCTAATG...     wt TEAD1;TEAD2;TEAD3;TEAD4

    #                             |               |
    #  A  [ 218  734  212 1132    0  265   29    0  621  280 ]
    #  C  [ 454  100  920   63   10   33 1132 1132  109  346 ]
    #  G  [ 147  398   21    0    4    0    2    1  206  138 ]
    #  T  [ 314   51   18   32 1132  867   24  291  511  367 ]

   browser()
   tbl.snp

    # target.gene chromosome      loc       snp shoulder genome
    #        AQP4      chr18 26864410 rs3763040     2000   hg38
    #        AQP4      chr18 26865469 rs3875089     2000   hg38  A->G (reverse)
    #        AQP4      chr18 26855623  rs335929     2000   hg38
    #        AQP4      chr18 26855854 rs3763043     2000   hg38  A->G (reverse?)
    #        AQP4      chr18 26850565 rs9951307     2000   hg38


} # identifyMotifs
#----------------------------------------------------------------------------------------------------
buildModel <- function()
{
   if(!exists("tbl.count"))
      load("aqp4.tbl.counts.5kbUpstream.RData", envir=.GlobalEnv)

   printf("motifs mapped to genes in motifDB: %d", length(intersect(toupper(tbl.counts$motifDB.tfs), rownames(mtx)))) # 101
   tfs.oldStyle <- unique(unlist(strsplit((tbl.counts$tfs), ";")))
   printf("tfs.oldStyle: %d", length(tfs.oldStyle)) # 714
   printf("motifs mapped to genes old style:  %d", length(intersect(tfs.oldStyle, rownames(mtx))))  # 466

   tfs <- intersect(rownames(mtx), tfs.oldStyle)
   solver.wt <- RandomForestSolver(mtx, targetGene="AQP4", candidateRegulators=tfs)
   model.wt <- run(solver.wt)
   tfs.strong <- rownames(subset(model.wt$edges, IncNodePurity > 2))

} # buildModel
#----------------------------------------------------------------------------------------------------
init.tv <- function(all.tracks=FALSE)
{
   tv <- TReNA.Viz(portRange=11011:11051)

    addBedTrackFromHostedFile(tv,
                             trackName="EncodeDHSclustered",
                             uri="http://pshannon.systemsbiology.net/annotations/dhsClusters_hg38.bed.gz",
                             index.uri="http://pshannon.systemsbiology.net/annotations/dhsClusters_hg38.bed.gz.tbi",
                             displayMode="SQUISHED",
                             color="blue")


   addBedTrackFromHostedFile(tv,
                            trackName="AQP4 snps",
                            uri="http://pshannon.systemsbiology.net/annotations/aqp4-chr18-snps.bed",
                            index.uri=NA,
                            displayMode="SQUISHED",
                            color="red")
   showRegion(tv, "chr18:26,865,340-26,865,658")
   if(all.tracks){
      displayAllTracks()
      } # all.tracks

   tv


} # init.tv
#----------------------------------------------------------------------------------------------------
runBasic <- function(chrom="chr18", start=26865462, end=26865867, min.motif.score=80L)
{
   chrom="chr18"; start=26865462; end=26865867; min.motif.score=80L  # short and simple, just 400bp
   chrom="chr18"; start=26849820; end=26871166   # 21kb
   tv <- init.tv()
   mm <- MotifMatcher(name="mm", genomeName="hg38", quiet=TRUE)
   tbl.regions.noSeq <- data.frame(chrom=chrom,
                                   start=start,   # min(tbl.snp$loc) - 100,
                                   end=end,       #max(tbl.snp$loc) + 5000,
                                   stringsAsFactors=FALSE)
   tbl.wt <- getSequence(mm, tbl.regions.noSeq)
   printf("findMatchesByChromosomalRegion, size: %d", 1 + end - start);
   x.wt <- findMatchesByChromosomalRegion(mm, tbl.regions.noSeq, pwmMatchMinimumAsPercentage=min.motif.score)

   tfs.oldStyle <- unique(unlist(strsplit((x.wt$tbl$tf), ";")))
   printf("tfs.oldStyle: %d", length(tfs.oldStyle)) # 714
   printf("motifs mapped to genes old style:  %d", length(intersect(tfs.oldStyle, rownames(mtx))))  # 466

   tfs <- intersect(rownames(mtx), tfs.oldStyle)
   solver.wt <- RandomForestSolver(mtx, targetGene="AQP4", candidateRegulators=tfs)
   printf("calling randomForestSolver");
   model.wt <- run(solver.wt)
   tfs.strong <- rownames(subset(model.wt$edges, IncNodePurity > 2))

   showRegion(tv, sprintf("%s:%d-%d", chrom, start, end))
   #addBedTrackFromDataFrame(tv, "all motifs", x.wt$tbl[, c(2,3,4,1, 6)], color="grey", displayMode="COLLAPSED")
   for(tf.strong in tfs.strong){
      trackName <- sprintf("%s motifs", tf.strong)
      #addBedTrackFromDataFrame(tv, trackName, x.wt$tbl[grep(tf.strong, x.wt$tbl$tf), c(2,3,4,1, 6)], color="green")
      }

     # expand x.wt$tbl (motifs & their tfs) to one row per motif/tf.  then remove all but tfs.strong
     # the resulting table should be ready for TReNA-Viz::geneRegulatoryModelToGraph

   printf("----- expanding tfs")
   tbl.trimmed <- subset(x.wt$tbl, nchar(tf) != 0)
   tfs.split <- strsplit(tbl.trimmed$tf, ";")
   length(tfs.split) # [1] 36929
   counts <- unlist(lapply(tfs.split, length))
   tfs.split.vec <- unlist(tfs.split)
   tbl.motifs <- expandRows(tbl.trimmed, counts, count.is.col=FALSE, drop=FALSE)
   checkEquals(length(tfs.split.vec), nrow(tbl.motifs))
   tbl.motifs$tf <- tfs.split.vec

   tbl.motifs <- subset(tbl.motifs, tf %in% tfs.strong)
   tbl.motifs$distance.from.tss <- tss - tbl.motifs$motifStart
   printf("----- expanding tfs, done")

   count <- nrow(model.wt$edges)
   tbl.model <- data.frame(tf=rownames(model.wt$edges),
                           randomForest=model.wt$edges$IncNodePurity,
                           pearson=model.wt$edges$gene.cor,
                           spearman=rep(0, count),
                           betaLasso=rep(0, count),
                           pcaMax=rep(0, count),
                           concordance=rep(0, count),
                           stringsAsFactors=FALSE)
   tbl.model <- subset(tbl.model, randomForest >= 2)
   printf("converting model to graph, %d tfs, %d motif sites", nrow(tbl.model), nrow(tbl.motifs))
   system.time(g <- geneRegulatoryModelToGraph(tv, "AQP4", tbl.model, tbl.motifs))
   printf("adding model layout");
   system.time(g.lo <- TReNA:::addGeneModelLayout(g))
   printf("addGraph")
   addGraph(tv, g.lo)
   loadStyle(tv, "style.js")
   browser()
   xyz <- 99

} # runBasic
#----------------------------------------------------------------------------------------------------
runWithAugmentedFootprints <- function()
{
   start.81kb <- 26824211
   end.81kb <- 26905511
   chrom <- "chr18"

   start.300bp <- 26865352
   end.300bp   <- 26865619


   start.1425bp <- 26864212
   end.1425bp <- 26865636

   start <- start.1425bp
   end <- end.1425bp

   tbl.regions <- displayAllTracks(chrom, start, end)
   dim(tbl.regions)  # 88938 x 5

   mm <- MotifMatcher(name="mm", genomeName="hg38", quiet=TRUE)
   tbl.regions.noSeq <- data.frame(chrom=chrom,
                                   start=start,   # min(tbl.snp$loc) - 100,
                                   end=end,       #max(tbl.snp$loc) + 5000,
                                   stringsAsFactors=FALSE)

   tbl.wt <- getSequence(mm, tbl.regions.noSeq)
   printf("findMatchesByChromosomalRegion, size: %d", 1 + end - start);
   x.wt <- findMatchesByChromosomalRegion(mm, tbl.regions.noSeq, pwmMatchMinimumAsPercentage=75)
   #colnames(tbl.regions) <- c("chrom", "start", "endpos",
   x.fp <- findMatchesByChromosomalRegion(mm, tbl.regions, pwmMatchMinimumAsPercentage=75)

} # runWithAugmentedFootprints
#----------------------------------------------------------------------------------------------------
tead1.motifs <- unique(subset(tbl.mg, tf.gene=="TEAD1")$motif)
