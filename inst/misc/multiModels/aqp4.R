library(TReNA)
library(RSQLite)
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
createModel <- function(target.gene, chrom, start, end)
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
   variants <- NA_character_
   #variants <- "rs3875089"

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
                                             pwmMatchPercentageThreshold=70L,
                                             geneInfoDatabase.uri=geneInfoDB,
                                             regionsSpec=regionsSpec,
                                             geneCenteredSpec=geneCenteredSpec,
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
test.createModel <- function()
{
  m <- createModel("AQP4", "chr18", 26865450, 26865480) # around tbl.snp[2,]
  head(m$model.fp$edges)

    #        IncNodePurity  gene.cor
    # TEAD1      39.052237 0.7605383
    # ATF7       19.679876 0.7163899
    # SMAD9      15.450478 0.6694617
    # SP3        13.706970 0.6837160
    # NFE2L2      9.666957 0.6886920
    # GLI2        9.514838 0.5750392

   m$candidates.fp$tbl[grep("TEAD", m$candidates.fp$tbl$tf),]
    #      chrom    start      end motifName length strand score1  score2   score3                      tf
    # 5851 chr18 26865890 26865897  MA0808.1      8      -      9 12.3265 6.51e-05 TEAD3;TEAD1;TEAD2;TEAD4

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
   x <- createModel(target.gene, tssUpstream, tssDownstream)

   chrom <- "chr18"
   start <- tss - 1000
   end   <- tss + 1000
   chromLoc <- sprintf("%s:%d-%d", chrom, start, end)

   tbl.snp$loc - tss # [1]  -1474   -415 -10261 -10030 -15319
   filterSpec <- list(genomeDB=genome.db.uri,
                      footprintDB=footprint.db.uri,
                      geneCenteredSpec=list(targetGene=target.gene,
                         tssUpstream=500,
                         tssDownstream=500),
                     regionsSpec=list())

   filter <- FootprintFilter(filterSpec$genomeDB,
                             filterSpec$footprintDB,
                             filterSpec$geneCenteredSpec,
                             filterSpec$regionsSpec,
                             quiet=TRUE)

   list.out <- getCandidates(filter)



   # rs3763040: 18:26864410 A/C/T
   # rs3875089: 18:26865469 A/G
   # rs335929:  18:26855623 A/C/G
   # rs3763043: 18:26855854 A/G
   # rs9951307: 18:26850565 A/G


browser()

result.2k <- list()
#for(r in 1:nrow(tbl.snp)){
   r <- 2
   gene <- tbl.snp$target.gene[r]
   result.2k[[r]] <- createModels(as.list(tbl.snp[r,]), mtx)
   printf("back from createModels, in main loop")
   xyz <- 99
#   }


   recipe <- list(genomeDB=genome.db.uri,
                  footprintDB=footprint.db.uri,
                  geneCenteredSpec=list(),
                  regionsSpec=c(chromLoc))

   filter <- FootprintFilter(recipe$genomeDB,
                             recipe$footprintDB,
                             recipe$geneCenteredSpec,
                             recipe$regionsSpec,
                             quiet=TRUE)

   list.out <- getCandidates(filter)
   tbl.fp <- list.out$tbl
   candidates <- list.out$tfs


} # explore.aqp4
#----------------------------------------------------------------------------------------------------

tead1.motifs <- unique(subset(tbl.mg, tf.gene=="TEAD1")$motif)
