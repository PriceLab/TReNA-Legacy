# Using the PostgreSQL databases on Whovian, grab a subset and create a SQLite database for test_TReNA
#----------------------------------------------------------------------------------------------------
library(RPostgreSQL)
library(RSQLite)

# Specify the database connection parameters and connect
driver <- PostgreSQL()
genome.dbname <- "hg38"
project.dbname <- "brain_wellington"
host <- "whovian"
genome.db <- dbConnect(driver, user = "trena", password = "trena", dbname = genome.dbname, host = host)
project.db <- dbConnect(driver, user = "trena", password = "trena", dbname = project.dbname, host = host)

# Specify the query parameters for genome/project dbs, then run the queries
target.gene <- "MEF2C"

# Note that we need the following tables for the FootprintFilter functionality:
# motifsgenes (genome.db; for a given set of motifs)
# gtf (genome.db; for a given gene and type)
# regions (project.db; for a given chromasome, start, and end positions)
# hits (project.db; for a given sent of locations from the regions table)
# 

genome.gtf.query <- sprintf("select * from gtf where gene_name ='%s' and moleculetype = 'gene'",target.gene)
tbl.genome.gtf <- dbGetQuery(genome.db, genome.gtf.query)

chrom <- tbl.genome.gtf$chr[1]
if(tbl.genome.gtf$strand[1] == "+") {
    tss <- tbl.genome.gtf$start[1]
} else {tss <- tbl.genome.gtf$end[1]
}

shoulder.size <- 10000 # Look at 10Kb upstream and downstream
    
project.regions.query <- sprintf("select * from regions where chrom = '%s' and start > %d and start < %d",
                                 chrom, tss - shoulder.size, tss + shoulder.size)
tbl.project.regions <- dbGetQuery(project.db, project.regions.query) # 60 rows

loc.set <- sprintf("('%s')", paste(tbl.project.regions$loc, collapse="','"))
project.hits.query <- sprintf("select * from hits where loc in %s", loc.set)
tbl.project.hits <- dbGetQuery(project.db, project.hits.query) # 90 rows

motifs <- unique(c(tbl.project.regions$name, tbl.project.hits$name))
collected.motifs <- sprintf("('%s')", paste(motifs, collapse="','"))
genome.motifsgenes.query <- sprintf("select * from motifsgenes where motif in %s", collected.motifs)
tbl.genome.motifsgenes <- dbGetQuery(genome.db, genome.motifsgenes.query)

# Close the PostgreSQL connections
dbDisconnect(genome.db)
dbDisconnect(project.db)

# Create the SQLite connections and write the tables
genome.con <- dbConnect(SQLite(), dbname = "genome.sub.db")
project.con <- dbConnect(SQLite(), dbname = "project.sub.db")

dbWriteTable(genome.con, "motifsgenes", tbl.genome.motifsgenes, overwrite = TRUE)
dbWriteTable(genome.con, "gtf", tbl.genome.gtf, overwrite = TRUE)
dbWriteTable(project.con, "regions", tbl.project.regions, overwrite = TRUE)
dbWriteTable(project.con, "hits", tbl.project.hits, overwrite = TRUE)

# Finish by closing the connection
dbDisconnect(genome.con); dbDisconnect(project.con)
