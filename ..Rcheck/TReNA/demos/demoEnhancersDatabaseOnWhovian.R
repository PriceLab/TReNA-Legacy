library(RPostgreSQL)
db <- dbConnect(PostgreSQL(), user= "trena", password="trena", dbname="enhancers.hg38", host="whovian")
dbListTables(db)  # [1] "enhancers"
dbGetQuery(db, "select count(*) from enhancers") # 3,435,042 rows
tbl <- dbGetQuery(db, "select * from enhancers where distance > 219000000")
as.data.frame(t(tbl))
#                              1               2               3               4               5               6
# chr1                      chr1            chr1            chr1            chr1            chr1            chr1
# start1                15341761        15684332        15684390        15684541        15764499        15764953
# end1                  15341761        15684332        15684390        15684541        15764499        15764953
# genename                 FHAD1         PLEKHM2         PLEKHM2         PLEKHM2          FBLIM1          FBLIM1
# transcriptname       FHAD1-004     PLEKHM2-001     PLEKHM2-002     PLEKHM2-004      FBLIM1-003      FBLIM1-004
# transcriptid   ENST00000444385 ENST00000375799 ENST00000375793 ENST00000462455 ENST00000441801 ENST00000332305
# geneid         ENSG00000142621 ENSG00000116786 ENSG00000116786 ENSG00000116786 ENSG00000162458 ENSG00000162458
# chr2                      chr1            chr1            chr1            chr1            chr1            chr1
# start2               234774253       234764253       234764253       234764253       234774253       234774253
# end2                 234784253       234789253       234789253       234789253       234784253       234784253
# distance             219432491       219079920       219079862       219079711       219009753       219009299
# celltype               GM12878            K562            K562            K562         GM12878         GM12878
# method                    Hi-C            Hi-C            Hi-C            Hi-C            Hi-C            Hi-C
# cor                          0               0               0               0               0               0
