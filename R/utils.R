#------------------------------------------------------------------------------------------------------------------------
# "chr10:118441047-118445892" -> list(chrom="chr10", start=118441047, end=118445892)
parseChromLocString <- function(chromLocString)
{
   chromLocString <- gsub(",", "", chromLocString);
   tokens.0 <- strsplit(chromLocString, ":", fixed=TRUE)[[1]]
   stopifnot(length(tokens.0) == 2)
   chrom <- tokens.0[1]
   if(!grepl("chr", chrom))
      chrom <- sprintf("chr%s", chrom)


   tokens.1 <- strsplit(tokens.0[2], "-")[[1]]
   stopifnot(length(tokens.1) == 2)
   start <- as.integer(tokens.1[1])
   end <- as.integer(tokens.1[2])

   return(list(chrom=chrom, start=start, end=end))

} # .parseChromLocString
#------------------------------------------------------------------------------------------------------------------------
# "postgres://whovian/gtf" >- list(brand="postgres", host="whovian", name="gtf"
parseDatabaseUri <- function(database.uri)
{
   topLevel.tokens <- strsplit(database.uri, "://")[[1]]
   database.brand <- topLevel.tokens[1]
   secondLevel.tokens <- strsplit(topLevel.tokens[2], "/(?=[^/]+$)", perl = TRUE)[[1]]
   host <- secondLevel.tokens[1]
   database.name <- secondLevel.tokens[2]

   list(brand=database.brand, host=host, name=database.name)

} # .parseDatabaseUri
#----------------------------------------------------------------------------------------------------
