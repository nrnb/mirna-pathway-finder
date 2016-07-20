## Set working directory to the location of this file.
## Change to reflect appropriate path on your computer
setwd("~/Sites/mirna-pathway-finder/wp-mir-table-builder")

#############################################################
## using some code from import_mirbase.R in mirtools by    ##
## James F. Reid <james.reid@ifom-ieo-campus.it>           ##
##                                                         ##
## DEPENDENCIES                                            ##
## sqlite3                                                 ##
##                                                         ##
## NOTES                                                   ##
## For first run, change downloadMirBase to TRUE           ##
## To update releases, change:                             ##
##         mirbaseVersion                                  ##
##         mirbaseData                                     ##
##                                                         ##
#############################################################

#############################################################
## Table of Contents:
##  1. variable names and function definitions
##  2. download data from mirbase.org
##  3. prepare sqlite database (table definitions)
##  4. extract extra information not available in database
##  5. import data into database
#############################################################

## Load dependencies
if (!require("pacman")) install.packages("pacman")
pacman::p_load("RCurl", "RSQLite", "dplyr", "tidyr", "httpuv", "httr", "devtools")

#############################################################
## variable names and paths
##
db                <- "mirbase"
mirbaseVersion    <- "21.0"
mirbaseDate       <- "06 Jun 2014"
mirbaseDomain     <- "mirbase.org"
mirbaseCurrent    <- "pub/mirbase/CURRENT"
dbFTP             <- sprintf("ftp://%s/%s/", mirbaseDomain, mirbaseCurrent)
dbFTPdata         <- sprintf("%sdatabase_files/", dbFTP)

## original data directory
dbRootPath      <- file.path("databases", db)
dbRootDataPath  <- file.path(dbRootPath, "data")

localdbFTP      <- file.path(dbRootPath, sub("ftp://", "", dbFTP))
localdbFTPdata  <- file.path(dbRootPath, sub("ftp://", "", dbFTPdata))

## template directory
dbFile       <- file.path(tempdir(), sprintf("%s.sqlite", db))

dbFilesToExclude <- c("mirna_target_links.txt.gz", "mirna_target_url.txt.gz", "wikipedia.txt.gz")
na_strings <- c("NA", "", "\\N", "\\")
quoting_characters <- c("\"")

## create directories if necc.
if (!file.exists(dbRootDataPath)) dir.create(dbRootDataPath, recursive = TRUE)
if (!file.exists(localdbFTPdata)) dir.create(localdbFTPdata, recursive=TRUE)

downloadMirBase <- FALSE

#############################################################
## functions
##

## download all ".txt.gz", ".str.gz" and ".sql" files from an ftp directory
downloadFtpPath <- function(url, destPath=tempdir()) {
  downloadError <- function(e) stop("Could not download file.\n", e)
  
  ## get directory contents
  fileList <- tryCatch(getURL(url, ftp.use.epsv=FALSE),
                       error=downloadError)
  fileList <- strsplit(fileList, "\r*\n")[[1]]
  ## remove directories
  fileList <- fileList[sapply(fileList, function(i) {
    substr(strsplit(i, " ")[[1]][1], 1, 1)}) != "d"]
  ## ".txt.gz" and ".sql" files only
  fileList <- Filter(function(file_name) {
    return(length(grep("\\.((txt\\.gz)|(str\\.gz)|(sql))$", file_name, perl = T)) == 1)
  }, fileList)
  
  ## names only
  fileList <- sapply(fileList, function(i) {
    tmp <- unlist(strsplit(i, " "))
    tmp[length(tmp)]})
  
  ## exclude files we don't want to import
  fileList <-setdiff(
    fileList,
    dbFilesToExclude)
  
  fileList <- sprintf("%s%s", url, fileList)
  
  ## download each file within directory
  for (ff in fileList) {
    filename <- basename(ff)
    tryCatch(download.file(url=ff,
                           destfile=file.path(destPath, filename),
                           mode="wb"), error=downloadError)
  }
}

## open table data file
readDataFile <- function(file, headers) {
  return(read.delim(file=gzfile(file),
                    header=FALSE,
                    col.names = headers,
                    quote = quoting_characters,
                    stringsAsFactors = FALSE,
                    na.strings = na_strings))
}

## write table data file
writeDataFile <- function(data, f) {
  con <- file(f, open="w")
  write.table(data, file=con, sep="\t", row.names=FALSE, col.names=FALSE)
  close(con)
}

# In a data frame, keep only the rows that have at least one non-NA value
filterOutNAOnlyRows <- function(df) {
  column_names <- colnames(df)
  
  query <- Reduce(function(acc, column_name) {
    subquery <- paste0("(!is.na(", column_name, "))")
    return(paste(c(acc, subquery), collapse=" | "))
  }, column_names, NULL)
  
  result <- df %>%
    filter_(query)
  
  return(result)
}

# The database schema may have fields with a non-null constraint.
# For a data frame that we are about to insert into the database,
# we need to filter out all rows with field values of NA if the
# 
filterOutNARowsWhenFieldNotnull <- function(df, notnull_field_names) {
  column_names <- intersect(
    colnames(df),
    notnull_field_names)
  
  query <- Reduce(function(acc, column_name) {
    subquery <- paste0("(!is.na(", column_name, "))")
    return(paste(c(acc, subquery), collapse=" & "))
  }, column_names, NULL)
  
  result <- df %>%
    filter_(query)
  
  return(result)
}

#############################################################
## download data from miRBase ftp site
##
if (downloadMirBase) {
  cat("Downloading data...\n")
  ## root files
  downloadFtpPath(url=dbFTP, destPath=localdbFTP)
  ## database files (MySQL data dumps)
  downloadFtpPath(url=dbFTPdata, destPath=localdbFTPdata)
}

#############################################################
## create extra tables
##

mirnaHairpinFile <- file.path(dbRootDataPath, "mirna_hairpin.txt.gz")
if (!file.exists(mirnaHairpinFile)) {
  cat("Creating extra 'mirna_hairpin' table.\n")
  ## add folding conformation of stem-loop sequence (not in db data)
  ## created with RNAfold program from the ViennaRNA suite.
  ## Hofacker IL, Stadler PF. Memory efficient folding algorithms for
  ## circular RNA secondary structures. Bioinformatics. 2006 May 15;
  ## 22(10):1172-6. PMID: 16452114
  miRNAstrFile <- file.path(localdbFTP, "miRNA.str.gz")
  con <- gzfile(miRNAstrFile)
  miRNAstr <- readLines(con)
  close(con)
  miRNAstr <- miRNAstr[which(miRNAstr != "")]
  ## scan each entry (six lines per sequence)
  strIndex <- grep(">", miRNAstr)
  if (!all(seq(1, length(miRNAstr), by=6) == strIndex)) {
    stop("Problem scanning ", miRNAstrFile)
  }
  tmp1 <- strsplit(miRNAstr[strIndex], " ")
  ## miRNA IDs
  rfIDs <- sub(">", "", unlist(lapply(tmp1, function(i) i[1])))
  ## minimum free energy
  rfMFE <- as.numeric(sub("\\(", "",
                          sub("\\)", "",
                              unlist(lapply(tmp1, function(i) i[2])))))
  ## mature/minor matches (this info is in mirna_mature table)
  ##tmp2 <- strsplit(miRNAstr[strIndex], "\\[")
  ##rfMMpos <- unlist(lapply(tmp2, function(i) {
  ##    paste(sub("\\]", "", sub(" ", "", i[2:length(i)])), collapse=", ")}))
  ## stem-loop folded sequence
  rfStemLoop <- sapply(strIndex, function(i) {
    paste(miRNAstr[(i+1):(i+5)], collapse="\n")})
  
  x <- data.frame('mirna_id' = rfIDs,
                  ##'mature_pos' = rfMMpos,
                  'hairpin' = rfStemLoop,
                  'mfe' = rfMFE,
                  check.names=FALSE, stringsAsFactors=FALSE)
  writeDataFile(x, mirnaHairpinFile)
}

mirnaClusterFile <- file.path(dbRootDataPath, "mirna_cluster.txt.gz")
if (!file.exists(mirnaClusterFile)) {
  cat("Creating extra 'mirna_cluster' table.\n")
  ## search for 'clustered' mirna in each species
  ## ie. other mirnas <10kb from any given mirna
  clusterWindow <- 10000
  
  ## read-in genomic coordinates
  coords <- readDataFile(file.path(localdbFTPdata,
                                   "mirna_chromosome_build.txt.gz"),
                         c('_id', 'xsome', 'contig_start',
                           'contig_end', 'strand'))

  ## read-in mirna ids
  mirna <- readDataFile(file.path(localdbFTPdata, "mirna.txt.gz"),
                        c('_id', 'mirna_acc', 'mirna_id', 'description',
                          'sequence', 'comment', 'auto_species'))
  
  mirnaCoords <- cbind(mirna[match(coords[, '_id'], mirna[, '_id']),
                             c(1, 3, 7)], coords[, 2:5])
  
  ## initialize mirna_cluster table and cluster id
  mirnaCluster <- data.frame()
  ## iterate over each organism (48 with genomic coordinates - v.14)
  for (orgID in unique(mirnaCoords$auto_species)) {
    cat("\n", orgID, "\n")
    ## list of miRNAs belonging to this species
    mirnaID <- mirnaCoords[mirnaCoords$auto_species == orgID, 'mirna_id']
    ## compute cluster window around each member
    
    ## extract chromosome name
    mirnaChr <- mirnaCoords[(mirnaCoords[, 'mirna_id'] %in% mirnaID),
                            c('mirna_id', 'xsome')]
    ## remove multiple mappings (mirbase.org does NOT do this...)
    mult <- unique(mirnaChr[duplicated(mirnaChr[, 1]), 1])
    mirnaChr <- mirnaChr[!(mirnaChr[, 'mirna_id'] %in% mult), ]
    
    ## skip loop if nothing is left
    if (nrow(mirnaChr) == 0) next
    
    ## append start and end of mirna
    mcMatch <- (mirnaCoords[, 'mirna_id'] %in% mirnaChr[, 'mirna_id'])
    mirnaChr <- cbind(mirnaChr,
                      mirnaCoords[mcMatch, 'contig_start'],
                      mirnaCoords[mcMatch, 'contig_end'])
    ## and start (s) and end (e) of cluster window
    ## (note: strand is not taken into account)
    mirnaChr <- cbind(mirnaChr,
                      c(mirnaChr[, 3] - clusterWindow),
                      c(mirnaChr[, 4] + clusterWindow))
    colnames(mirnaChr) <- c('mirna_id', 'chromosome',
                            'cs', 'ce', 's', 'e')
    ## order by chromosome then by start
    mirnaChr <- with(mirnaChr, mirnaChr[order(mirnaChr[, 'chromosome'],
                                              mirnaChr[, 's']), ])
    ## scan through each chromosome
    for (chr in unique(mirnaChr$chromosome)) {
      cat("\t", chr, ".")
      tmp <- mirnaChr[mirnaChr$chromosome == chr, ]
      if (nrow(tmp) == 1) next
      ## overlap matches
      tmpM <- sapply(tmp$cs, function(mir) {
        (mir > tmp$s & mir < tmp$e)})
      rownames(tmpM) <- colnames(tmpM) <- tmp[, 'mirna_id']
      ## collect results
      clMember <- lapply(1:nrow(tmpM), function(i) {
        rownames(tmpM)[tmpM[, i]]})
      clCluster <- lapply(1:nrow(tmpM), function(i) {
        rep(i, length(clMember[[i]]))})
      clId <- lapply(1:nrow(tmpM), function(i) {
        rep(rownames(tmpM)[i], length(clMember[[i]]))})
      mirnaCluster <- rbind(mirnaCluster,
                            cbind(unlist(clId),
                                  unlist(clMember),
                                  unlist(clCluster)))
    }
  }
  colnames(mirnaCluster) <- c('mirna_id', 'member', 'cluster')
  writeDataFile(mirnaCluster, mirnaClusterFile)
}

#############################################################
## modify original tables
##
cat("Modifying original data tables\n")

cat("Filling empty 'previous_mirna_id' fields with value from 'mirna_id'...\n")
mirna_file <- "mirna.txt.gz"
mirna_tbl <- tbl_df(readDataFile(file.path(localdbFTPdata, mirna_file),
                                 c("auto_mirna", "mirna_acc", "mirna_id", "previous_mirna_id", "description",
                                   "sequence", "comment", "auto_species"))
                    ) %>%
  mutate(previous_mirna_id = ifelse(is.na(previous_mirna_id), mirna_id, previous_mirna_id))
writeDataFile(mirna_tbl, file.path(localdbFTPdata, mirna_file))

cat("Filling empty 'previous_mature_id' fields with value from 'mature_name'...\n")
mirna_mature_file <- "mirna_mature.txt.gz"
mirna_mature_tbl <- tbl_df(readDataFile(file.path(localdbFTPdata, mirna_mature_file),
                                 c("auto_mature", "mature_name", "previous_mature_id", "mature_acc", "evidence",
                                   "experiment", "similarity"))
                           ) %>%
  mutate(previous_mature_id = ifelse(is.na(previous_mature_id), mature_name, previous_mature_id))
writeDataFile(mirna_mature_tbl, file.path(localdbFTPdata, mirna_mature_file))

#############################################################
## create sqlite db
##
cat("Create SQL database...\n")

## initialize database
unlink(dbFile)

# convert mysql dump to sqlite dump and run it to create tables
# using https://github.com/dumblob/mysql2sqlite
sqlite_dump_dest <- file.path(localdbFTPdata, "tables.sqlite")
if (!file.exists(sqlite_dump_dest)) {
  system(paste0("chmod +x ./mysql2sqlite.sh"))
  system(paste0("./mysql2sqlite.sh ", file.path(localdbFTPdata, "tables.sql"), " > ", sqlite_dump_dest))
}
if (!file.exists(dbFile)) {
  system(paste0("sqlite3 ", dbFile, " < ", file.path(localdbFTPdata, "tables.sqlite"))) 
}

# connect to sqlite database
drv <- dbDriver("SQLite")
con <- dbConnect(drv, dbname=dbFile)

dbTables <- intersect(
    db_list_tables(con),
    Map(function(x) {
      return(gsub("\\.txt.gz$", "", x))
    }, list.files(localdbFTPdata, pattern=".*\\.txt.gz$", include.dirs=F)))

## loop over each table
for (tt in dbTables) {
  column_names_expected <- dbListFields(con, tt)
  column_count_expected <- length(column_names_expected)
  
  db_table_metadata <- dbGetQuery(con, paste0("PRAGMA table_info(", tt, ");"))
  
  field_details_df <- tbl_df(list(name=column_names_expected, notnull = db_table_metadata$notnull)) %>%
    mutate(notnull = (notnull == 1))
  
  notnull_field_names <- (field_details_df %>%
                            dplyr::filter(notnull))$name
  
  ff_abs <- file.path(localdbFTPdata, sprintf("%s.txt.gz", tt))
  column_count_observed <- max(Filter(Negate(is.na), count.fields(gzfile(ff_abs),
                                            sep = "\t",
                                            quote = quoting_characters)))
  
  unexpected_column_count <- column_count_observed - column_count_expected
  
  if (unexpected_column_count != 0) {
    unexpected_columns_warning_message <- paste("SQL schema for", tt, "expected",
                                                column_count_expected, "columns, but",
                                                column_count_observed, "observed in data file")
    
    warning(unexpected_columns_warning_message)
  }
  
  column_names_observed <- c()

  if (unexpected_column_count > 0) {
    placeholders <- rep(NA, unexpected_column_count)
    column_names_observed <- c(column_names_expected, placeholders)
  } else {
    if (unexpected_column_count < 0) {
      column_names_observed <- column_names_expected[1:column_count_observed]
      warning(paste("Missing column observed:", column_names_expected[(column_count_observed + 1):column_count_expected], " "))
    } else {
      column_names_observed <- column_names_expected
    }
  }
  column_vector <- c(1:min(column_count_expected, column_count_observed))
  x <- tbl_df(readDataFile(ff_abs, column_names_observed)) %>%
    dplyr::select(column_vector) %>%
    filterOutNARowsWhenFieldNotnull(notnull_field_names)
  
  # insert data
  db_insert_into(con, tt, x)
}

# Get mappings between sequence:
# * names (stem loop and mature)
# * accessions (stem loop and mature)
# * ncbigene identifiers
# * hgnc identifiers
mirnas <- dbGetQuery(con, paste("SELECT mirna_acc, mirna_id, mature_name, mature_acc, ",
                                "db_link, db_id",
                                "FROM mirna ",
                                "INNER JOIN mirna_pre_mature ",
                                    "ON mirna.auto_mirna=mirna_pre_mature.auto_mirna ",
                                "INNER JOIN mirna_database_links ",
                                    "ON mirna.auto_mirna=mirna_database_links.auto_mirna ",
                                "INNER JOIN mirna_mature ",
                                    "ON mirna_pre_mature.auto_mature=mirna_mature.auto_mature ",
                                "WHERE (mirna_database_links.db_id='ENTREZGENE' ",
                                        "OR mirna_database_links.db_id='HGNC')")) %>%
  mutate(row = row_number()) %>%
  rename(mirbase = mirna_acc) %>%
  rename(mirbase.mature = mature_acc) %>%
  rename(stem_loop_name = mirna_id)

## clean and disconnect
dbGetQuery(con, "VACUUM")
dbDisconnect(con)

mappings_mirbase <- full_join(
  mirnas %>%
    filter(db_id == "HGNC") %>%
    spread(db_id, db_link) %>%
    mutate(row = NULL),
  mirnas %>%
    filter(db_id == "ENTREZGENE") %>%
    spread(db_id, db_link) %>%
    mutate(row = NULL),
  by=c("mirbase", "mirbase.mature", "mature_name", "stem_loop_name")) %>%
  rename(hgnc = HGNC) %>%
  rename(ncbigene = ENTREZGENE)

mappings_biomart <- tbl_df(read.table("gene-mir-map.txt",
                               na.strings = na_strings,
                                header=T,
                                sep="\t",
                                fill=T,
                                stringsAsFactors=F)) %>%
  na.omit() %>%
  rename(stem_loop_name = miRBase.ID.s.) %>%
  rename(ncbigene = EntrezGene.ID) %>%
  mutate(ncbigene=as.character(ncbigene)) %>%
  unique()

# in_mirbase_not_biomart <- setdiff(
#   mappings %>%
#     select(stem_loop_name, ncbigene),
#   mappings_ap
# )
# length(in_mirbase_not_biomart$ncbigene)
# 
# in_biomart_not_mirbase <- setdiff(
#   mappings_ap,
#   mappings %>%
#     select(stem_loop_name, ncbigene)
# )
# length(in_biomart_not_mirbase$ncbigene)
# 
# length(mappings_biomart$ncbigene)
# length(mappings_mirbase$ncbigene)
# 
# length(mappings_biomart$ncbigene) - length((mappings_biomart %>% unique())$ncbigene)
# length(mappings_mirbase$ncbigene) - length((mappings_mirbase %>% unique())$ncbigene)

mappings <- full_join(
  mappings_biomart %>%
    rename(ncbigene_biomart = ncbigene),
  mappings_mirbase %>%
    rename(ncbigene_mirbase = ncbigene),
  by=c("stem_loop_name")
) %>%
  mutate(ncbigene = ifelse(is.na(ncbigene_biomart), ncbigene_mirbase, ncbigene_biomart)) %>%
  select(c(stem_loop_name, ncbigene, mirbase, mature_name, mirbase.mature, hgnc)) %>%
  unique()

# length(mappings$ncbigene) - length((mappings %>% unique())$ncbigene)
# 
# mappings_mirbase %>%
#   filter(is.na(ncbigene))

# Given a list of pathways with genes:
# * Convert any miRNA genes to all possible mature sequence names
# * Get all the targets for the pathway genes as mature sequence names,
#   with the targets that each has on the pathway

# Given one or more mature sequence names, get:
# * All pathways with a matching miRNA on the pathway AND/OR a gene targeted by one or more of the inputs
# * For each pathway, include:
#   - any matching miRNA stem loop sequences on the pathway
#   - any targeting miRNAs, each with the genes (on the pathway) that they target

## Download daily gene set collection (WikiPathways data with gene products converted to Entrez Gene);
datestamp <- format(Sys.time(), "%Y%m%d")
wp_gmt_filepath <- paste0("wp_gmt_",datestamp,".txt")
if (!file.exists(wp_gmt_filepath)) {
  wp_gmt_request <- GET(url = "http://www.pathvisio.org/data/bots/gmt/wikipathways.gmt",
                        write_disk(wp_gmt_filepath, overwrite = TRUE))
  wp_gmt_text <- content(wp_gmt_request, "text")
  wp_gmt_datasource <- textConnection(wp_gmt_text)
} else {
  wp_gmt_datasource <- wp_gmt_filepath
}

wp_gmt_file_column_count <- max(count.fields(wp_gmt_filepath, sep = "\t"))
wp_gmt_file_non_gene_column_names <- c("percent_delim_field", "link")
wp_gmt_file_non_gene_column_count <- length(wp_gmt_file_non_gene_column_names)
wp_gmt_file_gene_column_count <- wp_gmt_file_column_count - wp_gmt_file_non_gene_column_count
wp_gmt_file_column_names = c(wp_gmt_file_non_gene_column_names, 1:wp_gmt_file_gene_column_count)
### load file
wp_gmt <- tbl_df(read.table(wp_gmt_datasource,
                               header=F,
                               sep="\t",
                               fill=T,
                               stringsAsFactors=F, 
                               col.names=wp_gmt_file_column_names)) %>%
  separate(col="percent_delim_field",
           into=c("pathway_name", "pathway_release", "wikipathways", "pathway_organism"),
           sep="%") %>%
  filter(pathway_organism == "Homo sapiens") %>%
  gather(gene_header,
         ncbigene,
         num_range("X", c(1:wp_gmt_file_gene_column_count)),
         na.rm=T) %>%
  mutate(ncbigene=as.character(ncbigene)) %>%
  mutate(gene_header = NULL)

### perform mapping
wp_gmt_mapped <- left_join(
  wp_gmt,
  mappings,
  by="ncbigene")

## targets: miRNA (targeter) to gene (target). Technically, the miRNA targets
##          one the products of the specified gene (the messenger RNA), but
##          it's common to say that the miRNA targets the gene.
targets <- left_join(
  tbl_df(read.table(
    "mir-gene-targets.txt",
    colClasses=c("target"="character"),
    na.strings = na_strings,
    header=T,
    sep="\t",
    fill=T,
    stringsAsFactors=F)) %>%
    na.omit() %>%
    unique() %>%
    rename(targeter_mature_name = miRNA) %>%
    rename(ncbigene = target),
  mappings %>%
    setNames(paste0('targeter_', names(.))),
  by=c("targeter_mature_name"))

### perform targeting
wp_gmt_mapped_targeted <- left_join(
    wp_gmt_mapped,
    targets,
    by="ncbigene") %>%
  filter(!is.na(mature_name) | !is.na(targeter_mature_name))

# WP1992_data <- wp_gmt_mapped_targeted %>%
#   filter(wikipathways == "WP1992")
# 
# WP1559_data <- wp_gmt_mapped_targeted %>%
#   filter(wikipathways == "WP1559")
# 
# both_target_and_targeter <- wp_gmt_mapped_targeted %>%
#   filter(!is.na(mature_name) & !is.na(targeter_mature_name))

# Some pathways have a targeter miRNA but have no corresponding
# targets on the pathway, e.g., hsa-mir-21 on pathway
# http://www.wikipathways.org/index.php?title=Pathway:WP1559&oldid=68890
# The column "shown_or_inferred_mature_name" ensures such a targeter
# will be included in targeter_count.
with_all_mature_names <- wp_gmt_mapped_targeted %>%
  mutate(shown_or_inferred_mature_name = ifelse(!is.na(mature_name), mature_name, targeter_mature_name))

write.table(with_all_mature_names,
            file=file.path(paste0("../wp-mir-table-hs-", datestamp, ".tsv")),
            sep="\t",
            row.names=FALSE,
            col.names=TRUE)