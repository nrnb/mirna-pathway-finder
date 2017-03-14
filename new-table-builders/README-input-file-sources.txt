################
## SOURCE INPUTS


* wikipathways-YYYYMMDD-gmt-Homo_sapiens -- downloaded via http://data.wikipathways.org/current/gmt/. Contains pathway names, links and gene sets unified to entrez gene.

* hsa_MTI_r#.#.xlsx -- downloaded from http://mirtarbase.mbc.nctu.edu.tw/php/download.php. Contains miRTarBase IDs for targeting miRNAs and Entrez for gene targets. Also includes strong/weak and pmid refs. PROCESSING: Release number added to filename.

* [IGNORE] mart_export_v##.csv -- exported from BioMART at http://www.ensembl.org/biomart/. See screenshot "mart_export.png". Filter for lincRMA, miRNA and processed_transcript; Attributes: entrez and mirbase ids. Export as CSV, unique results only. PROCESSING: Version number added to filename and extension set to CSV. Removed rows with blanks in Excel and renamed header rows: "entrez", "mirbase"

* [OPTIONAL] mirna_mature.txt -- downloaded from ftp://mirbase.org/pub/mirbase/CURRENT/database_files/ (log in as Guest). This file can be processed into a CSV to be used to prepare a _direct_ mapping from mirtarbase ids to pathway elements, i.e., via an additional MIMAT->Entrez mapping set by BridgeDb.  PROCESSING: Open in Excel. Filter for mirtarbase mature names that begin with "hsa"; copy mirtarbase and mimat id columns to new tab and export as mirna_mature.csv.


###################
## PROCESSED INPUTS

* mirtarbase-entrez-targets-***.csv
  1. open hsa_MTI.xlsx in Excel
  2. filter for strong support type: "Functional MTI" (i.e., not weak)
  3. copy columns to new tab called "strong": miRNA, Target Gene (Entrez)
  4. filter for both strong and weak support types: "Functional MTI" or "Functional MIT (Weak)"
  5. copy columns to new tab called "funcional": miRNA, Target Gene (Entrez)
  6. confirm paste lengths, remove rows with blanks and rename column header rows: "mirtarbase", "entrez targets"
  6. save each tab as: mirtarbase-entrez-targets-(strong|funcional).csv
  7. save hsa_MTI.xlsx

* mirt_targets_***_hash.csv/robj
  1. source mirt_targets-builder.R (shiny)

* mirt_targeting_***_hash.csv/robj
  1. source mirt-targeting-builder.R (shiny)
  

