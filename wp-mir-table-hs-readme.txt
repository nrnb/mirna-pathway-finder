### Description of wp-mir-table-hs.csv

## Script

annotatePathways.R starts with a tab-delimited gmt file from wikipathway.org, a mapping file for entrez gene to mirtarbase IDs, and a target interaction file for mirtarbase IDs to entrez genes, and then produces a table of pathways with columns of mirtarbase IDs that map directly and indirectly (via target interactions).

## Input data

* wp_gmt_YYYYMMDD.txt -- downloaded via script from http://wikipathways.org/index.php/Download_Pathways#GMT_Gene_Sets. Contains pathway names, links and gene sets unified to entrez gene.

* gene-mir-map.txt -- generated at biomart: latest release (r81) of human genome, filter for all genes with Entrez Gene annotation, annotate with mirtarbase IDs. Downloaded table and edited in Excel to keep complete, unique rows of Entrez Gene and mirtarbase ID columns.

* mir-gene-targets.txt -- from hsa_MTI.xlsx (http://mirtarbase.mbc.nctu.edu.tw/php/download.php). Deleted all columns except: mirtarbase IDs and Entrez Gene targets. Script reads and makes unique.

## Output

* wp-mir-table-hs.csv -- table of human pathways, links, and comma-separated lists of:
	* ghits -- Entrez Gene IDs of miRNAs on pathway
	* mhits -- mirtarbase IDs of miRNAs on pathway, paired 1:1 with ghits
	* mthits -- mirtarbase IDs of miRNAs that target genes on pathway

