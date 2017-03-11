## Direct Table
#
# index: Entrez ID
# col1: list of WPIDs 

#######
## LIBS

#install.packages("hash")
library(hash)

#########
## INPUTS

# gmt file
fname<-"inputs/wikipathways-20170210-gmt-Homo_sapiens.gmt"


############
## FUNCTIONS

processGMT <- function(fname) {
    tmp = readLines(fname)
    tmp2 = sapply(tmp,strsplit,"\t|%")
    # removing the first two columns 
    tmp3 = sapply(tmp2,'[',-1:-2)
    tmp4 = sapply(tmp3,'[',-2:-3)
    # remove any list items with length 0
    # These were blank lines in the original file
    tmp5 = tmp4[sapply(tmp4,length)>0]
    return(tmp5)
}

buildInvertedHash <- function(lst1){
    h1 = hash(keys=sapply(lst1,'[[',1), values=sapply(lst1,'[',-1))
    h2 = invert(h1)
    return(h2)
}

exportHash2Csv <- function(h1){
    entrez<-keys(h1)
    wpids<-paste(values(h1),sep=",")
    #remove list syntax cruft
    wpids=lapply(wpids,function(y) gsub('[c|(|)|\\"]','',y))
    c<-cbind(entrez,wpids)
    d<-as.data.frame(c)
    e <- data.frame(lapply(d, as.character), stringsAsFactors=FALSE)
    write.csv(e, file ="direct_entrez_wp_hash.csv",row.names=FALSE)
}

######
## RUN

wp_entrez_list<-processGMT(fname)

entrez_wp_hash<-buildInvertedHash(wp_entrez_list)

exportHash2Csv(entrez_wp_hash)
save(entrez_wp_hash, file="direct_entrez_wp_hash.robj")
