## WPID Annotation Table
#
# index: WPID
# col1: Title
# col2: Version
# col3: Species
# col4: URL

#######
## LIBS

#install.packages("hash")
library(hash)

#########
## INPUTS

# gmt file
fname<-"wikipathways-20170210-gmt-Homo_sapiens.gmt"


############
## FUNCTIONS

processGMT <- function(fname) {
    tmp = readLines(fname)
    tmp2 = sapply(tmp,strsplit,"\t|%")
    # keep the first five columns
    tmp3=lapply(tmp2,function(y) y[1:5])
    # remove any list items with length 0
    # These were blank lines in the original file
    tmp4 = tmp3[sapply(tmp3,length)>0]
    return(tmp4)
}

buildHash <- function(lst1){
    h1 = hash(keys=sapply(lst1,'[[',3), values=lapply(lst1,function(y) y[-3]))
    return(h1)
}

exportList2Csv <- function(l){
    d<-as.data.frame(l)
    e <- data.frame(lapply(d, as.character),row.names=c("title","version","wpid","species","url"), stringsAsFactors=FALSE)
    write.csv(t(e), file ="wp_annot_hash.csv",row.names=FALSE)
}

######
## RUN

wp_annot_list<-processGMT(fname)

wp_annot_hash<-buildHash(wp_annot_list)

exportList2Csv(wp_annot_list)
save(wp_annot_hash, file="wp_annot_hash.robj")
