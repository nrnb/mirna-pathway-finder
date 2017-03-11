## WPID Annotation Table
#
# index: WPID
# col1: Title
# col2: Version
# col3: Species
# col4: URL
# col5: Total gene count

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

addCounts<- function(tmp,cnt){
    i=1
    while(i<=length(tmp)){
        tmp[[i]][[6]]=cnt[[i]]
        i=i+1
    }
    return(tmp)
}

processGMT <- function(fname) {
    tmp = readLines(fname)
    tmp2 = sapply(tmp,strsplit,"\t|%")
    # count genes
    cnt = lapply(tmp2,function(y) length(y)-5)
    # keep the first five columns then add count column
    tmp3=lapply(tmp2,function(y) y[1:5])
    tmp4=addCounts(tmp3,cnt)
    # remove any list items with length 0
    # These were blank lines in the original file
    tmp5 = tmp4[sapply(tmp4,length)>0]
    return(tmp5)
}

buildHash <- function(lst1){
    h1 = hash(keys=sapply(lst1,'[[',3), values=lapply(lst1,function(y) y[-3]))
    return(h1)
}

exportList2Csv <- function(l){
    d<-as.data.frame(l)
    e<-data.frame(lapply(d, as.character),row.names=c("title","version","wpid","species","url","total gene count"), stringsAsFactors=FALSE)
    #reorder rows
    e<-e[c("wpid","title","version","species","url","total gene count"),]
    #write transposed df
    write.csv(t(e), file ="wp_annot_hash.csv",row.names=FALSE)
}

######
## RUN

wp_annot_list<-processGMT(fname)

wp_annot_hash<-buildHash(wp_annot_list)

exportList2Csv(wp_annot_list)
save(wp_annot_hash, file="wp_annot_hash.robj")
