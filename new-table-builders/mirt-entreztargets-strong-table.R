## mirTarBase Entrez Targets Strong Table
#
# index: mirtarbase
# col1: entrez targets


#######
## LIBS

#install.packages("hash")
library(hash)

#########
## INPUTS

# mirtarebase_entrez targets csv files
targets_strong<-"inputs/mirtarbase-entrez-targets-strong.csv"

############
## FUNCTIONS

## libs
source("libs/alleq.R")
source("libs/sunique.R")
source("libs/collapse.rows.R")
source("libs/vector2hashtable.R")

# build hash from row-collapsed dataframe
buildHash<-function(cdf){
    return(hash(keys=cdf[,1],values=lapply(cdf[,2],function(y) unlist(strsplit(y," // ")))))
}

exportHash2Csv <- function(h1){
    mirtarbase<-keys(h1)
    targets<-paste(values(h1),sep=",")
    #remove list syntax cruft
    targets=lapply(targets,function(y) gsub('[c|(|)|\\"]','',y))
    c<-cbind(mirtarbase,targets)
    d<-as.data.frame(c)
    e <- data.frame(lapply(d, as.character), stringsAsFactors=FALSE)
    write.csv(e, file ="mirt_entreztargets_strong_hash.csv",row.names=FALSE)
}


######
## RUN

mirt_entreztargets_strong<-read.csv(targets_strong,stringsAsFactors = FALSE)

mirt_entreztargets_strong_hash<-buildHash(collapse.rows(mirt_entreztargets_strong,1,2))

exportHash2Csv(mirt_entreztargets_strong_hash)
save(mirt_entreztargets_strong_hash, file="mirt_entreztargets_strong_hash.robj")
