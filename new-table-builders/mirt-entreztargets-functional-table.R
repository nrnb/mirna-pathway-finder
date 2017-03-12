## mirTarBase Entrez Targets Functional Table
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
targets_functional<-"inputs/mirtarbase-entrez-targets-functional.csv"

############
## FUNCTIONS

## libs
source("libs/alleq.R")
source("libs/sunique.R")
source("libs/collapse.rows.R")
source("libs/vector2hashtable.R")

# build hash from row-collapsed dataframe
buildHash<-function(cdf){
    return(hash(keys=cdf[,1],values=lapply(cdf[,2],function(y) strsplit(y," // "))))
}

exportHash2Csv <- function(h1){
    mirtarbase<-keys(h1)
    targets<-paste(values(h1),sep=",")
    #remove list syntax cruft
    targets=lapply(targets,function(y) gsub('[c|(|)|\\"]','',y))
    c<-cbind(mirtarbase,targets)
    d<-as.data.frame(c)
    e <- data.frame(lapply(d, as.character), stringsAsFactors=FALSE)
    write.csv(e, file ="mirt_entreztargets_functional_hash.csv",row.names=FALSE)
}


######
## RUN

mirt_entreztargets_functional<-read.csv(targets_functional,stringsAsFactors = FALSE)

mirt_entreztargets_functional_hash<-buildHash(collapse.rows(mirt_entreztargets_functional,1,2))

exportHash2Csv(mirt_entreztargets_functional_hash)
save(mirt_entreztargets_functional_hash, file="mirt_entreztargets_functional_hash.robj")
