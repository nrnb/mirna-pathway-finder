## WPID Annotation Table
#
# index: WPID
# col1: Title
# col2: Version
# col3: Species
# col4: URL
# col5: Total gene count
# col6: Total targeting mirna count

#######
## LIBS

#install.packages("hash")
library(hash)

library(plyr)

#########
## INPUTS

# gmt file
fname<-"inputs/wikipathways-20170210-gmt-Homo_sapiens.gmt"

# mirt_targeting_hash
load("outputs/hsa_mirt_targeting_strong_hash.robj")

##########
## OUTPUTS

outcsv<-"outputs/hsa_wp_annot_hash.csv"
outrobj<-"outputs/hsa_wp_annot_hash.robj"

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

# build hash from list of lists, using 3rd column as key
buildHash <- function(lst1){
    return(hash(keys=sapply(lst1,'[[',3), values=lapply(lst1,function(y) y[-3])))
}

# merges hashes and exports csv
mergeHashes<- function(h1,h2){
    l1=as.list(h1)
    l2=as.list(invert(h2))
    # count mirna per pathway
    c2=lapply(l2,function(y) length(y))
    exportLists2Csv(l1,c2)
    # gather unique set of keys
    keys=unique(c(names(l1),names(c2)))
    # build new combined list
    l3=setNames(mapply(c,l1[keys],c2[keys]),keys)
    #l4=do.call(rbind, lapply(names(l3), function(u) transform(l3[1][u], type=u)))
    return(hash(l3))
}

exportLists2Csv <- function(l1,c2){
    d1<-as.data.frame(l1)
    d2<-as.data.frame(c2)
    #e1<-data.frame(lapply(d1, as.character),row.names=c("title","version","species","url","total gene count"), stringsAsFactors=FALSE)
    #e2<-data.frame(lapply(d2, as.character),row.names=c("total mirt count"), stringsAsFactors=FALSE)
    # combind dfs
    e<-rbind.fill(d1,d2)
    f<-data.frame(lapply(e, as.character),stringsAsFactors=FALSE)
    ft<-as.data.frame(t(f))
    g<-transform(ft,wpid=colnames(f))
    g<-g[,c(7,1,2,3,4,5,6)]
    names(g)<-c("wpid","title","version","species","url","total gene count","total mirt count")
    #write transposed df
    write.csv(g, file =outcsv,row.names=FALSE)
}

######
## RUN

wp_gmt_list<-processGMT(fname)

wp_gmt_hash<-buildHash(wp_gmt_list)

wp_annot_hash<-mergeHashes(wp_gmt_hash,mirt_targeting_hash)

save(wp_annot_hash, file=outrobj)
