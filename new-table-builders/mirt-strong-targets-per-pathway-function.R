## mirTarBase Strong Targets Per Pathway Function
#
# Takes a list of miRTarBase IDs, e.g., query=c("hsa-miR-10a-5p","hsa-miR-10b-3p")
#
# Returns a named list of WPIDs and total unique counts of entrez targets per pathway

#######
## LIBS

#install.packages("hash")
library(hash)

library(foreach)

#########
## INPUTS

# entrez_wp_hash
load("direct_entrez_wp_hash.robj")

# mirt_entreztargets_strong_hash
load("mirt_entreztargets_strong_hash.robj")

############
## FUNCTIONS

lookupListInHash<-function(l, h){
    hits<-NULL
    #identify list items that match hash keys
    hit_vector<-unlist(lapply(l,function(y) !is.null(h[[y]])))
    hit_list<-l[hit_vector]
    hash_hits<-h[hit_list]
    if (!is.empty(hash_hits)){
        inv_hash_hits<-invert(hash_hits)
        hits<-unlist(lapply(values(inv_hash_hits), function(y) length(y)))
    }
    return(hits)
}

mirtStrongTargetsPerPathway<-function(mirtList){
    targetLists<-lookupListInHash(mirtList,mirt_entreztargets_strong_hash)
    wpLists<-lookupListInHash(names(targetLists),entrez_wp_hash)
    return(wpLists)
}

#########
## TESTS

# query=c("hsa-miR-10a-5p","hsa-miR-10b-3p")
# mirtStrongTargetsPerPathway(query)
