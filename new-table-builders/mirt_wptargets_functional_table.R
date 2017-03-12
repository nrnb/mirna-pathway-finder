## mirTarBase WikiPathways Targets Functional Table
#
# index: mirtarbase
# col1: WPID targets

#######
## LIBS

#install.packages("hash")
library(hash)

#########
## INPUTS

# entrez_wp hash
load("direct_entrez_wp_hash.robj")

# mirt_entreztargets hash
load("mirt_entreztargets_functional_hash.robj")

############
## FUNCTIONS

lookupEntWp<-function(targets){
    wptargets<-NULL
    #identify which targets match hash keys
    hit_vector<-unlist(lapply(targets,function(y) !is.null(entrez_wp_hash[[y]])))
    hit_targets<-targets[hit_vector]
    #subset hash with vector of hits/keys; error if vector includes non-keys!
    ent_wp_hits<-entrez_wp_hash[hit_targets]
    wp_ent_hits<-invert(ent_wp_hits)
    #flatten into wptargets vector
    if (!is.empty(wp_ent_hits)){
        wptargets<-keys(wp_ent_hits)
    }
    return(wptargets)
}

exportHash2Csv <- function(h1){
    mirtarbase<-keys(h1)
    wptargets<-paste(values(h1),sep=",")
    #remove list syntax cruft
    wptargets=lapply(wptargets,function(y) gsub('[c|(|)|\\"]','',y))
    c<-cbind(mirtarbase,wptargets)
    d<-as.data.frame(c)
    e <- data.frame(lapply(d, as.character), stringsAsFactors=FALSE)
    write.csv(e, file ="mirt_wptargets_functional_hash.csv",row.names=FALSE)
}

######
## RUN

all_list<-lapply(values(mirt_entreztargets_functional_hash), lookupEntWp)
non_null_list<-all_list[!sapply(all_list,is.null)]
mirt_wptargets_functional_hash<-hash(non_null_list)

exportHash2Csv(mirt_wptargets_functional_hash)
save(mirt_wptargets_functional_hash, file="mirt_wptargets_functional_hash.robj")
