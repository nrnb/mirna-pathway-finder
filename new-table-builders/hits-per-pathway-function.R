## Hits Per Pathways Function
#
# Parameters:
#   - query: List of identifiers, e.g., query=c("hsa-miR-10a-5p","hsa-miR-10b-3p")
#   - queryType: Type of identifiers. Perform ID mapping to match available types. 
#       queryType="mirtarbase"|"entrez"(default)
#   - hitType: Hash to query. Use "direct" with queryType="entrez" to lookup pathway
#       elements directly. Use "targeting" to identify which query identifiers target
#       elements on a given pathway. Use "targets" to identify the elements targeted by
#       the query set. 
#       hitType="targets"|"targeting"|"direct"(default)
#   - supportType: Strength of evidence supporting mirtarbase interactions.The "strong"
#       type is a subset of "all" functional evidence types. Ignored when hashType="direct".
#       supportType="strong"|"all"(default)
#   - returnType: Values to return in named list. Type "list" returns a named list of 
#       lists of entrez ids; type "count" returns a list of integers.
#       returnType="count"|"list"(default)
#
# Examples:
#   - hitsPerPathway(query) =>returns named list of lists of entrez query ids 
#       found per pathway. 
#   - hitsPerPathway(query,returnType="count") =>returns named list of counts of 
#       entrez query ids found per pathway. 
#   - hitsPerPathway(query,"mirtarbase","targeting") => returns named list of lists of
#       mirtarbase query ids that target each pathway based on all functional evidence.
#   - hitsPerPathway(query,"mirtarbase","targets","strong","count") => returns named list
#       of counts of entrez targets of mirtarbase query ids per pathway based on only
#       strong functional evidence.

#######
## LIBS

#install.packages("hash")
library(hash)


#########
## INPUTS

# entrez_wp_hash
load("outputs/hsa_direct_entrez_wp_hash.robj")

############
## FUNCTIONS

lookupListInHash<-function(l,h,r){
    hits<-NULL
    #identify list items that match hash keys
    hit_vector<-unlist(lapply(l,function(y) !is.null(h[[y]])))
    hit_list<-l[hit_vector]
    hash_hits<-h[hit_list]
    if (!is.empty(hash_hits)){
        inv_hash_hits<-invert(hash_hits)
        if(r=="count"){
            hits<-unlist(lapply(values(inv_hash_hits), function(y) length(y)))
        } else if (r=="list"){
            hits<-values(inv_hash_hits)
        } else {
            stop(cat("Return type",r,"is not supported!"))
        }
    }
    return(hits)
}

hitsPerPathway<-function(query,queryType="entrez",hitType="direct",supportType="all",returnType="list"){
    
    ### MIRTARBASE QUERY ###
    if(queryType=="mirtarbase"){
        if(hitType=="targets"){
            targetLists<-NULL
            if(supportType=="strong"){
                load("outputs/hsa_mirt_targets_strong_hash.robj")
            } else if (supportType=="all"){
                load("outputs/hsa_mirt_targets_functional_hash.robj")
            } else {
                stop(cat("Support type",supportType,"is not supported!"))
            }
            targetLists<-lookupListInHash(query,mirt_targets_hash,"count")
            hitList<-lookupListInHash(names(targetLists),entrez_wp_hash,returnType)
            message(cat("Aggregate",returnType,"of entrez ids targetted by query ids per pathway based on",supportType,"functional miRTarBase evidence:"))
            return(hitList)
        } else if (hitType=="targeting"){
            hitList<-NULL
            if(supportType=="strong"){
                load("outputs/hsa_mirt_targeting_strong_hash.robj")
            } else if (supportType=="all"){
                load("outputs/hsa_mirt_targeting_functional_hash.robj")
            } else {
                stop(cat("Support type",supportType,"is not supported!"))
            }
            hitList<-lookupListInHash(query,mirt_targeting_hash,returnType)
            message(cat("Aggregate",returnType,"of query ids targetting elements per pathway based on",supportType,"functional miRTarBase evidence:"))
            return(hitList)
        } else if (hitType=="direct"){
            stop("Finding direct hits of mirtarbase query ids is not supported yet.")
        }
        
    ### ENTREZ QUERY ###   
    } else if (queryType=="entrez"){
        if(hitType=="targets"){
            stop("Finding targets of entrez query ids is not supported yet.")
        } else if (hitType=="targeting"){
            stop("Finding entrez query ids that target other pathway elements is not supported yet.")
        } else if (hitType=="direct"){
            load("outputs/hsa_direct_entrez_wp_hash.robj")
            hitList<-lookupListInHash(query,entrez_wp_hash,returnType)
            message(cat("Aggregate",returnType,"of query ids represented per pathway:"))
            return(hitList)
        }
    } else {
        stop(cat("Query type",queryType,"is not supported!"))
    }

    stop("The function specified by the args you've provided is not supported!")
}

#########
## TESTS
# 
# query1=c("5294","6885","5728")
# query2=c("hsa-miR-10a-5p","hsa-miR-10b-3p")
# hitsPerPathway(query1)
# hitsPerPathway(query1, returnType="count")
# hitsPerPathway(query2,"mirtarbase","targeting")
# hitsPerPathway(query2,"mirtarbase","targets","strong","count")

## WRITE RESULTING COUNTS TO CSV
library(plyr) 
counts_to_csv <- function(listfordf){
    
    df <- list(list.element = listfordf)
    class(df) <- c("tbl_df", "data.frame")
    attr(df, "row.names") <- .set_row_names(length(listfordf))
    
    if (!is.null(names(listfordf))) {
        df$name <- names(listfordf)
    }
    write.csv(df, file = "~/Desktop/counts.csv")
}


## WRITE RESULTING LISTS TO CSV

list_to_csv <- function(l1){
    h1=hash(l1)
    wpid<-keys(h1)
    hits<-paste(values(h1),sep=",")
    #remove list syntax cruft
    hits=lapply(hits,function(y) gsub('[c|(|)|\\"]','',y))
    c<-cbind(wpid,hits)
    d<-as.data.frame(c)
    e <- data.frame(lapply(d, as.character), stringsAsFactors=FALSE)
    write.csv(e, file="~/Desktop/list.csv",row.names=FALSE)
}

# ## WRITE PATHWAY ANNOT FOR LIST NAMES

wp_annot_to_csv <- function(l1){
    load("outputs/hsa_wp_annot_hash.robj")
    l=names(l1)
    h=wp_annot_hash
    hit_vector<-unlist(lapply(l,function(y) !is.null(h[[y]])))
    hit_list<-l[hit_vector]
    hash_hits<-h[hit_list]
    m1=values(hash_hits)
    write.csv(t(m1), file="~/Desktop/wpannot.csv",row.names=TRUE)
}





