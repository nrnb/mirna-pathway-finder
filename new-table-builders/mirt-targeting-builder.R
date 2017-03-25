## mirTarBase Targeting Table/Hash Builder
#
# index: mirtarbase
# col1: WPID targets

#######
## LIBS

#install.packages("hash")
library(hash)

#install.packages("shiny")
library(shiny)

#########
## INPUTS

# entrez_wp hash
load("outputs/hsa_direct_entrez_wp_hash.robj")

# mirt_targets hashes
targets_functional<-"outputs/hsa_mirt_targets_functional_hash.robj"
targets_strong<-"outputs/hsa_mirt_targets_strong_hash.robj"

##########
## OUTPUTS

outcsv<-NULL
outcsv_strong<-"outputs/hsa_mirt_targeting_strong_hash.csv"
outcsv_functional<-"outputs/hsa_mirt_targeting_functional_hash.csv"
outrobj<-NULL
outrobj_strong<-"outputs/hsa_mirt_targeting_strong_hash.robj"
outrobj_functional<-"outputs/hsa_mirt_targeting_functional_hash.robj"

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

exportHash2Csv <- function(h1,outfile){
    mirtarbase<-keys(h1)
    wptargets<-paste(values(h1),sep=",")
    #remove list syntax cruft
    wptargets=lapply(wptargets,function(y) gsub('[c|(|)|\\"]','',y))
    c<-cbind(mirtarbase,wptargets)
    d<-as.data.frame(c)
    e <- data.frame(lapply(d, as.character), stringsAsFactors=FALSE)
    write.csv(e, file =outfile,row.names=FALSE)
}

######
## RUN

shinyApp(
    ui = fluidPage(
        titlePanel("miRTarBase Targeting Table Builder"),
        p("Builds a table of pathways with targeted elements per mirtarbase id."),
        br(),
        selectInput(
            "select", 
            label = h5("Which level of supporting evidence do you want to use:"),
            choices = list("Strong functional" = 1, "All functional" = 2), 
            selectize=TRUE
        ),
        actionButton("build","Build Table"),
        hr(),
        fluidRow(column(3, verbatimTextOutput("status")))
    ),
    server = function(input, output) {        
        observeEvent(input$build,{
            withProgress({
                setProgress(message = "Working...")
                choice<-input$select
                if(choice=="1"){
                    load(targets_strong)
                    outcsv<-outcsv_strong
                    outrobj<-outrobj_strong
                } else if (choice=="2"){
                    load(targets_functional)
                    outcsv<-outcsv_functional
                    outrobj<-outrobj_functional          
                }
                
                all_list<-lapply(values(mirt_targets_hash), lookupEntWp)
                non_null_list<-all_list[!sapply(all_list,is.null)]
                mirt_targeting_hash<-hash(non_null_list)
                
                setProgress(message = "Exporting csv...")
                exportHash2Csv(mirt_targeting_hash,outcsv)
                setProgress(message = "Exporting robj...")
                save(mirt_targeting_hash, file=outrobj)
                setProgress(message = "Done!")
                output$status<-renderText(paste("Generated files:","\n",outcsv,"\n",outrobj))
            })
        })
    }
)




