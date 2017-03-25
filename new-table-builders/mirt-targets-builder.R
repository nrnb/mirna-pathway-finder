## mirTarBase Targets Table/Hash Builder
#
# index: mirtarbase
# col1: entrez targets


#######
## LIBS

#install.packages("hash")
library(hash)

#install.packages("shiny")
library(shiny)

#########
## INPUTS

# mirtarebase_entrez targets csv files
targets<-NULL
targets_strong<-"inputs/hsa-mirtarbase-entrez-targets-strong.csv"
targets_functional<-"inputs/hsa-mirtarbase-entrez-targets-functional.csv"

##########
## OUTPUTS

outcsv<-NULL
outcsv_strong<-"outputs/hsa_mirt_targets_strong_hash.csv"
outcsv_functional<-"outputs/hsa_mirt_targets_functional_hash.csv"
outrobj<-NULL
outrobj_strong<-"outputs/hsa_mirt_targets_strong_hash.robj"
outrobj_functional<-"outputs/hsa_mirt_targets_functional_hash.robj"

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

exportHash2Csv <- function(h1,outfile){
    mirtarbase<-keys(h1)
    hits<-paste(values(h1),sep=",")
    #remove list syntax cruft
    hits=lapply(hits,function(y) gsub('[c|(|)|\\"]','',y))
    c<-cbind(mirtarbase,hits)
    d<-as.data.frame(c)
    e <- data.frame(lapply(d, as.character), stringsAsFactors=FALSE)
    write.csv(e, file=outfile,row.names=FALSE)
}


######
## RUN

shinyApp(
    ui = fluidPage(
        titlePanel("miRTarBase Targets Table Builder"),
        p("Builds a table of entrez targets per mirtarbase id."),
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
                    targets<-targets_strong
                    outcsv<-outcsv_strong
                    outrobj<-outrobj_strong
                } else if (choice=="2"){
                    targets<-targets_functional
                    outcsv<-outcsv_functional
                    outrobj<-outrobj_functional          
                }
                
                mirt_targets<-read.csv(targets,stringsAsFactors = FALSE)
                
                mirt_targets_hash<-buildHash(collapse.rows(mirt_targets,1,2))
                setProgress(message = "Exporting csv...")
                exportHash2Csv(mirt_targets_hash,outcsv)
                setProgress(message = "Exporting robj...")
                save(mirt_targets_hash, file=outrobj)
                setProgress(message = "Done!")
                output$status<-renderText(paste("Generated files:","\n",outcsv,"\n",outrobj))
            })
        })
    }
)
