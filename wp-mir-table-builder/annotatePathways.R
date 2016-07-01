### Annotate Human Pathways with miRNAs

## Download daily gene set collection;
datestamp=format(Sys.time(), "%Y%m%d")
file=paste0("wp_gmt_",datestamp,".txt")
download.file("http://www.pathvisio.org/data/bots/gmt/wikipathways.gmt",file,method="curl")

#read as lines to resolve tab-delim issue
lines=readLines(file)

#error: read.table does not read this file properly; must construct manually from readLines
#gmtHs=read.table(text=lines[sub],header=F,sep="\t",fill=T,stringsAsFactors=F)

## swap separators for accurate parsing -- all species
#linesAll=gsub("\\\t(http://wikipathways.org/instance/WP[0-9]+)\\\t","~\\1~",lines)
#linesAll2=gsub("\\\t",",",linesAll)
#linesAll2
#gmt=read.table(text=linesAll2,header=F,sep="~",fill=F,stringsAsFactors=F,col.names=c("name","link","genes"))
#gmt[204,]

## swap separators for accurate parsing -- human-only
sub=grep("Homo sapiens",substr(lines,1,500))
length(lines[sub])
## old gmt format: linesHs=gsub("\\(Homo sapiens\\)\\\t(http://wikipathways.org/instance/WP[0-9]+)\\\t","~\\1~",lines[sub])
linesHs=gsub("\\%WikiPathways_.*(http://www.wikipathways.org/instance/WP[0-9]+)_r[0-9]+\\\t","~\\1~",lines[sub])
linesHs2=gsub("\\\t",",",linesHs)
linesHs2
gmtHs=read.table(text=linesHs2,header=F,sep="~",fill=F,stringsAsFactors=F,col.names=c("name","link","genes"))
gmtHs[204,]

## query webservice for synonyms -- slow and incomplete!
#ls=strsplit(gmtHs[204,3],",")
#for(i in 1:length(ls[[1]])){
#    s=paste("wget -qO- http://webservice.bridgedb.org/Human/attributes/L/",ls[[1]][i],"?attrName=Synonyms",sep="")
#    r=system(s,intern=TRUE)
#    print(r)
#}

## load entrez gene to mirbase mapping file
map=read.table("gene-mir-map.txt",na.strings=c("", "NA"),header=T,sep="\t",fill=T,stringsAsFactors=F)
map=na.omit(map)
head(map)

## perform mapping
allghits=list()
allmhits=list()
for(x in 1:nrow(gmtHs)){
    print(x)
ls=strsplit(gmtHs[x,3],",")
ghits=list()
mhits=list()
for(i in 1:length(ls[[1]])){
    if(ls[[1]][i] %in% map$E){
        for(j in 1:length(map$E)){
            if(ls[[1]][i] == map$E[j]){
                ghits=append(ghits,ls[[1]][i])
                mhits=append(mhits,map$m[j])
            }
        }
    }
}
allghits=append(allghits,paste(ghits,collapse=","))
allmhits=append(allmhits,paste(mhits,collapse=","))

}

gmtHs[["ghits"]] = unlist(allghits)
gmtHs[["mhits"]] = unlist(allmhits)

gmtHs[204,]

## load mirbase to gene target file
targets=read.table("mir-gene-targets.txt",na.strings=c("", "NA"),header=T,sep="\t",fill=T,stringsAsFactors=F)
targets=na.omit(targets)
targets=unique(targets)
head(targets)

## perform targeting
#allgthits=list()
allmthits=list()
for(x in 1:nrow(gmtHs)){  
    print(x)
    ls=strsplit(gmtHs[x,3],",")
    #gthits=list()
    mthits=list()
    for(i in 1:length(ls[[1]])){
        if(ls[[1]][i] %in% targets$t){
            for(j in 1:length(targets$t)){
                if(ls[[1]][i] == targets$t[j]){
                    #gthits=append(gthits,ls[[1]][i])
                    mthits=append(mthits,targets$m[j])
                }
            }
        }
    }
    #allgthits=append(allgthits,paste(unique(gthits),collapse=","))
    allmthits=append(allmthits,paste(unique(mthits),collapse=","))
    
}

#gmtHs[["gthits"]] = unlist(allgthits)
gmtHs[["mthits"]] = unlist(allmthits)

gmtHs[204,]
lapply(gmtHs,typeof)

## write out new table
write.csv(gmtHs,file="wp-mir-table-hs.csv",row.names=F)

