pair.epi.retrieval <- function(query, key){
    require(stringr) #str_split_fixed
    require(genbankr) #readGenBank
    require(rentrez) #set_entrez_key
    query<- read.csv(query, header = FALSE, col.names=c("id","source","recipient")) #load data 
    set_entrez_key(key) #set API key
    print(paste(nrow(query), "Transmission cluster(s) to process"))
    print(paste("Using entrez key ", key, sep=""))
    ids <<- character() #return vector to store ids
    for(id in 1:nrow(query)){ #loop clusters
        cluster <- query[id,1] #read cluster ID
        source <- gsub("'", '', (as.character(query[id,2]))) #read source id
        recipient <- gsub("'", '', (as.character(query[id,3]))) #read recipient id
        print(paste("Queuing the Los Alamos db for cluster ",cluster,", source patient ",source," & recipient patient ",recipient,sep=""))
        system(paste("load_hiv cluster ", cluster, sep="")) #launch python "load_hiv"
        info <- read.delim(paste("cluster_",cluster,"_clinical.tsv",sep=""), header = FALSE, stringsAsFactors =FALSE) #load clinical data
        info.gb <- read.delim(paste("cluster_",cluster,"_accessions.tsv",sep=""), header = TRUE, stringsAsFactors =FALSE) #load accession data
        system(paste("rm cluster_", cluster, "_clinical.tsv ", "cluster_",cluster,"_accessions.tsv ",  sep="")) #remove raw output
        pos<- as.data.frame(str_split_fixed(info.gb$pos, ":", 2), stringsAsFactors=FALSE) #split HBX2 coordinates in columns
        pos[] <- lapply(pos, as.numeric)
        names(pos)<- c("pos_ini", "pos_end") 
        info.gb<- cbind(info.gb, pos) 
        info <- as.data.frame(apply(info, 1, unlist), stringsAsFactors = FALSE) #format epidemiological info
        names(info) <- gsub(" ", "", info[1,], fixed = TRUE) 
        info<- info[-1,] 
        rownames(info)=NULL
        names(info) <- c('PatientId',names(info)[2:ncol(info)])
        source.col <- names(info[,1:2])[grepl(source,info[,1:2])] #match source entry with Patientid
        source.row <- as.numeric(rownames(info)[grepl(source,info[,source.col])]) #retrieve source row
        recipient.col <- names(info[,1:2])[grepl(recipient,info[,1:2])] #match recipient entry with Patientid
        recipient.row <- as.numeric(rownames(info)[grepl(recipient,info[,recipient.col])]) #retrieve recipient row
        info$States <- "NA" #format epidemiological info
        info[source.row,'States'] <- "source"
        info[recipient.row,'States'] <- "recipient"
        info$ClusterId <- c(rep(cluster, nrow(info)))
        source <- info[source.row, source.col] 
        recipient <- info[recipient.row, recipient.col]
        info.gb <- info.gb[which(info.gb$patient_id %in% c(source, recipient)),]
        info.gb$sampling.dates <- NA
        accessions <- info.gb$accession_id
        for(entrada in accessions) { #retrieve calendar dates from GenBank accessions 
            try.gba <- try(gba<- readGenBank(GBAccession(entrada), partial=TRUE), silent = T) 
            if(!class(try.gba) == "try-error"){
                gba<- readGenBank(GBAccession(entrada), partial=TRUE)
                if(!is.null(gba@sources@elementMetadata@listData$collection_date)){
                    info.gb$sampling.dates[which(info.gb$accession_id==entrada)] <- gba@sources@elementMetadata@listData$collection_date
                } else info.gb$sampling.dates[which(info.gb$accession_id==entrada)] <- NA
            }
        }
        #Create directory and save data 
        if(nrow(info)>2) {cluster.name<- paste(as.character((unique(info$ClusterName))),".",noquote(info[recipient.row,"PatientCode"]),sep="")} else{cluster.name<-as.character((unique(info$ClusterName)))}
        if(nchar(id)>20){cluster.name <-paste("id", cluster,gsub("^.*?\\.", ".", cluster.name), sep="")}
        ids <<- c(ids, cluster.name)
        write.table(ids, 'ids.txt', quote=F, row.names = F, col.names=F)
        system(paste("mkdir ", cluster.name,sep="")) #create dir to save data. Takes into account if pair or transmission cluster.
        write.csv(info, paste(cluster.name, "/info.csv", sep=""), quote=T, row.names = F)
        write.csv(info.gb, paste(cluster.name, "/info.gb.csv", sep=""), quote=T, row.names = F)
    }
}
