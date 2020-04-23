pair.epi.retrieval <- function(query, key) {
    ####################################
    #load libraries
    ####################################
    require(stringr) #str_split_fixed
    require(genbankr) #readGenBank
    require(rentrez) #set_entrez_key
    ####################################
    #load file, set  values for function arguments
    ####################################
    query <- read.csv(query,
                      header = FALSE,
                      col.names = c("id", "source", "recipient"))
    print(paste(nrow(query), "Transmission cluster(s) to process"))
    set_entrez_key(key) #set API key
    print(paste("Using entrez key ", key, sep = ""))
    ####################################
    #Query clusters to LANLdb
    ####################################
    ids <- character()
    for (id in 1:nrow(query)) {
        ####################################
        #get source/recipient ids from file
        ####################################
        cluster <- query[id, 1]
        source <- gsub("'", '', as.character(query[id, 2]))
        recipient <- gsub("'", '', as.character(query[id, 3]))
        print(
            paste(
                "Queuing the Los Alamos db for cluster ",
                cluster,
                ", source patient ",
                source,
                " & recipient patient ",
                recipient,
                sep = ""
            )
        )
        ####################################
        #launch python script
        ####################################
        system(paste("load_hiv cluster ", cluster, sep = ""))
        ####################################
        #load output of python script
        ####################################
        info <- read.delim(
            paste("cluster_",
                  cluster, "_clinical.tsv", sep = ""),
            header = FALSE,
            stringsAsFactors = FALSE
        )
        info.gb <- read.delim(
            paste("cluster_",
                  cluster, "_accessions.tsv", sep = ""),
            header = TRUE,
            stringsAsFactors = FALSE
        )
        ####################################
        #polish outputs
        ####################################
        pos <- as.data.frame(str_split_fixed(info.gb$pos, ":", 2),
                             stringsAsFactors = FALSE)
        pos[] <- lapply(pos, as.numeric)
        names(pos) <- c("pos_ini", "pos_end")
        info.gb <- cbind(info.gb, pos)
        info <-
            as.data.frame(apply(info, 1, unlist), stringsAsFactors = FALSE) #format epidemiological info
        names(info) <- gsub(" ", "", info[1, ], fixed = TRUE)
        info <- info[-1, ]
        rownames(info) <- NULL
        names(info) <- c('PatientId', names(info)[2:ncol(info)])
        ####################################
        #retrieve info from outputs
        ####################################
        source.col <- names(info[, 1:2])[grepl(source, info[, 1:2])]
        source.row <- as.numeric(rownames(info)[grepl(source, info[, source.col])])
        recipient.col <-
            names(info[, 1:2])[grepl(recipient, info[, 1:2])]
        recipient.row <- as.numeric(rownames(info)[grepl(recipient, info[, recipient.col])])
        info$States <- "NA"
        info[source.row, 'States'] <- "source"
        info[recipient.row, 'States'] <- "recipient"
        info$ClusterId <- c(rep(cluster, nrow(info)))
        source <- info[source.row, source.col]
        recipient <- info[recipient.row, recipient.col]
        info.gb <- info.gb[which(info.gb$patient_id
                                 %in% c(source, recipient)), ]
        info.gb$sampling.dates <- NA
        accessions <- info.gb$accession_id
        ####################################
        #Query genbank for additional info
        ####################################
        for (entrada in accessions) {
            try.gba <- try(gba <- readGenBank(GBAccession(entrada),
                                              partial = TRUE),
                           silent = T)
            if (!class(try.gba) == "try-error") {
                gba <- readGenBank(GBAccession(entrada), partial = TRUE)
                if (!is.null(gba@sources@elementMetadata@listData$collection_date)) {
                    info.gb$sampling.dates[which(info.gb$accession_id == entrada)] <-
                        gba@sources@elementMetadata@listData$collection_date
                } else
                    info.gb$sampling.dates[which(info.gb$accession_id == entrada)] <- NA
            }
        }
        ####################################
        #save data
        ####################################
        if (nrow(info) > 2) {
            cluster.name <-
                paste(as.character((unique(
                    info$ClusterName
                ))), ".", noquote(info[recipient.row, "PatientCode"]), sep = "")
        } else{
            cluster.name <- as.character((unique(info$ClusterName)))
        }
        if (nchar(id) > 20)
            cluster.name <-
            paste("id", cluster, gsub("^.*?\\.", ".", cluster.name), sep = "")
        write.table(
            ids,
            'ids.txt',
            quote = F,
            row.names = F,
            col.names = F
        )
        system(paste("mkdir ", cluster.name, sep = ""))
        write.csv(
            info,
            paste(cluster.name, "/info.csv", sep = ""),
            quote = T,
            row.names = F
        )
        write.csv(
            info.gb,
            paste(cluster.name, "/info.gb.csv", sep = ""),
            quote = T,
            row.names = F
        )
        system(
            paste(
                "rm cluster_",
                cluster,
                "_clinical.tsv ",
                "cluster_",
                cluster,
                "_accessions.tsv ",
                sep = ""
            )
        )
    }
    return(ids)
}
