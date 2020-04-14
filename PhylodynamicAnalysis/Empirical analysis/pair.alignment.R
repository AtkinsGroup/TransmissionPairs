pair.alignment <- function(ids, max.n, ref.gap.tolerance, gap.threshold, sets){
    require(ape) #as.matrix.DNAbin read.GenBank
    require(Biostrings) #readDNAStringSet
    require(R.utils) #seqToIntervals
    require(pegas) #nuc.div
    if(missing(max.n)) max.n = 300 #default value 
    if(missing(ref.gap.tolerance)) ref.gap.tolerance = 99 #default value 
    if(missing(gap.threshold)) gap.threshold = 1 #default value 
    if(missing(sets)) sets = "set.csv" #default value
    HBX2<- as.matrix.DNAbin(read.GenBank('K03455'))
    for(id in ids){
        cluster.name <- id
        print(paste("Processing ", cluster.name, sep=""))
        #read data
        info.gb <- read.csv(paste(cluster.name, "/", "info.gb.csv", sep=""), stringsAsFactors = F)
        info.gb<- info.gb[,c("row_id", "patient_id", "accession_id", #df to save the sets
                             "seq_length", "pos","pos_ini","pos_end")]
        info <- read.csv(paste(cluster.name, "/", "info.csv", sep=""),stringsAsFactors = F)
        cluster <- unique(info$ClusterId) #get cluster ID
        source <- na.omit(info$PatientId[info$States=='source']) #get source
        recipient <- na.omit(info$PatientId[info$States=='recipient']) #get recipient
        #load sets
        set.id<- read.csv(paste(id, "/", sets, sep=""), row.names = NULL, stringsAsFactors = F) #load sets.csv
        #Add colums to sets
        for(set in 1:nrow(set.id)){ #load every sets
            if(set.id$counts[set] > 300){
                set.id$alignment.l[set] <- NA
                set.id$gap.percentage[set] <- NA
                next
            }
            system(paste("mkdir ", id, "/", id, ".",set.id$ref[set], sep="")) #create dir to save data
            accessions <- unlist(strsplit(set.id$accessions[set], ",")) #retrieve accessions for set
            accessions.source <- accessions[which(accessions %in% as.character(info.gb$accession_id[which(info.gb$patient_id ==source)]))]
            accessions.recipient <- accessions[which(accessions %in% as.character(info.gb$accession_id[which(info.gb$patient_id ==recipient)]))]
            accessions<- c(accessions.source,accessions.recipient)
            accessions.gb <- read.GenBank(accessions)
            accessions.source <- paste("source.", accessions.source, sep="")
            accessions.recipient <- paste("recipient.", accessions.recipient, sep="")
            accessions<- c(accessions.source,accessions.recipient)
            names(accessions.gb) <- accessions
            #alignment: align HXB2 region with genomic subtype refs
            outgroup.HBX2 <-  HBX2[labels(HBX2), set.id$min[set]:set.id$max[set]] #HXB2 fragment
            write.dna(outgroup.HBX2, paste(id, "/", id, ".",set.id$ref[set], "/outgroup.HBX2",sep=""), format = "fasta", colsep = "") ##export HXB2 fragment in fasta
            outgroup <- read.FASTA(paste(id, "/outgroup.fasta",sep=""))
            write.dna(outgroup, paste(id, "/", id, ".",set.id$ref[set], "/outgroup.fasta",sep=""), format = "fasta", colsep = "") ##export subtype specific ref genome sequences in fasta
            system(paste("cat ", id, "/", id, ".",set.id$ref[set], "/outgroup.HBX2 ",id, "/", id, ".",set.id$ref[set], "/outgroup.fasta > ",id, "/", id, ".",set.id$ref[set], "/output", sep="")) 
            system(paste("muscle -in ",id, "/", id, ".",set.id$ref[set], "/output -out ",id, "/", id, ".",set.id$ref[set], "/outgroup.HBX2", sep="")) #aligning
            #alignment: load alignment and transform to matrix
            outgroup.gb <- read.dna(paste(id, "/", id, ".",set.id$ref[set], "/outgroup.HBX2", sep=""),format="fasta", as.matrix=TRUE)
            if(set.id$subtype[set] == "B") outgroup.gb <- outgroup.gb["K03455",] else outgroup.gb <- outgroup.gb[labels(outgroup),] #outgroup alignment
            coordinates <- readDNAStringSet(paste(id, "/", id, ".",set.id$ref[set], "/outgroup.HBX2", sep=""))
            coordinates <- data.frame(names(coordinates), paste(coordinates))
            coordinates<- coordinates[match(c("K03455", labels(outgroup)),coordinates$names.coordinates.),] #reroder
            coordinates <- strsplit(as.character(coordinates$paste.coordinates.), "")
            #alignment: use HBX2 ref. Transform alignment to df or non-gap intervals
            coordinates<- (coordinates[[1]] != "-")*(coordinates[[1]] != "-")
            coordinates.limits <- which(coordinates==1)
            coordinates.limits<- as.data.frame(seqToIntervals(coordinates.limits))
            #alignment: select parts of the ref alignment to be used
            if(nrow(coordinates.limits)>1){ 
                coordinates.limits$range <- coordinates.limits$to - coordinates.limits$from
                coordinates.limits$gap <- NA
                for(i in 1:(nrow(coordinates.limits)-1)) coordinates.limits[i,4]=coordinates.limits[i+1,1]-coordinates.limits[i,2] 
                coordinates.limits$gap[nrow(coordinates.limits)] <- ref.gap.tolerance+1
                #get info for longer interval
                bottom<- coordinates.limits[which.max(coordinates.limits$range),1]
                top<- coordinates.limits[which.max(coordinates.limits$range),2]
                #use HBX2 coordinates for those intervals where: gaps < ref.gap.tolerance
                if(which.max(coordinates.limits$range) != nrow(coordinates.limits)){
                    for(i in (which.max(coordinates.limits$range)):nrow(coordinates.limits)){
                        if(coordinates.limits[i,4] <= ref.gap.tolerance) top = coordinates.limits[i+1,2] else break
                    }     
                }
                if(which.max(coordinates.limits$range) > 1){
                    for(i in (which.max(coordinates.limits$range)-1):1){
                        if(coordinates.limits[i,4] <= ref.gap.tolerance) bottom = coordinates.limits[i,1] else break
                    }
                }
                coordinates <- c(bottom, top)
            } else coordinates<- c(min(which(coordinates == "1")),max(which(coordinates == "1")))
            #alignment: save regions of references to be used
            if(set.id$subtype[set] == "B") outgroup.gb <- outgroup.gb["K03455", coordinates[1]:coordinates[2]] else outgroup.gb <- outgroup.gb[labels(outgroup), coordinates[1]:coordinates[2]] #outgroup alignment
            outgroup.gb <- del.gaps(outgroup.gb)
            write.dna(accessions.gb, paste(id, "/", id, ".",set.id$ref[set], "/ingroup.fasta",sep=""), format = "fasta", colsep = "") ##export in fasta
            write.dna(outgroup.gb, paste(id, "/", id, ".",set.id$ref[set], "/outgroup.fasta",sep=""), format = "fasta", colsep = "") ##export in fasta
            #alignment: include ingroup and outgroup
            system(paste("cat ", id, "/", id, ".",set.id$ref[set], "/outgroup.fasta ",id, "/", id, ".",set.id$ref[set], "/ingroup.fasta > ",id, "/", id, ".",set.id$ref[set], "/output", sep="")) 
            system(paste("muscle -in ",id, "/", id, ".",set.id$ref[set], "/output -out ",id, "/", id, ".",set.id$ref[set], "/output.fasta", sep="")) #alingning with ref dataset
            dat.dna <- read.dna(paste(id, "/", id, ".",set.id$ref[set], "/output.fasta", sep=""), format="fasta", as.matrix=TRUE) 
            if(set.id$subtype[set] == "B") dat<- dat.dna[c('K03455',accessions),] else dat <- dat.dna[c(labels(outgroup)[1],accessions),] #subset and order
            #BETA testing: limit the alignment to the coordinates of the ref
            coordinates <- readDNAStringSet(paste(id, "/", id, ".",set.id$ref[set], "/output.fasta", sep=""))
            coordinates <- data.frame(names(coordinates), paste(coordinates))
            if(set.id$subtype[set] == "B")  coordinates <- coordinates[coordinates$names.coordinates.=="K03455",]  else coordinates <- coordinates[coordinates$names.coordinates.==labels(outgroup)[1],]
            coordinates <- strsplit(as.character(coordinates$paste.coordinates.), "")
            coordinates<- (coordinates[[1]] != "-")*(coordinates[[1]] != "-")
            coordinates.limits <- which(coordinates==1)
            coordinates.limits<- as.data.frame(seqToIntervals(coordinates.limits))
            #alignment: select parts of the ref alignment to be used
            if(nrow(coordinates.limits)>1){ 
                coordinates.limits$range <- coordinates.limits$to - coordinates.limits$from
                coordinates.limits$gap <- NA
                for(i in 1:(nrow(coordinates.limits)-1)) coordinates.limits[i,4]=coordinates.limits[i+1,1]-coordinates.limits[i,2] 
                coordinates.limits$gap[nrow(coordinates.limits)] <- ref.gap.tolerance+1
                #get info for longer interval
                bottom<- coordinates.limits[which.max(coordinates.limits$range),1]
                top<- coordinates.limits[which.max(coordinates.limits$range),2]
                #use HBX2 coordinates for those intervals where: gaps < ref.gap.tolerance
                if(which.max(coordinates.limits$range) != nrow(coordinates.limits)){
                    for(i in (which.max(coordinates.limits$range)):nrow(coordinates.limits)){
                        if(coordinates.limits[i,4] <= ref.gap.tolerance) top = coordinates.limits[i+1,2] else break
                    }     
                }
                if(which.max(coordinates.limits$range) > 1){
                    for(i in (which.max(coordinates.limits$range)-1):1){
                        if(coordinates.limits[i,4] <= ref.gap.tolerance) bottom = coordinates.limits[i,1] else break
                    }
                }
                coordinates <- c(bottom, top)
            } else coordinates<- c(min(which(coordinates == "1")),max(which(coordinates == "1")))
            dat <- dat[,coordinates[1]:coordinates[2]]
            #remove columns with gaps
            ##Import/read/open a fasta file to a df ##Biostrings.
            dat <- del.colgapsonly(dat, threshold = gap.threshold) #rm gaps cols
            #alignment:save alignment and info
            set.id$alignment.l[set] <- dim(dat)[2]
            set.id$gap.percentage[set] <- round((dim(dat)[2]-(set.id$max[set]-set.id$min[set]))*100/(set.id$max[set]-set.id$min[set]),2)
            write.FASTA(as.list(dat), paste(id, "/", id, ".",set.id$ref[set], "/", id, ".",set.id$ref[set], ".fasta", sep=""))
            write.nexus.data(as.list(dat), paste(id, "/", id, ".",set.id$ref[set], "/", id, ".",set.id$ref[set], ".nex", sep=""), interleaved = FALSE, format = "dna", charsperline=0)
            system(paste("rm ", id, "/", id, ".",set.id$ref[set], "/outgroup.fasta ", id, "/", id, ".",set.id$ref[set], "/output", sep=""))
            #create df with individual states (0, source, 1, recipient)
            df.source<- data.frame('accession'=accessions.source, 'state'=rep("source", length(accessions.source)))
            df.recipient<- data.frame('accession'=accessions.recipient, 'state'=rep("rec", length(accessions.recipient)))
            clades<- rbind(df.source, df.recipient)
            write.table(clades, paste(id, "/", id, ".",set.id$ref[set], "/", id, ".",set.id$ref[set], ".clades",sep=""),col.names = FALSE, row.names = FALSE, quote = FALSE)
            if(set.id$subtype[set] == "B") df.outgroup<- data.frame('id'="K03455", 'state'="?",stringsAsFactors = FALSE) else df.outgroup<- data.frame('id'=labels(outgroup)[1], 'state'="?",stringsAsFactors = FALSE)
            df.source<- data.frame('id'=accessions.source, 'state'=rep(0, length(accessions.source)))
            df.recipient<- data.frame('id'=accessions.recipient, 'state'=rep(1, length(accessions.recipient)))
            states <- rbind(df.outgroup, df.source, df.recipient)
            write.table(states, paste(id, "/", id, ".",set.id$ref[set], "/", id, ".",set.id$ref[set], ".states",sep=""),col.names = FALSE, row.names = FALSE, quote = FALSE)
            #calculate nei distances
            dat <- read.dna(paste(id, "/", id, ".",set.id$ref[set], "/", id, ".",set.id$ref[set],".fasta",sep=""), format = "fasta", as.matrix = TRUE)
            states <- read.table(paste(id, "/", id, ".",set.id$ref[set], "/", id, ".",set.id$ref[set],".states",sep=""), header=FALSE, col.names=c('accession', 'state'), stringsAsFactors = F)
            dat <- dat[states$accession[states$state!="?"],]
            dat.source <- dat[states$accession[states$state=="0"],]
            dat.rec <- dat[states$accession[states$state=="1"],]
            set.id$pair.nuc.div[set] <- nuc.div(dat)
            set.id$source.nuc.div[set] <- nuc.div(dat.source)
            set.id$rec.nuc.div[set] <- nuc.div(dat.rec)
        }
        write.csv(set.id, paste(id, "/", sets, sep=""), row.names = FALSE)
    }
}
