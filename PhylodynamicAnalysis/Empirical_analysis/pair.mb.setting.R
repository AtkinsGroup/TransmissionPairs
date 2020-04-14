pair.mb.setting <- function(ids, ngen, samplefreq, printfreq, diagnfreq, burninfrac, sets){
        require(ape) #as.matrix.DNAbin read.GenBank
        if(missing(ngen)) ngen = as.integer(10000000) #default value for 
        if(missing(samplefreq)) samplefreq = as.integer(ngen*0.0001) #default value for 
        if(missing(printfreq)) printfreq = as.integer(ngen*0.01) #default value for 
        if(missing(diagnfreq)) diagnfreq = as.integer(ngen*0.1) #default value for 
        if(missing(burninfrac)) burninfrac = 0.5 #default value for 
        if(missing(sets)) sets = "set.csv" #default value
        ngen = as.integer(ngen)
        samplefreq = as.integer(samplefreq)
        printfreq = as.integer(printfreq)
        diagnfreq = as.integer(diagnfreq)
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
                                set.id$MCMC.ngen[set] <- NA
                                next
                        }
                     dat <- read.dna(paste(id, "/", id, ".",set.id$ref[set], "/", id, ".",set.id$ref[set], ".fasta", sep=""),format="fasta", as.matrix=TRUE)
                     #save/update info to sets
                     set.id$alignment.l[set] <- dim(dat)[2]
                     set.id$gap.percentage[set] <- round((dim(dat)[2]-(set.id$max[set]-set.id$min[set]))*100/(set.id$max[set]-set.id$min[set]),2)
                     set.id$MCMC.ngen[set] <- as.integer(ngen)
                     write.nexus.data(as.list(dat), paste(id, "/", id, ".",set.id$ref[set], "/", id, ".",set.id$ref[set], ".nex", sep=""), interleaved = FALSE, format = "dna", charsperline=0)
                     #create files for mrbayes with ancestral state reconstruction
                     fileConn<- file(paste(id, "/", id, ".",set.id$ref[set], "/", id, ".",set.id$ref[set], ".mb.nex", sep=""))
                     writeLines(c(
                             "#NEXUS",
                             "BEGIN DATA;",
                             paste("DIMENSIONS NTAX=", dim(dat)[1]," NCHAR=", dim(dat)[2]+1, ";",sep=""),
                             paste("FORMAT DATATYPE=mixed(DNA:1-", dim(dat)[2],",Standard:", dim(dat)[2]+1,") MISSING=? GAP=- INTERLEAVE=YES;", sep=""),
                             "MATRIX"
                     ),fileConn)
                     close(fileConn)
                     system(paste("awk '/MATRIX/{flag=1;next}/;/{flag=0}flag' ", id, "/", id, ".",set.id$ref[set], "/", id, ".",set.id$ref[set], ".nex >> ", id, "/", id, ".",set.id$ref[set], "/", id, ".",set.id$ref[set], ".mb.nex", sep=""))
                     system(paste("echo >> ", id, "/", id, ".",set.id$ref[set], "/", id, ".",set.id$ref[set], ".mb.nex; cat ",  id, "/", id, ".",set.id$ref[set], "/", id, ".",set.id$ref[set], ".states >> ", id, "/", id, ".",set.id$ref[set], "/", id, ".",set.id$ref[set], ".mb.nex ", sep=""))
                     write("\n;\nend;",file=paste(id, "/", id, ".",set.id$ref[set], "/", id, ".",set.id$ref[set], ".mb.nex", sep=""),append=TRUE)
                     fileConn<-file(paste(id, "/", id, ".",set.id$ref[set], "/mrbayes.nex", sep=""))
                     #writeLines(c("begin mrbayes;","set autoclose=yes nowarn=yes;", paste("execute ", id, ".",set.id$ref[set], ".mb.nex;", sep=""), #server version
                     writeLines(c("begin mrbayes;","set autoclose=yes nowarn=yes;", paste("execute ", id, "/", id, ".",set.id$ref[set], "/", id, ".",set.id$ref[set], ".mb.nex;", sep=""),
                                     "outgroup 1;", #to better edit the outgroup choice
                                     paste("charset DNA = 1-",dim(dat)[2], ";", sep=""),
                                     paste("charset Transmission = ",dim(dat)[2]+1,";",sep=""),
                                     "partition ancstates = 2: DNA, Transmission;",
                                     "set partition = ancstates;",
                                     "lset applyto=(1) nst=6 rates=invgamma;",
                                     "unlink statefreq = (all) revmat = (all) pinvar = (all) shape = (all);",
                                     "prset ratepr=variable;",
                                     paste("constraint ingroup = 2 - ", dim(dat)[1], ";", sep=""),
                                     "prset topologypr = constraints(ingroup);",
                                     "report applyto=(2) ancstates=yes;",
                                     paste("mcmc nchains=2 ngen=",as.integer(ngen)," samplefreq=",samplefreq," printfreq=",printfreq," diagnfreq=",diagnfreq," burninfrac=",burninfrac,";", sep=""),
                                     "sump;",
                                     "end;"), fileConn)
                        close(fileConn)
                }
                write.csv(set.id, paste(id, "/", sets, sep=""), row.names = FALSE)
        }
}
