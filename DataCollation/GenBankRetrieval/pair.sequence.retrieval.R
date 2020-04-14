pair.sequence.retrieval <- function(ids, min.overlap, max.length.diff, min.seqs){
    require(DescTools) #Overlap()
    require(ape) #read.GenBank
    require(stringr) #str_split
    '%!in%' <- function(x, y) ! ('%in%'(x, y)) #define operator
    if(missing(min.overlap)) min.overlap = 70 
    if(missing(max.length.diff)) max.length.diff = 10 
    if(missing(min.seqs)) min.seqs = 5 
    print(paste(length(ids), "Transmission cluster(s) to process"))
    for(id in ids){     #loop every transmission pair 
        cluster.name <- id
        print(paste("Processing ", cluster.name, sep="")) #load pair data
        info.gb <- read.csv(paste(cluster.name, "/", "info.gb.csv", sep=""), stringsAsFactors = F) #read data
        info.gb<- info.gb[,c("row_id", "patient_id", "accession_id", 
                                        "seq_length", "pos","pos_ini","pos_end", "subtype")]
        info <- read.csv(paste(cluster.name, "/", "info.csv", sep=""),stringsAsFactors = F)
        cluster <- unique(info$ClusterId) 
        source <- na.omit(info$PatientId[info$States=='source']) 
        recipient <- na.omit(info$PatientId[info$States=='recipient']) 
        #create data.frame to save sets
        set.id <- data.frame( "cluster_name"=character(), "cluster_id"=integer(), "donor.id"=integer(), "recipient.id"=integer(),
                              "subtype"=character(), "ref" = character(), "counts" = integer(), "source.n"=integer(), "recipient.n"=integer(), "min" = integer(), "max" = integer(),
                              "length" = integer(), "accessions" = character(), stringsAsFactors = FALSE) 
        while(nrow(info.gb) >= min.seqs){ #set genomic region 
            newSet<- character() #vector to save hits that will be included
            toRemove <- character() #vector to save hits that will be dropped
            accession.min <- numeric() 
            accession.max <- numeric() 
            minSeq_length <- as.integer(min(info.gb$seq_length)) #identify shortest sequence
            ref<- as.character(info.gb$accession_id[which(info.gb$seq_length == minSeq_length)])[1] 
            ref.ini <-info.gb$pos_ini[which(info.gb$accession_id==ref)] 
            ref.end <- info.gb$pos_end[which(info.gb$accession_id==ref)] 
            for(accession in as.character(info.gb$accession_id)){ #evaluate overlaps
                accession.ini <- info.gb$pos_ini[which(info.gb$accession_id==accession)] 
                accession.end<- info.gb$pos_end[which(info.gb$accession_id==accession)] 
                if(round(Overlap(accession.ini:accession.end, ref.ini:ref.end)/length(ref.ini:ref.end)*100, 0) >=min.overlap){ 
                    newSet <- c(newSet, accession)
                } 
                #drop single.region accessions
                if(round(Overlap(accession.ini:accession.end, ref.ini:ref.end)/length(ref.ini:ref.end)*100, 0) >=min.overlap & 
                   length(accession.ini:accession.end) - length(ref.ini:ref.end) <= round(max.length.diff/100 * length(ref.ini:ref.end), 0)){
                    toRemove <- c(toRemove, accession)
                    accession.min <- c(accession.min, info.gb$pos_ini[which(info.gb$accession_id==accession)])
                    accession.max <- c(accession.max, info.gb$pos_end[which(info.gb$accession_id==accession)])
                }
            } 
            #validate the set (Minimum number of sequences criteria) and return data
            states <- info.gb$patient_id[which(info.gb$accession_id %in% newSet)]
            if(dim(table(states)) == 2){ 
                if(as.numeric(table(states)[1]) >= min.seqs & as.numeric(table(states)[2]) >= min.seqs){
                    set.subtype <- as.character(unique(info.gb$subtype))
                    if(length(set.subtype)>1){
                        if(any(unlist(lapply(lapply(set.subtype, function(x) unlist(strsplit(x, ""))), function(x) length(x))) == 1)){
                            set.subtype <- set.subtype[(unlist(lapply(lapply(set.subtype, function(x) unlist(strsplit(x, ""))), function(x) length(x))) == 1)]
                            if(length(set.subtype)>1) set.subtype <- sort(set.subtype)[1]
                        } else set.subtype <- sort(set.subtype)[1]
                    }
                    set.id.row <- data.frame(
                        "cluster_name"=cluster.name, "cluster_id"=cluster, "donor.id"=source, "recipient.id"=recipient,
                        "subtype"=set.subtype,
                        "ref" = gsub(':','.', as.character(info.gb$pos[which(info.gb$accession_id==ref)])),
                        "counts" = length(newSet), "source.n"=sum(states==source), "recipient.n"=sum(states==recipient),
                        "min" = min(accession.min), "max" = max(accession.max),
                        "length"= max(accession.max) - min(accession.min),
                        "accessions"= paste(newSet, collapse=","), 
                        "accessions.source"= paste(info.gb$accession_id[which(info.gb$patient_id %in% source)], collapse=","), 
                        "accessions.rec"= paste(info.gb$accession_id[which(info.gb$patient_id %in% source)], collapse=","),
                        stringsAsFactors = FALSE)
                    set.id.row[] <- lapply(set.id.row, as.character)
                    set.id <- rbind(set.id.row,set.id)
                }
            }
            info.gb <- info.gb[!info.gb$accession_id %in% toRemove,] #update the pool of sequences to be evaluated
        }
        set.id <- set.id[!duplicated(set.id)] #format results
        set.id$counts <- as.numeric(set.id$counts)
        set.id$min<- as.numeric(set.id$min)
        set.id$max<- as.numeric(set.id$max)
        set.id$length<- as.numeric(set.id$length)
        ##add epidemiological info to set
        set.id[,c('source.dfi', 'source.dfs', 'source.y', 'source.fiebig', 'recipient.dfi', 'recipient.dfs', 'recipient.y', 'recipient.fiebig')] <- list(NA)
        for(set in 1:nrow(set.id)){
            info.gb <- read.csv(paste(cluster.name, "/", "info.gb.csv", sep=""), stringsAsFactors = F)
            accessions <- unlist(strsplit(set.id$accessions[set], ",")) #retrieve accessions for set
            accessions.source <- accessions[which(accessions %in% as.character(info.gb$accession_id[which(info.gb$patient_id ==source)]))]
            accessions.recipient <- accessions[which(accessions %in% as.character(info.gb$accession_id[which(info.gb$patient_id==recipient)]))]
            info.gb.source <- info.gb[which(info.gb$accession_id %in% accessions.source),]
            info.gb.recipient <- info.gb[which(info.gb$accession_id %in% accessions.recipient),]
            set.id$source.y[[set]] <- trimws(paste((unique(info.gb.source$sampling_year)), collapse=" "))
            source.fiebig<- character()
            for (i in 1:length(unlist(unique(info.gb.recipient$fiebig_stage)))){
                source.fiebig <- c(source.fiebig,tail(unlist(str_split(unlist(unique(info.gb.recipient$fiebig_stage))[i], " ")),1))
            }
            set.id$source.fiebig[[set]] <- trimws(paste(source.fiebig, collapse=" "))
            set.id$source.dfi[[set]] <- trimws(paste((unique(info.gb.source$days_from_infection)), collapse=" "))
            set.id$source.dfs[[set]] <- trimws(paste((unique(info.gb.source$days_from_seroconversion)), collapse=" "))
            set.id$recipient.y[[set]] <- trimws(paste((unique(info.gb.recipient$sampling_year)), collapse=" "))
            recipient.fiebig<- character()
            for (i in 1:length(unlist(unique(info.gb.recipient$fiebig_stage)))){
                recipient.fiebig <- c(recipient.fiebig,tail(unlist(str_split(unlist(unique(info.gb.recipient$fiebig_stage))[i], " ")),1))
            }
            set.id$recipient.fiebig[[set]] <- trimws(paste(recipient.fiebig, collapse=" "))
            set.id$recipient.dfi[[set]] <- trimws(paste((unique(info.gb.recipient$days_from_infection)), collapse=" "))
            set.id$recipient.dfs[[set]] <- trimws(paste((unique(info.gb.recipient$days_from_seroconversion)), collapse=" "))
        }
        write.csv(set.id, paste(id, "/sets.csv", sep=""), row.names = FALSE) #save sets
        #retrieve fasta sequences
        source.gb <- read.GenBank(as.character(info.gb$accession_id[which(info.gb$patient_id == source)]))
        names(source.gb) <- paste('source.', names(source.gb),sep="")
        recipient.gb <- read.GenBank(as.character(info.gb$accession_id[which(info.gb$patient_id == recipient)]))
        names(recipient.gb) <- paste('recipient.', names(recipient.gb),sep="")
        #retrieve subtype transmission pair specific fasta sequences.
        df.subtype <- data.frame("subtype"=c("A","A","A","A1","A1","A1","A2","A2","A2","B","B","B","B","C","C","C","C","D","D","D","D","F1","F1","F1","F1","F2","F2","F2","F2","G","G","G","G","H","H","H","H","J","J","J","K","K","01_AE","01_AE","01_AE","02_AG","02_AG","02_AG","03_AB","04_cpx","04_cpx","04_cpx","05_DF","05_DF","05_DF","06_cpx","06_cpx","06_cpx","07_BC","07_BC","07_BC","08_BC","08_BC","09_cpx","09_cpx","09_cpx","09_cpx","10_CD","10_CD","10_CD","11_cpx","11_cpx","11_cpx","12_BF","12_BF","12_BF","13_cpx","13_cpx","13_cpx","14_BG","14_BG","14_BG","15_01B","15_01B","15_01B","16_A2D","16_A2D","17_BF","17_BF","17_BF","18_cpx","18_cpx","18_cpx","19_cpx","19_cpx","19_cpx","20_BG","21_A2D","21_A2D","21_A2D","22_01A1","22_01A1","23_BG","23_BG","24_BG","24_BG","24_BG","25_cpx","25_cpx","25_cpx","26_AU","26_AU","26_AU","27_cpx","27_cpx","28_BF","28_BF","28_BF","29_BF","29_BF","29_BF","31_BC","31_BC","31_BC","32_06A1","33_01B","33_01B","33_01B","34_01B","35_AD","35_AD","35_AD","36_cpx","36_cpx","37_cpx","37_cpx","38_BF1","38_BF1","38_BF1","39_BF","39_BF","39_BF","40_BF","40_BF","40_BF","42_BF","43_02G","43_02G","43_02G","44_BF","45_cpx","45_cpx","45_cpx","46_BF","46_BF","46_BF","47_BF","47_BF","49_cpx","49_cpx","49_cpx"),
                                 "accession"=c("DQ676872","AB253421","AB253429","DQ676872","AB253421","AB253429","AF286238","GU201516","AF286237","K03455","AY423387","AY173951","AY331295","U52953","U46016","AF067155","AY772699","K03454","AY371157","AY253311","U88824","AF077336","AF005494","AF075703","AJ249238","AY371158","AJ249236","AJ249237","AF377956","AF084936","AF061641","U88826","AY612637","AF190127","AF190128","AF005496","FJ711703","EF614151","GU237072","AF082394","AJ249235","AJ249239","GQ477441","GU564221","U54771","AY271690","AB485636","L39106","AF193276","AF049337","AF119820","AF119819","AF076998","AF193253","AY227107","AF064699","AY535659","AB286851","EF368372","EF368370","AF286230","HM067748","AY008715","AJ866553","AY093605","AY093603","AY093607","AF289548","AF289549","AF289550","AF492624","AF492623","AJ291718","AF408629","AF408630","AF385936","DQ845388","DQ845387","AF460972","AF450096","AF450097","GU230137","DQ354120","AF516184","AF530576","AY945736","AF286239","EU581825","EU581827","EU581828","AF377959","AY586541","AY894993","AY588971","AY588970","AY894994","AY586545","AY945737","AF457051","AF457072","AY371159","GQ229529","AY900571","AY900572","AY900574","AY900575","FJ670526","EU693240","EU697906","EU697908","FM877780","FM877782","FM877777","AJ404325","AM851091","DQ085872","DQ085873","DQ085874","DQ085876","AY771590","DQ085871","EF091932","AY727526","AY727527","AY535660","AB547464","DQ366659","DQ366662","EF165541","EF158043","EF158040","EF158041","EF087995","EF087994","EF116594","AF377957","FJ213781","FJ213782","FJ213780","EU735534","EU735536","EU735535","EU735538","EU735540","EU735539","EU170155","EU697904","EU697907","EU697909","FJ358521","FN392874","FN392876","FN392877","DQ358801","DQ358802","HM026456","GQ372987","FJ670529","HQ385477","HQ385479","HQ385478"), stringsAsFactors = FALSE)
        subtype <- as.character(unique(info.gb$subtype))
        if(length(subtype)>1) subtype <- sort(subtype)[1]
        if(subtype %!in% sort(unique(df.subtype$subtype))) subtype <- unlist(strsplit(unique(unlist(strsplit(subtype, "[^A-Z]+"))), ""))[1]
        outgroup<- read.GenBank(df.subtype$accession[which(subtype == df.subtype$subtype)])
        ingroup<-c(source.gb, recipient.gb) #save fasta files
        write.dna(ingroup, paste(id, "/ingroup.fasta",sep=""), format = "fasta", colsep = "") ##export in fasta
        write.dna(outgroup, paste(id, "/outgroup.fasta",sep=""), format = "fasta", colsep = "") ##export in fasta
    }
}
