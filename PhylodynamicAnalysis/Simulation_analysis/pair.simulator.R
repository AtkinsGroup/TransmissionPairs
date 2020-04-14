pair.simulator <- function(pair.id, vts, simulation, transmission, sampling.source, sampling.rec, index, n.source, n.rec, vt.source, vt.rec, tau, Ne0, K, t50, evo.r, gene.ref, subs.model, subs.param, HBX2, acute, recent, chronic) {
    pot.sampler <- function(n, pots) 1+table(factor(sample(pots, n-pots, rep = T), levels = 1:pots))
    require(ape) #as.matrix.DNAbin read.GenBank
    require(phylotools)# to load phylip data 
    require(seqRFLP) #df to fasta
    require(adegenet) #fasta2DNAbin
    require(pegas) #nuc.div
    #set default values
    if(missing(vts)) vts <- "/Users/julian/Documents/Tools/VirusTreeSimulator-master/out/artifacts/VirusTreeSimulator_jar/VirusTreeSimulator.jar"
    if(missing(pair.id)) pair.id <- 'unknown'
    if(missing(simulation)) simulation <- 1
    if(missing(transmission)) transmission <- 30  #default value for transmission 
    if(missing(sampling.source)) sampling.source <- 30  #default value for sampling of the source
    if(missing(sampling.rec)) sampling.rec <- 30  #default value for sampling of the recipient
    if(missing(acute)) acute<- c(13:90)
    if(missing(recent)) recent<- c(91:180)
    if(missing(chronic)) chronic<- c(181:7300)
    if(!is.numeric(transmission)){
      if(length(unlist(strsplit(transmission, " ")))==1){
        if(transmission=='chronic') transmission<- sample(chronic,1) else {
          if(transmission=='acute') transmission <- sample(acute,1) else {
            if(transmission=='recent') transmission <- sample(recent,1) else {
            transmission <- unique(na.omit(as.numeric(unlist(strsplit(transmission, "[^0-9]+")))))
            if(transmission==max(chronic)) transmission <- max(chronic) else transmission <- sample(transmission:max(chronic),1)
            }
          }
        }
      }
      else transmission <- sample(as.numeric(unlist(strsplit(transmission, " ")))[1]:as.numeric(unlist(strsplit(transmission, " ")))[2],1)
    }
    if(transmission < 13) transmission <- 13  #adjuts for  early estimations
    if(!is.numeric(sampling.source)){
      sampling.source<- sample(as.numeric(unlist(strsplit(sampling.source, " ")))[1]:as.numeric(unlist(strsplit(sampling.source, " ")))[2],1)
    }
    if(!is.numeric(sampling.rec)){
      sampling.rec<- sample(as.numeric(unlist(strsplit(sampling.rec, " ")))[1]:as.numeric(unlist(strsplit(sampling.rec, " ")))[2],1)
    }  
    if(sampling.rec < 13) sampling.rec <- 13  #adjuts for  early estimations
    if(sampling.source<0 & abs(sampling.source) > transmission - 13) sampling.source <- (transmission-13)*-1 #adjuts for  early estimations when sampling source before transmission
    transmission <- transmission/365
    sampling.source <- sampling.source/365
    sampling.rec <- sampling.rec/365
    if(transmission<91/365) transmission.char = "acute" else transmission.char = "chronic"  #default value for transmission 
    if(sampling.source<91/365) sampling.source.char="early" else sampling.source.char="late"
    if(sampling.rec<91/365) sampling.rec.char="early" else sampling.rec.char="late"
    if(missing(n.source)) n.source <- 5  #default value for sampled sequences
    if(missing(n.rec)) n.rec <- 5  #default value for sampled sequences
    if(missing(vt.source)) vt.source <- 1 #default value for particles at transmission
    if(missing(index)) index <- 365*3
    index <- index/365
    if(missing(vt.rec)) vt.rec <- 1  #default value for F/T at transmission
    if(missing(tau)) tau <- 1.8 #default value for log. dem. pop. generation time parameter
    tau <- tau/365
    if(missing(Ne0)) Ne0 <- 1 #default value for log. dem. pop. Ne0 at time of Infection/Transmission
    if(missing(K)) K <- 300 #default value for log. dem. pop. final Ne (carrying capacity)
    if(missing(t50)) t50 <- -2 #default value for log. dem. pop. -t50 parameter (VirusTreeSimulator Input)
    if(missing(evo.r)) evo.r <- 0.01148 #default value for evolutionary rate
    if(missing(gene.ref)) gene.ref <- c(6758, 7757) #default value for gene.length
    if(missing(subs.model)) subs.model <- "GTR" #default value for substitution.model
    if(missing(subs.param)) subs.param <- c(0.48,4,0.18,0.41,0.17,0.20,0.22,3.01,5.59,1.21,1.39,8.25,1.0) #default value Env-GTR. gamma,cat,inv,fa,fc,fg,ft,rates)
    if(missing(HBX2)) HBX2<- as.matrix.DNAbin(read.GenBank('K03455'))
    #calculate log. dem. pop. parameters for VTS:
    C= (K*Ne0)/(K-Ne0) #Katie's formulation to keep K (carrying capacity) constant.
    r=(log(K/C))*(1/-t50) #calculate log. dem. pop. rate (VTS Input)
    N0<- tau * Ne0 #calculate N0. Ne0 at time of Infection/Transmission * generation time (VTS Input)
     #Create Transmission/Sampling files for VTS
    samples.source <- pot.sampler(n.source,vt.source)
    samples.rec <- pot.sampler(n.rec,vt.rec)
    transmitting.pots <- character()
    sink("inf_test.csv")
    cat("IDREC,IDTR,TIME_TR\nindex,NA,0\n")
    for(i in 1:length(samples.source)) cat(paste("source_pot",i,",index,",index,"\n", sep=""))
    for(i in 1:length(samples.rec)){
        j<- sample(1:length(samples.source),1) #transmitting pot
        transmitting.pots <- c(transmitting.pots,j)
        cat(paste("rec_pot",i,",source_pot",j,",",transmission+index,"\n",sep=""))
    }
    sink()
    sink("samp_test.csv")
    cat("IDPOP,TIME_SEQ,SEQ_COUNT\n")
    for(i in 1:length(samples.source)) cat(paste("source_pot",i,",",sampling.source+transmission+index,",",samples.source[[i]],"\n",sep=""))
    for(i in 1:length(samples.rec)) cat(paste("rec_pot",i,",",sampling.rec+transmission+index,",",samples.rec[[i]],"\n",sep=""))
    sink()
    #run VTS and load simulated trees
    stdout<- system(paste("java -jar ",vts," -demoModel Logistic -N0 ",N0," -growthRate ",r," -t50 ",t50," inf_test.csv samp_test.csv simLog", sep=""), intern = TRUE)
    stdout<- gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.", "", paste(stdout,collapse = ""), perl=T)
    tr <- read.nexus("simLogID_index_simple.nex")
    tr.d<- read.nexus("simLogID_index_detailed.nex")
    system("rm simLogID_index_detailed.nex simLogID_index_simple.nex inf_test.csv samp_test.csv")
    tips <- tr$tip.label
    states<-character()
    states.binary<-character()
    for(tip in tips){
        if(any(unlist(strsplit(tip, "_")) %in% "rec")){states<- c(states, "rec")
        states.binary<-c(states.binary,"1")} else{
            if(any(unlist(strsplit(tip, "_")) %in% "source")){states<- c(states, "source")
            states.binary<-c(states.binary,"0")} 
        }
    }
    clades<- data.frame("id"=tips, "state"=states, stringsAsFactors = FALSE)
    states<-data.frame("id"=tips, "state"=states.binary, stringsAsFactors = FALSE)
    t.p.vectors <- topology.class(tr=as.multiPhylo(tr), states=states, outgroup = FALSE, p.ancestral.logs=FALSE, logs = FALSE)
    sim.lineages.source <- t.p.vectors[['vector.lineages.source']] 
    sim.lineages.recipient <- t.p.vectors[['vector.lineages.recipient']]
    sim.topology <- t.p.vectors[['vector.topology']]
    sim.p.ancestral.pars <- as.numeric(t.p.vectors[['vector.p.ancestral.pars']])
    sim.p.ancestral.ml.ER <- as.numeric(t.p.vectors[['vector.p.ancestral.ml.ER']])
    sim.tr=tr
    #set the trees to simulate sequences at transmission:
    singletons <- as.numeric(names(which(unlist(lapply(listDD(tr.d), length) ==1))))
    singletons <- singletons[which(singletons>min(singletons))]
    singleton.count.source <- 0
    singleton.count.rec <- 0
    singleton.names.source <- character()
    singleton.names.rec <- character()
    while (length(singletons) >0) {
        Singleton <- min(singletons)
        numDescendants <- getDescendants(tr.d, Singleton) #list descendants
        chartips <- tr.d$tip.label[c(numDescendants)] #get names of descendants
        chartips <- chartips[!is.na(chartips)] #remove internal nodes, only keep tips
        if (as.numeric(distRoot(tr.d, tips = Singleton))<transmission+index-1/365){
            chartips<- chartips[!chartips %in% grep(paste0("rec", collapse = "|"), chartips, value = T)]
            singleton.count.source <- unique(as.numeric(gsub("[^0-9]", "", sub("^([^_]+_){2}([^_]+).*", "\\2", chartips)))) #retrieve name of infected source pot
            singleton.names.source <- c(singleton.names.source, paste("tv.source.", singleton.count.source, sep=""))
            tr.d <- bind.tip(tr.d, paste("tv.source.", singleton.count.source, sep=""), 0.00001, Singleton)
            singletons <- as.numeric(names(which(unlist(lapply(listDD(tr.d), length) ==1))))
            singletons <- singletons[which(singletons>min(singletons))]
        } else{
            chartips<- chartips[!chartips %in% grep(paste0("source", collapse = "|"), chartips, value = T)]
            singleton.count.rec <- unique(as.numeric(gsub("[^0-9]", "", sub("^([^_]+_){2}([^_]+).*", "\\2", chartips)))) #retrieve name of infected rec pot
            singleton.names.rec <- c(singleton.names.rec, paste("tv.rec.", singleton.count.rec, sep=""))
            tr.d <- bind.tip(tr.d, paste("tv.rec.", singleton.count.rec, sep=""), 0.00001, Singleton)
            singletons <- as.numeric(names(which(unlist(lapply(listDD(tr.d), length) ==1))))
            singletons <- singletons[which(singletons>min(singletons))]
        }
    }
    tr.d<- collapse.singles(tr.d)
    #Set the parameters for SeqGen
    dna <-  HBX2[,gene.ref[1]:gene.ref[2]] #HXB2 fragment
    gene.length <- length(dna)
    dna <- toupper(as.alignment(dna)$seq)
    if(subs.model=='GTR') model <- paste("seq-gen -m",subs.model," -n1 -s",evo.r," -a",subs.param[1]," -g",subs.param[2]," -i",subs.param[3]," -f",noquote(paste(subs.param[4:7],collapse=',')),
                                         " -r",noquote(paste(subs.param[8:13],collapse=','))," -k1 -wa -or < sim.nwk > sim.phy", sep="")
    if(subs.model=='HKY')  model <- paste("seq-gen -m",subs.model," -n1 -s",evo.r," -a",subs.param[1]," -g",subs.param[2]," -i",subs.param[3]," -f",noquote(paste(subs.param[4:7],collapse=',')),
                                          " -t",subs.param[8]," -k1 -wa -or < sim.nwk > sim.phy", sep="")
    #simulate sequences with seq-gen
    fileConn<- file("sim.nwk")
    writeLines(c(paste("1 ", nchar(dna), sep=""),
                 paste("HXB2 ", dna, sep=""),
                 "1"),fileConn)
    close(fileConn)
    write.tree(tr.d, "sim.tre")
    system("cat sim.tre >> sim.nwk && rm sim.tre")
    system(model)
    dat<- sequences.sim <- read.table('sim.phy', skip=1, header = F, col.names = c('id', 'seq'))
    sequences.sim$simulation <- simulation
    #save ALL simulated sequences
    if(!file.exists("sequences.sim.csv")) write.table(sequences.sim, 'sequences.sim.csv', row.names = F, sep = ",") else write.table(append=T, sequences.sim, 'sequences.sim.csv', row.names = F, col.names = F, sep = ",")
    #subset simulated sequences
    dat.tv.source <- dat[which(dat$id %in% singleton.names.source),]
    dat.tv.rec <- dat[which(dat$id %in% singleton.names.rec),]
    dat <- dat[which(dat$id %in% clades$id),]
    #calculate genetic distances
    dataframe2fas(dat.tv.source, file = "sim.fas") ##first column the seq' names, second column DNA sequences.
    dat.tv.source<- read.dna("sim.fas", format = "fasta", as.matrix = TRUE)
    dataframe2fas(dat.tv.rec, file = "sim.fas") ##first column the seq' names, second column DNA sequences.
    dat.tv.rec<- read.dna("sim.fas", format = "fasta", as.matrix = TRUE)
    dataframe2fas(dat, file = "sim.fas") ##first column the seq' names, second column DNA sequences.
    dat<- read.dna("sim.fas", format = "fasta", as.matrix = TRUE)
    dat.source <- dat[clades$id[clades$state=='source'],]
    dat.rec <- dat[clades$id[clades$state=='rec'],]
    tv.source.haplotypes <- dim(haplotype(dat.tv.source))[[1]]
    if(dim(dat.tv.source)[1] > 1){
        tv.source.nuc.div <- mean(dist.dna(dat.tv.source, "raw"))
        tv.source.nuc.div.max <- max(dist.dna(dat.tv.source, "raw"))
        tv.source.nuc.div.min <- min(dist.dna(dat.tv.source, "raw"))
    } else tv.source.nuc.div<- tv.source.nuc.div.max <- tv.source.nuc.div.min <- NA
    tv.rec.haplotypes <- dim(haplotype(dat.tv.rec))[[1]]
    if(dim(dat.tv.rec)[1] > 1){
        tv.rec.nuc.div <- mean(dist.dna(dat.tv.rec, "raw"))
        tv.rec.nuc.div.max <- max(dist.dna(dat.tv.rec, "raw"))
        tv.rec.nuc.div.min <- min(dist.dna(dat.tv.rec, "raw"))
    } else tv.rec.nuc.div <- tv.rec.nuc.div.max <- tv.rec.nuc.div.min <- NA
    source.nuc.div <- mean(dist.dna(dat.source, "raw"))
    source.haplotypes <- dim(haplotype(dat.source))[[1]]
    rec.nuc.div <- mean(dist.dna(dat.rec, "raw"))
    rec.haplotypes <- dim(haplotype(dat.rec))[[1]]
    pair.nuc.div <- mean(dist.dna(dat, "raw"))
    pair.haplotypes <- dim(haplotype(dat))[[1]]
    #ML tree
    dat<- read.table('sim.phy', skip=1, header = F, col.names = c('id', 'seq'))
    dat <- dat[which(dat$id %in% clades$id),]        
    outgroup<- data.frame("id"="HXB2", "seq"=dna)
    dat<- rbind(outgroup, dat)
    dataframe2fas(dat, file = "sim.fas") ##first column the seq' names, second column DNA sequences.
    dat<- fasta2DNAbin("sim.fas")
    #create files with pair information
    tips <- labels(dat)
    states<-character()
    states.binary<-character()
    for(tip in tips){
        if(any(unlist(strsplit(tip, "_")) %in% "rec")){states<- c(states, "rec")
        states.binary<-c(states.binary,"1")} else{
            if(any(unlist(strsplit(tip, "_")) %in% "source")){states<- c(states, "source")
            states.binary<-c(states.binary,"0")} else{
                if(any(unlist(strsplit(tip, "_")) %in% "HXB2")){
                    states.binary<-c(states.binary,"?")}
            }
        }
    }
    clades<- data.frame("id"=tips[-1], "state"=states, stringsAsFactors = FALSE)
    states<-data.frame("id"=tips, "state"=states.binary, stringsAsFactors = FALSE)
    write.table(states, paste("states"),col.names = FALSE, row.names = FALSE, quote = FALSE)
    write.table(clades, paste("clades"),col.names = FALSE, row.names = FALSE, quote = FALSE)
    #run iq-tree and load tree
    system("iqtree -s sim.fas")
    trees<- read.tree("sim.fas.treefile")
    ## set zero-length branches to be 1/1000000 total tree length
    trees$edge.length[trees$edge.length==0]<-max(nodeHeights(trees))*1e-6
    system("rm sim.fas.log sim.fas.model.gz sim.fas.uniqueseq.phy sim.fas.bionj sim.fas.ckp.gz sim.fas.mldist sim.fas.treefile sim.fas.iqtree sim.nwk sim.phy sim.fas clades states")
    t.p.vectors <- topology.class(tr=as.multiPhylo(trees), states=states, outgroup =TRUE, p.ancestral.logs=FALSE, logs=FALSE)
    #save simulated sequences and data
    result<- data.frame("cluster_name"=pair.id, "simulation"=simulation,  "id"=as.numeric(paste(vt.source, vt.rec, sep="")),
                        "transmission"=round(transmission*365,1), "index"=index*365, "sampling.source"=round(sampling.source*365,1), "sampling.rec"=round(sampling.rec*365,1),
                        "transmission.char"=transmission.char, "sampling.source.char"=sampling.source.char, "sampling.rec.char"=sampling.rec.char,
                        "transmission.num"=transmission, "sampling.source.num"=sampling.source, "sampling.rec.num"=sampling.rec,
                        "div.time"=round((sampling.source+sampling.rec)*365,1),
                        "n.source"=n.source, "n.rec"=n.rec, 'vt.source'=vt.source, 'pots.source'=paste(samples.source, collapse=","), 'vt.rec'=vt.rec,
                        'pots.rec'=paste(samples.rec, collapse=","), 'transmitting.pots'=paste(transmitting.pots, collapse=","), 'tau'=tau, 'Ne0'= Ne0, 'K'= K, "N0"=N0, "r"=r, 't50'=t50,
                        'gene.ref'=paste(gene.ref, collapse="."), 'evolutionary.r'=evo.r, 'gene.length'=gene.length,
                        "stdout"=stdout,
                        "sim.topology"=sim.topology, "sim.lineages.source"=sim.lineages.source, "sim.lineages.recipient"=sim.lineages.recipient,
                        "sim.p.ancestral.pars"=sim.p.ancestral.pars, "sim.p.ancestral.ml.ER"=sim.p.ancestral.ml.ER, "sim.TL"=max(distTips(sim.tr)),
                        "sim.SBL"=sum(sim.tr$edge.length), "sim.SPD"=sum(distTips(sim.tr)),
                        "sim.MBL"=mean(sim.tr$edge.length), "sim.MPD"=mean(distTips(sim.tr)),
                        "sim.VBL"=var(sim.tr$edge.length), "sim.VPD"=var(distTips(sim.tr)), 
                        "l"=dim(dat)[2], 'tv.source.haplotypes'=tv.source.haplotypes,
                        'tv.source.nuc.div'=tv.source.nuc.div, 'tv.source.nuc.div.min'=tv.source.nuc.div.max, 'tv.source.nuc.div.max'=tv.source.nuc.div.max, 
                        'tv.rec.haplotypes'=tv.rec.haplotypes, 'tv.rec.nuc.div'=tv.rec.nuc.div, 'tv.rec.nuc.div.min'=tv.rec.nuc.div.min, 'tv.rec.nuc.div.max'=tv.rec.nuc.div.max,
                        'source.nuc.div'=source.nuc.div, 'rec.nuc.div'=rec.nuc.div, 'Nei.nuc.div'=pair.nuc.div,
                        'source.haplotypes'=source.haplotypes, 'rec.haplotypes'=rec.haplotypes,'pair.haplotypes'=pair.haplotypes,
                        "topology"=t.p.vectors[['vector.topology']],
                        "lineages.source"=t.p.vectors[['vector.lineages.source']],
                        "lineages.rec"=t.p.vectors[['vector.lineages.recipient']],
                        "direction.par"=t.p.vectors[['vector.pars.direction']],
                        "direction.par.p"=t.p.vectors[['vector.p.ancestral.pars']],
                        "direction.ml.ER"=t.p.vectors[['vector.ml.ER.direction']],
                        "direction.ml.ER.p"=t.p.vectors[['vector.p.ancestral.ml.ER']],
                        stringsAsFactors = FALSE)
        if(!file.exists('results.sim.csv')) write.table(result,'results.sim.csv', row.names = F, sep=",") else write.table(append=T, result,'results.sim.csv', row.names = F, sep=",", col.names = F)
        write.tree(sim.tr, file = "trees.sim.tre", append = TRUE, digits = 10, tree.names = FALSE)
        write.tree(trees, file = "trees.tre", append = TRUE, digits = 10, tree.names = FALSE)
}   
