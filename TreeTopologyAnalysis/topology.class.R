topology.class <- function(tr, states, outgroup, p.ancestral.logs, logs) {
    require(ape)
    require(phytools) #reroot
    require(adephylo) #listDD
    require(dplyr) #pipes
    require(tibble) #column_to_rownames
    require(phangorn)
    #useful mini functions (operator inverse of %in%)
    '%!in%' <- function(x, y) ! ('%in%'(x, y))
    #create states variables for individuals
    if(missing(states)) states = states
    names(states) <- c("accession", "state")
    states.source = 0
    states.recipient = 1
    if(missing(outgroup)) outgroup = TRUE
    if(missing(p.ancestral.logs)) p.ancestral.logs = TRUE
    if(missing(logs)) logs = logs
    if(!exists("simulation")) simulation <- 1
    #vectors to save data
    vector.lineages.source <- integer()
    vector.lineages.recipient <- integer()
    vector.topology <- character()
    vector.p.ancestral.pars <- numeric()
    vector.p.ancestral.ml.ER <- numeric()
    vector.p.ancestral.bayes <- numeric()
    vector.pars.direction <- character()
    vector.ml.ER.direction <- character()
    vector.bayes.direction <- character()
    vector.tree.length <- numeric() 
    vector.tree.SBL <- numeric() 
    vector.tree.SPD <- numeric() 
    vector.tree.MBL <- numeric() 
    vector.tree.MPD <- numeric() 
    vector.tree.VBL <- numeric()
    vector.tree.VPD <- numeric()
    if(outgroup) {
        states <- states[2:nrow(states),] #remove HXB2
        rownames(states) <- NULL
        states$state <- factor(states$state)
    }
    for (i in 1:length(tr)){
        if(outgroup) {
            troot <- reroot(tr[[i]], 1) #phytools #root to outgroup
            troot<- drop.tip(troot, 1) #remove HXB2
        } else troot <- tr[[i]] 
        if(i %in% c(1,seq(length(tr)*0.1,length(tr),length(tr)*0.1))) print(paste("reading tree number ",i,sep=""))
        #read tree
        taxa.source <- as.character(states$accession[which(states$state == states.source)]) #vector of source names
        taxa.rec <- as.character(states$accession[which(states$state == states.recipient)]) #vector of recipients names
        rootNode<- Ntip(troot)+1 #get node number for root
        nodes <- 1:troot$Nnode+Ntip(troot) #create vector of internal nodes
        nodes.eval <- nodes #create vector of internal nodes that will be iteratively updated
        lineages.recipient<- 0 #create vector of lineage counts (in recipient) to be iteratively updated
        lineages.source<- 0 #create ector of lineage counts (in source) to be iteratively updated 
        tips<- 1:Ntip(troot) #create vector of tips
        DD <-listDD(troot, nameBy=c("label","number")) #require(adephylo) #List direct descendants for all nodes of a tree
        for(node in nodes){ #Loop internal nodes
            if(node %in% nodes.eval){ #Check if to be evaluated
                Children <- as.integer(unlist(DD[as.character(node)])) #retrieve children of node (direct descendants)
                for(child in Children){ #evaluate each child separately
                    numDescendants <- getDescendants(troot, child) #list descendants
                    numnodes <- numDescendants[which(numDescendants %in% nodes)] #keep internal nodes, remove tips
                    charDescendants <- troot$tip.label[c(numDescendants)] #get names of descendants
                    chartips <- charDescendants[!is.na(charDescendants)] #remove internal nodes, only keep tips
                    taxtips <- as.character(states$state[which(states$accession %in% chartips)]) #get states (source/recipient) of tips
                    if(child %in% tips){ #the child is a tip
                        if(troot$tip.label[child] %in% taxa.rec){lineages.recipient <- lineages.recipient + 1} else {lineages.source <- lineages.source + 1}
                    }
                    else{
                        if(length(unique(taxtips)) == 1){ #if all descedants are the same
                            nodes.eval <- nodes.eval[! nodes.eval %in% c(child, numnodes)] #remove nodes (they are monophyletic) from nodes.eval
                            if((unique(taxtips))==states.recipient) { lineages.recipient <- lineages.recipient + 1 } else {lineages.source <- lineages.source + 1} 
                        }
                    } 
                }
                rm(Children, child, numDescendants, numnodes, charDescendants, chartips, taxtips)
            }
        }
        if(lineages.recipient==1 & lineages.source==1) topology="MM" else { #MM, PM, PP?
            if (lineages.recipient>=2 & lineages.source==1 | lineages.recipient==1 & lineages.source>=2)  topology="PM" else topology = "PP"
        }
        #direction
        states.pars <- states  %>% column_to_rownames("accession") %>% as.matrix %>% as.phyDat(type="USER", levels=c("0","1")) %>% ancestral.pars(troot,., type="MPR")
        p.ancestral.pars <- states.pars[[Ntip(troot)+1]][1]
        states.ml.ER <- setNames(factor(states$state),factor(states$accession)) %>% ace(.,troot,method="ML", model="ER",type="discrete") 
        p.ancestral.ml.ER <- states.ml.ER$lik.anc[1]
        if(p.ancestral.pars >= 0.9){pars.direction<-"consistent"}else{
            if(p.ancestral.pars <= 0.1){pars.direction<-"inconsistent"}else{
                pars.direction<-"equivocal"
            }
        }
        if(p.ancestral.ml.ER >= 0.9){ml.ER.direction<-"consistent"}else{
            if(p.ancestral.ml.ER <= 0.1){ml.ER.direction<-"inconsistent"}else{
                ml.ER.direction<-"equivocal"
            }
        }
        #append data to vectors
        vector.lineages.source <- c(vector.lineages.source, lineages.source)
        vector.lineages.recipient <- c(vector.lineages.recipient, lineages.recipient)
        vector.topology <- c(vector.topology, topology)
        vector.p.ancestral.pars <- c(vector.p.ancestral.pars, p.ancestral.pars) 
        vector.p.ancestral.ml.ER <- c(vector.p.ancestral.ml.ER, p.ancestral.ml.ER) 
        vector.pars.direction <- c(vector.pars.direction, pars.direction) 
        vector.ml.ER.direction <- c(vector.ml.ER.direction, ml.ER.direction) 
        
        #p.BI.signal #this required a corresponding loaded log file
        if(p.ancestral.logs){
            p.ancestral.bayes <- logs[i, 19]
            if(p.ancestral.bayes >= 0.9){bayes.direction<-"consistent"}else{
                if(p.ancestral.bayes <= 0.1){bayes.direction<-"inconsistent"}else{
                    bayes.direction<-"equivocal"
                }
            }
            vector.p.ancestral.bayes <- c(vector.p.ancestral.bayes, p.ancestral.bayes) 
            vector.bayes.direction <- c(vector.bayes.direction, bayes.direction)
        }
        vector.tree.length <- c(vector.tree.length, max(distRoot(troot)))
        vector.tree.SBL <- c(vector.tree.SBL,sum(troot$edge.length)) 
        vector.tree.SPD <- c(vector.tree.SPD,sum(distTips(troot))) 
        vector.tree.MBL <- c(vector.tree.MBL, mean(troot$edge.length)) 
        vector.tree.MPD <- c(vector.tree.MPD, mean(distTips(troot))) 
        vector.tree.VBL <- c(vector.tree.VBL, var(troot$edge.length))
        vector.tree.VPD <- c(vector.tree.VPD, var(distTips(troot)))
    }
    if(p.ancestral.logs){
        logs[, "p(0).pars"] <- vector.p.ancestral.pars
        logs[,"p(0).ml.ER"] <- vector.p.ancestral.ml.ER
        logs[, "simulation"] <- simulation
    }
    return(list('vector.lineages.source'=vector.lineages.source, 'vector.lineages.recipient'=vector.lineages.recipient, 
                'vector.topology'=vector.topology, 'vector.p.ancestral.pars'=vector.p.ancestral.pars,
                'vector.p.ancestral.ml.ER'=vector.p.ancestral.ml.ER, 'vector.p.ancestral.bayes'=vector.p.ancestral.bayes,
                'vector.pars.direction'=vector.pars.direction,'vector.ml.ER.direction'=vector.ml.ER.direction,
                'vector.bayes.direction'=vector.bayes.direction, 'vector.tree.length'=vector.tree.length, 
                'vector.tree.SBL'=vector.tree.SBL, 'vector.tree.SPD'=vector.tree.SPD, 
                'vector.tree.MBL'=vector.tree.MBL, 'vector.tree.MPD'=vector.tree.MPD, 
                'vector.tree.VBL'=vector.tree.VBL, 'vector.tree.VPD'=vector.tree.VPD, 'logs'=logs))
    
}

