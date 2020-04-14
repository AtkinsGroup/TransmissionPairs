pair.founder.p.R <- function(sims,empirical){
    require(dplyr)
    vector.is.empty <- function(x) return(length(x) ==0 )
    Mode <- function(x) {
        ux <- unique(x)
        ux[which.max(tabulate(match(x, ux)))]
    }
    gm_mean = function(x, na.rm=TRUE){
        exp(sum(log(x+1), na.rm=na.rm) / length(x))-1
    }
    mysims <- read.csv(sims, stringsAsFactors = FALSE)
    empirical <- read.csv(empirical,stringsAsFactors = FALSE)
    mysims$tv.rec.haplotypes <- as.factor(mysims$tv.rec.haplotypes)
    mysims$tv.source.haplotypes <- as.factor(mysims$tv.source.haplotypes)
    empirical <- empirical[which(empirical$cluster_name %in% unique(mysims$cluster_name)),]
    names(mysims)[which(names(mysims)=="cluster_name")] <- "pair"
    names(empirical)[which(names(empirical)=="cluster_name")] <- "pair"
    #probabilities
    for(loopedVt in unique(mysims$vt.rec)){
        df.a <- mysims %>% filter(vt.rec %in% c(loopedVt) & transmission.char=='acute')
        df.c <- mysims %>% filter(vt.rec %in% c(loopedVt) & transmission.char=='chronic')
        df<- rbind(df.a, df.c)
        df$tv.rec.haplotypes <- factor(df$tv.rec.haplotypes)
        df.p <- data.frame('pair'=character(), 'transmission.char'=character(), 'risk'=character(),
                           'transmission'=numeric(), 'sampling.source'=numeric(), 'sampling.rec'=numeric(),
                           'n.source'=numeric(),'n.rec'=numeric(),
                           'MM'=numeric(), 'pMM'=numeric(), 'pMM1'=numeric(), 'pMM2'=numeric(), 
                           'PM'=numeric(), 'pPM'=numeric(), 'pPM1'=numeric(), 'pPM2'=numeric(), 
                           'PP'=numeric(), 'pPP'=numeric(), 'pPP1'=numeric(), 'pPP2'=numeric())
        for(pair in sort(as.character(unique(df$pair)))){ #edit this properly when AD17 is fine
            myvector <- numeric()
            n_obs <- nrow(df[which(df$pair==pair),])
            risk <- paste(as.character(unique(df[which(df$pair==pair),]$risk)), collapse="/")
            transmission.char <- paste(as.character(unique(df[which(df$pair==pair),]$transmission.char)), collapse="/")
            if(transmission.char=="acute/chronic") transmission.char <- "acute"
            transmission <- median(df[which(df$pair==pair),]$transmission)
            sampling.source <- median(df[which(df$pair==pair),]$sampling.source)
            sampling.rec <- median(df[which(df$pair==pair),]$sampling.rec)
            n.source <- median(df[which(df$pair==pair),]$n.source)
            n.rec <- median(df[which(df$pair==pair),]$n.rec)
            for(topology in c('MM', 'PM', 'PP')){
                
                df %>% filter(pair=={{pair}}) %>% group_by(topology) %>%
                    group_by(tv.source.haplotypes, add=TRUE) %>% group_by(tv.rec.haplotypes, add=TRUE) %>%
                    summarise(n=n()) %>% group_by(topology) %>% mutate(prop=n*100/sum(n)) -> mydata
                mydata %>% filter(topology=={{topology}}) -> mydata
                mydata$tv.rec.haplotypes <- as.numeric(mydata$tv.rec.haplotypes)
                mydata$tv.source.haplotypes <- as.numeric(mydata$tv.source.haplotypes)
                p11<- mydata %>% filter(tv.source.haplotypes==1 & tv.rec.haplotypes==1) %>% .['prop'] %>% mutate(prop=prop*0.73) 
                p21<- mydata %>% filter(tv.source.haplotypes==2 & tv.rec.haplotypes==1) %>% .['prop'] %>% mutate(prop=prop*0.27) 
                p12<- mydata %>% filter(tv.source.haplotypes==1 & tv.rec.haplotypes>1) %>% .['prop'] %>% mutate(prop=prop*0.73) 
                p22<- mydata %>% filter(tv.source.haplotypes==2 & tv.rec.haplotypes>1) %>% .['prop'] %>% mutate(prop=prop*0.27) 
                p11<- ifelse(nrow(p11)>0, sum(p11), 0)
                p21<- ifelse(nrow(p21)>0, sum(p21), 0)
                p12<- ifelse(nrow(p12)>0, sum(p12), 0)
                p22<- ifelse(nrow(p22)>0, sum(p22), 0)
                #p(T|Fr=1)
                #p(T|Fr=2)
                pTr1 <- (p11+p21)/(p11+p21+p12+p22)
                pTr2 <- (p12+p22)/(p11+p21+p12+p22)
                myvector<-c(myvector, sum(mydata$n), sum(mydata$n)*100/n_obs, pTr1, pTr2)
            }
            myvector<- c(pair, transmission.char, risk, transmission, 
                         sampling.source, sampling.rec,
                         n.source,n.rec, myvector)
            names(myvector) <- c('pair','transmission.char', 'risk',
                                 'transmission', 'sampling.source', 'sampling.rec',
                                 'n.source','n.rec',
                                 'MM','pMM','pMM1', 'pMM2',
                                 'PM', 'pPM','pPM1', 'pPM2', 
                                 'PP', 'pPP', 'pPP1', 'pPP2')
            myvector <- t(as.data.frame(myvector))
            rownames(myvector)<- NULL
            myvector<- as.data.frame(myvector)
            df.p <- rbind(df.p, myvector)
            rm(p11, p12, p21, p22, pair, transmission, 
               sampling.source, sampling.rec, pTr1, pTr2,
               n.source,n.rec, transmission.char, myvector, topology)
        }
        df.p$pMM1<- as.numeric(as.character(df.p$pMM1))
        df.p$pPM1<- as.numeric(as.character(df.p$pPM1))
        df.p$pPP1<- as.numeric(as.character(df.p$pPP1))
        df.p<- merge(empirical, df.p, by="pair")
        assign(paste("df.p.", loopedVt, sep = ""), df.p)
    }
    #particles
    bestmodel <- data.frame('status'=character(), 'risk'=character(), 'pair'=character(), 'model1'=numeric(),'model2'=numeric(),'model3'=numeric(), 'model4'=numeric(), 'model5'=numeric(), 'model6'=numeric(),
                            'model7'=numeric(), 'model8'=numeric(), 'model9'=numeric(), 'model10'=numeric(),'model11'=numeric(), 'model12'=numeric(),
                            'D'=character(), "M1"=character(), "M2"=character(), "M3"=character(), "M4"=character(), "M5"=character(), "M6"=character(), "M7"=character(), "M8"=character(),
                            "M9"=character(), "M10"=character(), "M11"=character(), "M12"=character())
    for(i in unique(mysims$pair)){
        simulations<- data.frame('model1'=if(vector.is.empty(as.numeric(as.character(unlist(df.p.1[which(df.p.1$pair==i),c('MM','PM','PP')]))))) c(0,0,0) else as.numeric(as.character(unlist(df.p.1[which(df.p.1$pair==i),c('MM','PM','PP')]))),
                                 'model2'=if(vector.is.empty(as.numeric(as.character(unlist(df.p.2[which(df.p.2$pair==i),c('MM','PM','PP')]))))) c(0,0,0) else as.numeric(as.character(unlist(df.p.2[which(df.p.2$pair==i),c('MM','PM','PP')]))),
                                 'model3'=if(vector.is.empty(as.numeric(as.character(unlist(df.p.3[which(df.p.3$pair==i),c('MM','PM','PP')]))))) c(0,0,0) else as.numeric(as.character(unlist(df.p.3[which(df.p.3$pair==i),c('MM','PM','PP')]))),
                                 'model4'=if(vector.is.empty(as.numeric(as.character(unlist(df.p.4[which(df.p.4$pair==i),c('MM','PM','PP')]))))) c(0,0,0) else as.numeric(as.character(unlist(df.p.4[which(df.p.4$pair==i),c('MM','PM','PP')]))),
                                 'model5'=if(vector.is.empty(as.numeric(as.character(unlist(df.p.5[which(df.p.5$pair==i),c('MM','PM','PP')]))))) c(0,0,0) else as.numeric(as.character(unlist(df.p.5[which(df.p.5$pair==i),c('MM','PM','PP')]))),
                                 'model6'=if(vector.is.empty(as.numeric(as.character(unlist(df.p.6[which(df.p.6$pair==i),c('MM','PM','PP')]))))) c(0,0,0) else as.numeric(as.character(unlist(df.p.6[which(df.p.6$pair==i),c('MM','PM','PP')]))),
                                 'model7'=if(vector.is.empty(as.numeric(as.character(unlist(df.p.7[which(df.p.7$pair==i),c('MM','PM','PP')]))))) c(0,0,0) else as.numeric(as.character(unlist(df.p.7[which(df.p.7$pair==i),c('MM','PM','PP')]))),
                                 'model8'=if(vector.is.empty(as.numeric(as.character(unlist(df.p.8[which(df.p.8$pair==i),c('MM','PM','PP')]))))) c(0,0,0) else as.numeric(as.character(unlist(df.p.8[which(df.p.8$pair==i),c('MM','PM','PP')]))),
                                 'model9'=if(vector.is.empty(as.numeric(as.character(unlist(df.p.9[which(df.p.9$pair==i),c('MM','PM','PP')]))))) c(0,0,0) else as.numeric(as.character(unlist(df.p.9[which(df.p.9$pair==i),c('MM','PM','PP')]))),
                                 'model10'=if(vector.is.empty(as.numeric(as.character(unlist(df.p.10[which(df.p.10$pair==i),c('MM','PM','PP')]))))) c(0,0,0) else as.numeric(as.character(unlist(df.p.10[which(df.p.10$pair==i),c('MM','PM','PP')]))),
                                 'model11'=if(vector.is.empty(as.numeric(as.character(unlist(df.p.11[which(df.p.11$pair==i),c('MM','PM','PP')]))))) c(0,0,0) else as.numeric(as.character(unlist(df.p.11[which(df.p.11$pair==i),c('MM','PM','PP')]))),
                                 'model12'=if(vector.is.empty(as.numeric(as.character(unlist(df.p.12[which(df.p.12$pair==i),c('MM','PM','PP')]))))) c(0,0,0) else as.numeric(as.character(unlist(df.p.12[which(df.p.12$pair==i),c('MM','PM','PP')]))))
        mydata<- round(as.numeric(as.character(unlist(df.p.1[which(df.p.1$pair==i),c('MM.p','PM.p','PP.p')]))))
        raw.simulations <- simulations
        for(j in 1:ncol(simulations)){
            if(sum(simulations[,j]) != 0) simulations[,j] <- simulations[,j] + 1
            
        }
        myrawdata <- mydata
        mydata <- mydata + 1
        likelihood1 <- if (class(try(dmultinom(mydata, prob=as.numeric(simulations[,1])), silent = TRUE)) == "try-error") NA else log(dmultinom(mydata, prob=as.numeric(simulations[,1])))
        likelihood2 <- if (class(try(dmultinom(mydata, prob=as.numeric(simulations[,2])), silent = TRUE)) == "try-error") NA else log(dmultinom(mydata, prob=as.numeric(simulations[,2])))
        likelihood3 <- if (class(try(dmultinom(mydata, prob=as.numeric(simulations[,3])), silent = TRUE)) == "try-error") NA else log(dmultinom(mydata, prob=as.numeric(simulations[,3])))
        likelihood4 <- if (class(try(dmultinom(mydata, prob=as.numeric(simulations[,4])), silent = TRUE)) == "try-error") NA else log(dmultinom(mydata, prob=as.numeric(simulations[,4])))
        likelihood5 <- if (class(try(dmultinom(mydata, prob=as.numeric(simulations[,5])), silent = TRUE)) == "try-error") NA else log(dmultinom(mydata, prob=as.numeric(simulations[,5])))
        likelihood6 <- if (class(try(dmultinom(mydata, prob=as.numeric(simulations[,6])), silent = TRUE)) == "try-error") NA else log(dmultinom(mydata, prob=as.numeric(simulations[,6])))
        likelihood7 <- if (class(try(dmultinom(mydata, prob=as.numeric(simulations[,7])), silent = TRUE)) == "try-error") NA else log(dmultinom(mydata, prob=as.numeric(simulations[,7])))
        likelihood8 <- if (class(try(dmultinom(mydata, prob=as.numeric(simulations[,8])), silent = TRUE)) == "try-error") NA else log(dmultinom(mydata, prob=as.numeric(simulations[,8])))
        likelihood9 <- if (class(try(dmultinom(mydata, prob=as.numeric(simulations[,9])), silent = TRUE)) == "try-error") NA else log(dmultinom(mydata, prob=as.numeric(simulations[,9])))
        likelihood10 <- if (class(try(dmultinom(mydata, prob=as.numeric(simulations[,10])), silent = TRUE)) == "try-error") NA else log(dmultinom(mydata, prob=as.numeric(simulations[,10])))
        likelihood11 <- if (class(try(dmultinom(mydata, prob=as.numeric(simulations[,11])), silent = TRUE)) == "try-error") NA else log(dmultinom(mydata, prob=as.numeric(simulations[,11])))
        likelihood12 <- if (class(try(dmultinom(mydata, prob=as.numeric(simulations[,12])), silent = TRUE)) == "try-error") NA else log(dmultinom(mydata, prob=as.numeric(simulations[,12])))
        bestmodel <- rbind(bestmodel, data.frame('status'=as.character(df.p.1[which(df.p.1$pair==i),c("transmission.char")]), 
                                                 'risk'=as.character(df.p.1[which(df.p.1$pair==i),c("risk")]), 
                                                 'pair'=i, 
                                                 'model1'=likelihood1,'model2'=likelihood2,'model3'=likelihood3, 'model4'=likelihood4,'model5'=likelihood5,'model6'=likelihood6,
                                                 'model7'=likelihood7,'model8'=likelihood8,'model9'=likelihood9, 'model10'=likelihood10,'model11'=likelihood11,'model12'=likelihood12,
                                                 'D'=paste(myrawdata, collapse = " "),
                                                 "M1"=paste(raw.simulations[,1], collapse = " "),
                                                 "M2"=paste(raw.simulations[,2], collapse = " "),
                                                 "M3"=paste(raw.simulations[,3], collapse = " "),
                                                 "M4"=paste(raw.simulations[,4], collapse = " "),
                                                 "M5"=paste(raw.simulations[,5], collapse = " "),
                                                 "M6"=paste(raw.simulations[,6], collapse = " "),
                                                 "M7"=paste(raw.simulations[,7], collapse = " "),
                                                 "M8"=paste(raw.simulations[,8], collapse = " "),
                                                 "M9"=paste(raw.simulations[,9], collapse = " "),
                                                 "M10"=paste(raw.simulations[,10], collapse = " "),
                                                 "M11"=paste(raw.simulations[,11], collapse = " "),
                                                 "M12"=paste(raw.simulations[,12], collapse = " ")))
    }
    #ci
    model<- character()
    for(i in unique(mysims$pair)){
        paste(min(as.numeric(na.omit(as.numeric(unlist(strsplit(names(bestmodel[which(bestmodel$pair==i),c(names(bestmodel[4:15]))])[which(bestmodel[which(bestmodel$pair==i), c(names(bestmodel[4:15]))]==max(bestmodel[which(bestmodel$pair==i),c(names(bestmodel[4:15]))],na.rm = TRUE))], "[^0-9]+")))))))
        model <- c(model,     paste(min(as.numeric(na.omit(as.numeric(unlist(strsplit(names(bestmodel[which(bestmodel$pair==i),c(names(bestmodel[4:15]))])[which(bestmodel[which(bestmodel$pair==i), c(names(bestmodel[4:15]))]==max(bestmodel[which(bestmodel$pair==i),c(names(bestmodel[4:15]))],na.rm = TRUE))], "[^0-9]+")))))))
        )
    }
    model.ci.u <- character()
    model.ci.l <- character()
    for(i in unique(mysims$pair)){
        best <- paste(min(as.numeric(na.omit(as.numeric(unlist(strsplit(names(bestmodel[which(bestmodel$pair==i),c(names(bestmodel[4:15]))])[which(bestmodel[which(bestmodel$pair==i), c(names(bestmodel[4:15]))]==max(bestmodel[which(bestmodel$pair==i),c(names(bestmodel[4:15]))],na.rm = TRUE))], "[^0-9]+")))))))
        cl.l <- max(bestmodel[which(bestmodel$pair==i),c(names(bestmodel[4:15]))],na.rm = TRUE) - 1.92
        cl.u <- max(bestmodel[which(bestmodel$pair==i),c(names(bestmodel[4:15]))],na.rm = TRUE) + 1.92
        model.l <- as.numeric(na.omit(as.numeric(unlist(strsplit(names(bestmodel[which(bestmodel$pair==i),c(names(bestmodel[4:15]))])[which(bestmodel[which(bestmodel$pair==i), c(names(bestmodel[4:15]))]>cl.l & bestmodel[which(bestmodel$pair==i), c(names(bestmodel[4:15]))]<cl.u)], "[^0-9]+")))))
        model.ci.u <- c(model.ci.u, max(model.l))
        model.ci.l <- c(model.ci.l, min(model.l))
    }
    bestmodel$bestmodel <- model
    bestmodel$model.ci.u <- model.ci.u
    bestmodel$model.ci.l <- model.ci.l
    #best
    Ds <- character()
    pFr1s <- numeric()
    for(i in unique(mysims$pair)){
        D <- paste("p", as.character(df.p.1$topology[df.p.1$pair==i]), 1, sep="")
        if(length(unlist(strsplit(bestmodel$bestmodel[bestmodel$pair==i], " ")))==1){
            if(bestmodel$bestmodel[bestmodel$pair==i]==1){
                pFr1 <- df.p.1[,D][df.p.1$pair==i]
            }
            if(bestmodel$bestmodel[bestmodel$pair==i]==2){
                pFr1 <- df.p.2[,D][df.p.2$pair==i]
            }
            if(bestmodel$bestmodel[bestmodel$pair==i]==3){
                pFr1 <- df.p.3[,D][df.p.3$pair==i]
            }
            if(bestmodel$bestmodel[bestmodel$pair==i]==4){
                pFr1 <- df.p.4[,D][df.p.4$pair==i]
            }
            if(bestmodel$bestmodel[bestmodel$pair==i]==5){
                pFr1 <- df.p.5[,D][df.p.5$pair==i]
            }
            if(bestmodel$bestmodel[bestmodel$pair==i]==6){
                pFr1 <- df.p.6[,D][df.p.6$pair==i]
            }
            if(bestmodel$bestmodel[bestmodel$pair==i]==7){
                pFr1 <- df.p.7[,D][df.p.7$pair==i]
            }
            if(bestmodel$bestmodel[bestmodel$pair==i]==8){
                pFr1 <- df.p.8[,D][df.p.8$pair==i]
            }
            if(bestmodel$bestmodel[bestmodel$pair==i]==9){
                pFr1 <- df.p.9[,D][df.p.9$pair==i]
            }
            if(bestmodel$bestmodel[bestmodel$pair==i]==10){
                pFr1 <- df.p.10[,D][df.p.10$pair==i]
            }
            if(bestmodel$bestmodel[bestmodel$pair==i]==11){
                pFr1 <- df.p.11[,D][df.p.11$pair==i]
            }
            if(bestmodel$bestmodel[bestmodel$pair==i]==12){
                pFr1 <- df.p.12[,D][df.p.12$pair==i]
            }
            
        } else pFr1 <- NA
        Ds <- c(Ds, as.character(df.p.1$topology[df.p.1$pair==i]))
        pFr1s <- c(pFr1s, pFr1)
    }
    Ds <- character()
    pFr1s <- numeric()
    pFr1Fulls <- numeric()
    for(i in unique(mysims$pair)){
        D <- paste("p", as.character(df.p.1$topology[df.p.1$pair==i]), 1, sep="")
        MM<-empirical$MM.p[empirical$pair==i]/100
        PM<-empirical$PM.p[empirical$pair==i]/100
        PP<-empirical$PP.p[empirical$pair==i]/100
        if(length(unlist(strsplit(bestmodel$bestmodel[bestmodel$pair==i], " ")))==1){
            if(bestmodel$bestmodel[bestmodel$pair==i]==1){
                pFr1Full <- sum(as.vector(na.omit(c(df.p.1[,"pMM1"][df.p.1$pair==i]*MM, df.p.1[,"pPM1"][df.p.1$pair==i]*PM, df.p.1[,"pPP1"][df.p.1$pair==i]*PP))))
                pFr1 <- df.p.1[,D][df.p.1$pair==i]
                
            }
            if(bestmodel$bestmodel[bestmodel$pair==i]==2){
                pFr1Full <- sum(as.vector(na.omit(c(df.p.2[,"pMM1"][df.p.2$pair==i]*MM , df.p.2[,"pPM1"][df.p.2$pair==i]*PM , df.p.2[,"pPP1"][df.p.2$pair==i]*PP))))
                pFr1 <- df.p.2[,D][df.p.2$pair==i]
            }
            if(bestmodel$bestmodel[bestmodel$pair==i]==3){
                pFr1Full <- sum(as.vector(na.omit(c(df.p.3[,"pMM1"][df.p.3$pair==i]*MM , df.p.3[,"pPM1"][df.p.3$pair==i]*PM , df.p.3[,"pPP1"][df.p.3$pair==i]*PP))))
                pFr1 <- df.p.3[,D][df.p.3$pair==i]
                
            }
            if(bestmodel$bestmodel[bestmodel$pair==i]==4){
                pFr1Full <- sum(as.vector(na.omit(c(df.p.4[,"pMM1"][df.p.4$pair==i]*MM , df.p.4[,"pPM1"][df.p.4$pair==i]*PM , df.p.4[,"pPP1"][df.p.4$pair==i]*PP))))
                pFr1 <- df.p.4[,D][df.p.4$pair==i]
                
            }
            if(bestmodel$bestmodel[bestmodel$pair==i]==5){
                pFr1Full <- sum(as.vector(na.omit(c(df.p.5[,"pMM1"][df.p.5$pair==i]*MM , df.p.5[,"pPM1"][df.p.5$pair==i]*PM , df.p.5[,"pPP1"][df.p.5$pair==i]*PP))))
                pFr1 <- df.p.5[,D][df.p.5$pair==i]
            }
            if(bestmodel$bestmodel[bestmodel$pair==i]==6){
                pFr1Full <- sum(as.vector(na.omit(c(df.p.6[,"pMM1"][df.p.6$pair==i]*MM , df.p.6[,"pPM1"][df.p.6$pair==i]*PM , df.p.6[,"pPP1"][df.p.6$pair==i]*PP))))
                pFr1 <- df.p.6[,D][df.p.6$pair==i]
            }
            if(bestmodel$bestmodel[bestmodel$pair==i]==7){
                pFr1Full <- sum(as.vector(na.omit(c(df.p.7[,"pMM1"][df.p.7$pair==i]*MM , df.p.7[,"pPM1"][df.p.7$pair==i]*PM , df.p.7[,"pPP1"][df.p.7$pair==i]*PP))))
                pFr1 <- df.p.7[,D][df.p.7$pair==i]
            }
            if(bestmodel$bestmodel[bestmodel$pair==i]==8){
                pFr1Full <- sum(as.vector(na.omit(c(df.p.8[,"pMM1"][df.p.8$pair==i]*MM , df.p.8[,"pPM1"][df.p.8$pair==i]*PM , df.p.8[,"pPP1"][df.p.8$pair==i]*PP))))
                pFr1 <- df.p.8[,D][df.p.8$pair==i]
            }
            if(bestmodel$bestmodel[bestmodel$pair==i]==9){
                pFr1Full <- sum(as.vector(na.omit(c(df.p.9[,"pMM1"][df.p.9$pair==i]*MM , df.p.9[,"pPM1"][df.p.9$pair==i]*PM , df.p.9[,"pPP1"][df.p.9$pair==i]*PP))))
                pFr1 <- df.p.9[,D][df.p.9$pair==i]
            }
            if(bestmodel$bestmodel[bestmodel$pair==i]==10){
                pFr1Full <- sum(as.vector(na.omit(c(df.p.10[,"pMM1"][df.p.10$pair==i]*MM , df.p.10[,"pPM1"][df.p.10$pair==i]*PM , df.p.10[,"pPP1"][df.p.10$pair==i]*PP))))
                pFr1 <- df.p.10[,D][df.p.10$pair==i]
            }
            if(bestmodel$bestmodel[bestmodel$pair==i]==11){
                pFr1Full <- sum(as.vector(na.omit(c(df.p.11[,"pMM1"][df.p.11$pair==i]*MM , df.p.11[,"pPM1"][df.p.11$pair==i]*PM , df.p.11[,"pPP1"][df.p.11$pair==i]*PP))))
                pFr1 <- df.p.11[,D][df.p.11$pair==i]
            }
            if(bestmodel$bestmodel[bestmodel$pair==i]==12){
                pFr1Full <- sum(as.vector(na.omit(c(df.p.12[,"pMM1"][df.p.12$pair==i]*MM , df.p.12[,"pPM1"][df.p.12$pair==i]*PM , df.p.12[,"pPP1"][df.p.12$pair==i]*PP))))
                pFr1 <- df.p.12[,D][df.p.12$pair==i]
            }
            
        } else pFr1 <- NA
        Ds <- c(Ds, as.character(df.p.1$topology[df.p.1$pair==i]))
        pFr1s <- c(pFr1s, pFr1)
        pFr1Fulls <- c(pFr1Fulls, pFr1Full)
        
    }
    bestmodel$topologyD <- Ds    
    bestmodel$pFr1 <- pFr1s
    bestmodel$pFr1Full <- pFr1Fulls
    #founders
    Mode <- function(x) {
        ux <- unique(x)
        ux[which.max(tabulate(match(x, ux)))]
    }
    singleHaplotypeProp <- meanNhaplotypes <- modeNhaplotypes <- maxNhaplotypes <- meanNhaplotypes.adjusted <- numeric()
    for(i in unique(bestmodel$pair)){
        bestmodel.num <- bestmodel$bestmodel[bestmodel$pair==i]
        mypair <- as.numeric(mysims$tv.rec.haplotypes[mysims$pair==i & mysims$vt.rec==bestmodel.num])
        singleHaplotypeProp <- c(singleHaplotypeProp, sum(mypair==1)*100/length(mypair))
        meanNhaplotypes <- c(meanNhaplotypes, mean(mypair))
        modeNhaplotypes <- c(modeNhaplotypes, Mode(mypair))
        maxNhaplotypes <- c(maxNhaplotypes, max(mypair))
        mypair <- mypair[mypair!=1]
        if(vector.is.empty(mypair)) meanNhaplotypes.adjusted <- c(meanNhaplotypes.adjusted, 0) else  meanNhaplotypes.adjusted <- c(meanNhaplotypes.adjusted, mean(mypair))
    }
    bestmodel$singleHaplotypeProp <- singleHaplotypeProp
    bestmodel$meanNhaplotypes <- meanNhaplotypes
    bestmodel$modeNhaplotypes <- modeNhaplotypes
    bestmodel$maxNhaplotypes <- maxNhaplotypes
    bestmodel$meanNhaplotypes.adjusted <- meanNhaplotypes.adjusted
    write.csv(bestmodel, 'bestmodel.csv', row.names = F)
    print(paste("probability of one founder strains is ", gm_mean(bestmodel$pFr1Full), sep=""))
}



