pair.mb.summary <- function(ids, sets) {
    ####################################
    #load libraries
    ####################################
    require(ape)
    require(mcmcse) #ESS calculations
    ####################################
    #set values for function arguments
    ####################################
    if (missing(sets))
        sets = "set.csv"
    ####################################
    #loop every cluster in ids
    ####################################
    for (id in ids) {
        ####################################
        #get source/recipient ids from files
        ####################################
        cluster.name <- id
        print(paste("Processing ", cluster.name, sep = ""))
        #read data
        info.gb <-
            read.csv(paste(cluster.name, "/", "info.gb.csv", sep = ""),
                     stringsAsFactors = F)
        info.gb <-
            info.gb[, c(
                "row_id",
                "patient_id",
                "accession_id",
                "seq_length",
                "pos",
                "pos_ini",
                "pos_end"
            )]
        info <-
            read.csv(paste(cluster.name, "/", "info.csv", sep = ""),
                     stringsAsFactors = F)
        cluster <- unique(info$ClusterId) 
        source <-
            na.omit(info$PatientId[info$States == 'source']) 
        recipient <-
            na.omit(info$PatientId[info$States == 'recipient']) 
        ####################################
        #load datasets
        ####################################
        set.id <-
            read.csv(
                paste(id, "/", sets, sep = ""),
                row.names = NULL,
                stringsAsFactors = F
            ) 
        ####################################
        #loop datasets--ignore large datasets
        ####################################
        for (set in 1:nrow(set.id)) {
            if (set.id$counts[set] > 300) {
                set.id$MCMC.AvgStdDev[set] <- NA
                set.id$minESS[set] <- NA
                next
            }
            ####################################
            #list files to be processed
            ####################################
            t.files <-
                list.files(
                    path = paste(id, "/", id, ".", set.id$ref[set], sep = ""),
                    pattern = "\\.t$",
                    full.names = TRUE
                )
            p.files <-
                list.files(
                    path = paste(id, "/", id, ".", set.id$ref[set], sep = ""),
                    pattern = "\\.p$",
                    full.names = TRUE
                )
            ####################################
            #list relevant param
            ####################################
            ngen <- as.integer(set.id$MCMC.ngen[set])
            mcmc <-
                read.table(
                    paste(
                        id,
                        "/",
                        id,
                        ".",
                        set.id$ref[set],
                        "/",
                        id,
                        ".",
                        set.id$ref[set],
                        ".mb.nex.mcmc",
                        sep = ""
                    ),
                    header = TRUE,
                    comment.char = "["
                )
            AvgStdDev <- mcmc[nrow(mcmc), "AvgStdDev.s."]
            ####################################
            #check analysis integrity
            ####################################
            if (is.nan(AvgStdDev)) {
                set.id$MCMC.AvgStdDev[set]  <- NA
                set.id$minESS[set]  <- NA
                set.id$trees[set]  <- NA
                set.id$trees.TL[set]  <- NA
                set.id$trees.SBL[set]  <- NA
                set.id$trees.SPD[set]  <- NA
                set.id$trees.MBL[set]  <- NA
                set.id$trees.MPD[set]  <- NA
                set.id$trees.VBL[set]  <- NA
                set.id$trees.VPD[set]  <- NA
                set.id$MM[set]  <- NA
                set.id$PM[set]  <- NA
                set.id$PP[set]  <- NA
                set.id$topology[set]  <- NA
                set.id$MM.p[set]  <- NA
                set.id$PM.p[set]  <- NA
                set.id$PP.p[set]  <- NA
                set.id$confidence.95.t[set]  <- NA
                set.id$lineages.median[set]  <- NA
                set.id$lineages.min[set]  <- NA
                set.id$lineages.max[set]  <- NA
                set.id$ft.median[set]  <- NA
                set.id$ft.min[set]  <- NA
                set.id$ft.max[set]  <- NA
                set.id$Consistent.par[set]  <- NA
                set.id$inconsistent.par[set]  <- NA
                set.id$equivocal.par[set]  <- NA
                set.id$consistent.par.p[set]  <- NA
                set.id$inconsistent.par.p[set]  <- NA
                set.id$equivocal.par.p[set]  <- NA
                set.id$confidence.par.90.s[set]  <- NA
                set.id$Consistent.ml[set]  <- NA
                set.id$inconsistent.ml[set]  <- NA
                set.id$equivocal.ml[set]  <- NA
                set.id$consistent.ml.p[set]  <- NA
                set.id$inconsistent.ml.p[set]  <- NA
                set.id$equivocal.ml.p[set]  <- NA
                set.id$confidence.ml.90.s[set]  <- NA
            } else{
                ####################################
                #process log files
                ####################################
                system(
                    paste(
                        "head -2 ",
                        p.files[1],
                        " > ",
                        id,
                        "/",
                        id,
                        ".",
                        set.id$ref[set],
                        "/",
                        id,
                        ".",
                        set.id$ref[set],
                        ".COMB.log; tail -n +",
                        as.integer(ngen / (ngen * 0.0001)) / 2 + 4,
                        " -q ",
                        id,
                        "/",
                        id,
                        ".",
                        set.id$ref[set],
                        "/*.p >> ",
                        id,
                        "/",
                        id,
                        ".",
                        set.id$ref[set],
                        "/",
                        id,
                        ".",
                        set.id$ref[set],
                        ".COMB.log",
                        sep = ""
                    )
                )
                resampling <-
                    as.integer(seq(
                        as.integer(ngen * 0.0001),
                        ngen,
                        as.integer(ngen * 0.0001)
                    ))
                resampling.names <- paste("gen_", resampling, sep = "")
                logs <-
                    read.delim(
                        paste(
                            id,
                            "/",
                            id,
                            ".",
                            set.id$ref[set],
                            "/",
                            id,
                            ".",
                            set.id$ref[set],
                            ".COMB.log",
                            sep = ""
                        ),
                        comment.char = "["
                    )
                logs$Gen <- as.integer((1:nrow(logs)) * (ngen * 0.0001))
                logs$Joint <- logs$LnL + logs$LnPr
                logs <- logs[which(logs$Gen %in% resampling), ]
                ESS <-
                    as.integer(min(as.vector(na.omit(
                        ess(logs[, 2:ncol(logs)])
                    ))))
                ####################################
                #process tree files
                ####################################
                system(
                    paste(
                        "sed -i '' -- 's/gen./gen_/g' ",
                        id,
                        "/",
                        id,
                        ".",
                        set.id$ref[set],
                        "/*.t",
                        sep = ""
                    )
                ) 
                system(
                    paste(
                        "head -",
                        6 + set.id$counts[set],
                        " ",
                        t.files[1],
                        " > ",
                        id,
                        "/",
                        id,
                        ".",
                        set.id$ref[set],
                        "/",
                        id,
                        ".",
                        set.id$ref[set],
                        ".COMB.trees; awk 'FNR >= ",
                        as.integer(ngen / (ngen * 0.0001)) / 2 + 8 + set.id$counts[set],
                        " && FNR <= ",
                        as.integer(ngen / (ngen * 0.0001)) / 2 + 7 + set.id$counts[set] + ngen /
                            (ngen * 0.0001) / 2,
                        "' ",
                        id,
                        "/",
                        id,
                        ".",
                        set.id$ref[set],
                        "/*.t >> ",
                        id,
                        "/",
                        id,
                        ".",
                        set.id$ref[set],
                        "/",
                        id,
                        ".",
                        set.id$ref[set],
                        ".COMB.trees ; echo \"end;\" >> ",
                        id,
                        "/",
                        id,
                        ".",
                        set.id$ref[set],
                        "/",
                        id,
                        ".",
                        set.id$ref[set],
                        ".COMB.trees",
                        sep = ""
                    )
                )
                trees <-
                    read.nexus(
                        paste(
                            id,
                            "/",
                            id,
                            ".",
                            set.id$ref[set],
                            "/",
                            id,
                            ".",
                            set.id$ref[set],
                            ".COMB.trees",
                            sep = ""
                        )
                    )
                rm(p.files, t.files)
                names(trees) <-
                    paste("gen_", as.integer((1:length(
                        trees
                    )) * (ngen * 0.0001)), sep = "")
                trees <-
                    trees[which(names(trees) %in% resampling.names)]
                write.tree(
                    trees,
                    file = paste(
                        id,
                        "/",
                        id,
                        ".",
                        set.id$ref[set],
                        "/",
                        id,
                        ".",
                        set.id$ref[set],
                        ".COMB.newick.trees",
                        sep = ""
                    )
                ) 
                ####################################
                #load df with transmission pair identity
                ####################################
                states <-
                    read.table(
                        paste(
                            id,
                            "/",
                            id,
                            ".",
                            set.id$ref[set],
                            "/",
                            id,
                            ".",
                            set.id$ref[set],
                            ".states",
                            sep = ""
                        ),
                        header = FALSE,
                        col.names = c('accession', 'state')
                    )
                ####################################
                #calculate topology class
                ####################################
                topological.pattern.list <-
                    topology.class(
                        tr = trees,
                        states = states,
                        outgroup = TRUE,
                        p.ancestral = TRUE,
                        logs = logs
                    )
                ####################################
                #save calculation to logs
                ####################################
                print("Saving tree metrics to logs")
                logs$lineages.source <-
                    vector.lineages.source <-
                    topological.pattern.list[['vector.lineages.source']]
                logs$lineages.recipient <-
                    vector.lineages.recipient <-
                    topological.pattern.list[['vector.lineages.recipient']]
                logs$topology <-
                    vector.topology <- topological.pattern.list[['vector.topology']]
                logs$par.direction.p <-
                    vector.par.direction.p  <-
                    topological.pattern.list[['vector.p.ancestral.pars']]
                logs$par.direction <-
                    vector.par.direction  <-
                    topological.pattern.list[['vector.pars.direction']]
                logs$ml.direction.p <-
                    vector.ml.direction.p <-
                    topological.pattern.list[['vector.p.ancestral.ml.ER']]
                logs$ml.direction <-
                    vector.ml.direction <-
                    topological.pattern.list[['vector.ml.ER.direction']]
                logs$tree.length <-
                    vector.tree.length <-
                    topological.pattern.list[['vector.tree.length']]
                logs$tree.SBL <-
                    vector.tree.SBL <- topological.pattern.list[['vector.tree.SBL']]
                logs$tree.SPD <-
                    vector.tree.SPD <- topological.pattern.list[['vector.tree.SPD']]
                logs$tree.MBL <-
                    vector.tree.MBL <- topological.pattern.list[['vector.tree.MBL']]
                logs$tree.MPD <-
                    vector.tree.MPD <- topological.pattern.list[['vector.tree.MPD']]
                logs$tree.VBL <-
                    vector.tree.VBL <- topological.pattern.list[['vector.tree.VBL']]
                logs$tree.VPD <-
                    vector.tree.VPD <- topological.pattern.list[['vector.tree.VPD']]
                write.table(
                    logs,
                    paste(
                        id,
                        "/",
                        id,
                        ".",
                        set.id$ref[set],
                        "/",
                        id,
                        ".",
                        set.id$ref[set],
                        ".ANN.log",
                        sep = ""
                    ),
                    sep = "\t",
                    row.names = FALSE,
                    quote = FALSE
                )
                ####################################
                #summarize findings
                ####################################
                topology <- c(
                    'MM' = length(which(vector.topology == "MM")),
                    'PM' = length(which(vector.topology == "PM")),
                    'PP' = length(which(vector.topology == "PP"))
                )
                topology.p <- round(topology / nrow(logs) * 100, 2)
                if (as.numeric(topology.p[names(which(topology.p == max(topology.p)))]) >= 95) {
                    confidence.t = "pass"
                } else{
                    confidence.t = "fail"
                }
                pass.par.direction <-
                    vector.par.direction[vector.topology == names(which(topology == max(topology)))]
                pass.ml.direction <-
                    vector.ml.direction[vector.topology == names(which(topology == max(topology)))]
                vector.lineages.recipient <-
                    vector.lineages.recipient[vector.topology == names(which(topology == max(topology)))]
                vector.lineages.source <-
                    vector.lineages.source[vector.topology == names(which(topology == max(topology)))]
                direction.par <-
                    c(
                        'consistent' = length(which(
                            pass.par.direction == "consistent"
                        )),
                        'inconsistent' = length(which(
                            pass.par.direction == "inconsistent"
                        )),
                        'equivocal' = length(which(
                            pass.par.direction == "equivocal"
                        ))
                    )
                direction.ml <-
                    c(
                        'consistent' = length(which(
                            pass.ml.direction == "consistent"
                        )),
                        'inconsistent' = length(which(
                            pass.ml.direction == "inconsistent"
                        )),
                        'equivocal' = length(which(
                            pass.ml.direction == "equivocal"
                        ))
                    )
                direction.par.p <-
                    round(direction.par / sum(direction.par) * 100, 2)
                direction.ml.p <-
                    round(direction.ml / sum(direction.ml) * 100, 2)
                if (as.numeric(direction.par.p[names(which(direction.par.p ==
                                                           max(direction.par.p)))]) >= 90) {
                    confidence.par.d = "pass"
                } else{
                    confidence.par.d = "fail"
                }
                if (as.numeric(direction.ml.p[names(which(direction.ml.p ==
                                                          max(direction.ml.p)))]) >= 90) {
                    confidence.ml.d = "pass"
                } else{
                    confidence.ml.d = "fail"
                }
                ####################################
                #save calculations to datasets
                ####################################
                print("Appending information to set")
                set.id$MCMC.AvgStdDev[set] <- round(AvgStdDev, 5)
                set.id$minESS[set] <- round(ESS, 2)
                set.id$trees[set] <- length(trees)
                set.id$trees.TL[set] <-
                    paste(
                        paste(as.vector(
                            summary(vector.tree.length)
                        ), collapse = " "),
                        var(vector.tree.length),
                        collapse = " "
                    )
                set.id$trees.SBL[set] <-
                    paste(paste(as.vector(
                        summary(vector.tree.SBL)
                    ), collapse = " "),
                    var(vector.tree.SBL),
                    collapse = " ")
                set.id$trees.SPD[set] <-
                    paste(paste(as.vector(
                        summary(vector.tree.SPD)
                    ), collapse = " "),
                    var(vector.tree.SPD),
                    collapse = " ")
                set.id$trees.MBL[set] <-
                    paste(paste(as.vector(
                        summary(vector.tree.MBL)
                    ), collapse = " "),
                    var(vector.tree.MBL),
                    collapse = " ")
                set.id$trees.MPD[set] <-
                    paste(paste(as.vector(
                        summary(vector.tree.MPD)
                    ), collapse = " "),
                    var(vector.tree.MPD),
                    collapse = " ")
                set.id$trees.VBL[set] <-
                    paste(paste(as.vector(
                        summary(vector.tree.VBL)
                    ), collapse = " "),
                    var(vector.tree.VBL),
                    collapse = " ")
                set.id$trees.VPD[set] <-
                    paste(paste(as.vector(
                        summary(vector.tree.VPD)
                    ), collapse = " "),
                    var(vector.tree.VPD),
                    collapse = " ")
                set.id$MM[set] <- topology[1]
                set.id$PM[set] <- topology[2]
                set.id$PP[set] <- topology[3]
                set.id$topology[set] <-
                    names(topology)[which(topology == max(topology))]
                set.id$MM.p[set] <- topology.p[1]
                set.id$PM.p[set] <- topology.p[2]
                set.id$PP.p[set] <- topology.p[3]
                set.id$confidence.95.t[set] <- confidence.t
                set.id$lineages.median[set] <-
                    median(vector.lineages.source)
                set.id$lineages.min[set] <-
                    min(vector.lineages.source)
                set.id$lineages.max[set] <-
                    max(vector.lineages.source)
                set.id$ft.median[set] <-
                    median(vector.lineages.recipient)
                set.id$ft.min[set] <- min(vector.lineages.recipient)
                set.id$ft.max[set] <- max(vector.lineages.recipient)
                set.id$Consistent.par[set] <- direction.par[1]
                set.id$inconsistent.par[set] <- direction.par[2]
                set.id$equivocal.par[set] <- direction.par[3]
                set.id$consistent.par.p[set] <- direction.par.p[1]
                set.id$inconsistent.par.p[set] <- direction.par.p[2]
                set.id$equivocal.par.p[set] <- direction.par.p[3]
                set.id$confidence.par.90.s[set] <- confidence.par.d
                set.id$Consistent.ml[set] <- direction.ml[1]
                set.id$inconsistent.ml[set] <- direction.ml[2]
                set.id$equivocal.ml[set] <- direction.ml[3]
                set.id$consistent.ml.p[set] <- direction.ml.p[1]
                set.id$inconsistent.ml.p[set] <- direction.ml.p[2]
                set.id$equivocal.ml.p[set] <- direction.ml.p[3]
                set.id$confidence.ml.90.s[set] <- confidence.ml.d
                rm(
                    vector.lineages.recipient,
                    vector.lineages.source,
                    vector.ml.direction,
                    vector.par.direction,
                    vector.topology,
                    vector.tree.length,
                    vector.tree.MBL,
                    vector.tree.MPD,
                    vector.tree.SBL,
                    vector.tree.SPD,
                    vector.tree.VBL,
                    vector.tree.VPD
                )
                
            }
        }
        ####################################
        #save data
        ####################################
        write.csv(set.id, paste(id, "/", sets, sep = ""), row.names = FALSE)
    }
    
}
