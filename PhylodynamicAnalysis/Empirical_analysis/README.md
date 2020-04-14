## *pair.alignment* <br/><br/>

Aligns sequence data and create files for downstream analysis.

If there are multiple genetic region (rows) in the input file, each one is analyzed separately.

**Outputs:**

- Genetic regions input file (e.g. either set.csv or sets.csv) is updated with the alignment length, the gap percentage (compared to the expected length based on HBX2 coordinates) and calculations of Nei (1987) nucleotide diversity. 
- Directory with LANLdb transmission cluster names and HBX2 coordinates of the aligned genetic region:
  - ingroup.fasta (viral genetic sequences in fasta format)
  - outgroup.HBX2 (multiple sequence alignment of outgroup and HBX2)
  - output.fasta (multiple sequence alignment of ingroup/outgroup in fasta format before trimming)
  - *pairName.HBX2coordinates*.clades (source/recipient identity)
  - *pairName.HBX2coordinates*.states (binary source/recipient identity)
  - *pairName.HBX2coordinates*.fasta (multiple sequence alignment of ingroup/outgroup in fasta format)
  - *pairName.HBX2coordinates*.nex (multiple sequence alignment of ingroup/outgroup in nexus format)

**Requirements:** 

- Transmission cluster data stored as by *pair.epi.retrieval*.
- [MrBayes](https://nbisweden.github.io/MrBayes/download.html) as mb in your local bin directory. 
- R dependencies: *ape*, *Biostrings*, *R.utils*

**Arguments:**

- ids: character vector with LANLdb transmission cluster names.
- max.n: *Default=300*. Ignore genetic region if the number of sequences is higher than this value.
- ref.gap.tolerance: *Default=99*. Maximum value for an extended gap. Extended gaps higher than this value are trimmed from the alignment.
- gap.threshold: *Default=1.0*. the largest gap proportion to delete a column from the alignment
- sets: *Default="set.csv"*. A file name containing the genomic region(s) to be analysed. 

**Usage example:**

`pair.alignment(ids=ids, max.n=300, ref.gap.tolerance=99, gap.threshold=1, sets="set.csv")` <br/><br/>

## *pair.mb.setting* <br/><br/>

Sets files for MrBayes

**Outputs:**

- mrbayes.nex (instructions to run mr.bayes)
- *pairName.HBX2coordinates*.mb.nex (sequences in nexus file with source/recipient identity for ancestral state reconstruction)
- Genetic regions input file (e.g. either set.csv or sets.csv) is updated with the number of programmed MCMC generations. 
- If you run *mrbayes.nex* (e.g. `mb mrbayes.nex`) you will obtain the following mrbayes output files:
  - *pairName.HBX2coordinates*.mb.nex.mcmc (mcmc convergence diagnostics)
  - *pairName.HBX2coordinates*.mb.nex.run*#*.p (sampled values of the parameters over the runs)
  - *pairName.HBX2coordinates*.mb.nex.pstat (parameter summaries)
  - *pairName.HBX2coordinates*.mb.nex.run*#*.t (sampled trees over the runs)
  - *pairName.HBX2coordinates*.mb.nex.ckp (checkpoint files)
  - *pairName.HBX2coordinates*.mb.nex.lstat (estimated marginal likelihoods for runs)

**Requirements:**


- Transmission cluster data stored as by *pair.epi.retrieval* function.
- R Dependencies: ape


**Arguments:**

- query: query file
- ngen: *Default=10000000*.
- samplefreq: *Default=ngen x 0.0001*.
- printfreq: *Default=ngen x 0.01*.
- diagnfreq: *Default=ngen x 0.1*.
- burninfrac: *Default=0.5*.
- sets: *Default="set.csv"*. A file name containing the genomic region(s) to be analysed. 

**Usage example:**

`pair.mb.setting(ids=ids, ngen=10000000, samplefreq=1000, printfreq=100000, diagnfreq=1000000, burninfrac=0.5, sets="set.csv")`<br/><br/>

## *pair.mb.summary* <br/><br/> 

Process log and tree files.

**Outputs:**

- Genetic regions input file (e.g. either set.csv or sets.csv) is updated with minimun ESS, standard deviation of split frequencies, frequency summaries for tree metrics, the frequency and proportion of every topology class, 95% confidence level evaluation for topology class, frequency summaries for the ancestral state reconstructions and 90% confidence level evaluation for ancestral state reconstructions.
- *pairName.HBX2coordinates*.COMB.log (combined log files)
- *pairName.HBX2coordinates*.ANN.log (combined log files with calculations of *topology.class* function)
- *pairName.HBX2coordinates*.COMB.trees (combined tree files in nexus format)
- *pairName.HBX2coordinates*.newick.trees (combined tree files in newick format)

**Requirements:**

- Transmission cluster data stored as by pair.retrieval function.
- MrBayes output as specified by *pair.mb.setting* function.
- R Dependencies: ape, mcmcse.
- *topology.class* function.

**Arguments:**

- ids: character vector with LANLdb transmission cluster names.
- sets: *Default="set.csv"*. A file name containing the genomic region(s) to be analysed. 

**Usage example:**

`pair.mb.summary(ids=ids, sets="set.csv")`<br/><br/>


