## *pair.sequence.retrieval* <br/><br/>

Uses info.gb and info.gb.csv to subset data from transmission pairs into genomic regions based on HXB2 coordinates. Subtype specific reference sequences are informed by Los Alamos National Laboratory HIV sequence database (LANLdb). Viral genetic sequences are retrieved from GenBank. 

Evaluation starts with the shortest sequence and a set (genetic region) is defined using its HXB2 coordinates. All sequences are evaluated and overlapping sequences (*Argument min.overlap*) are added to the set. Sequences may be dropped from further evaluation based on length (*Argument max.length.diff*). The set is returned if it meets the minimum number of sequences requirement (*Argument min.seqs*). Evaluation continues until there are less than (*Argument min.seqs*) sequences to be evaluated. 

**Outputs:**

- sets.csv (viral genomic regions with genbank accessions and epidemiological information).
- ingroup.fasta (viral genetic sequences in fasta format).
- outgroup.fasta (subtype specific reference viral genetic sequences from LANLdb).

**Requirements:**

- Transmission cluster data stored as by *pair.epi.retrieval*.
- R dependencies: *genbankr*, *rentrez*, *stringr*.

**Arguments:**

- ids: character vector with LANLdb transmission cluster names.
- min.overlap: *Default=70*. Sequences with a minumum overlap of [70%] are considered to be the same region.
- max.length.diff: *Default=10*. Overlapping sequence that are shorter than then current genetic - region length + [10%] are dropped from the analyses, i.e. these sequences are not long enough to cover any other genetic region.
- Argument min.seqs: *Default=5*. Minimum number of sequences to delineate a set.

**Usage example:**

`pair.sequence.retrieval(ids=ids, min.overlap=70, max.length.diff=10, min.seqs=5)`<br/><br/>

## *pair.sets.filter* <br/><br/>

[Optional] Subset the longest genomic region. If there is any indication of multiple time points per individual, *attempts* to reduce the number of genetic sequences to one single time point per individual (giving preferent to the closest to the time of infection) using the epidemiological data from Los Alamos National Laboratory HIV sequence database (LANLdb) and GenBank.

This is an *optional automated filter* and it is constrained by the availability of information associated to every genetic sequence. If you opt out from this step, manually subset the line from sets.csv (genomic region) that interest you the most and saved it as set.csv.

**Outputs:**

- set.csv (longest viral genomic region with genbank accessions and epidemiological information).

**Requirements:**

- Transmission cluster data stored as by *pair.epi.retrieval*.
- R dependencies: *stringr*.

**Arguments:**

- ids: character vector with LANLdb transmission cluster names.
- Argument min.seqs: *Default=5*. Minimum number of sequences to delineate the set.

**Usage example:**

`pair.sets.filter(ids=ids, min.seqs=5)`<br/><br/>