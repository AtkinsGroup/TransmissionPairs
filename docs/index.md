# Code and data for reproducing: 
  Villabona-Arenas *et al.* Number of HIV-1 founder variants is determined by the recency of the source partner infection

The functions are separated into three directories:

1. Data Collation
2. Tree Topology Classification
3. Phylodynamic Analysis

**1. Data Collation.**<br/><br/>

Here we list two functions:

- *pair.epi.retrieval*, which use the python script alamos-extract to help retrieve the metadata from queried HIV transmission pairs that are indexed in Los Alamos Database.
- *pair.sequence.retrieval* helps, in a downstrain way, to retrieve the genetic sequences from queried HIV transmission pairs that are indexed in Los Alamos Database.

**2. Tree Topology Classification.** <br/><br/>

Here we list the function ```topology.class``` which determines the tree topology class (i.e. MM, PM or PP). ```topology.class``` is sourced during the phylodynamic analysis

**3. Phylodynamic Analysis.** <br/><br/>

Here, we separate the functions in two parts:

- Empirical analysis: Here we list the functions to (a) align data, (b) run MrBayes and (c) summarize run results. 
- Simulation analysis: Here we list the function to (a) simulate trees using the metadata from the transmission pairs and (b) to compute the likelihoods by combining both the empirical and simulation outcomes.
