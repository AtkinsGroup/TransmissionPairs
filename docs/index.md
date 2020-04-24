# Code and data for reproducing: 
  Villabona-Arenas *et al.* Number of HIV-1 founder variants is determined by the recency of the source partner infection

The functions are separated in three directories:

1. Data Collation
2. Tree Topology Classification
3. Phylodynamic Analysis

**1. Data Collation.**<br/><br/>

Here we list two functions:

- *pair.epi.retrieval*, which use the python script alamos-extract to help retrieve the metadata from queried HIV transmission pairs that are indexed in Los Alamos Database.
- *pair.sequence.retrieval* helps, in a downstrain way, to retrieve the genetic sequences from queried HIV transmission pairs that are indexed in Los Alamos Database.

**2. Tree Topology Classification.** <br/><br/>

Here we list the function *topology.class* which is sourced during the phylodynamic analysis and allow to use determine the tree topology class.

**3. Phylodynamic Analysis.** <br/><br/>

Here, we separate the functions in two parts:

- Empirical analysis: Here we list the functions to align data, run mrbayes and summarize run results. 
- Simulation analsysis: Here we list the function to simulate trees using the metadata from the transmission pairs and to compute probabilities by combining both, the empirical and simulation outcomes.
