
## This repository supplies the code and data to generate the results in Villabona-Arenas *et al.* (2020) Number of HIV-1 founder variants is determined by the recency of the source partner infection

## To use the code please cite the following:
  - The DOI of the code version: XXX
  - The original publication: XXX 
  

### Repository overview


## Running information

As per figure 1B from the manuscript, we have our methods structured in four stages: Data collation, Data pre-processing, Phylogenetic analysis and Model Calibration.

-	Data collation: here we retrieve epidemiological data from los LANL HIV and GenBank databases, as well as genetic sequences from the latter, using HIV transmission cluster (and partners) ids as queries. In order to help retrieving this, we write two functions which are saved to the directory Data_Collation. The function *pair.epi.retrieval* queries the HIV transmission pair input file, collates the epidemiological data and stores it in individual directories. The function *pair.sequence.retrieval* allows the user to retrieve the genetic sequences from GenBank and store them in the corresponding directories.
-	Data pre-processing: here we either build sequence alignments or build the epidemiological timelines. For the alignment, In order to run the alignments we have written the function *pair.alignment.R* and saved it to Phylodynamic Analysis/empirical_analysis. In order to build the timelines, most of the epidemiological data is organized in sub-sets by *pair.epi.retrieval* and *pair.sets.filter* (Analysis/empirical_analysis) allows the user to filter (if possible) the sequences based on collection dates. The function *pair.simulator* in Phylodynamic Analysis/simulation_analysis builds the timeline for any input subset before running the simulation.
-	Phylogenetic analysis: here we infer phylogenetic trees using either empirical or simulated genetic data from the HIV transmission pairs. For the empirical analysis we wrote the functions *pair.mb.setting.R* and *pair.mb.summary.R* and saved them to Phylodynamic Analysis/empirical_analysis. These functions set the mrbayes running files and summarize the empirical outputs, respectively. For the simulations, the same function that build the timelines, *pair.simulator*, produces an output with the phylogenetic findings.
-	Model Calibration. Here we combine the empirical and simulated outcomes and calculate the probability of one founder strain. We wrote a function for this purpose which is saved to PhylodynamicAnalysis/Simulation_analysis as *pair.founder.p.R*


## Directory organisation

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
