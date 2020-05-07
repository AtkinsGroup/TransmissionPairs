
## This repository supplies the code and data to generate the results in Villabona-Arenas *et al.* (2020) Number of HIV-1 founder variants is determined by the recency of the source partner infection

## To use the code please cite the following:
  - The DOI of the code version: XXX
  - The original publication: XXX 
  

## Repository overview and organisation

Our methods are structured in four stages: Data collation, Data pre-processing, Phylogenetic analysis and Model Calibration (please refer to Fig. 1B in the accompanying publication). These methods are organised into the following directories:

1. Data Collation (containing Data collation and Data pre-processing steps)
2. Tree Topology Classification (containing phylogenetic analysis steps)
3. Phylodynamic Analysis (containing model calibration steps)
4. Data (containing output files from both Data Collation and Tree Topology Classification / Phylodynamic analysis)


**1. Data Collation.**<br/><br/>

First we retrieve epidemiological data from los LANL HIV and GenBank databases, as well as genetic sequences from the latter, using HIV transmission cluster / partner IDs. For this, we provide two R functions that can be called:

```pair.epi.retrieval``` which use the python script alamos-extract to help retrieve the metadata from queried HIV transmission pairs that are indexed in Los Alamos Database

```pair.sequence.retrieval``` allows the user to retrieve the genetic sequences indexed in the Los Alamos Datacase from GenBank and store them in the corresponding directories.  

Next, we pre-process the data by building sequence alignments and by building the epidemiological timelines used in the analysis. 

To generate the alignments, we provide an R function:

```pair.alignment``` and saved it to Phylodynamic Analysis/empirical_analysis. -- KEA:WHAT DOES THIS MEAN? 

For the timelines, 

most of the epidemiological data is organized in sub-sets by *pair.epi.retrieval* and *pair.sets.filter* (Analysis/empirical_analysis) allows the user to filter (if possible) the sequences based on collection dates. --KEA:WHAT DOES THIS MEAN? 

The function *pair.simulator* in Phylodynamic Analysis/simulation_analysis builds the timeline for any input subset before running the simulation. -KEA:WHAT DOES THIS MEAN? 

**2. Tree Topology Classification.** <br/><br/>

Here we list the function ```topology.class``` which determines the tree topology class (i.e. MM, PM or PP). ```topology.class``` is sourced during the phylodynamic analysis

**3. Phylodynamic Analysis.** <br/><br/>

Here, we separate the functions in two parts:

- Empirical analysis: Here we list the functions to (a) align data, (b) run MrBayes and (c) summarize run results. 
- Simulation analysis: Here we list the function to (a) simulate trees using the metadata from the transmission pairs and (b) to compute the likelihoods by combining both the empirical and simulation outcomes.

**4. Data.**

This directory contains: 
- *S1_SITable_EpiGeneticData.csv*: output from the Data Collation step (both automated and manual retrieval) detailing the 112 transmission pairs used in the analysis (their epidemiological data and all metadata, publication details and ethics information) 
- *S2_SITable_AnalysisData.csv*: output from the Tree Topology Classification and Phylodynamics analysis steps detailing the results used to generate the publication figures
- *S3_SITable_ColumnNamesKey.csv*: The column headings for data tables S1 and S2 together with the publication figures numbers where these columns are plotted
- *S4_Alignments.zip*:  aligned sequence data for the 112 transmission pairs used in the analysis as .fasta files with files names as the "LANLdb_cluster_name" listed in .csv files

### Running information



-	Phylogenetic analysis: here we infer phylogenetic trees using either empirical or simulated genetic data from the HIV transmission pairs. For the empirical analysis we wrote the functions *pair.mb.setting.R* and *pair.mb.summary.R* and saved them to Phylodynamic Analysis/empirical_analysis. These functions set the mrbayes running files and summarize the empirical outputs, respectively. For the simulations, the same function that build the timelines, *pair.simulator*, produces an output with the phylogenetic findings.
-	Model Calibration. Here we combine the empirical and simulated outcomes and calculate the probability of one founder strain. We wrote a function for this purpose which is saved to PhylodynamicAnalysis/Simulation_analysis as *pair.founder.p.R*


