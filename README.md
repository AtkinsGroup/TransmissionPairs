# HIV Transmission Pairs

This is the code and data repository for the HIV Transmission Pairs project.

This README includes:
- How to cite the code and paper
- An overview of repository structure
- Installation instructions
- Example code for running analysis


## Please remember to cite the following:
- The original publication: Villabona-Arenas CJ, Hall M, Lythgoe KA, Gaffney SG, Regoes RR, Hué S, Atkins KE (2020) "**Number of HIV-1 founder variants is determined by the recency of the source partner infection**", *Science*
- And, if citing the code or data specifically, the DOI as specified in the Science publication
- Please refer to CITATION file for up-to-date references

## Repository overview and organisation
Our methods are structured in four stages: Data collation, Data pre-processing, Phylogenetic analysis and Model Calibration (please refer to Fig. 1B in the accompanying publication). These methods are organised into the following directories:
1. Data Collation
2. Phylodynamic Analysis (data pre-processing, phylogenetic analysis and model calibration)
3. Tree Topology Classification (inferring phylogenetic tree topology class)


We also provide all the data used to plot the figures, and some minimal example data for running the analysis from scratch:

4. Data
5. SH_ZM221 (an example HIV-1 transmission pair used in the worked example below)

More detailed documentation is provided in the README files within the sub-directories. 


### 1. Data Collation.
First we retrieve epidemiological data from LANL HIV and GenBank databases, as well as genetic sequences from the latter, using HIV transmission cluster / partner IDs. For this, we provide two R functions that can be called:
`pair.epi.retrieval` which uses the python script `alamos-extract` that retrieves the metadata from HIV transmission pairs that are indexed in Los Alamos Database. 

`pair.sequence.retrieval` which retrieves the genetic sequences indexed in the Los Alamos Database from GenBank and creates subdirectories for each pair.

 
### 2. Phylodynamic Analysis.
*Empirical analysis* consists of (a) aligning the sequence data collated, (b) writing the  MrBayes running files that will generate the phylogenetic trees of the empirical data and (c) summarizing MrBayes outputs.
The function `pair.sets.filter ` selects the sequences to be used in the analysis. The function `pair.alignment` then aligns the sequence data.
`pair.mb.setting` generates a .nex file to run MrBayes using the sequence data.
`pair.mb.summary` generates a summary file of the phylogenetic analysis completed in MrBayes. 
 
*Simulation analysis* consists of (a) simulating phylogenetic trees for a transmission pair using the epidemiological timelines for that pair for a range of transmission models, (b) computing the likelihood associated with each transmission model by combining both the empirical and simulation outcomes for a transmission pair (Model calibration), (c) calculating the probability of one founder strain for each transmission pair’s maximum likelihood model, and (d) calculating the relative risk that one founder strain initiates infection in acute pair transmissions than in chronic pair transmissions (across all pairs)
 
pair.simulator generates the transmission timelines for a transmission pair based on their epidemiological data, simulates sequences using VirusTreeSimulator (xxxx) and for each simulation, calculates a maximum likelihood tree for the pair; using these ML trees, we then calculate the number of topology classifications (i.e. PP, PM, MM) for each pair across all the simulations.  
 
`pair.founder.p.R` uses maximum likelihood to calculate the best fit model for each transmission pair and then, for this model, calculates the probability that one variant was transmitted. Finally, we calculate the relative risk of acute transmissions vs. chronic transmissions initiating infections with one founder strain across pairs.
 
 
### 3. Tree Topology Classification.
`topology.class` automatically determines the tree topology class (i.e. MM, PM or PP) and is called during the Phylodynamic analysis. 


### 4. Data.
This directory contains:
S1_SITable_EpiGeneticData.csv: output from the Data Collation step (both automated and manual retrieval) detailing the 112 transmission pairs used in the analysis (their epidemiological data and all metadata, publication details and ethics information)
S2_SITable_AnalysisData.csv: output from the Tree Topology Classification and Phylodynamics analysis steps detailing the results used to generate the publication figures
S3_SITable_ColumnNamesKey.csv: The column headings for data tables S1 and S2 together with the publication figures numbers where these columns are plotted
S4_Alignments.zip: aligned sequence data for the 112 transmission pairs used in the analysis as .fasta files with files names as the "LANLdb_cluster_name" listed in .csv files


## Installation

There are two options for installation.

1. **Manual installation**: The included R scripts use CRAN and Bioconductor packages listed in `renv.lock`, that can be installed using the `renv` package. Additional tools required are:
- `mrbayes`: https://nbisweden.github.io/MrBayes/index.html
- `VirusTreeSimulator`: https://github.com/PangeaHIV/VirusTreeSimulator
- `Seq-Gen`: https://github.com/rambaut/Seq-Gen/releases/tag/1.3.4
- `IQ-TREE`: http://www.iqtree.org/

2. **Docker**: All code, data, dependencies, and RStudio are included in a Docker image available on Docker Hub (https://hub.docker.com/r/atkinsgroup/transmissionpairs). With Docker [installed](https://docs.docker.com/get-docker/), you can launch RStudio in the TransmissionPairs container with a single line:
```bash
docker run -it -e PASSWORD=YOURPASSWORD -p 8787:8787 atkinsgroup/transmissionpairs
```
- Modify the above line to specify a password, changing 'YOURPASSWORD' to a new value.
- RStudio will be available in your browser at `localhost:8787`, with username `rstudio` and password as above.
- Repository code is in the `/transmissionpairs` directory.
- To make directories on your host machine available to the container, e.g. to save output, use the `-v` argument to bind a directory, e.g. `-v /some/local/path:/transmissionpairs/out`
- To build and run the container using `Singularity` instead of Docker, run the following:
```bash
# Build the Singularity image (you only need to do this once)
singularity build transmissionpairs.sif docker://atkinsgroup/transmissionpairs
# Launch the container
PASSWORD=YOURPASSWORD singularity exec transmissionpairs.sif rserver --auth-none=0  --auth-pam-helper-path=pam-helper
```
  - as with the Docker installation, modify the command with suitable password.
  - RStudio will be available on port 8787 with your username and the chosen password.
  - For more information on setup, refer to the Rocker instructions for Singularity, as the container is built on a Rocker container: https://www.rocker-project.org/use/singularity/


## Running the analysis: example code
 
```R
setwd(<filepath to TransmissionPairs>)
#request a key in https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/, query.txt example is specified 
source("DataCollation/LANLRetrieval/pair.epi.retrieval.R")
ids <- pair.epi.retrieval(query = "DataCollation/LANLRetrieval/example_query.csv", key = "myAPIkeyNumber")

#write .fasta file with sequences of ‘ids’ to be used in the analysis as sets.csv
source("DataCollation/GenBankRetrieval/pair.sequence.retrieval.R")
pair.sequence.retrieval(ids=ids, min.overlap=70, max.length.diff=10, min.seqs=5)

#update sets.csv to take only the sequences closest to the recipient infection time
source("DataCollation/GenBankRetrieval/pair.sets.filter.R")
pair.sets.filter(ids=ids, min.seqs=5)

# align the sequence to be used in the analysis
source("PhylodynamicAnalysis/Empirical_analysis/pair.alignment.R")
pair.alignment(ids=ids, max.n=300, ref.gap.tolerance=99, gap.threshold=1, sets="set.csv")

# read in the source files for the topology analysis
source("TreeTopologyClassification/topology.class.R")

# generate mrbayes.nex for MrBayes
source("PhylodynamicAnalysis/Empirical_analysis/pair.mb.setting.R")
pair.mb.setting(ids=ids, ngen=10000000, samplefreq=1000, printfreq=100000, diagnfreq=1000000, burninfrac=0.5, sets="set.csv")

# Now run MrBayes via R using:
system("mb SH_ZM221/SH_ZM221.6639.7478/mrbayes.nex")


# Note that example results are attached as a .zip 
# Unzip MrBayes results via R
unzip("SH_ZM221/SH_ZM221.6639.7478/SH_ZM221.6639.7478.zip", exdir = "SH_ZM221/SH_ZM221.6639.7478/")

# Generate phylogenetic findings for empirical tree and overwrite set.csv with results
source("PhylodynamicAnalysis/Empirical_analysis/pair.mb.summary.R")
pair.mb.summary(ids=ids, sets="set.csv")

# Simulates trees and summarises results for a specific epidemiological timeline
pair.simulator(pair.id="mypair", vts="/VirusTreeSimulator.jar", simulation=1, transmission=30, sampling.source=30, sampling.rec=30, index=1095, n.source=5, n.rec=5, vt.source=1, vt.rec=1, tau=1.8, Ne0=1, K=300, t50=-2, evo.r=0.01148, gene.ref=c(6758, 7757), subs.model="GTR", subs.param=c(0.48,4,0.18,0.41,0.17,0.20,0.22,3.01,5.59,1.21,1.39,8.25,1.0)S, acute=c(13:90), recent=c(91:180), chronic=c(181:7300))

# Calibrate the simulations to the empirical data (and save the best model to a file. Uses results from all the pairs, example simulation files are provided
source("PhylodynamicAnalysis/Simulation_analysis/pair.founder.p.R")
unzip("PhylodynamicAnalysis/Simulation_analysis/example_simulations.csv.zip", exdir="PhylodynamicAnalysis/Simulation_analysis/")
pair.founder.p(simulations="PhylodynamicAnalysis/Simulation_analysis/example_simulations.csv", empirical="PhylodynamicAnalysis/Empirical_analysis/example_empirical.csv")
```
