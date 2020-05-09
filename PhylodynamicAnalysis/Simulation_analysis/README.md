## *pair.simulator* <br/><br/>

Simulates phylogenetic trees using epidemiological information from a transmission pair.

**Outputs:**

- results.sim.csv (simulation information)
- sequences.sim.csv (simulated viral sequences)
- trees.sim.tre (simulated transmission tree)
- trees.tre (reconstructed trees)

**Requirements:**

- [VirusTreeSimulator](https://github.com/PangeaHIV/VirusTreeSimulator).
- [Seq-Geb](https://github.com/rambaut/Seq-Gen/releases/tag/1.3.4) as seq-gen in your local bin directory. 
- [IQ-TREE](http://www.iqtree.org/) as iqtree in your local bin directory. 
- *topology.class.R* function.
- R Dependencies: adegenet, ape, pegas, phylotools, seqRFLP

**Arguments:**

- pair.id: *Default=unknown*. A character string as identifier for the simulations.
- vts: full path to VirusTreeSimulator.jar
- simulation: *Default=1*. An integer  as identifier for the simulation.
- transmission: *Default=30*. An integer informing days from infection of the transmitter at transmission OR a character string indicating a range c("30 60"), a truncated range c(">365") or full uncertainty ("acute", "recent" or "chronic").
- sampling.source: *Default=30*. An integer informing sampling of the transmitter (in days) after transmission OR a character string indicating a range c("30 60").
- sampling.rec: *Default=30*. An integer  informing sampling of the recipient (in days) after infection OR a character string indicating a range c("30 60").
- index: *Default=3*. An integer  informing days from infection of the index case at transmission to the transmitter partner from the transmission pair.
- n.source: *Default=5*. An integer informing the number of sequences sampled from the transmitter.
- n.rec: *Default=5*. An integer informing the number of sequences sampled from the recipient.
- vt.source: *Default=1*. An integer informing the number of founding particles in the transmitter.
- vt.rec: *Default=1*. An integer informing the number of founding particles in the recipient.
- tau: *Default=1.8*. viral generation time in days for the intra-host logistic population model.
- Ne0: *Default=1*. Initial effective population size for the intra-host logistic population model.
- K: *Default=300*. Carrying capacity for the intra-host logistic population model.
- t50: *Default=-2*. Backwards time to half the carrying capacity of the viral population.
- evo.r: *Default=0.01148*. Viral substitution rate for the simulated sequences
- gene.ref: *Default=c(6758, 7757)*. An integer vector with two HBX2 coordinates. Specify the ancestral sequence for Seq-Gen. 
- subs.model: *Default='GTR'*. The model of substitution for Seq-Gen (either GTR or HKY).
- subs.param: *Default=c(0.48,4,0.18,0.41,0.17,0.20,0.22,3.01,5.59,1.21,1.39,8.25,1.0)*. A numeric vector with parameters for Seq-Gen. For GTR: alpha, categories, proportion of invariants, A frequency, C frequency, G frequency, substitution rates. For HKY, only two substitution rates are needed (if more are supplied, only the first two are used).
- HBX2: object with HBX2 as a matrix DNAbin. 
- acute: *Default=c(13:90)* An integer vector indicating the uncertainty around the transmission time for a chronic transmitter.
- recent: *Default=c(91:180)* An integer vector indicating the uncertainty around the transmission time for a recent transmitter.
- chronic: *Default=c(181:7300)* An integer vector indicating the uncertainty around the transmission time for an acute transmitter.

**Usage example:**

`source('topology.class.R')` 

`HBX2<- as.matrix.DNAbin(read.GenBank('K03455'))`

`pair.simulator(pair.id="mypair", vts="/VirusTreeSimulator.jar", simulation=1, transmission=30, sampling.source=30, sampling.rec=30, index=1095,` `n.source=5, n.rec=5, vt.source=1, vt.rec=1, tau=1.8, Ne0=1, K=300, t50=-2, evo.r=0.01148, gene.ref=c(6758, 7757), subs.model="GTR",` `subs.param=c(0.48,4,0.18,0.41,0.17,0.20,0.22,3.01,5.59,1.21,1.39,8.25,1.0), HBX2=HBX2, acute=c(13:90), recent=c(91:180), chronic=c(181:7300))`<br/><br/>


## *pair.founder.p* <br/><br/>

Calculates probability of one founder strain by combining information from the empirical and simulated data.

**Outputs:**

- bestmodel.csv (includes likelihood calculations for partcile models, the best-fit model selection, the probability of one founder strain and mode of founder variants)

**Requirements:**

- R Dependencies: adegenet, ape, pegas, phylotools, seqRFLP

**Arguments:**

- simulations: combined output of multiple simulated pairs using *pair.simulator* function.
- empirical: combined output of multiple empirical pairs using *pair.mb.summary* function (e.g. merged set.csv files).


**Usage example:**

`pair.founder.p(sims=simulations, empirical=empirical)`

