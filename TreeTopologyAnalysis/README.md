## *topology.class* <br/><br/>

Reads a transmission pair tree or a set of trees and computes the topology class (monophyletic-monophyletic [MM], paraphyletic-monophyletic [PM] or paraphyletic-polyphyletic [PP]) and phylogenetic metrics. 

**Outputs:**

- A list of vectors:
  - vector.lineages.source (Number of Monophyletic groups for the source partner)
  - vector.lineages.recipient (Number of Monophyletic groups for the recipient partner)
  - vector.topology (Topology class: MM, PM or PP)
  - vector.p.ancestral.pars (Ancentral state reconstruction using parsimony)
  - vector.p.ancestral.ml.ER (Ancentral state reconstruction using an equal-rates model and Maximum Likelihood [ML])
  - vector.p.ancestral.bayes (Ancentral state reconstruction using the full hierarchical Bayesian approach [BI])
  - vector.pars.direction (Inferred transmission history using parsimnoy: consistent, inconsistent or equivocal)
  - vector.ml.ER.direction (Inferred transmission history using ML: consistent, inconsistent or equivocal)
  - vector.bayes.direction (Inferred transmission history using BI: consistent, inconsistent or equivocal)
  - vector.tree.length (tree length)
  - vector.tree.SBL (phylogenetic metric: sum of branch lengths)
  - vector.tree.SPD (phylogenetic metric: sum of phylogenetic distances)
  - vector.tree.MBL (phylogenetic metric: mean branch length)
  - vector.tree.MPD (phylogenetic metric: mean pairwise distance)
  - vector.tree.VBL (phylogenetic metric: variation of of branch lengths)
  - vector.tree.VPD (phylogenetic metric: variation of pairwise distances)

**Requirements:**

- a data frame with tips and character states. The first column correspond to the tips labels and have to match the tip labels of the trees. The second column informs the states in the transmission pair, either 0 (source) or 1 (recipient). If there is an outgroup (one tip only) it should be listed in the first row (state=?).
- (optional) a log file from a bayesian analyses using MrBayes that included Ancestral State Reconstruction. Column 19 have to correspond with the probability of the source being the ancestral state for the ingroup.
- R Dependencies: adephylo, ape, dplyr, phangorn, phytools, tibble 

**Arguments:**

- tr: phylogenetic tree(s) in (multi)"phylo" format.
- states: a data frame with tips and character states. 
- outgroup: logical value indicating whether there is an outgroup. This should be the first enumerated tip in the object tr and the first row in the element states.
- p.ancestral.logs: logical value indicating whether there was an Ancentral state reconstruction using the full hierarchical Bayesian approach
- logs: if p.ancestral.logs=TRUE then the name of the logs object to append the information.

**Usage example:**

`topology.class(tr=trees, states=states, outgroup=TRUE, p.ancestral.logs=TRUE, logs=logs)` 

