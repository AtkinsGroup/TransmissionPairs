alamos-extract
==============


Extract cluster and patient info from Los Alamos HIV Database.


## Installation


Install in develop mode, directly from GitHub:

```bash
pip install -e git+https://github.com/sggaffney/alamos-extract.git#egg=alamos-extract
```

This should put the command line tool 'load_hiv' in your path.

## Example

```bash
load_hiv cluster 684
```

output:
```
Cluster: SM_cluster_ABC
Description: Source patient and 3 partners (A, B, C) with MSM transmission.
6 patients: 0558, 0559, RecA, RecB, RecC, SourceABC
117 accessions.
Clinical data written to cluster_684_clinical.tsv
Accession data written to cluster_684_accessions.tsv
```
