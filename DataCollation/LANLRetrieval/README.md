## *pair.epi.retrieval* <br/><br/>

Retrieve epidemiological tables from transmission clusters (by pairs) stored in Los Alamos National Laboratory HIV sequence database (LANLdb) and retrieve calendar dates from GenBank Accessions.

**Outputs:**

- Directory with LANLdb transmission cluster names (if name is too long (>20 characters) clusterID is used instead).
- info.csv (summary of epidemiological info).
- info.gb.csv (epidemiological info by GenBank Accession).
- ids.txt and object "ids" (lists all queried transmission clusters using the transmission name).

**Requirements:**

- [Entrez API key](https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/)
- [alamos-extract](https://github.com/AtkinsGroup/alamos-extract) installed in your local Python distribution
- R Dependencies: *genbankr*, *rentrez*, *stringr*.

**Arguments:**

- query: query file
- key: Entrez API key

**Usage example:**

`pair.epi.retrieval("query.csv", "myAPIkeyNumber")` 

**Query file specs:**

Three identifiers are required per row in a csv file (without header): LANLdbClusterID, LANLdbPatienID(Source), LANLdbPatienID(Recipient)\
Example:<br/>
`757,31107,31093`\
`238,4534,4535` <br/><br/>
