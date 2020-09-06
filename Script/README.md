# Helper Scripts for Using Database NITREc

Yu Wan

Release: 1 August 2020; latest update: 5 September 2020

<br/>

Nine scripts are developed to facilitate use of database NITREcMut.



## Functional classification of scripts

**Gene screen**

- `screenGenes.pbs` and its configuration profile `screenGenes.config`
- `compileBLAST.py` compiles outputs of `screenGenes.pbs` or users' BLASTn commands

**Sequence manipulation**

- `extractRegionFromContig.py`
- `rmSeqAnnot.py`
- `splitMultiFASTA.py`

**Protein-level mutation identification**

- `translateDNA.py`: Translate coding sequences into protein sequences, given a codon table.
- `missenseFinder.py`: Identify amino acid substitutions in query protein sequences against their most similar reference protein sequences, assuming an identical length of all sequences.

- `findKnownMutations.py`: Identify mutations of interest in the output of `missenseFinder.py`.