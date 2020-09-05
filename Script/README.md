# Helper Scripts for Using Database NITREcMut

_Yu Wan_

_1 August 2020_



Eight scripts are developed to facilitate use of database NITREcMut.



## Functional classification of scripts

**Gene screen**

- `screenGenes.pbs` and its configuration profile `screenGenes.config`
- `compileBLAST.py` compiles outputs of `screenGenes.pbs` or users' BLASTn commands

**Sequence maniulation**

- `extractRegionFromContig.py`
- `rmSeqAnnot.py`
- `splitMultiFASTA.py`

**Protein-level mutation identification**

- `translateDNA.py`
- `missenseFinder.py`