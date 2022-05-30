# Helper Scripts for Using Database NITREc

Yu Wan

Release: 1 Aug 2020; latest update: 29 May 2022

<br/>

This directory comprises scripts developed to facilitate use of database NITREc. The scripts can also be used for other analyses.

**Table of contents**

1. [Installation](#sec1)
2. [Functional classification of scripts](#sec2)
3. [Gene detection based on nucleotide sequences](#sec3)
    - [Step 1.  Create a FASTA file of query sequences and a list of subject genome names](#step1)
    - [Step 2. Configure the detection job](#step2)
    - [Step 3. Use `searchGenes.pbs` to detect genes from FASTA files](#step3)
    - [Step 4. Compile output TSV files into a table](#step4)

<br/>

## 1. <a name = "sec1">Installation</a>

```bash
git clone https://github.com/wanyuac/NITREc.git
```
<br/>

## <a name = "sec2">2. Functional classification of scripts</a>

**Gene detection (DNA)**

- `searchGenes.pbs` and its configuration profile `searchGenes.config`
- `compileBLAST.py` compiles outputs of `screenGenes.pbs` or users' BLASTn commands

**Sequence manipulation**

- `extractRegionFromContig.py`
- `rmSeqAnnot.py`
- `filterMultiFASTA.py`
- `rmProteinsByLength.py`

**Protein-level mutation identification**

- `translateDNA.py`: Translate coding sequences into protein sequences, given a codon table.
- `missenseFinder.py`: Identify amino acid substitutions in query protein sequences against their most similar reference protein sequences, assuming an identical length of all sequences.

- `findKnownMutations.py`: Identify mutations of interest in the output of `missenseFinder.py`.

**Prediction of nitrofurantoin susceptibility**

- `scoreHitsNITR.R`: An R script offering a score function for prediction based on genetic presence-absence and alterations.

<br/>

## <a name = "sec3">3. Gene detection based on nucleotide sequences</a>

###  <a name = "step1">Step 1.  Create a FASTA file of query sequences and a list of subject genome names</a>

A query is the DNA sequence to be searched against a genome assembly (the subject). FASTA files in `NITREc/Seq/Nucl/` are examples of query sequences. The FASTA file of query sequences is an essential input for script `searchGenes.pbs`.

Each subject genome (for example, `isolate_01`) is represented by a FASTA file with the genome name in the base of its filename (for example, `isolate_01.fna`). An example of the list of subject genome names (say, `subject_genomes.txt`):

```text
isolate_01
isolate_02
isolate_03
isolate_04
isolate_05
```

Note that `subjects_genomes.txt` must be a plain-text file where each line is ended by the Unix/Linux style newline characters `\n`. Otherwise, isolate names may not be read correctly by `searchGenes.pbs` and cause problems to its output. For Windows/Mac OS users, the newline character can be specified in a text editor [with the utility such as "Convert > Line Break Style > To UNIX (LF only)"] or with the `dos2unix` program on Linux:

```bash
dos2unix subject_genomes.txt
```



### <a name = "step2">Step 2. Configure the detection job</a>

Copy `searchGenes.config` to your working directory (for instance, `~/work) and configure its variable statement in a plain-text editor. Assuming the `NITREc` repository has been downloaded to `~/bin`:

```bash
cp ~/bin/NITREc/Script/searchGenes.config ~/work/searchGenes_demo.config  # It's not mandatory to rename the configuration file.
```

The configuration file `searchGenes.config` is a mandatory bash script to be executed by `searchGenes.pbs`. Therefore, any changes to this file should observe the bash syntax. An example of the populated configuration file:

```bash
# Configuration of parameters for script searchGenes.pbs
# Copyright (C) 2020-2021 Yu Wan <wanyuac@126.com>
# Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
# Creation: 28 August 2020; the latest modification: 16 September 2021

# Environmental constants (optional)
# Leave MODULE and CONDA_ENV empty if programs makeblastdb and blastn are accessible through environment variable $PATH.
MODULE='anaconda3/personal'  # Name of module Conda in your system
CONDA_ENV='python3'  # Name of Conda environment to be loaded
CONDA_ACTIVATE_METHOD='source'  # Or 'conda', depending on whether 'source activate' or 'conda activate' is used by the system
BLASTN_DIR=""  # The directory where program blastn is stored (not needed if both MODULE and CONDA_ENV are specified)

# Input data (mandatory)
GENOMES="$HOME/work/subjects_genomes.txt"  # Path to a newline-delimited text file of sample-genome names (basename of filenames of FASTA files). Every line must be followed by a Linux newline character '\n', including the last line.
DB_DIR="blast_db"  # Directory for output BLAST databases (can be deleted by users afterwards)
SUBJECT_DIR="$HOME/work/assemblies"  # Directory of input subject sequences (Genome assemblies in FASTA format)
SUBJECT_SUFFIX='fna'  # Filename extension for subject sequences. Alternative values can be 'fasta', 'fa', etc.
QUERIES="$HOME/work/queries.fna"  # Path to the FASTA file of query sequences.
OUT_DIR="results"  # Output directory

# BLAST arguments
TASK='megablast'  # Or 'blastn' for short query DNA sequences
MIN_IDENTITY='70'  # Minimum nucleotide identity required for every hit
MAX_TARGET_SEQS='10'  # Maximum number of hits returned per query sequence
MAX_EVALUE='1e-5'  # Maximum acceptable e-value
```



### <a name = "step3">Step 3. Use `searchGenes.pbs` to detect genes from FASTA files<a/>

**Approach 1: PBS mode**

On a high-performance computing cluster where the PBS job scheduler is enabled:

```bash
cd ~/work
qsub -v c="searchGenes_demo.config" ~/bin/NITREc/Script/searchGenes.pbs > searchGenes.log
```

**Approach 2: bash mode**

Detect genes without submitting a PBS job:

```bash
cd ~/work
bash ~/bin/NITREc/Script/searchGenes.pbs searchGenes_demo.config > searchGenes.log
```

or

```bash
cd ~/work
chmod 750 ~/bin/NITREc/Script/searchGenes.pbs
~/bin/NITREc/Script/searchGenes.pbs searchGenes_demo.config > searchGenes.log
```



### <a name = "step4">Step 4. Compile output TSV files into a table</a>

Python script `compileBLAST.py` has been developed to compile output tab-delimited (TSV) files into a large TSV file, which can be imported into R, Excel, or LibreOffice Calc as a data frame or table for further analysis. This script requires Python v3 and a compatible [Biopython](https://biopython.org/) package.

```bash
cd ~/work
python compileBLAST.py --input results/*__megablast.tsv --delimiter '__' --genes 'gene1,gene2,gene3' --output demo --codon_table 11 --add_sample_name > compile_blast.log
```

Alternatively, a plain-text file of a list of gene names (one gene per line) can be given to the `--genes` argument using the syntax `--genes File:$file_path`. For example:

```bash
python compileBLAST.py --input results/*__megablast.tsv --delimiter '__' --genes 'File:gene_names.txt' --output demo --codon_table 11 --add_sample_name > compile_blast.log
```

The gene-name list can be created from the query FASTA file using the command that follows:

```bash
cat queries.fna | grep '>' | cut -d ' ' -f 1 | sed 's/>//g' > gene_querys.txt
```

