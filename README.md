# A Database of Chromosomal Genetic Alterations Related to Nitrofurantoin Resistance in _Escherichia coli_

Yu Wan

Creation: 20 May 2020; latest update: 11 March 2022.

<br/>

This database consists of nucleotide and amino acid sequences of chromosomal genes involving in nitrofurantoin resistance in *E. coli*.

**Citation**

Yu Wan, Ewurabena Mills, Rhoda C.Y. Leung, Ana Vieira, Elita Jauneikaite, Xiangyun Zhi, Nicholas J Croucher, Neil Woodford, Matthew J. Ellington, Shiranee Sriskandan. **Alterations in chromosomal genes *nfsA*, *nfsB*, and *ribE* are associated with nitrofurantoin resistance in *Escherichia coli* from the United Kingdom**. *Microbial Genomics* 2021; 7. DOI: [10.1099/mgen.0.000702](https://doi.org/10.1099/mgen.0.000702).

<!-- Yu Wan, Ewurabena Mills, Rhoda C.Y. Leung, Ana Vieira, Elita Jauneikaite, Xiangyun Zhi, Nicholas J Croucher, Neil Woodford, Matthew J. Ellington, Shiranee Sriskandan. Diverse Genetic Determinants of Nitrofurantoin Resistance in UK *Escherichia coli*. *bioRxiv* 2021.05.27.446087; doi: https://doi.org/10.1101/2021.05.27.446087. -->



## 1. Format of sequence headers

NITREcMut database uses what we call ISPA (**I**D, **s**ource, **p**roduct, and **a**dditional note) format. Specifically, the header of every sequence in the NITREcMut database is comprised of two standard domains, namely, sequence ID and sequence annotation (Both domains constitute a sequence description), that are accessible through an object of BioPython's `SeqRecord` class:

```fasta
>{Sequence ID} {Sequence annotation}
```



Specifically, we use the sequence ID for the allele name, while the domain of sequence annotation consists of bar-delimited fields describing gene or cluster names, sequence source, product, etc., as shown below:

```fasta
>[Allele name] [Gene or cluster name]|[Genome name]|[NCBI nucleotide accession or contig name]|[Coordinates]|[Coding strand for the current product]|[Locus tag]|[NCBI protein accession]|[Product name];[Additional information]
```



Notes about these fields:

- Allele name: An index is appended to the allele name with a full-stop delimiter when allele names clash in a database. For example, `nfsA_1.1` distinguishes allele `nfsA_1` from the other allele of the same name.
- Gene or cluster name: In addition to conventional or pre-defined gene names, users may cluster alleles in various ways. For example, `CD-HIT-EST` is widely used to group alleles based on their nucleotide identity. Duplicated gene or cluster names can be distinguished through the same way as allele names. For instance, `ribE.1` represents a different gene to `ribE.2`, although they originally share the same name `ribE` in reference _E. coli_ genomes.
- Genome name: The isolate or strain name retrieved from an NCBI nucleotide record.
- NCBI nucleotide accession or contig name: A contig name is used when the allele sequence is not stored in the NCBI nucleotide database.
- Coordinates: Start and end positions of the allele in the forward strand (which is arbitrarily designated by a genome assembler, and in a SPAdes assembly graph, the forward strand is indicated by the plus sign at the end of node names) of a contig are separated by a dash and recorded in the same format as that in a GenBank file — namely, the start position is always less than the end position. For example, a CDS may have coordinates `608139-608792`.
- Coding strand of the current product: On which strand is the coding sequence. '+' stands for the forward strand of the NCBI record accessed through the nucleotide accession or of a source contig, and '-' stands for the reverse complementary strand of the same record or contig.
- Locus tag: A locus tag in a GenBank file if the tag is available. Otherwise, this filed equals 'NA'.
- NCBI protein accession: Accession of the current product in the NCBI Protein database if the accession is available.
- Product name: Usually it is the name of the product protein.
- Additional information, such as mutation information, can be appended to the sequence description using a semicolon as a separator.
- 'NA' is a space holder for empty fields. For example, it is used when the protein accession number or locus tag has not been assigned.



Two examples of legit headers of nucleotide and protein sequences:

```fasta
>ribE ribE|EC958|NZ_HG941718.1|463949-464419|+|EC958_RS02180|WP_001021161.1|6,7-dimethyl-8-ribityllumazine synthase

>RibE ribE|IN01|13|54543-55013|+|NA|6,7-dimethyl-8-ribityllumazine synthase

>nfsA nfsA|ATCC25922|NZ_CP009072.1|4377122-4377844|-|DR76_RS21800|WP_000189167.1|Nitroreductase NfsA

>NfsA nfsA|IN07|7|141755-142477|-|NA|NA|Nitroreductase NfsA
```



The design of sequence headers aims to be compatible to the [ResFinder](https://cge.cbs.dtu.dk/services/ResFinder/) database and keep database curation simple. Although I am familiar with the [SRST2-compatible format](https://github.com/katholt/srst2), my experience of creating its [ARGannot_r2.fasta](https://github.com/katholt/srst2/blob/master/data/ARGannot_r2.fasta) database suggests that this stringent format requires much effort for expanding the database.



### A note on performing sequence clustering for confirmation of database non-redundancy

Sometimes users may want to run `cd-hit-est` or `cd-hit` on their curated gene databases to check for non-redundancy. It is easier to inspect clusters if sequence IDs are isolate names rather than allele names (which maybe the same in each multi-FASTA file). Assuming sequence headers follow the ISPA format, then users can perform the following two steps to reformat sequence headers.

Suppose the database to be clustered is comprised of three sequences:

```fasta
>nfsA nfsA|UMN026|NC_011751.1|1074103-1074825|+|ECUMN_RS06175|WP_000189140.1|Nitroreductase NfsA
ATGACGCCAACCATTGAACTTATTTGTGGCCATCGCTCCATTCGCCATTTCACTGATGAACCCATTTCCG...

>nfsA nfsA|BR02|NZ_CP035320.1|1056121-1056843|+|EK474_RS05475|WP_000189159.1|Nitroreductase NfsA
ATGACGCCAACCATTGAACTTATTTGTGGCCATCGCTCCATTCGCCATTTCACTGATGAACCCATTTCCG...

>nfsA nfsA|ATCC25922|NZ_CP009072.1|4377122-4377844|-|DR76_RS21800|WP_000189167.1|Nitroreductase NfsA
ATGACGCCAACCATTGAACTTATTTGTGGCCATCGCTCCATTCGCCATTTCACTGATGAACCCATTTCCG...
```

- Step 1: Substitute "nfsA nfsA|" (for the corresponding protein database, use "NfsA nfsA|") with empty characters "" in a text editor using function "Find and Replace". The outcome will be:

```fasta
>UMN026|NC_011751.1|1074103-1074825|+|ECUMN_RS06175|WP_000189140.1|Nitroreductase NfsA
ATGACGCCAACCATTGAACTTATTTGTGGCCATCGCTCCATTCGCCATTTCACTGATGAACCCATTTCCG...

>BR02|NZ_CP035320.1|1056121-1056843|+|EK474_RS05475|WP_000189159.1|Nitroreductase NfsA
ATGACGCCAACCATTGAACTTATTTGTGGCCATCGCTCCATTCGCCATTTCACTGATGAACCCATTTCCG...

>ATCC25922|NZ_CP009072.1|4377122-4377844|-|DR76_RS21800|WP_000189167.1|Nitroreductase NfsA
ATGACGCCAACCATTGAACTTATTTGTGGCCATCGCTCCATTCGCCATTTCACTGATGAACCCATTTCCG...
```

- Step 2: Substitute "|" with white spaces " " in the same way as step 1, generating a FASTA file as follows:

```fasta
>UMN026 NC_011751.1 1074103-1074825 + ECUMN_RS06175 WP_000189140.1 Nitroreductase NfsA
ATGACGCCAACCATTGAACTTATTTGTGGCCATCGCTCCATTCGCCATTTCACTGATGAACCCATTTCCG...

>BR02 NZ_CP035320.1 1056121-1056843 + EK474_RS05475 WP_000189159.1 Nitroreductase NfsA
ATGACGCCAACCATTGAACTTATTTGTGGCCATCGCTCCATTCGCCATTTCACTGATGAACCCATTTCCG...

>ATCC25922 NZ_CP009072.1 4377122-4377844 - DR76_RS21800 WP_000189167.1 Nitroreductase NfsA
ATGACGCCAACCATTGAACTTATTTGTGGCCATCGCTCCATTCGCCATTTCACTGATGAACCCATTTCCG...
```

Since `cd-hit-est` and `cd-hit` read sequence headers till the first space when parameter `-d` equals zero, only isolate names will appear in the cluster file, thus solving the problem of sequence clustering.



## 2. Filename convention

The name of every FASTA file (`.fna` for nucleotide sequences and `.faa` for protein sequences) of this database consists of two fields: gene or cluster name, and an indicator for nitrofurantoin susceptible ('S') or resistant ('R') alleles. Double underscores separate these two fields. For example, `nfsA__R.fna` stores alleles of the gene _nfsA_ in nitrofurantoin resistant _E. coli_ (NITREc).

Gene or cluster names are kept in sequence headers for convenience of concatenating sequence files for some software or analyses, despite presence of these names in the names of FASTA files.



## 3. Content of nucleotide and protein sequences

For each gene in this database, some nucleotide sequences (in subdirectory `Nucl`) do not have their corresponding protein sequences (matched by sequence descriptions) in subdirectory `Prot`, and vice versa, because of synonymous mutations and sometimes literature does not report nucleotide sequences that encode protein mutants. Nonetheless, the protein database (`Prot`) is often adequate for functional prediction of proteins.