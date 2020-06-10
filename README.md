# A Database of Chromosomal Genetic Alterations Related to Nitrofurantoin Resistance in _Escherichia coli_

_Yu Wan_

Creation: 20 May 2020; latest update: 10 June 2020.

<br/>

This database consists of nucleotide and amino acid sequences of chromosomal genes involving in nitrofurantoin resistance in *E. coli*.



## 1. Format of sequence headers

The header of every sequence in the NITREcMut database is comprised of two standard domains, namely, sequence ID and sequence description, that are accessible through an object of BioPython's `SeqRecord` class:

```fasta
>{Sequence ID} {Sequence description}
```

Specifically, the sequence ID is the allele name, while the domain of sequence description consists of bar-delimited fields describing gene or cluster names, sequence source, product, etc., as shown below:

```fasta
>[Allele name] [Gene or cluster name]|[Genome name]|[NCBI nucleotide accession or contig name]|[Coding strand]|[Coordinates]|[Coordinate strand]|[Locus tag]|[NCBI protein accession]|[Product name];[Additional information]
```



Notes about these fields:

- Allele name: An index is appended to the allele name with a full-stop delimiter when allele names clash in a database. For example, `nfsA_1.1` distinguishes allele `nfsA_1` from the other allele of the same name.
- Gene or cluster name: In addition to conventional or pre-defined gene names, users may cluster alleles in various ways. For example, `CD-HIT-EST` is widely used to group alleles based on their nucleotide identity. Duplicated gene or cluster names can be distinguished through the same way as allele names. For instance, `ribE.1` represents a different gene to `ribE.2`, although they originally share the same name `ribE` in reference _E. coli_ genomes.
- Genome name: the isolate or strain name retrieved from an NCBI nucleotide record.
- NCBI nucleotide accession or contig name: a contig name is used when the allele sequence is not stored in the NCBI nucleotide database.
- Coding strand: on which strand is the coding sequence. '+' stands for the forward strand of the NCBI record accessed through the nucleotide accession or of a source contig, and '-' stands for the reverse complementary strand of the same record or contig.
- Coordinates: start and end positions of the allele in the forward strand, recorded in the same order as that in a GenBank file, and are separated by a dash. For example, `608139-608792`.
- Coordinate strand: which strand do the coordinates refer to. For sequences from an NCBI nucleotide record, the coordinates always refer to the forward strand ('+'), whereas sequences extracted from the reverse strand ('-') of an unpublished contig is usually noted by coordinates on the same strand. 
- Product name: usually it is the name of the product protein.
- Additional information, such as mutation information, can be appended to the sequence description using a semicolon as a separator.
- 'NA' is a space holder for empty fields. For example, it is used when the protein accession number or locus tag has not been assigned.



Two examples of legit headers of nucleotide and protein sequences:

```fasta
>ribE ribE|EC958|NZ_HG941718.1|+|463949-464419|+|EC958_RS02180|WP_001021161.1|6,7-dimethyl-8-ribityllumazine synthase

>RibE ribE|IN01|13|+|54543-55013|+|NA|6,7-dimethyl-8-ribityllumazine synthase

>nfsA nfsA|ATCC25922|NZ_CP009072.1|-|4377122-4377844|+|DR76_RS21800|WP_000189167.1|Nitroreductase NfsA

>NfsA nfsA|IN07|7|-|141755-142477|-|NA|NA|Nitroreductase NfsA
```



The design of sequence headers aims to be compatible to the [ResFinder](https://cge.cbs.dtu.dk/services/ResFinder/) database and keep database curation simple. Although I am familiar with the [SRST2-compatible format](https://github.com/katholt/srst2), my experience of creating its [ARGannot_r2.fasta](https://github.com/katholt/srst2/blob/master/data/ARGannot_r2.fasta) database suggests that this stringent format requires much effort for expanding the database.



## 2. Filename convention

The name of every FASTA file (`.fna` for nucleotide sequences and `.faa` for protein sequences) of this database consists of two fields: gene or cluster name, and an indicator for nitrofurantoin susceptible ('S') or resistant ('R') alleles. Double underscores separate these two fields. For example, `nfsA__R.fna` stores alleles of the gene _nfsA_ in nitrofurantoin resistant _E. coli_ (NITREc).

Gene or cluster names are kept in sequence headers for convenience of concatenating sequence files for some software or analyses, despite presence of these names in the names of FASTA files.

