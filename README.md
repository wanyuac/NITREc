# A Nucleotide database of Chromosomal Mutations Related to Nitrofurantoin Resistance in _Escherichia coli_

_Yu Wan_

Creation: 20 May 2020; latest update: 20 May 2020.



## 1. Format of sequence headers

The header of every sequence in the NITREcMut database is comprised of two standard domains, namely, sequence ID and sequence description, that are accessible through BioPython's `SeqRecord` object:

```fasta
>{Sequence ID} {Sequence description}
```

Specifically, the sequence ID is the allele name, while the domain of sequence description consists of bar-delimited fields describing gene or cluster names, sequence source, product, etc., as shown below:

```fasta
>[Allele name] [Gene or cluster name]|[NCBI nucleotide accession or contig name]|[Strand]|[Coordinates][Locus tag][NCBI protein accession]|[Product name]
```



Notes about these fields:

- Allele name: An index is appended to the allele name with a full-stop delimiter when allele names clash in a database. For example, `nfsA_1.1` distinguishes allele `nfsA_1` from the other allele of the same name.
- Gene or cluster name: In addition to conventional or pre-defined gene names, users may cluster alleles in various ways. For example, `CD-HIT-EST` is widely used to group alleles based on their nucleotide identity. Duplicated gene or cluster names can be distinguished through the same way as allele names. For instance, `ribE.1` represents a different gene to `ribE.2`, although they originally share the same name `ribE` in reference _E. coli_ genomes.
- NCBI nucleotide accession or contig name: a contig name is used when the allele sequence is not stored in the NCBI nucleotide database. The contig name may be preceded by a genome name, separated by a hash sign.
- Strand: '+' stands for the forward strand of the NCBI record accessed through the nucleotide accession or of a source contig, and '-' stands for the reverse complementary strand of the same record or contig.
- Coordinates: start and end positions of the allele in the forward strand, recorded in the same format as that in a GenBank file. For example, `608139..608792`.
- Product name: usually it is the name of the product protein.
- 'NA' is a space holder for empty fields. For example, it is used when the protein accession number is unavailable.



Two examples of legit sequence headers:

```fasta
>ribE ribE|NZ_HG941718.1|+|463949..464419|EC958_RS02180|WP_001021161.1|6,7-dimethyl-8-ribityllumazine synthase

>>ribE ribE|IN01#13|+|54543..55013|NA|6,7-dimethyl-8-ribityllumazine synthase
```