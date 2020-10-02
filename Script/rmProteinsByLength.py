#!/usr/bin/env python

"""
Remove protein sequences that are shorter than (<=) a given cutoff and also remove corresponding nucleotide sequences.
This script is useful for creating a multi-sequence alignment for estimating dN/dS ratios. Note that this sequence does not transfer
annotations from sequence headers.

Parameters:
    in_p: input FASTA file of protein sequences
    in_n: input FASTA file of nucleotide sequences in which sequence IDs match to those in input protein sequences
    out_pp: output FASTA file of protein sequences passed the length filter
    out_pf: output FASTA file of protein sequences failed the length filter
    out_np: output FASTA file of nucleotide sequences passed the length filter
    out_nf: output FASTA file of nucleotide sequences failed the length filter
    aa: minimum length of protein sequences

Command line:
    python rmProteinsByLength.py --in_p NfsA.faa --in_n nfsA.fna --out_pp NfsA_l.faa --out_pf NfsA_s.faa --out_np nfsA_l.fna --out_nf nfsA_s.fna --aa 240 1> NfsA_lengths.tsv

Copyright (C) 2020 Yu Wan <wanyuac@126.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
Publication: 2 Oct 2020; the latest modification: 2 Oct 2020
"""

import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from argparse import ArgumentParser


def parse_arguments():
    parser = ArgumentParser(description = "Remove protein and corresponding nucleotide sequences based on a minimum protein length",\
             epilog = "This is a helper script for database NITREc.")
    parser.add_argument("--in_p", type = str, required = True, help = "Input FASTA file of protein sequences")
    parser.add_argument("--in_n", type = str, required = True, help = "Input FASTA file of nucleotide sequences")
    parser.add_argument("--out_pp", type = str, required = True, help = "Output FASTA file of protein sequences passed the length filter")
    parser.add_argument("--out_pf", type = str, required = True, help = "Output FASTA file of protein sequences failed the length filter")
    parser.add_argument("--out_np", type = str, required = True, help = "Output FASTA file of nucleotide sequences passed the length filter")
    parser.add_argument("--out_nf", type = str, required = True, help = "Output FASTA file of nucleotide sequences failed the length filter")
    parser.add_argument("--aa", type = int, required = True, help = "A minimum protein length as the filter")
    
    return parser.parse_args()


def main():
    args = parse_arguments()

    # Sanity check
    if (not os.path.exists(args.in_p)) or (not os.path.exists(args.in_n)) or (args.aa <= 0):
        print("Argument error: please check.", file = sys.stderr)
        sys.exit(1)

    # Import sequences
    ps = SeqIO.to_dict(SeqIO.parse(args.in_p, "fasta"))
    ns = SeqIO.to_dict(SeqIO.parse(args.in_n, "fasta"))
    cdss = list()  # A list of CDSs

    # Create a list of CDS objects
    for i, pt in ps.items():
        if i in ns.keys():
            cdss.append(CDS(i, str(ns[i].seq), str(pt.seq)))
        else:
            print("Warning: protein sequence %s is not found in the input FASTA file of nucleotide sequences.", file = sys.stderr)

    # Filter and save sequences
    f_pp = open(args.out_pp, "w")
    f_pf = open(args.out_pf, "w")
    f_np = open(args.out_np, "w")
    f_nf = open(args.out_nf, "w")

    for c in cdss:
        i = c.seqid
        nt_record = ">%s\n%s\n" % (i, c.dna)
        pt_record = ">%s\n%s\n" % (i, c.prot)
        pt_len = c.len_aa
        if (pt_len < args.aa):  # Fail the length filter
            f_pf.write(pt_record)
            f_nf.write(nt_record)
            print("%s\t%i\tFail" % (i, pt_len), file = sys.stdout)
        else:  # Pass the filter
            f_pp.write(pt_record)
            f_np.write(nt_record)
            print("%s\t%i\tPass" % (i, pt_len), file = sys.stdout)

    f_pp.close()
    f_pf.close()
    f_np.close()
    f_nf.close()

    return


class CDS:
    """
    Managing a coding sequence and its corresponding protein sequence.
    Use command 'import CDS from rmProteinsByLength' to import this class into other scripts.
    """

    def __init__(self, seqid, nucl, prot):
        self.__id = seqid
        self.__dna = nucl
        self.__prot = prot

        return

    @property
    def seqid(self):
        return self.__id

    @property
    def dna(self):
        return self.__dna

    @property
    def prot(self):
        return self.__prot
    
    @property
    def len_aa(self):
        return len(self.__prot)

    @property
    def len_bp(self):
        return len(self.__dna)


if __name__ == "__main__":
    main()
