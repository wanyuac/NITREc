#!/usr/bin/env python
"""
This script is derived from script extractSeqFromMultiFASTA.py and it aims to assign sequences
from an input multi-FASTA file into one or two multi-FASTA files based on a list of sequence
IDs (default) or descriptions. Specifically, this script considers the format of sequence
headers in FASTA files as follows:

>[Sequence description]

Sometimes, the description is treated as a combination of a sequence ID and a sequence annotation,
separated by a while space. Namely, [Sequence description] = [Sequence ID] [Sequence annotation].
For instance, >nfsA_1 nfsA|isolate_1 is comprised of sequence ID "nfsA_1" and annotation "nfsA|isolate_1".

Usage:
    python splitMultiFASTA.py -i all_sequences.fna -l ids.txt -o1 in_list.fna -o2 out_list.fna  # Two output files
    python splitMultiFASTA.py -i all_sequences.fna -l ids.txt -o1 in_list.fna  # One output file

Copyright (C) 2020 Yu Wan <wanyuac@126.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
Publication: 27 July 2020; last modification: 24 July 2022
"""

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from argparse import ArgumentParser


def parse_arguments():
    parser = ArgumentParser(description = "Split a multi-FASTA file into two files")
    parser.add_argument("-i", type = str, required = True, help = "Input multi-FASTA file")
    parser.add_argument("-l", type = str, required = True, help = "Input list of sequence IDs or descriptions for inclusion. One ID/description per line.")
    parser.add_argument("-o1", type = str, required = True, help = "Output FASTA file of sequences whose IDs comprise the ID list.")
    parser.add_argument("-o2", type = str, required = False, default = None, help = "Output FASTA file of sequences whose IDs are not on the ID list. (Optional)")
    parser.add_argument("-d", action = "store_true", help = "Use the whole sequence description rather than ID for separation")
    parser.add_argument("-r", action = "store_true", help = "Remove sequence ID from sequence description: useful for CD-HIT and is disabled when -d is not turned on")
    parser.add_argument("-c", action = "store_true", help = "Concise sequence headers: removing the sequence annotation from each sequence description and only leaving the sequence ID. This argument overrides '-r'.")

    return parser.parse_args()


def main():
    args = parse_arguments()
    items_inc = read_inc_list(args.l)
    
    out_1 = open(args.o1, "w")
    if args.o2 != None:  # Save sequences that are not listed in a second output file
        out_2 = open(args.o2, "w")
        for seq in SeqIO.parse(args.i, "fasta"):  # read the input FASTA file from stdin
            if args.d:
                item = seq.description  # Sequence description = sequence ID + white space + sequence annotation
                if args.r:
                    seq.description = seq.description[len(seq.id) + 1 : ]
            else:
                item = seq.id  # Not include the sequence annotation
            
            # Assign the current sequence
            if item in items_inc:
                if args.c:
                    write_seq(seq.seq, seq.id, out_1)
                else:
                    write_seq(seq.seq, seq.description, out_1)
            else:
                if args.c:
                    write_seq(seq.seq, seq.id, out_2)
                else:
                    write_seq(seq.seq, seq.description, out_2)
        out_2.close()
    else:  # Only produce output file 1 and discard sequences whose IDs are not listed
        for seq in SeqIO.parse(args.i, "fasta"):
            if args.d:
                item = seq.description
                if args.r:
                    seq.description = seq.description[len(seq.id) + 1 : ]  # Drop sequence ID from the description
            else:
                item = seq.id  # Otherwise, the whole sequence description is copied for both the header comparison and the new sequence.
            
            # Determine whether the current sequence should be written into the output file
            if item in items_inc:
                if args.c:
                    write_seq(seq.seq, seq.id, out_1)
                else:
                    write_seq(seq.seq, seq.description, out_1)
    out_1.close()

    return


def read_inc_list(list_file):
    # Import items to be included for file out_1
    with open(list_file, "r") as f:
        ids = f.read().splitlines()

    return ids


def write_seq(seq, descr, fasta_handle):
    print(">" + descr, file = fasta_handle)
    print(seq, file = fasta_handle)

    return


if __name__ == "__main__":
    main()
    