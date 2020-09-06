"""
This script identifies amino acid substitutions against reference protein sequences.
Assumptions of this script for input sequences: (1) identical length, (2) high homology, (3) distinct sequence IDs.

Example command:
    python missenseFinder.py -q queries.faa -r NfsA__S.faa -o NfsA

Principal output files:
    *_hits.tsv: BLASTp output in format 6
    *_subs.tsv: Amino acid substitutions
    *_num.tsv: Number of amino acid substitutions per query sequence
    *.log: Log file

Dependencies:
    1. Programs blastn and makeblastdb must be accessible through the environmental variable $PATH.
    2. Python 3 and package BioPython.

Copyright (C) 2020 Yu Wan <wanyuac@126.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
Creation: 5 Sep 2020
"""

from argparse import ArgumentParser
from Bio import Seq
from Bio import SeqRecord
from Bio import SeqIO
import os
import sys
import logging
import subprocess


def parse_arguments():
    parser = ArgumentParser(description= "Identify amino acid substitutions in query sequences compared to their most closely related references.",\
        epilog = "This script is a helper for using database NITREc.")
    parser.add_argument("-q", "--queries", dest = "queries", type = str, required = True, help = "A multi-FASTA file of query protein sequences")
    parser.add_argument("-r", "--refs", dest = "refs", type = str, required = True, help = "A multi-FASTA file of reference protein sequences")
    parser.add_argument("-o", "--out", dest = "out", type = str, required = False, default = "protein", help = "Prefix of output filenames")
    parser.add_argument("-n", "--new_db", dest = "new_db", action = "store_true", help = "Override previous BLAST database, if it has been created.")

    return parser.parse_args()


def main():
    args = parse_arguments()

    # Set up a logger
    logging.basicConfig(filename = args.out + ".log", level = logging.INFO)  # Level INFO covers WARNING, ERROR, and CRITICAL levels.
    logger = logging.getLogger()

    # Import protein sequences
    logger.info("Read query and reference sequences from %s and %s." % (args.queries, args.refs))
    qs = SeqIO.to_dict(SeqIO.parse(args.queries, "fasta"))  # Key: sequence IDs; values: SeqRecord objects
    rs = SeqIO.to_dict(SeqIO.parse(args.refs, "fasta"))

    # Create a reference BLAST database
    db_name = os.path.basename(args.refs)
    if args.new_db and os.path.exists(db_name + "psq"):  # For simplicity, .phr and pin files are not checked, although this is not stringent.
        logger.info("Skip overriding existing BLAST database.")
    else:
        logger.info("Create BLAST database %s." % db_name)
        cmd = ["makeblastdb", "-in", args.refs, "-out", db_name, "-dbtype", "prot"]  # The database will appear under the current working directory.
        try:
            mkdb = subprocess.run(cmd, check = True, text = True, capture_output = True)  # Both stdout and stderr are redirected to subprocess.PIPE.
            logger.info(mkdb.stdout)
        except subprocess.CalledProcessError:
            logger.error("Failed to create a BLAST database.")
            logger.error(mkdb.stderr)

    # Run BLASTp to identify the best match per query sequence
    logger.info("Run BLASTp search.")
    blastp_out = args.out + "_hits.tsv"
    cmd = ["blastp", "-query", args.queries, "-db", db_name, "-out", blastp_out, "-task", "blastp-fast", "-evalue", "1e-5",\
           "-max_target_seqs", "1", "-outfmt", "6 qseqid sseqid mismatch gapopen length qlen slen qstart qend sstart send"]
    try:
        blastp = subprocess.run(cmd, check = True, text = True, capture_output = True)  # Normal completion does not return any message.
    except subprocess.CalledProcessError:
        logger.error("Failed to run BLASTp")
        logger.error(blastp.stderr)
    top_hits = import_blastp_output(blastp_out)

    # Identify mutations
    logger.info("Identify mutations in each query sequence against its most similar reference sequence.")
    find_mutations(qs, rs, top_hits, args.out + "_subs.tsv", args.out + "_num.tsv")

    return


def import_blastp_output(tsv):
    """ Read the first two columns in the output TSV file """
    top_hits = dict()  # query : reference
    with open(tsv, "r") as f:
        lines = f.read().splitlines()
    
    for hit in lines:
        fields = hit.split("\t")
        top_hits[fields[0]] = fields[1]

    return top_hits


def find_mutations(queries, refs, top_hits, sub_out, num_out):
    """ Identify amino acid substitutions in each query sequence """
    f_sub = open(sub_out, "w")
    f_num = open(num_out, "w")
    print("\t".join(["Query", "Reference", "Ref", "Pos", "Alt"]), file = f_sub)  # Ref and Alt: reference and alternative amino acid
    print("\t".join(["Query", "Reference", "Substitutions"]), file = f_num)

    for q_name, q_seq in queries.items():
        if q_name in top_hits.keys():  # The current query sequence has a hit.
            r_name = top_hits[q_name]  # Get the name of the best matched reference sequence
            r = str(refs[r_name].seq)  # Reference sequence
            q = str(q_seq.seq)  # Query sequence
            aa_num = min(len(q), len(r))  # Although it has been assumed that len(q) = len(r), I use the min function for safety.
            sub_count = 0
            for i in range(0, aa_num):  # Reference of this code block: mutations_final_v2.py by Ana Vieira.
                aa_q = q[i]  # Only access this element once in order to save time
                aa_r = r[i]
                if aa_q != aa_r:  # A mutation is found.
                    print("\t".join([q_name, r_name, aa_r, str(i + 1), aa_q]), file = f_sub)
                    sub_count += 1
            print("\t".join([q_name, r_name, str(sub_count)]), file = f_num)
        else:  # The current query sequence does not have any hit in the reference database, probably due to a lack of homology.
            print("\t".join([q_name, "NA", "NA"]), file = f_num)  # No reference sequence is available for the current query sequence.

    f_sub.close()
    f_num.close()

    return


if __name__ == "__main__":
    main()
