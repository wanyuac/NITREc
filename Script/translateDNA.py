#!/usr/bin/env python3

"""
Translate DNA sequences into protein sequences using a given codon table. The translation terminates at the
first stop codon.

Command:
    translateDNA.py -c [codon table number] -k

Example commands:
cat seq.fna | python translateDNA.py > seq.faa  # Print error messages on screen
cat seq.fna | python translateDNA.py -c 11 -k 1>seq.faa 2>seq.err
cat seq.fna | python translateDNA.py -c 10 | python rmSeqDescr.py > seq.faa

Default codon table number is 11 (for prokaryotic genes).

Copyright (C) 2020 Yu Wan <wanyuac@126.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
Publication: 15 June 2020; the latest modification: 23 May 2021
"""

import sys
import warnings
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Data.CodonTable import TranslationError
from Bio import BiopythonWarning
from argparse import ArgumentParser


def parse_arguments():
    parser = ArgumentParser("Convert nucleotide sequences into protein sequences")
    parser.add_argument("-c", type = int, required = False, default = 11, action = "store", help = "Codon table (Default: 11)")
    parser.add_argument("-k", action = "store_true", help = "Use sequence IDs as protein names without capitalise the first letter (Default: off)")
    parser.add_argument("-s", action = "store_true", help = "Check each sequence for a valid start codon (Default: off)")
    parser.add_argument("-f", action = "store_true", help = "Force translation even if partial codons are found. (Default: do not translate)")
    return parser.parse_args()


def main():
    # Initialisation
    args = parse_arguments()
    codon_tab = args.c
    print("Codon table: %i" % codon_tab, file = sys.stderr)
    if args.s:
        print("Warning: alternative start codon will be translated to \"M\" when args.s = True.", file = sys.stderr)

    # Go through sequences from the stdin
    for rec in SeqIO.parse(sys.stdin, "fasta"):
        # Process the sequence description
        seqid = rec.id
        seqid_len = len(seqid)
        descr = "" if seqid_len == len(rec.description) else rec.description[(seqid_len + 1) : ]  # Drop seqid from rec.description
        seqid_new = seqid if args.k else geneName_to_proteinName(seqid)
        
        # Translate the DNA sequence into a protein sequence
        with warnings.catch_warnings():
            """
            https://raw.githubusercontent.com/biopython/biopython/master/Bio/Seq.py:
            elif n % 3 != 0:
                warnings.warn(
                "Partial codon, len(sequence) not a multiple of three. "
                "Explicitly trim the sequence or add trailing N before "
                "translation. This may become an error in future.",
                BiopythonWarning,)
            The following code aims to catch this BiopythonWarning in order to deal with partial codons.
            """
            warnings.simplefilter("error", BiopythonWarning)  # Turn warnings into exceptions
            try:
                # Set cds = True to check error: 'First codon is not a start codon'. Exception: Bio.Data.CodonTable.TranslationError.
                # Note that the start codon will always become "M" in this setting, which is not desirable.
                # Set to_stop = False to print "*" that represents stop codons. Otherwise, the asterisk is not printed.
                rec_prot = rec.translate(table = codon_tab, id = seqid_new, description = descr, to_stop = True, cds = args.s, stop_symbol = "*")
                trans_succ = True
            except KeyError:
                print("Warning: sequence \"%s\" cannot be translated." % rec.description, file = sys.stderr)
                #rec_prot = SeqRecord(Seq(''), id = '', name = '', description = '')
                trans_succ = False
            except TranslationError:
                print("First codon '%s' of %s is not a start codon. Skip translation of this sequence." % (str(rec.seq)[ : 3], seqid), file = sys.stderr)
                trans_succ = False
            except BiopythonWarning:
                if args.f:  # Force translation
                    print("Warning: partial codon is found in sequence %s. Translate this sequence as per the argument '-f'." % seqid, file = sys.stderr)
                    """ It must be rerun as the result from the last translation command will be discarded before assigning it to rec_prot when the warning is captured.
                    Otherwise, an error of \"local variable 'rec_prot' referenced before assignment\" or duplicated record occur. """
                    rec_prot = rec.translate(table = codon_tab, id = seqid_new, description = descr, to_stop = True, cds = False, stop_symbol = "*")
                    trans_succ = True
                else:
                    print("Warning: partial codon is found in sequence %s. Skip translation of this sequence." % seqid, file = sys.stderr)
                    trans_succ = False

        # Print protein sequences
        if trans_succ:
            if descr == "":
                print(">" + rec_prot.id, file = sys.stdout)
            else:
                print(">%s %s" % (rec_prot.id, rec_prot.description), file = sys.stdout)
            print(rec_prot.seq, file = sys.stdout)
    return


def geneName_to_proteinName(g):
    # Converting a gene name to protein name through capitalising the first letter of the gene name
    p = g[0].upper() + g[1 : ]
    return p


if __name__ == "__main__":
    main()
