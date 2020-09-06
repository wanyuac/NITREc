"""
This script finds out amino acid substitutions of interest, given a table of known substitutions (mutations).
Specifically, it takes as input *_subs.tsv from the output of script missenseFinder.py and annotates mutations
based on the table of known mutations. It creates a TSV file as output. I decided not to merge functions of
this script into missenseFinder.py for flexibility.

The table of known mutations should be a header-free TSV file with three columns: position, reference amino acid,
and alternative (substituted) amino acid, which is similar to a VCF file.

Example command:
    python findKnownMutations.py -i NfsA_subs.tsv -t NfsA_mutations.tsv -o NfsA_subs_annot.tsv

Copyright (C) 2020 Yu Wan <wanyuac@126.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
Creation: 6 Sep 2020
"""

from argparse import ArgumentParser
from collections import defaultdict

# Define annotation code of mutation types
CODE_KNOWN = "1"
CODE_NOVEL = "0"


def parse_arguments():
    parser = ArgumentParser(description = "Identify known mutations in a mutation table",\
        epilog = "This script is a helper for using database NITREc.")
    parser.add_argument("-i", dest = "i", required = True, type = str, help = "Input table of amino acid substitutions (output of missenseFinder.py)")
    parser.add_argument("-t", dest = "t", required = True, type = str, help = "TSV-format, header-free table of known mutations (columns: Pos, Ref, Alt")
    parser.add_argument("-o", dest = "o", required = True, type = str, help = "Name of output TSV file")

    return parser.parse_args()


def main():
    args = parse_arguments()
    subs = import_mutations(args.i)
    refs = import_known_mutations(args.t)
    annots = open(args.o, "w")
    print("\t".join(["Query", "Reference", "Ref", "Pos", "Alt", "Known"]), file = annots)

    for line in subs:
        q_name, r_name, aa_r, pos, aa_q = line.split("\t")
        if pos in refs.keys():
            record = refs[pos]
            if aa_r in record.keys():
                if aa_q in record[aa_r]:
                    print("\t".join([q_name, r_name, aa_r, pos, aa_q, CODE_KNOWN]), file = annots)
                else:
                    print("\t".join([q_name, r_name, aa_r, pos, aa_q, CODE_NOVEL]), file = annots)
            else:
                print("\t".join([q_name, r_name, aa_r, pos, aa_q, CODE_NOVEL]), file = annots)
        else:  # Not a known mutation
            print("\t".join([q_name, r_name, aa_r, pos, aa_q, CODE_NOVEL]), file = annots)  # 0: a novel mutation

    annots.close()

    return


def import_mutations(f):
    """ Import amino acid substitutions in query proteins """
    with open(f, "r") as f_in:
        lines = f_in.read().splitlines()
        lines = lines[1 : ]  # Skip the header line
    
    return lines


def import_known_mutations(f):
    """
    Read the table of known mutations into a two-dimentional dictionary having the structure:
    {Position : {Reference : [alt_1, alt_2, ...]}}
    Note that although it is often the case that one position has only one reference amino acid,
    at the species level the position may not be conserved.
    """
    refs = defaultdict(dict)
    with open(f, "r") as f_in:
        lines = f_in.read().splitlines()

    for line in lines:
        pos, ref, alt = line.split("\t")
        if pos in refs.keys():
            if ref in refs[pos].keys():
                refs[pos][ref].append(alt)
            else:
                refs[pos][ref] = [alt]
        else:  # A new position is encountered: initialise an element
            refs[pos] = {ref : [alt]}

    return refs


if __name__ == "__main__":
    main()
