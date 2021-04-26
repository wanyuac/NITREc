#!/usr/bin/env python

"""
Remove sequence annotation from sequence headers (descriptions) in a FASTA file.

For example, input sequence
>allele_1 Annotation

becomes
>allele_1

after using this script.

Example command: cat seq.fna | python rmSeqAnnot.py 1> seq_out.fna

Copyright (C) 2020 Yu Wan <wanyuac@126.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
Publication: 15 June 2020; the latest modification: 3 October 2020
"""

import sys

def main():
    for line in sys.stdin:
        line = line.rstrip()
        if line.startswith(">"):
            try:
                line, _ = line.split(" ", 1)  # Split the header by the first space character.
            except ValueError:  # not enough values to unpack (expected 2, got 1)
                print("Warning: Line %s is directly copied to the output as it does not have the annotation field." % line, file = sys.stderr)
        print(line, file = sys.stdout)

    return

if __name__ == "__main__":
    main()
