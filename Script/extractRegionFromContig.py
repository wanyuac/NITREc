"""
This script extracts regions in contigs from several FASTA files following genomic coordinates specified.

Arguments
    -i: A tab-delimited input file (TSV) of a five/six-column table: path to each FASTA file, genome name,
        name of target contig (sequence ID), start position, end position, and optionally, a name for the
        target region (e.g., a gene). The table must not contain column names. One row per target region.
    -o: Filename of the output multiFASTA file. Every sequence header in this file follows one of the formats:
        '>[genome name].[region name] [description]', when the region name is specified;
        '>[genome name] [description]', when the region name is not specified.

Example command:
    python extractRegionFromContig.py -i coords.tsv -o targets.fna

Notes
    1. Since this script does not validate start and end coordinates for each target region, please ensure
       these coordinates do not exceed contig boundaries.
    2. Following the NCBI convention, sequence coordinates are always refer to the positive strand.
    2. Python version 2 is not supported by this script.
    3. I do not know whether this script works properly if start = end.
    4. This script is derived from extractNuclRegionFromFASTA.py in repository github.com/wanyuac/BINF_toolkit
       with adaptations for use cases of database NITREc. This script is an improved version of extractNuclRegionFromFASTA.py.

Copyright (C) 2020 Yu Wan <wanyuac@126.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
Creation: 3 Aug 2020; the latest update: 19 May 2022.
"""

import sys
from argparse import ArgumentParser
from collections import namedtuple
from Bio import SeqIO
from Bio.Seq import Seq
# from Bio.Alphabet import generic_dna  # Is no longer supported by BioPython.


def parse_args():
    parser = ArgumentParser(description = "Extract regions from contigs in several FASTA files")
    parser.add_argument("-i", "--input", dest = "i", type = str, required = True,\
                        help = "A header-free TSV of five or six columns: File path, Genome name, Contig name, Start, End, Region name (optional)")
    parser.add_argument("-o", "--output", dest = "o", type = str, required = False, default = "target_regions.fna",\
                        help = "Name of the output FASTA file")
    return parser.parse_args()


def main():
    args = parse_args()
    regions = import_regions(args.i)
    fasta_out = open(args.o, "w")

    # Loop through FASTA files
    for fasta_in, region_list in regions.items():
        contigs = SeqIO.to_dict(SeqIO.parse(fasta_in, "fasta"))  # Import the whole FASTA file
        
        # Go through target regions specified for the current FASTA file
        for region in region_list:         
            if region.contig in contigs.keys():  # The following code block assumes that start != end.
                seq_id = region.name + "." + region.genome if region.name != None else region.genome
                descr_fields = [region.name, region.contig, str(region.start) + "-" + str(region.end)]

                # Create a Seq object so the reverse complementary sequence can be easily determined.
                contig = contigs[region.contig]  # Contig is a SeqRecord object.
                seq = str(contig.seq)  # Convert a Seq object into a string
                if region.end > region.start:  # The sequence to be extracted is on the positive strand of the contig.
                    descr_fields.append("+")
                    new_seq = seq[region.start - 1 : region.end]
                else:  # Extract a sequence from the reverse complementary strand of the contig
                    descr_fields.append("-")
                    # new_seq = str(Seq(seq[region.end - 1 : region.start], generic_dna).reverse_complement())
                    new_seq = str(Seq(seq[region.end - 1 : region.start]).reverse_complement())

                # Write the extracted sequence
                print(">%s %s" % (seq_id, "|".join(descr_fields)), file = fasta_out)
                print(new_seq, file = fasta_out)  # No line wrap applies here for ease of sequence alignment.
            else:
                print("Warning: contig %s does not exist in FASTA file %s." % (region.contig, fasta_in))
    fasta_out.close()
    return


def import_regions(tsv):
    """
    This function groups regions by FASTA files. Therefore, this script allows regions extracted from the same
    genome stored in different FASTA files.
    """
    Region = namedtuple("Region", ["genome", "contig", "start", "end", "name"])
    regions = dict()  # A dictionary of lists
    with open(tsv, "r") as f:
        lines = f.read().splitlines()
        for line in lines:
            fields = line.split("\t")
            n = len(fields)
            if n == 6:  # The region name is provided in the current line.
                file_path = fields[0]
                new_region = Region(genome = fields[1], contig = fields[2], start = int(fields[3]), end = int(fields[4]), name = fields[5])
            elif n == 5:
                file_path = fields[0]
                new_region = Region(genome = fields[1], contig = fields[2], start = int(fields[3]), end = int(fields[4]), name = None)
            else:
                print("Error: row '%s' must have five or six fields." % line)
                sys.exit(1)
            
            # Push the new region into the dictionary of target regions
            if file_path in regions.keys():
                regions[file_path] += [new_region]  # Multiple regions from one FASTA file
            else:
                regions[file_path] = [new_region]
    return regions


if __name__ == "__main__":
    main()
