"""
This script extracts regions in contigs from several FASTA files following genomic coordinates specified.

Arguments
    -i: A tab-delimited input file (TSV) of a five-column table: Genome name, Path to each FASTA file,
        Name of target contig (sequence ID), Start position, and End position. No column name should be
        contained. One row per target region.
    -n: A name (such as gene name or feature name) constantly added to the first field of output sequence
        description. Default: None.
    -o: Filename of the output FASTA file

Example command line:
    python extractRegionFromContig.py -i targets.tsv -n ribE -o targets.fna
    
Notes
    1. Since this script does not validate start and end coordinates for each target region, please ensure
       these coordinates do not exceed contig boundaries.
    2. Following the NCBI convention, sequence coordinates are always refer to the positive strand.
    2. Python version 2 is not supported by this script.
    3. I do not know whether this script works properly if start = end.
    4. This script is derived from extractNuclRegionFromFASTA.py in repository github.com/wanyuac/BINF_toolkit
       with adaptations for use cases of database NITREcMut. This script is an improved version of script
       extractNuclRegionFromFASTA.py.

Copyright (C) 2020 Yu Wan <wanyuac@126.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
Creation: 03 Aug 2020
"""

from argparse import ArgumentParser
from collections import namedtuple
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna


def parse_args():
    parser = ArgumentParser(description = "Extract regions from contigs in several FASTA files")
    parser.add_argument("-i", type = str, required = True, help = "A no-header TSV of five columns: Genome name, File path, Contig name, Start, End")
    parser.add_argument("-n", type = str, required = False, default = None, help = "A name constantly added to the first field of output sequence")
    parser.add_argument("-o", type = str, required = False, default = "target_regions.fna", help = "Name of the output FASTA file")

    return parser.parse_args()


def main():
    args = parse_args()
    regions = import_regions(args.i)
    add_name = not args.n == None
    fasta_out = open(args.o, "w")

    # Loop through FASTA files
    for fasta_in, region_list in regions.items():
        contigs = SeqIO.to_dict(SeqIO.parse(fasta_in, "fasta"))  # Import the whole FASTA file
        
        # Go through target regions specified for the current FASTA file
        for region in region_list:         
            if region.contig in contigs.keys():  # The following code block assumes that start != end.
                # Create invariable fields of sequence description
                descr_fields = [region.genome, region.contig, str(region.start) + "-" + str(region.end)]
                if add_name:
                    descr_fields = [args.n] + descr_fields
                
                # Create a Seq object
                contig = contigs[region.contig]  # Contig is a SeqRecord object.
                seq = str(contig.seq)  # Convert a Seq object into a string
                if region.end > region.start:  # The sequence to be extracted is on the positive strand of the contig.
                    descr_fields.append("+")
                    new_seq = seq[region.start - 1 : region.end]
                else:  # Extract a sequence from the reverse complementary strand of the contig
                    descr_fields.append("-")
                    new_seq = str(Seq(seq[region.end - 1 : region.start], generic_dna).reverse_complement())

                # Write the extracted sequence
                print(">" + "|".join(descr_fields), file = fasta_out)
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
    Region = namedtuple("Region", ["genome", "contig", "start", "end"])
    regions = dict()  # A dictionary of lists
    with open(tsv, "r") as f:
        lines = f.read().splitlines()
        for line in lines:
            genome, file_path, contig, start, end = line.split("\t")
            new_region = Region(genome = genome, contig = contig, start = int(start), end = int(end))
            if file_path in regions.keys():
                regions[file_path] += [new_region]
            else:
                regions[file_path] = [new_region]

    return regions


if __name__ == "__main__":
    main()
