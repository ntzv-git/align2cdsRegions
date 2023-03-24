#! /bin/python3
# -*- coding: utf-8 -*-

"""
The program adds to the input alignments, the region where they match, and the cds information nearest to them. These
regions can be a 5' flanking region (5UTR), a 3' flanking region (3UTR), a CDS region (CDS), a small intergenic
region (UTR), a intercistronic region (ICR), an overlapping CDS region (OVL), or a region located elsewhere on the
sequence (OTHER).
The program starts by parsing the Gene Features File (GFF) and genome file (FASTA) which respectively contain the CDS
and the length informations of the subject. It then transforms the subject sequence into a region dictionary with the
region coordinates (start, end) as key and the subject region name (CDS, UTR, etc...) as value.
For each alignment, the program returns the CDS nearest to the query, using the coordinates of the subject center. Thus,
for each alignment, a nearest CDS is added, unless no CDS is annotated on the subject strand on which the query matches.
An OTHER annotation is indicated when the region where the query matches is neither a CDS, nor a UTR, nor an ICR, nor
an OVL. The query matches may be far from an annotated CDS or may be on a CDS located on the other strand.
If verbose, the program returns in the log file the total sizes of each region (on both strands).
Note that in the input alignment file, the subject end position must be greater than the subject start position.


Author  : Emmanuel Clostres
Mail    : emmanuel.clostres@univ-rennes.fr
Python  : v3.8
Version : 1.0
Date    : 24/03/2023
"""

import getopt
import gzip
import os
import sys
import time

from datetime import datetime
from Bio import SeqIO


# Default parameters
INPUT = None
OUTPUT = None
GFF = None
FASTA = None
SSEQID = None
SSTART = None
SEND = None
STRAND = None
FS = "\t"
HEADER = False
VERBOSE = False
FORCE = False
UTR5_SIZE = 20
UTR3_SIZE = 150
ICR_SIZE = 5
LOG_FILE = f"align2cdsRegions.log"
REGION_CONVERT = {
    "CDS": "CDS", "5UTR": "5UTR", "3UTR": "3UTR", "UTR": "UTR", "OTHER": "OTHER", "ICR": "ICR", "OVL": "OVL",
    "CDS-CDS": "CDS-CDS", "5UTR-3UTR": "5UTR-3UTR", "3UTR-5UTR": "5UTR-3UTR",
    "CDS-5UTR": "CDS-5UTR", "5UTR-CDS": "CDS-5UTR", "CDS-3UTR": "CDS-3UTR", "3UTR-CDS": "CDS-3UTR",
    "CDS-UTR": "CDS-UTR", "UTR-CDS": "CDS-UTR", "CDS-OVL": "CDS-OVL", "OVL-CDS": "CDS-OVL",
    "CDS-ICR": "CDS-ICR", "ICR-CDS": "CDS-ICR", "5UTR-OTHER": "5UTR-OTHER", "OTHER-5UTR": "5UTR-OTHER",
    "3UTR-OTHER": "3UTR-OTHER", "OTHER-3UTR": "3UTR-OTHER", "CDS-OTHER": "CDS-OTHER", "OTHER-CDS": "CDS-OTHER"}


class GFFasDict:
    """
    Converts a GFF file into an object with CDS and region informations.
    """
    def __init__(self, gff_path: str, fasta_path: str):
        list_gfflines = parse_file_lines(gff_path)
        dict_seqid_len = self.__parse_fasta(fasta_path)
        self.dict_seq_data = self.__parse_gff(list_gfflines, dict_seqid_len)
        self.dict_region_locannot, self.dict_region_cdslines = self.__parse_region()
        self.region_sizes = self.__genome_region_sizes()

    @staticmethod
    def __parse_fasta(file_path: str) -> dict:
        """
        Fetches all sequence length from a given fasta file.
        :param file_path: Path to the fasta file
        :return: Dictionary with sequence id as keys and sequence length as values
        """
        dict_seqid_len = dict()
        with open_file(file_path, "rt") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                dict_seqid_len[record.id] = len(record.seq)

        # for record in SeqIO.parse(file_path, "fasta"):
        #     dict_seqid_len[record.id] = len(record.seq)
        return dict_seqid_len

    @staticmethod
    def __parse_gff(list_gfflines: list, dict_seqid_len: dict) -> dict:
        """
        Parses all the lines of the GFF file to get for each sequence id, its length and the CDS coordinates.
        It returns a dictionary formated as :
        dict_seq_data = {'seqid': {'len': seqlength, '+': {'CDS':{(start, end): cds_line, ...}, ...}, '-': ...}, ...}
        :param list_gfflines: List of all lines in the GFF file
        :param dict_seqid_len: Dictionary with sequence id as keys and sequence length as values
        :return: Dictionary with for each sequence id, its length and the CDS coordinates
        """
        # Parses sequence_id and sequence_length using FASTA
        dict_seq_data = dict()
        for seqid, seqlen in dict_seqid_len.items():
            dict_seq_data[seqid] = {"+": dict(), "-": dict(), "len": seqlen}
            for strand in ["+", "-"]:
                dict_seq_data[seqid][strand] = {"CDS": dict()}

        # Parses all feature_type data using GFF
        for line in list_gfflines:
            # Pass commented lines
            if line[0] == "#":
                continue
            else:
                # Parse the feature information
                line_split = line[:-1].split(sep="\t")
                seqid = line_split[0]
                feature_type = line_split[2]
                feature_start = line_split[3]
                feature_end = line_split[4]
                feature_strand = line_split[6]

                if feature_type == "CDS":
                    dict_seq_data[seqid][feature_strand][feature_type][(int(feature_start), int(feature_end))] = line

        return dict_seq_data

    def __parse_region(self) -> (dict, dict):
        """
        Considering all sequence ids and both strands, this method returns two dictionaries : the first one with
        coordinates (start, end) as keys and corresponding regions as values ; and the second one with coordinates
        (start, end) as keys and corresponding CDS informations as values.
        :return: Return a dictionary of regions coordinates and a dictionary of CDS coordinates
        """
        # Initialization
        dict_region_locannot = dict()  # e.g. {seq_id: {'+': {(region_start, region_end): region, ...}, '-': ...}
        dict_region_cdsloclines = dict()  # e.g. {seq_id: {'+': {(cds-region_start, cds-region_end): cds_line,...}, ...}

        # Iteration on each sequence id
        for seq_id in self.dict_seq_data.keys():
            dict_region_locannot[seq_id] = {"+": dict(), "-": dict()}
            dict_region_cdsloclines[seq_id] = {"+": dict(), "-": dict()}
            seq_len = self.dict_seq_data[seq_id]["len"]

            # Iteration on each strand
            for strand in ["+", "-"]:
                dict_region_locannot[seq_id][strand], dict_region_cdsloclines[seq_id][strand] = \
                    self.__parse_region_cds(self.dict_seq_data[seq_id], seq_len, strand)

        return dict_region_locannot, dict_region_cdsloclines

    @staticmethod
    def __parse_region_cds(dict_contig_cdsloc: dict, contig_len: int, strand: chr) -> (dict, dict):
        """
        Parse the CDS dict to returns a region dict {(start, stop) : region, ...} and the corresponding cds
        dict {(start, pos): cds_line, ...} for a given sequence and a given strand.
      - CDS   : coding sequence
      - 5UTR  : 5' flanking region
      - 3UTR  : 3' flanking region
      - UTR   : short intergenic region
      - ICR   : intercistronic region
      - OVL   : overlap region of 2 CDS
      - OTHER : other region on the genome (in all other cases; result not interpretable)

        :param dict_contig_cdsloc: Dictionary of all CDS coordinates in a given sequence id
        :param contig_len: Length of the sequence
        :param strand: Strand of the targeted sequence
        :return: Returns a dictionary of regions and a dictionary of CDS coordinates (cut by region if overlapping)
        """
        # Initialization
        region_dict = dict()  # dictionary storing all region positions
        cds_idx_dict = dict()  # dictionary storing all cds lines
        region_start = 1
        region_end = contig_len
        startend_line_dict = dict_contig_cdsloc[strand]["CDS"]
        startend_list = list(startend_line_dict.keys())

        if strand == "+":
            left_ext = UTR5_SIZE
            right_ext = UTR3_SIZE
            chr_left = "5"
            chr_right = "3"
        elif strand == "-":
            left_ext = UTR3_SIZE
            right_ext = UTR5_SIZE
            chr_left = "3"
            chr_right = "5"
        else:
            raise ValueError("The strand can only be worth '+' or '-'.")

        # Processing
        # If there is no CDS on the sequence
        if len(startend_list) == 0:
            region_dict[(region_start, region_end)] = "OTHER"
            return region_dict, cds_idx_dict
        else:
            startend_list.sort(key=lambda x: x[0])

        # Fetches the first CDS
        cds_start, cds_end = startend_list[0]
        if (cds_start - region_start) > left_ext:
            region_dict[(region_start, (cds_start - left_ext - 1))] = "OTHER"
            region_dict[((cds_start - left_ext), (cds_start - 1))] = chr_left + "UTR"
        elif (cds_start - region_start) > 0:
            region_dict[(region_start, (cds_start - 1))] = "UTR"
        prev_regioncds_start = cds_start
        prev_regioncds_end = cds_end
        prev_cds_start, prev_cds_end = startend_list[0]

        # Searching for and adding all CDS
        for startend in startend_list[1:]:
            cds_start, cds_end = startend

            # Overlap region
            if (cds_start - (prev_regioncds_end + 1)) < 0:
                cds_idx_dict[(prev_regioncds_start, (prev_regioncds_end - (prev_regioncds_end - cds_start) - 1))] = \
                    startend_line_dict[(prev_cds_start, prev_cds_end)]
                region_dict[(prev_regioncds_start, (prev_regioncds_end - (prev_regioncds_end - cds_start) - 1))] = "CDS"
                region_dict[((prev_regioncds_end - (prev_regioncds_end - cds_start)), prev_regioncds_end)] = "OVL"
                cds_start = (prev_regioncds_end + 1)  # For prev_cds_start to take the value '(prev_cds_stop + 1)'

            # No intercistronic region
            elif (cds_start - (prev_regioncds_end + 1)) == 0:
                cds_idx_dict[(prev_regioncds_start, prev_regioncds_end)] = \
                    startend_line_dict[(prev_cds_start, prev_cds_end)]
                region_dict[(prev_regioncds_start, prev_regioncds_end)] = "CDS"

            # Intercistronic region
            elif (cds_start - (prev_regioncds_end + 1)) <= ICR_SIZE:
                cds_idx_dict[(prev_regioncds_start, prev_regioncds_end)] = startend_line_dict[
                    (prev_cds_start, prev_cds_end)]
                region_dict[(prev_regioncds_start, prev_regioncds_end)] = "CDS"
                region_dict[((prev_regioncds_end + 1), cds_start - 1)] = "ICR"

            # Short intergenic region
            elif (cds_start - (prev_regioncds_end + 1)) <= (left_ext + right_ext):
                cds_idx_dict[(prev_regioncds_start, prev_regioncds_end)] = startend_line_dict[
                    (prev_cds_start, prev_cds_end)]
                region_dict[(prev_regioncds_start, prev_regioncds_end)] = "CDS"
                region_dict[((prev_regioncds_end + 1), (cds_start - 1))] = "UTR"

            # Large intergenic region
            else:
                # elif (cds_start - (prev_regioncds_end + 1)) >= (left_ext + right_ext):
                cds_idx_dict[(prev_regioncds_start, prev_regioncds_end)] = startend_line_dict[
                    (prev_cds_start, prev_cds_end)]
                region_dict[(prev_regioncds_start, prev_regioncds_end)] = "CDS"
                region_dict[((prev_regioncds_end + 1), (prev_regioncds_end + right_ext))] = chr_right + "UTR"
                region_dict[((prev_regioncds_end + right_ext + 1), (cds_start - left_ext - 1))] = "OTHER"
                region_dict[((cds_start - left_ext), (cds_start - 1))] = chr_left + "UTR"

            prev_regioncds_start = cds_start
            prev_regioncds_end = cds_end
            prev_cds_start, prev_cds_end = startend

        # Addition of the last CDS
        cds_idx_dict[(prev_regioncds_start, prev_regioncds_end)] = startend_line_dict[(prev_cds_start, prev_cds_end)]
        region_dict[(prev_regioncds_start, prev_regioncds_end)] = "CDS"
        if (region_end - prev_regioncds_end) > right_ext:
            region_dict[((prev_regioncds_end + 1), (prev_regioncds_end + right_ext))] = chr_right + "UTR"
            region_dict[((prev_regioncds_end + right_ext + 1), region_end)] = "OTHER"
        elif (region_end - prev_regioncds_end) > 0:
            region_dict[((prev_regioncds_end + 1), region_end)] = "UTR"

        return region_dict, cds_idx_dict

    def __genome_region_sizes(self) -> str:
        """
        Returns the sum of the sizes (in nucleotide) of each region in all contigs in the targeted genome.
        Sizes format is 5UTR;3UTR;UTR;CDS;OVL;ICR;OTHER
        :Return: String with total size of each region in the genome
        """
        region_prop = {"5UTR": 0, "3UTR": 0, "UTR": 0, "CDS": 0, "OVL": 0, "ICR": 0, "OTHER": 0}

        for seq_id in self.dict_region_locannot:
            for strand in ['+', '-']:
                for position, region in self.dict_region_locannot[seq_id][strand].items():
                    start = position[0]
                    stop = position[1]
                    region_prop[region] += (stop - start + 1)
                    if region not in ["5UTR", "3UTR", "UTR", "CDS", "OVL", "ICR", "OTHER"]:
                        raise ValueError(f"Unknown region '{region}'")

        p_5utr = region_prop["5UTR"]
        p_3utr = region_prop["3UTR"]
        p_utr = region_prop["UTR"]
        p_cds = region_prop["CDS"]
        p_ovl = region_prop["OVL"]
        p_icr = region_prop["ICR"]
        p_other = region_prop["OTHER"]
        return f"5UTR={p_5utr};3UTR={p_3utr};UTR={p_utr};CDS={p_cds};OVL={p_ovl};ICR={p_icr};OTHER={p_other}"

    def get_nearest_cds(self, seq_id: str, strand: chr, target_position: int) -> str:
        """
        Returns the line of the nearest cds from the query_position and on the targeted strand.
        :param seq_id: Subject sequence id
        :param strand: Subject strand
        :param target_position: Position of the query on the subject to consider
        :return: Nearest CDS informations from the target position
        """
        # Initialization
        min_distance = 1_000_000_000  # arbitrary distance value
        nearest_cds_line = None

        for start_end, cds_line in self.dict_region_cdslines[seq_id][strand].items():
            cds_start, cds_end = start_end

            # Searches for the nearest CDS
            distance = min(abs(cds_start - target_position), abs(cds_end - target_position))
            if min(cds_start, cds_end) < target_position < max(cds_start, cds_end):
                return cds_line
            elif distance < min_distance:
                min_distance = distance
                nearest_cds_line = cds_line

        return nearest_cds_line


def get_region(target_start: int, target_end: int, region_dict: dict) -> str:
    """
    Parses the region dictionary and returns in which region(s) the start and stop positions of the target match.
    :param target_start: Start position of the target on the subject sequence
    :param target_end: Stop position of the target on the subject sequence
    :param region_dict: Region coordinates dictionary of a given .gff file
    :return: Return the annoted region(s) of this coordinates (separated by a '-', if multiple targeted regions)
    """
    annotation = ""

    for position, region in region_dict.items():
        start = position[0]
        stop = position[1]

        if (start - 1) < target_start < (stop + 1) and (start - 1) < target_end < (stop + 1):
            return region
        elif (start - 1) < target_start < (stop + 1):
            annotation = region + "-" + annotation
        elif (start - 1) < target_end < (stop + 1):
            annotation += region

    if not annotation:
        raise ValueError(f"No region found for target position ({target_start}, {target_end})")

    if annotation[-1] == "-":
        annotation = annotation[:-1]

    return REGION_CONVERT[annotation]


def main():
    # Initialization
    start_time = time.time()
    input_lines = parse_file_lines(INPUT, has_header=HEADER)
    input_line_number = len(input_lines)
    output_lines = []
    not_nearest_cds = 0

    # Display parameters in log file
    if VERBOSE:
        log_file = open_file(LOG_FILE, "at")
        log_file.write(f"#{'-' * 50}\n")
        time_now = datetime.now().strftime("%Y-%m-%d at %H:%M:%S")
        log_file.write(f"#The region annotation from the GFF files was added the {time_now}\n\n")
        log_file.write(f"#Parameters\n"
                       f"Input            : {INPUT}\n"
                       f"Output           : {OUTPUT}\n"
                       f"Subject GFF      : {GFF}\n"
                       f"Subject fasta    : {FASTA}\n"
                       f"Subject seq id   : column {SSEQID}\n"
                       f"Subject start    : column {SSTART}\n"
                       f"Subject end      : column {SEND}\n"
                       f"Strand           : column {STRAND}\n"
                       f"5'UTR size       : {UTR5_SIZE} nt\n"
                       f"3'UTR size       : {UTR3_SIZE} nt\n"
                       f"ICR size         : {ICR_SIZE} nt\n"
                       f"FS               : '{FS}'\n"
                       f"Header           : {HEADER}\n\n")
        log_file.close()

    # Parses the GFF and FASTA files to get the predicted regions
    gffasdict = GFFasDict(GFF, FASTA)

    # Fetches all region sizes (CDS, 5UTR...) in all sequences of the genome (on both strands)
    region_sizes = gffasdict.region_sizes

    # Parses header
    if HEADER:
        header = f"{get_header(INPUT)[:-1]}{FS}region{FS}cds_start{FS}cds_end{FS}cds_id\n"
        output_lines.append(header)

    # Processing : Iteration for each alignment in INPUT
    for line in input_lines:
        # Gets some alignment information
        line_split = line[:-1].split(sep=FS)
        sseqid = line_split[SSEQID - 1]
        target_start = int(line_split[SSTART - 1])
        target_end = int(line_split[SEND - 1])
        strand = line_split[STRAND - 1]

        # Searches for the nearest cds
        target_center = (target_start + target_end) // 2
        nearest_cds = gffasdict.get_nearest_cds(sseqid, strand, target_center)

        # Searches for regions corresponding to the start and stop positions on the subject
        region_dict = gffasdict.dict_region_locannot[sseqid][strand]
        region = get_region(target_start, target_end, region_dict)

        # Write the result in the output file OUTPUT_STD
        if not nearest_cds:
            cds_start = cds_stop = ""
            cds_id = ""
            not_nearest_cds += 1
        else:
            _, _, _, cds_start, cds_stop, _, _, _, cds_attributes = nearest_cds[:-1].split(sep="\t")
            cds_id = parse_cds_attributes(cds_attributes)

        output_line = f"{line[:-1]}{FS}{region}{FS}{cds_start}{FS}{cds_stop}{FS}{cds_id}\n"
        output_lines.append(output_line)

    # Ending : writes the output
    if OUTPUT:
        output_file = open_file(OUTPUT, 'wt')
        for line in output_lines:
            output_file.write(line)
        output_file.close()
    else:
        for line in output_lines:
            print(line, end="")

    if VERBOSE:
        log_file = open_file(LOG_FILE, "at")
        log_file.write(f"#Results\n"
                       f"CDS found for    : {input_line_number - not_nearest_cds}/{input_line_number} alignments\n"
                       f"Region sizes     : {region_sizes}\n"
                       f"Execution time   : {round(time.time() - start_time, 2)} sec\n\n")
        log_file.close()


def parse_cds_attributes(cds_attributes: str) -> str:
    """
    Parses the GFF line attributes to fetches the CDS ID.
    :param cds_attributes: Attributes of the line (column 9 of the GFF)
    :retrun: Return the CDS id of the attributes
    """
    attributes_as_list = cds_attributes.split(";")
    return attributes_as_list[0].split("=")[1]


def get_header(file_path: str) -> str:
    """
    Returns the first line (header) of a file.
    :param file_path: Path to the file
    :return: File header
    """
    file = open_file(file_path, "rt")
    header = file.readline()
    file.close()
    return header


def open_file(file: str, mode="rt"):
    """
    Opens a file according to its extension (gzip or not).
    :param file: File to open
    :param mode: Character/String that specifies how the file is opened
    :return: I/O object
    """
    if file.split(".")[-1] == "gz":
        return gzip.open(file, mode)
    else:
        return open(file, mode)


def parse_file_lines(file_path: str, has_header=False) -> list:
    """
    Fetches all lines from a given file.
    :param file_path: Path to the file
    :param has_header: Boolean : skip the first line if the file has header
    :return: All lines in a list
    """
    file = open_file(file_path, "rt")
    file_lines = file.readlines()
    file.close()
    if has_header:
        return file_lines[1:]
    else:
        return file_lines


def usage() -> None:
    """ Display the usage of the program in the terminal. """
    print("\nUsage :\n"
          "  python3 align2cdsRegions.py [arguments]\n"
          "\n"
          "Mandatory arguments :\n"
          "  -i, --input                 path to the input alignment file\n"
          "  -g, --gff                   path to the gene features file of the subject ('.gz' file allowed)\n"
          "  -f, --fasta                 path to the fasta sequences file of the subject  ('.gz' file allowed)\n"
          "  -s, --sseqid                [int] column number of the subject sequence id\n"
          "  -a, --sstart                [int] column number of start position in subject (sstart must be lower than "
          "send)\n"
          "  -e, --send                  [int] column number of end position in subject (send must be greater than "
          "sstart)\n"
          "  -t, --strand                [int] column number of the targeted strand on subject\n"
          "\n"
          "Optional arguments :\n"
          "  -5, --utr5_size             [int] size of the 5' flanking region of the CDS to consider (default is 20 "
          "nt)\n"
          "  -3, --utr3_size             [int] size of the 3' flanking region of the CDS to consider (default is 150 "
          "nt)\n"
          "  -r, --icr_size              [int] size of the intercistronic region between CDS (default is 5 nt)\n"
          "  -o, --output                path to write the output\n"
          "  -d, --delimiter             field separator of the input file (default is '\\t')\n"
          "  -l, --has_header            indicates that the input file has a header\n"
          "  -v, --verbose               write in a log file the program parameters, the number of nearest CDS found "
          "and the total sizes of each region (on both strands)\n"
          "  -F, --force                 delete the output file if it exists\n"
          "  -h, --help                  how to use the program\n"
          "\n"
          "Example:\n"
          "  python3 align2cdsRegions.py -i Example/align_MH-DSM.tsv -o align_MH-DSM_regions.tsv -g "
          "Example/MH-DSM.gff.gz -f Example/MH-DSM.fna.gz -s 2 -a 5 -e 6 -t 8 -l -v -F\n")


if __name__ == "__main__":
    # Fetches input arguments
    try:
        opts, _ = getopt.getopt(sys.argv[1:],
                                'i:o:g:f:d:s:a:e:t:5:3:r:lvFh',
                                ['input=', 'output=', 'gff=', 'fasta=', 'delimiter=', 'sseqid=', 'sstart=', 'send=',
                                 'strand=', 'utr5_size=', 'utr3_size=', 'icr_size=', 'has_header', 'verbose', 'force',
                                 'help'])
    except getopt.GetoptError as err:
        print_error(f"{str(err)[0].upper() + str(err)[1:]}")
        usage()
        sys.exit(2)

    for option, arg in opts:
        if option == '-i' or option == '--input':
            INPUT = arg
        elif option == '-o' or option == '--output':
            OUTPUT = arg
        elif option == '-g' or option == '--gff':
            GFF = arg
        elif option == '-f' or option == '--fasta':
            FASTA = arg
        elif option == '-d' or option == '--delimiter':
            FS = arg
        elif option == '-s' or option == '--sseqid':
            SSEQID = int(arg)
        elif option == '-a' or option == '--sstart':
            SSTART = int(arg)
        elif option == '-e' or option == '--send':
            SEND = int(arg)
        elif option == '-t' or option == '--strand':
            STRAND = int(arg)
        elif option == '-5' or option == '--utr5_size':
            UTR5_SIZE = int(arg)
        elif option == '-3' or option == '--utr3_size':
            UTR3_SIZE = int(arg)
        elif option == '-r' or option == '--icr_size':
            ICR_SIZE = int(arg)
        elif option == '-l' or option == '--has_header':
            HEADER = True
        elif option == '-v' or option == '--verbose':
            VERBOSE = True
        elif option == '-F' or option == '--force':
            FORCE = True
        elif option == '-h' or option == '--help':
            usage()
            sys.exit(0)

    # Checks the number of arguments and the presence of mandatory arguments
    for arg in [INPUT, GFF, FASTA, SSEQID, SSTART, SEND, STRAND]:
        if (not arg) or (len(opts) == 0):
            print("\033[31mError : One or more mandatories arguments are missing\033[0m")
            usage()
            sys.exit(2)

    # Checks if the output files already exist
    if OUTPUT and os.path.exists(OUTPUT):
        if FORCE:
            os.remove(OUTPUT)
        else:
            raise FileExistsError(f"File '{OUTPUT}' already exist.")

    # Run the program
    main()
