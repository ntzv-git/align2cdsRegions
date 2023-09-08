#! /bin/python3
# -*- coding: utf-8 -*-

"""
The program adds to the input alignments (e.g. blast), the region where the query match on the subject, as well as the distance,
coordinates and ID of the cds nearest to them. These regions can be a 5' flanking region (5FLR), a 3' flanking region
(3FLR), a CDS region (CDS), a small inter-feature region (SIR), an overlapping region (OVL), or a region located
elsewhere on the sequence (OTHER).
The program starts by parsing the Gene Features File (GFF) and genome file (FASTA) which respectively contain the CDS
and the length informations of the subject. It then transforms the subject sequence into a region dictionary with the
region coordinates (start, end) as key and the subject region name (CDS, FLR, OVL etc...) as value.
For each alignment, the program returns the CDS nearest to the query, using the coordinates of the subject center
(aligned query). Thus, for each alignment, a nearest CDS is added, except if no CDS is annotated on the subject strand
on which the query aligns.
An OTHER annotation is indicated when the region where the query matches is neither a CDS, a FLR, a SIR nor an OVL. In
this case, the query alignment may be far from an annotated CDS or may be on a CDS located on the other strand.
If verbose, the program returns all the parameters and total sizes of each region (on both strands) in the log file.

Note :
 - The calculation of the distance to the nearest CDS is optimized for query lengths shorter than those of the CDSs. In
 fact, the distance is the absolute value between the aligned center of query and the CDS start and stop positions. The
 "cds_dist" column can also be set to 0 if the query center is inside the annotated CDS.
 - The annotated region corresponds to the rounded center of the aligned query.
 - The tool has only been tested for the "cds" and "gene" feature types, but it can work for all GFF column 3 features
 (e.g. "tRNA", "rRNA" etc.)
 - If the input alignment file contains several genomes from different organisms, the best way to perform the analysis
 (in terms of performance) is to respectively pool all the FASTA and GFF files from all the organisms concerned into a
 single file (rather than running the tool one by one for each organism).

Author  : Emmanuel Clostres
Mail    : emmanuel.clostres@univ-rennes.fr
Python  : v3.8
Date    : 07/09/2023

New :
- Replacement of 5' and 3'UTR by 5' and 3'FLR (flanking region) respectively ; and UTR by SIR (short inter-feature
region)
- Removes intercistronic region (ICR) from possible annotated regions
- Adds the distance between the query center and the nearest CDS in a new column of the output
- Adds the possibility to perform analyses either on CDS, genes or other feature types (column 3 in the GFF)
- Improvement of region annotation
- The annotated region corresponds to the rounded center of the aligned query
- Replacement of the example by E. coli K12
- Readme update: example and program workflow diagram added
"""

import getopt
import gzip
import os
import sys
import time

from Bio import SeqIO

# Version
VERSION = "1.2"

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
FLR5_SIZE = 20
FLR3_SIZE = 150
FEATURE_TYPE = "CDS"
LOG_FILE = f"align2cdsRegions.log"


class GFFasDict:
    """
    Converts a GFF file into an object with feature (CDS/gene) and region informations.
    """
    def __init__(self, gff_path: str, fasta_path: str):
        list_gfflines = parse_file_lines(gff_path)
        dict_seqid_len = self.__parse_fasta(fasta_path)
        self.dict_features = self.__parse_gff(list_gfflines, dict_seqid_len)
        self.dict_locregions = self.__parse_region(self.dict_features)
        self.region_sizes = self.__genome_region_sizes(self.dict_locregions)

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

    def __parse_gff(self, list_gfflines: list, dict_seqid_len: dict) -> dict:
        """
        Parses all the lines of the GFF file to get for each sequence id, its length and the feature coordinates.
        :param list_gfflines: List of all lines in the GFF file
        :param dict_seqid_len: Dictionary with sequence id as keys and sequence length as values
        :return: Dictionary with for each sequence id, its length and the feature coordinates
        """
        # Parses sequence_id and sequence_length using FASTA
        dict_seq_data = dict()  # e.g. {'seqid': {'len': seqlength, '+': {(start, end): feature_id, ...}, '-': ...}...}
        for seqid, seqlen in dict_seqid_len.items():
            dict_seq_data[seqid] = {"+": dict(), "-": dict(), "len": seqlen}
            for strand in ["+", "-"]:
                dict_seq_data[seqid][strand] = dict()

        # Parses all feature_type data using GFF
        for line in list_gfflines:
            # Pass commented lines
            if line[0] == "#":
                continue
            else:
                # Parse the feature information
                line_split = line[:-1].split(sep="\t")
                seqid = line_split[0]
                feature_type = line_split[2].upper()
                feature_start = line_split[3]
                feature_end = line_split[4]
                feature_strand = line_split[6]

                if feature_type == FEATURE_TYPE:
                    start = min((int(feature_start), int(feature_end)))
                    end = max((int(feature_start), int(feature_end)))
                    dict_seq_data[seqid][feature_strand][(start, end)] = self.parse_featureid(line[:-1])

        return dict_seq_data

    @staticmethod
    def parse_featureid(gffline: str) -> str:
        """
        Parses the GFF line to get the feature id.
        :param gffline: Line of the feature in the GFF
        :retrun: Return the feature id
        """
        attribute = gffline.split(sep="\t")[8]
        return attribute.split(";")[0].split("=")[1]

    def __parse_region(self, dict_features: dict) -> dict:
        """
        Considering all sequence ids and both strands, this method returns a dictionnary with coordinates (start, end)
        as keys and corresponding regions as values.
        :param dict_features: Dictionary with for each sequence id, its length and the feature coordinates
        :return: Return a dictionary of regions coordinates and a dictionary of feature coordinates
        """
        # Initialization
        dict_locregions = dict()  # e.g. {seq_id: {'+': {(region_start, region_end): region, ...}, '-': ...}

        # Iteration on each sequence id
        for seq_id in dict_features.keys():
            dict_locregions[seq_id] = {"+": dict(), "-": dict()}
            seq_len = dict_features[seq_id]["len"]

            # Iteration on each strand
            for strand in ["+", "-"]:
                # Group overlapping features into a list
                list_startend = list(set(dict_features[seq_id][strand].keys()))
                list_relations = self.__get_relations(list_startend)

                # Process on all overlapping feature groups
                list_ovf = list()  # ovf : overlapping features
                dict_blocks = dict()
                for block in list_relations:
                    list_ovf += block
                    # Search for all positions where the features overlaps to each other
                    list_intersects = self.__interval_intersections(block)

                    # Merge all the feature intersections (removing of overlapping redundancies)
                    list_intersects = self.__intersection_merge(list_intersects)

                    # Convert block to regions
                    block_union = self.__interval_union(block)
                    dict_blocks[block_union] = self.__ovf2regions(list_intersects, block_union)

                # Search for and add non-overlapping features to the dictionary
                list_nof = list(set(list_startend) - set(list_ovf))  # nof : non-overlapping features
                for feature in list_nof:
                    dict_blocks[feature] = {feature: FEATURE_TYPE}

                # Annotates regions flanking blocks
                dict_locregions[seq_id][strand] = self.__blocks2regions(dict_blocks, seq_len, strand)

        return dict_locregions

    def __search_pairwise_relations(self, intervals: [tuple]) -> list:
        """
        Search for pairwise overlapping relationships between intervals in a list.
        :param intervals: List of tuple (start, end)
        :return: List of relations between the tuples
        """
        # Sort intervals
        intervals.sort(key=lambda x: x[1])
        intervals.sort(key=lambda x: x[0])

        len_inter = len(intervals)
        relations = list()  # e.g. [[a, b], [a, d], [b, c], [c, f]...] with a, b ... f are tuples
        for i in range(len_inter):
            for j in range((i + 1), len_inter):
                if self.__two_interval_intersect(intervals[i], intervals[j]):
                    relations.append([intervals[i], intervals[j]])
        return relations

    def __all_children(self, parent: tuple, graph: dict) -> set:
        """
        Recursive function that searches in a graph for all the children of a given parent.
        :param parent: Tuple of (start, end) position
        :param graph: Dictionary with parent as key and children in a list as values (directed acyclic graph)
        :return: Set of all childs of the parent
        """
        if parent not in graph:
            return set()
        # return set(graph[parent] + [j for i in graph[parent] for j in self.__all_children(i, graph)])
        ll = list()
        for i in graph[parent]:
            for j in self.__all_children(i, graph):
                ll.append(j)
        return set(graph[parent] + ll)

    def __get_relations(self, intervals: [tuple]) -> list:
        """
        Extract and group interval that overlap to each other.
        :param intervals: List of intervals
        :return: List of grouped overlapping intervals
        """
        # Define relationships between overlatpping intervals
        pairwise_relations = self.__search_pairwise_relations(intervals)

        # Make a directed acyclic graph based on relationships
        overlap_graph = {}
        for parent, child in pairwise_relations:
            overlap_graph.setdefault(parent, set()).add(child)
        for parent, childs in overlap_graph.items():
            if parent in childs:
                childs.discard(parent)
            overlap_graph[parent] = list(childs)

        # Search for all roots of the graph (parents who are not children of another)
        childs = set(child for child_list in overlap_graph.values() for child in child_list)
        graph_roots = set(overlap_graph) - childs

        # Group overlapping features into a dict
        dict_relations = {parent: self.__all_children(parent, overlap_graph) for parent in graph_roots}

        return [[parent] + list(childs) for parent, childs in dict_relations.items()]

    @staticmethod
    def __ovf2regions(intervals: [tuple], startend: tuple) -> dict:
        """
        Converts a list of intervals into a dictionary of CDS-OVL regions with position as keys and annotated regions
        as values.
        :param intervals: List of non overlapping intervals
        :param startend: Tuple of the start and the end position of the interval list
        :return: Dictionnary of annotated regions
        """
        intervals.sort(key=lambda x: x[1])
        intervals.sort(key=lambda x: x[0])
        dict_regions = dict()  # e.g.  {(1, 4): 'CDS', (5, 9): 'OVL', (10, 17): 'CDS', ...}
        start, end = startend
        cursor = min(start, intervals[0][0])

        # Process on first interval
        if cursor == intervals[0][0]:
            dict_regions[intervals[0]] = "OVL"
        else:
            dict_regions[(start, (intervals[0][0] - 1))] = FEATURE_TYPE
            dict_regions[intervals[0]] = "OVL"
        cursor = (intervals[0][1] + 1)
        # Process on next intervals
        for interval in intervals[1:]:
            if interval[0] != cursor:
                dict_regions[(cursor, (interval[0] - 1))] = FEATURE_TYPE
            dict_regions[interval] = "OVL"
            cursor = (interval[1] + 1)
        # Add last region
        if cursor <= end:
            dict_regions[(cursor, end)] = FEATURE_TYPE

        return dict_regions

    @staticmethod
    def __blocks2regions(dict_blocks: dict, seq_len: int, strand: chr) -> dict:
        """
        Add FLR, SIR and OTHER region using dict_block keys and according the strand and 5' and 3' FLR sizes.
        :param dict_blocks: Dictionnay with block positions as keys and block content as value.
        :param seq_len: Length of of the sequence
        :param strand: Sequence strand
        :return: Dictionnary with positions as keys and annotated regions as values.
        """
        if strand == "+":
            left_ext = FLR5_SIZE
            right_ext = FLR3_SIZE
            chr_left = "5"
            chr_right = "3"
        elif strand == "-":
            left_ext = FLR3_SIZE
            right_ext = FLR5_SIZE
            chr_left = "3"
            chr_right = "5"
        else:
            raise ValueError("The strand can only be worth '+' or '-'.")

        intervals = list(dict_blocks.keys())
        intervals.sort(key=lambda x: x[1])
        intervals.sort(key=lambda x: x[0])
        dict_regions = dict()  # e.g.  {(1, 5): 'OTHER', (6, 8): '5FLR', (9, 25): 'CDS', (26, 31): '3FLR', ...}
        cursor = 1
        if len(intervals) == 0:
            return dict_regions

        # Process on first interval
        if (intervals[0][0] - cursor) < left_ext:
            # SIR
            dict_regions[(cursor, (intervals[0][0] - 1))] = "SIR"
        elif (intervals[0][0] - cursor) == left_ext:
            # FLR
            dict_regions[(cursor, (intervals[0][0] - 1))] = chr_left + "FLR"
        else:
            # OTHER-FLR
            dict_regions[(cursor, ((intervals[0][0] - 1) - left_ext))] = "OTHER"
            dict_regions[((intervals[0][0] - left_ext), (intervals[0][0] - 1))] = chr_left + "FLR"
        dict_regions.update(dict_blocks[intervals[0]])
        cursor = (intervals[0][1] + 1)

        # Process on next intervals
        for interval in intervals[1:]:
            # No inter-feature region
            if cursor == interval[0]:
                pass
            # Short inter-feature region (SIR)
            elif (interval[0] - cursor) < (left_ext + right_ext):
                dict_regions[(cursor, (interval[0] - 1))] = "SIR"
            # Large inter-feature region
            else:
                if (interval[0] - cursor) == (left_ext + right_ext):
                    # FLR-FLR
                    dict_regions[(cursor, (cursor + right_ext - 1))] = chr_right + "FLR"
                    dict_regions[((interval[0] - left_ext), (interval[0] - 1))] = chr_left + "FLR"
                else:
                    # elif (interval[0] - cursor) > (left_ext + right_ext):
                    # FLR-OTHER-FLR
                    dict_regions[(cursor, (cursor + right_ext - 1))] = chr_right + "FLR"
                    dict_regions[((cursor + right_ext), (interval[0] - left_ext - 1))] = "OTHER"
                    dict_regions[((interval[0] - left_ext), (interval[0] - 1))] = chr_left + "FLR"
            dict_regions.update(dict_blocks[interval])
            cursor = (interval[1] + 1)

        # Add last region
        if (seq_len - (cursor - 1)) < right_ext:
            # SIR
            dict_regions[(cursor, seq_len)] = "SIR"
        elif (seq_len - (cursor - 1)) == right_ext:
            # FLR
            dict_regions[(cursor, (cursor + right_ext - 1))] = chr_right + "FLR"
        else:
            # FLR-OTHER
            dict_regions[(cursor, (cursor + right_ext - 1))] = chr_right + "FLR"
            dict_regions[((cursor + right_ext), seq_len)] = "OTHER"

        return dict_regions

    @staticmethod
    def __two_interval_intersect(interval1: tuple, interval2: tuple) -> tuple:
        """
        Returns the intersection of two intervals.
        :param interval1: Tuple of two values
        :param interval2: Tuple of two values
        :return: Intersection of 2 intervals (as tuple) / None if no interval
        """
        # If no intersection exists
        if interval2[0] > interval1[1] or interval2[1] < interval1[0]:
            return None
        else:
            return max(interval1[0], interval2[0]), min(interval1[1], interval2[1])

    @staticmethod
    def __interval_union(intervals: [tuple]) -> tuple:
        """
        Returns the union of a list of intervals.
        :param intervals: Interval list
        :return: Interval union (as tuple)
        """
        union = []
        union += [i[0] for i in intervals]
        union += [i[1] for i in intervals]
        return min(union), max(union)

    def __interval_intersections(self, intervals: [tuple]) -> list:
        """
        Searches for all positions where intervals overlaps to each other.
        :param intervals: List of tuples with 2 values (start, end)
        :return: List of all intersects
        """
        intersects = list()
        len_inter = len(intervals)
        for i in range(len_inter):
            for j in range((i + 1), len_inter):
                intersect = self.__two_interval_intersect(intervals[i], intervals[j])
                if intersect:
                    intersects.append(intersect)
        return intersects

    def __intersection_merge(self, intervals: list) -> list:
        """
        Groups and perform union on all intersections of a interval list.
        :param intervals: List of intervals
        :return: Merged intervals
        """
        # Group overlapping intervals into a list
        list_relations = self.__get_relations(intervals)
        list_unions = [self.__interval_union(relations) for relations in list_relations]

        # Add non overlapping features
        list_ovf = []
        for relations in list_relations:
            list_ovf += relations
        list_nof = list(set(intervals) - set(list_ovf))

        return list_unions + list_nof

    @staticmethod
    def __genome_region_sizes(dict_locregions: dict) -> str:
        """
        Returns the sum of the sizes (in nucleotide) of each region in all contigs in the targeted genome.
        Sizes format is 5FLR;3FLR;SIR;CDS;OVL;OTHER
        :param dict_locregions: Dictionnary with positions as keys and annotated regions as values
        :Return: String with total size of each region in the genome
        """
        region_prop = {"5FLR": 0, "3FLR": 0, "SIR": 0, FEATURE_TYPE: 0, "OVL": 0, "OTHER": 0}

        for seq_id in dict_locregions:
            for strand in ['+', '-']:
                for position, region in dict_locregions[seq_id][strand].items():
                    start = position[0]
                    stop = position[1]
                    region_prop[region] += (stop - start + 1)
                    if region not in ["5FLR", "3FLR", "SIR", FEATURE_TYPE, "OVL", "OTHER"]:
                        raise ValueError(f"Unknown region '{region}'")

        p_5flr = region_prop["5FLR"]
        p_3flr = region_prop["3FLR"]
        p_sir = region_prop["SIR"]
        p_cds = region_prop[FEATURE_TYPE]
        p_ovl = region_prop["OVL"]
        p_other = region_prop["OTHER"]
        return f"5FLR={p_5flr};3FLR={p_3flr};SIR={p_sir};{FEATURE_TYPE}={p_cds};OVL={p_ovl};OTHER={p_other}"

    def get_nearest_feature(self, seq_id: str, strand: chr, target_position: int) -> (str, float):
        """
        Returns the nearest feature data from the aligned position and on the aligned strand.
        :param seq_id: Subject sequence id
        :param strand: Subject strand
        :param target_position: Center of the query on the subject to consider
        :return: Feature data nearest to target position and its corresponding distance
        """
        if len(self.dict_features[seq_id][strand]) == 0:
            # If there is no feature on this sequence
            return None, None
        else:
            # Searches for the nearest feature (CDS/gene)
            min_distance = 1_000_000_000_000  # Arbitrary distance value to start the search
            nearest_feature = None

            for start_end, feature_id in self.dict_features[seq_id][strand].items():
                start, end = start_end
                if start < target_position < end:
                    # If target inside the feature (CDS/gene)
                    return (start, end, feature_id), 0.0
                distance = min(abs(start - target_position), abs(end - target_position))
                if distance < min_distance:
                    min_distance = distance
                    nearest_feature = (start, end, feature_id)
        return nearest_feature, min_distance

    def get_region(self, seq_id: str, strand: chr, target_position: int) -> str:
        """
        Parses the region dictionary and returns in which region the target match.
        :param seq_id: Subject sequence id
        :param strand: Subject strand
        :param target_position: Center of the query on the subject to consider
        :return: Return the annoted region of this coordinate
        """
        for start_end, region in self.dict_locregions[seq_id][strand].items():
            start, end = start_end
            if start <= target_position <= end:
                return region
        if len(self.dict_locregions[seq_id][strand]) == 0:
            # If no feature on this strand
            return ""
        else:
            raise ValueError(f"No region found for '{target_position}'")


def main():
    # Initialization
    start_time = time.time()
    input_lines = parse_file_lines(INPUT, has_header=HEADER)
    input_line_number = len(input_lines)
    output_lines = []
    not_nearest_feature = 0

    # Display parameters in log file
    if VERBOSE:
        log_file = open_file(LOG_FILE, "at")
        log_file.write(f"#{'-' * 50}\n")
        time_now = time.strftime("%Y-%m-%d at %H:%M:%S")
        log_file.write(f"#Log of align2cdsRegions v{VERSION} on {time_now}\n"
                       f"\n"
                       f"#Parameters\n"
                       f"Input             : {INPUT}\n"
                       f"Output            : {OUTPUT}\n"
                       f"Subject GFF       : {GFF}\n"
                       f"Subject fasta     : {FASTA}\n"
                       f"Subject seq id    : column {SSEQID}\n"
                       f"Subject start     : column {SSTART}\n"
                       f"Subject end       : column {SEND}\n"
                       f"Strand            : column {STRAND}\n"
                       f"5'FLR size        : {FLR5_SIZE} nt\n"
                       f"3'FLR size        : {FLR3_SIZE} nt\n"
                       f"Feature type      : {FEATURE_TYPE.lower()}\n"
                       f"FS                : '{FS}'\n"
                       f"Header            : {HEADER}\n\n")
        log_file.close()

    # Parses the GFF and FASTA files to get the predicted regions
    gffasdict = GFFasDict(GFF, FASTA)

    # Fetches all region sizes (CDS, 5FLR...) in all sequences of the genome (on both strands)
    region_sizes = gffasdict.region_sizes

    # Parses header
    if HEADER:
        header = (f"{get_header(INPUT)[:-1]}{FS}region{FS}{FEATURE_TYPE.lower()}_dist{FS}{FEATURE_TYPE.lower()}_start"
                  f"{FS}{FEATURE_TYPE.lower()}_end{FS}{FEATURE_TYPE.lower()}_id\n")
        output_lines.append(header)

    # Processing : Iteration for each alignment in INPUT
    for line in input_lines:
        # Gets some alignment information
        line_split = line[:-1].split(sep=FS)
        sseqid = line_split[SSEQID - 1]
        target_start = int(line_split[SSTART - 1])
        target_end = int(line_split[SEND - 1])
        strand = line_split[STRAND - 1]

        # Searches for the nearest feature (cds/gene) and its distance
        target_center = (target_start + target_end) / 2
        nearest_feature, nearest_distance = gffasdict.get_nearest_feature(sseqid, strand, target_center)

        # Searches for regions corresponding to the start and stop positions on the subject
        region = gffasdict.get_region(sseqid, strand, round(target_center, 0))

        # Write the result in the output file OUTPUT_STD
        if not nearest_feature:
            f_start = f_end = f_id = nearest_distance = ""
            not_nearest_feature += 1
        else:
            f_start, f_end, f_id = nearest_feature

        output_line = f"{line[:-1]}{FS}{region}{FS}{nearest_distance}{FS}{f_start}{FS}{f_end}{FS}{f_id}\n"
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
                       f"Feature found for : {input_line_number - not_nearest_feature}/{input_line_number} alignments\n"
                       f"Region sizes      : {region_sizes}\n"
                       f"Execution time    : {round(time.time() - start_time, 2)} sec\n\n")
        log_file.close()


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
    print(f"\nUsage of align2cdsRegions v{VERSION}\n"
          f"  python3 align2cdsRegions.py [arguments]\n"
          f"\n"
          f"Mandatory arguments :\n"
          f"  -i, --input                 [str] path to the input alignment file ('.gz' file allowed)\n"
          f"  -g, --gff                   [str] path to the gene features file of the subject ('.gz' file allowed)\n"
          f"  -f, --fasta                 [str] path to the fasta sequences file of the subject ('.gz' file allowed)\n"
          f"  -s, --sseqid                [int] column number of the subject sequence id\n"
          f"  -a, --sstart                [int] column number of start position in subject (sstart must be lower than "
          f"send)\n"
          f"  -e, --send                  [int] column number of end position in subject (send must be greater than "
          f"sstart)\n"
          f"  -t, --strand                [int] column number of the targeted strand on subject (column values must "
          f"contain '+' or '-' characters)\n"
          f"\n"
          f"Optional arguments :\n"
          f"  -5, --5flr_size             [int] size of the 5' flanking region sequence of the CDS to consider "
          f"(default is 20 nt)\n"
          f"  -3, --3flr_size             [int] size of the 3' flanking region sequence of the CDS to consider "
          f"(default is 150 nt)\n"
          f"  -x, --feature_type          [str] type of the feature to search (column 3) in the GFF (default is "
          f"'cds')\n"
          f"  -o, --output                [str] path to write the output\n"
          f"  -d, --delimiter             [chr] field separator of the input file (default is '\\t')\n"
          f"  -l, --has_header            indicates that the input file has a first line header (the program will "
          f"report it in the output file)\n"
          f"  -v, --verbose               write in a log file the program parameters, the number of nearest CDS found "
          f"and the total sizes of each region (on both strands)\n"
          f"  -F, --force                 delete the output file if it exists\n"
          f"  -h, --help                  how to use the program\n"
          f"\n"
          f"Example:\n"
          f"  python3 align2cdsRegions.py -i Example/align_MH-DSM.tsv -o align_MH-DSM_regions.tsv -g "
          f"Example/MH-DSM.gff.gz -f Example/MH-DSM.fna.gz -s 2 -a 5 -e 6 -t 8 -l -v -F\n")


if __name__ == "__main__":
    # Fetches input arguments
    try:
        opts, _ = getopt.getopt(sys.argv[1:],
                                'i:o:g:f:d:s:a:e:t:5:3:x:lvFh',
                                ['input=', 'output=', 'gff=', 'fasta=', 'delimiter=', 'sseqid=', 'sstart=',
                                 'send=', 'strand=', '5flr_size=', '3flr_size=', 'feature_type=', 'has_header',
                                 'verbose', 'force', 'help'])
    except getopt.GetoptError as err:
        print(f"\033[31mError : {str(err)[0].upper() + str(err)[1:]}\033[0m")
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
        elif option == '-5' or option == '--5flr_size':
            FLR5_SIZE = int(arg)
        elif option == '-3' or option == '--3flr_size':
            FLR3_SIZE = int(arg)
        elif option == '-x' or option == '--feature_type':
            FEATURE_TYPE = arg.upper()
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
