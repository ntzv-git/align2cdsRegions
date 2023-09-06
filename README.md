# align2cdsRegions

<pre>
Author  : Emmanuel Clostres
Mail    : emmanuel.clostres@univ-rennes.fr
Python  : v3.8+
Version : 1.2
</pre>

## Description

The program adds to the input alignments, the region where the query match on the subject, as well as the distance, coordinates and ID of the cds nearest to them. These regions can be a 5' flanking region (5FLR), a 3' flanking region (3FLR), a CDS region (CDS), a small inter-feature region (SIR), an overlapping region (OVL), or a region located elsewhere on the sequence (OTHER).

The program starts by parsing the Gene Features File (GFF) and genome file (FASTA) which respectively contain the CDS and the length informations of the subject. It then transforms the subject sequence into a region dictionary with the region coordinates (start, end) as key and the subject region name (CDS, FLR, OVL etc...) as value.

For each alignment, the program returns the CDS nearest to the query, using the coordinates of the subject center (aligned query). Thus, for each alignment, a nearest CDS is added, except if no CDS is annotated on the subject strand on which the query aligns.

An OTHER annotation is indicated when the region where the query matches is neither a CDS, a FLR, a SIR nor an OVL. In this case, the query alignment may be far from an annotated CDS or may be on a CDS located on the other strand.

If verbose, the program returns all the parameters and total sizes of each region (on both strands) in the log file.

### Notes

The calculation of the distance to the nearest CDS is optimized for query lengths shorter than those of the CDSs. In fact, the distance is the absolute value between the aligned center of query and the CDS start and stop positions. The "cds_dist" column can also be set to 0 if the query center is inside the annotated CDS.

The annotated region corresponds to the rounded center of the aligned query.

The tool has only been tested for the "cds" and "gene" feature types, but it can work for all GFF column 3 features (e.g. "tRNA", "rRNA" etc.).

If the input alignment file contains several genomes from different organisms, the best way to perform the analysis  (in terms of performance) is to respectively pool all the FASTA and GFF files from all the organisms concerned into a  single file (rather than running the tool one by one for each organism).

## Installation

```bash
git clone https://github.com/ntzv-git/align2cdsRegions.git
```

### Requirements

- python3.8

### Package dependencies

- getopt
- gzip
- os
- sys
- time
- SeqIO (from Bio package)

## Usage

```
python3 align2cdsRegions.py [arguments]

Mandatory arguments :
  -i, --input                 [str] path to the input alignment file ('.gz' file allowed)
  -g, --gff                   [str] path to the gene features file of the subject ('.gz' file allowed)
  -f, --fasta                 [str] path to the fasta sequences file of the subject ('.gz' file allowed)
  -s, --sseqid                [int] column number of the subject sequence id
  -a, --sstart                [int] column number of start position in subject (sstart must be lower than send)
  -e, --send                  [int] column number of end position in subject (send must be greater than sstart)
  -t, --strand                [int] column number of the targeted strand on subject (column values must contain '+' or '-' characters)

Optional arguments :
  -5, --5flr_size             [int] size of the 5' flanking region sequence of the CDS to consider (default is 20 nt)
  -3, --3flr_size             [int] size of the 3' flanking region sequence of the CDS to consider (default is 150 nt)
  -x, --feature_type          [str] type of the feature to search (column 3) in the GFF (default is 'cds')
  -o, --output                [str] path to write the output
  -d, --delimiter             [chr] field separator of the input file (default is '\t')
  -l, --has_header            indicates that the input file has a first line header (the program will report it in the output file)
  -v, --verbose               write in a log file the program parameters, the number of nearest CDS found and the total sizes of each region (on both strands)
  -F, --force                 delete the output file if it exists
  -h, --help                  how to use the program
```
  
### Example

In this example we want to know where the small non-coding RNA (sRNA) "sRNA-MIMAT0000177" hybridizes in the genome of _E. coli_ K12 ([GCF_000005845.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000005845.2/)).
We first perfom a quick analyse by aligned (e.g. blast) our sRNA with the genome of _E. coli_ K12 (without forgiven to add the subject strand information in the output tabular).
We then use align2cdsRegions.py to get the nearness of the sARN to a CDS/gene using the following command :

``` bash
python3 align2cdsRegions.py -i Example/align_coli-K12.tsv -o align_coli-K12_regions.tsv -g Example/ecoli-K12.gff.gz -f Example/ecoli-K12.fna.gz -s 2 -a 9 -e 10 -t 13 -l -v -F
```

Concerning the parameters : "-i align_coli-K12.tsv" is the input file (corresponding to a modified blast output file in which the "plus" and "minus" strands have been replaced by "+" and "-" respectively) ; "-o align_coli-K12_regions.tsv" is the output file ; "-s 2" indicates that the target sequence id is in column 2 ; "-a 9" indicates that the start of the aligned position of the query is in column 9 ; "-e 10" indicates that the end of the aligned position of the query is in column 10 ; "-t 13" indicates that the aligned strand is in column 13 ; "-l" indicates that the input alignment file has a first line header ; "-F" indicates to delete the output file if it already exists.

## License

Free and unrestricted use.
