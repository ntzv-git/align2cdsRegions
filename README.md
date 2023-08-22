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

Notes :
 - in the input alignment file, the subject end position must be greater than the subject start position.
 - the calculation of the distance to the nearest CDS is optimized for query lengths shorter than those of the CDSs. In fact, the distance is the absolute value between the aligned center of query and the CDS start and stop positions. The "cds_dist" column can also be set to 0 if the query center is inside the annotated CDS, or -1 if no CDS is annotated for the subject in the GFF.
 - if the input alignment file contains several genomes from different organisms, the best way to perform the analysis  (in terms of performance) is to respectively pool all the FASTA and GFF files from all the organisms concerned into a  single file (rather than running the tool one by one for each organism).

## Dependencies

- getopt
- gzip
- os
- sys
- time
- datetime (from datetime pakage)
- SeqIO (from Bio pakage)

## Usage

```
python3 align2cdsRegions.py [arguments]

Mandatory arguments :
  -i, --input                 path to the input alignment file ('.gz' file allowed)
  -g, --gff                   path to the gene features file of the subject ('.gz' file allowed)
  -f, --fasta                 path to the fasta sequences file of the subject ('.gz' file allowed)
  -s, --sseqid                [int] column number of the subject sequence id
  -a, --sstart                [int] column number of start position in subject (sstart must be lower than send)
  -e, --send                  [int] column number of end position in subject (send must be greater than sstart)
  -t, --strand                [int] column number of the targeted strand on subject

Optional arguments :
  -5, --5flr_size             [int] size of the 5' flanking region sequence of the CDS to consider (default is 20 nt)
  -3, --3flr_size             [int] size of the 3' flanking region sequence of the CDS to consider (default is 150 nt)
  -o, --output                path to write the output
  -d, --delimiter             field separator of the input file (default is '\t')
  -l, --has_header            indicates that the input file has a header
  -v, --verbose               write in a log file the program parameters, the number of nearest CDS found and the total sizes of each region (on both strands)
  -F, --force                 delete the output file if it exists
  -h, --help                  how to use the program
```
  
## Example

``` bash
python3 align2cdsRegions.py -i Example/align_MH-DSM.tsv -o align_MH-DSM_regions.tsv -g Example/MH-DSM.gff.gz -f Example/MH-DSM.fna.gz -s 2 -a 5 -e 6 -t 8 -l -v -F
```

## License

Free and unrestricted use.
