# align2cdsRegions

<pre>
Author  : Emmanuel Clostres
Mail    : emmanuel.clostres@univ-rennes.fr
Python  : v3.8+
Version : 1.0 
</pre>

## Description
The program adds to the input alignments, the region where they match, and the cds information nearest to them. These regions can be a 5' flanking region (5UTR), a 3' flanking region (3UTR), a CDS region (CDS), a small intergenic region (UTR), a intercistronic region (ICR), an overlapping CDS region (OVL), or a region located elsewhere on the sequence (OTHER).

The program starts by parsing the Gene Features File (GFF) and genome file (FASTA) which respectively contain the CDS and the length informations of the subject. It then transforms the subject sequence into a region dictionary with the region coordinates (start, end) as key and the subject region name (CDS, UTR, etc...) as value.

For each alignment, the program returns the CDS nearest to the query, using the coordinates of the subject center. Thus, for each alignment, a nearest CDS is added, except if no CDS is annotated on the subject strand on which the query matches.

An OTHER annotation is indicated when the region where the query matches is neither a CDS, nor a UTR, nor an ICR, nor an OVL. The query matches may be far from an annotated CDS or may be on a CDS located on the other strand.

If verbose, the program returns in a log file the total sizes of each region (on both strands).

Note that in the input alignment file, the subject end position must be greater than the subject start position.

## Dependencies
- getopt
- gzip
- os
- sys
- time
- datetime
- SeqIO

## Usage
```
python3 align2cdsRegions.py [arguments]

Mandatory arguments :
  -i, --input                 path to the input alignment file
  -g, --gff                   path to the gene features file of the subject ('.gz' file allowed)
  -f, --fasta                 path to the fasta sequences file of the subject  ('.gz' file allowed)
  -s, --sseqid                [int] column number of the subject sequence id
  -a, --sstart                [int] column number of start position in subject (sstart must be lower than send)
  -e, --send                  [int] column number of end position in subject (send must be greater than sstart)
  -t, --strand                [int] column number of the targeted strand on subject

Optional arguments :
  -5, --utr5_size             [int] size of the 5' flanking region of the CDS to consider (default is 20 nt)
  -3, --utr3_size             [int] size of the 3' flanking region of the CDS to consider (default is 150 nt)
  -r, --icr_size              [int] size of the intercistronic region between CDS (default is 5 nt)
  -o, --output                path to write the output
  -d, --delimiter             field separator of the input file (default is '\t')
  -l, --has_header            indicates that the input file has a header
  -v, --verbose               write in a log file the program parameters, the number of nearest CDS found and the total sizes of each region (on both strands)
  -F, --force                 delete the output file if it exists
  -h, --help                  how to use the program
```
  
## Example
`python3 align2cdsRegions.py -i Example/align_MH-DSM.tsv -o align_MH-DSM_regions.tsv -g Example/MH-DSM.gff.gz -f Example/MH-DSM.fna.gz -s 2 -a 5 -e 6 -t 8 -l -v -F
`

## License
Free and unrestricted use.
