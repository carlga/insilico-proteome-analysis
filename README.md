<img src="https://img.shields.io/badge/Language-Python-yellow.svg" style="zoom:100%;" /> <!--- <img src="https://visitor-badge.glitch.me/badge?page_id=carlga.insilico-proteome-analysis&right_color=red&left_text=Total%20Visits" alt="visitor badge"/> -->
<img src="https://github.com/carlga/insilico-proteome-analysis/blob/main/pics/ms_hist.png" width=125 align="right">

# *In silico* proteome analysis

> This repository contains a collection of Python scripts that can be used to computationally predict
> a reference proteome for a given genome sequence and optimize enzymatic digestion for MS-based approaches.
> The code was developed as a practice exercise and is provided as is for those learning and/or seeking for inspiration.
> Kindly drop a :star: if this is helpful!

![infographic](./pics/insilico-proteome-analysis.png)


## 1. *Ab initio* protein-coding gene prediction

`01_ORFfinder.py` is an Open Reading Frame (ORF) finder for prokariotic genomic sequences.
Each input file must be fasta formatted and at least one filename is required for the program to run.
A single path to a directory can be used as input if batch mode is enabled (`-b` option). By default, 
ORFs are identified in all six reading frames and the minimum length threshold is set to 300 nucleotides 
(99 amino acids). `STDOUT` is used by default for all ORF predictions, but an output file can be 
specified with `-o`. 

For demonstration purposes and to keep things simple we work with a small example sequence:

```
$ cat data/genome.fa
>ENA|HZ245980|HZ245980.1
ATGCCCCCCTACACCGTGGTGTACTTCCCCGTGAGAGGCAGATGCGCCGCCCTGAGAATGCTGCTGGCC
GACCAGGGCCAGAGCTGGAAGGAGGAGGTGGTGACCGTGGAGACCTGGCAGGAGGGCAGCCTGAAGGCC
AGCTGCCTGTACGGCCAGCTGCCCAAGTTCCAGGACGGCGACCTGACCCTGTACCAGAGCAACACCATC
CTGAGACACCTGGGCAGAACCCTGGGCCTGTACGGCAAGGACCAGCAGGAGGCCGCCCTGGTGGACATG
GTGAACGACGGCGTGGAGGACCTGAGATGCAAGTACATCAGCCTGATCTACACCAACTACGAGGCCGGC
AAGGACGACTACGTGAAGGCCCTGCCCGGCCAGCTGAAGCCCTTCGAGACCCTGCTGAGCCAGAACCAG
GGCGGCAAGACCTTCATCGTGGGCGACCAGATCAGCTTCGCCGACTACAACCTGCTGGACCTGCTGCTG
ATCCACGAGGTGCTGGCCCCCGGCTGCCTGGACGCCTTCCCCCTGCTGAGCGCCTACGTGGGCAGACTG
AGCGCCAGACCCAAGCTGAAGGCCTTCCTGGCCAGCCCCGAGTACGTGAACCTGCCCATCAACGGCAAC
GGCAAGCAGTAG
$
$ python3 01_ORFfinder.py data/genome.fa
>ENA|HZ245980|HZ245980.1 frame=+1 protein=1 loc=0:632 length=210
MPPYTVVYFPVRGRCAALRMLLADQGQSWKEEVVTVETWQEGSLKASCLYGQLPKFQDGDLTLYQSNTIL
RHLGRTLGLYGKDQQEAALVDMVNDGVEDLRCKYISLIYTNYEAGKDDYVKALPGQLKPFETLLSQNQGG
KTFIVGDQISFADYNLLDLLLIHEVLAPGCLDAFPLLSAYVGRLSARPKLKAFLASPEYVNLPINGNGKQ
>ENA|HZ245980|HZ245980.1 frame=-2 protein=2 loc=313:8 length=101
MYLHLRSSTPSFTMSTRAASCWSLPYRPRVLPRCLRMVLLWYRVRSPSWNLGSWPYRQLAFRLPSCQVST
VTTSSFQLWPWSASSILRAAHLPLTGKYTTV
$
$ python3 01_ORFfinder.py data/genome.fa -o data/proteome.fa
```

The minimum length for ORF detection can be changed by passing the `-ml` argument. For instance,
lowering the 300 nucleotides limit to 30 adds two more proteins to the output for this example.

```
$ python3 01_ORFfinder.py data/genome.fa -ml 30
>ENA|HZ245980|HZ245980.1 frame=+1 protein=1 loc=0:632 length=210
MPPYTVVYFPVRGRCAALRMLLADQGQSWKEEVVTVETWQEGSLKASCLYGQLPKFQDGDLTLYQSNTIL
RHLGRTLGLYGKDQQEAALVDMVNDGVEDLRCKYISLIYTNYEAGKDDYVKALPGQLKPFETLLSQNQGG
KTFIVGDQISFADYNLLDLLLIHEVLAPGCLDAFPLLSAYVGRLSARPKLKAFLASPEYVNLPINGNGKQ
>ENA|HZ245980|HZ245980.1 frame=-2 protein=2 loc=610:536 length=24
MGRFTYSGLARKAFSLGLALSLPT
>ENA|HZ245980|HZ245980.1 frame=-2 protein=3 loc=430:353 length=25
MKVLPPWFWLSRVSKGFSWPGRAFT
>ENA|HZ245980|HZ245980.1 frame=-2 protein=4 loc=313:8 length=101
MYLHLRSSTPSFTMSTRAASCWSLPYRPRVLPRCLRMVLLWYRVRSPSWNLGSWPYRQLAFRLPSCQVST
VTTSSFQLWPWSASSILRAAHLPLTGKYTTV
```

The tool also allows cdna output, frame selection and filtering of overlapping ORFs by 
passing the `-t`, `-f` and `-fov` arguments, respectively.

```
$ python3 01_ORFfinder.py data/genome.fa -ml 30 -f -2
>ENA|HZ245980|HZ245980.1 frame=-2 protein=1 loc=610:536 length=24
MGRFTYSGLARKAFSLGLALSLPT
>ENA|HZ245980|HZ245980.1 frame=-2 protein=2 loc=430:353 length=25
MKVLPPWFWLSRVSKGFSWPGRAFT
>ENA|HZ245980|HZ245980.1 frame=-2 protein=3 loc=313:8 length=101
MYLHLRSSTPSFTMSTRAASCWSLPYRPRVLPRCLRMVLLWYRVRSPSWNLGSWPYRQLAFRLPSCQVST
VTTSSFQLWPWSASSILRAAHLPLTGKYTTV
$
$ python3 01_ORFfinder.py data/genome.fa -ml 30 -fov length
>ENA|HZ245980|HZ245980.1 frame=+1 protein=1 loc=0:632 length=210
MPPYTVVYFPVRGRCAALRMLLADQGQSWKEEVVTVETWQEGSLKASCLYGQLPKFQDGDLTLYQSNTIL
RHLGRTLGLYGKDQQEAALVDMVNDGVEDLRCKYISLIYTNYEAGKDDYVKALPGQLKPFETLLSQNQGG
KTFIVGDQISFADYNLLDLLLIHEVLAPGCLDAFPLLSAYVGRLSARPKLKAFLASPEYVNLPINGNGKQ
$
$ python3 01_ORFfinder.py data/genome.fa -ml 30 -fov length -t cdna
>ENA|HZ245980|HZ245980.1 frame=+1 cdna=1 loc=0:632 length=633
ATGCCCCCCTACACCGTGGTGTACTTCCCCGTGAGAGGCAGATGCGCCGCCCTGAGAATGCTGCTGGCCG
ACCAGGGCCAGAGCTGGAAGGAGGAGGTGGTGACCGTGGAGACCTGGCAGGAGGGCAGCCTGAAGGCCAG
CTGCCTGTACGGCCAGCTGCCCAAGTTCCAGGACGGCGACCTGACCCTGTACCAGAGCAACACCATCCTG
AGACACCTGGGCAGAACCCTGGGCCTGTACGGCAAGGACCAGCAGGAGGCCGCCCTGGTGGACATGGTGA
ACGACGGCGTGGAGGACCTGAGATGCAAGTACATCAGCCTGATCTACACCAACTACGAGGCCGGCAAGGA
CGACTACGTGAAGGCCCTGCCCGGCCAGCTGAAGCCCTTCGAGACCCTGCTGAGCCAGAACCAGGGCGGC
AAGACCTTCATCGTGGGCGACCAGATCAGCTTCGCCGACTACAACCTGCTGGACCTGCTGCTGATCCACG
AGGTGCTGGCCCCCGGCTGCCTGGACGCCTTCCCCCTGCTGAGCGCCTACGTGGGCAGACTGAGCGCCAG
ACCCAAGCTGAAGGCCTTCCTGGCCAGCCCCGAGTACGTGAACCTGCCCATCAACGGCAACGGCAAGCAG
TAG
```

For additional usage information the help option `-h` can be passed as an argument.

```
$ python3 01_ORFfinder.py -h
usage: 01_ORFfinder.py [-h] [-b] [-o OUTPUT] [-f {+1,+2,+3,-1,-2,-3} [{+1,+2,+3,-1,-2,-3} ...]] [-ml MINIMUM_LENGTH]
                       [-a] [-fov {length,cai}] [-t {protein,cdna}] [-fmt FORMAT]
                       input [input ...]

Search open reading frames (ORFs) in a sequence.

positional arguments:
  input                 Enter file/s to analyze (for path to directory enable batch mode option.)

optional arguments:
  -h, --help            show this help message and exit
  -b, --batch           Enables batch mode for processing all '.fasta' files in a directory. Path to directory must be
                        specified as input.
  -o OUTPUT, --output OUTPUT
                        Allows output to specified file.
  -f {+1,+2,+3,-1,-2,-3} [{+1,+2,+3,-1,-2,-3} ...], --frames {+1,+2,+3,-1,-2,-3} [{+1,+2,+3,-1,-2,-3} ...]
                        Allows a customized selection of reading frames for ORF search. Default is all reading frames.
  -ml MINIMUM_LENGTH, --minimum_length MINIMUM_LENGTH
                        Allows a customized minimum length threshold for ORF detection. Default is 300 nucleotides.
  -a, --ambiguities     Enables the detection of ORFs containing ambiguous nucleotides.
  -fov {length,cai}, --filter_overlaps {length,cai}
                        Allows filtering overlapping ORFs based on the selected criteria ('length' or 'cai').
  -t {protein,cdna}, --type {protein,cdna}
                        Allows a custom ORF output type. Default is in protein sequences.
  -fmt FORMAT, --format FORMAT
                        Allows a customized line length for the ORF output. Default is 70 characters per line.
```

