# VarMatch
matching equivalent genetic variants with different complex representations

# Prerequisite
- GCC 4.7 or later for c++11 support

# Installation
**Quick Install Instruction:**
You can build VarMatch from source. 
```
git clone https://github.com/medvedevgroup/varmatch.git
cd varmatch
make all
```

# Usage
### Quick Usage:

*compare two vcf files to match variants*

```
./vm -b baseline.vcf -q query.vcf -g chr1.fa -o out
```
- `-b` baseline vcf file
- `-q` query vcf file
- `-g` genome fasta file
- `-o` output file prefix, default value is `out`

### Detail Usage

```
./vm  -g <file> -b <file> -q <file> [-o <string>] [-t <int>] [-u <0|1>]
     [-m <0|1>] [-s <0|1|2|3>] [-h]
```

Where:

   `-g` <file>,  `--genome_sequence` <file>
     (required)  genome sequence FASTA filename

   `-b` <file>,  `--baseline` <file>
     (required)  baseline variant VCF filename

   `-q` <file>,  `--query` <file>
     (required)  query variant VCF filename

   `-o` <string>,  `--output_prefix` <string>
     output filename prefix, default is "out"

   `-t` <int>,  `--thread_num` <int>
     number of threads, default is the number of available cores.

     If larger than number of available cores or less than 1, automatically
     set to default value

   `-u` <0|1>,  `--score_unit` <0|1>
     scoring function/score unit: (Default: 0)

     0 : the score that a VCF entry contributes is 1.

     1 : the score that a VCF entry contributes is the edit distance
     between the new allele and the reference one.


   `-m` <0|1>,  `--match_mode` <0|1>
     matching mode: (Default: 0)

     0 : a set of query entries match a set of baseline entries if, for
     each entry, we can select one of the alleles such that the inferred
     sequences are identical

     1 : a set of query entries match a set of baseline entries if there
     exist a phasing of each set such that the two inferred haplotypes from
     the query are equal to the two inferred haplotypes from the
     baseline.


   `-s` <0|1|2|3>,  `--score_scheme` <0|1|2|3>
     scoring scheme: (Default: 0)

     0 : find two subsets of non-overlapping equivalent variants such that
     the score of the matched variants is maximized (Default)

     1 : find two subsets of non-overlapping equivalent variants such that
     the score of the chosen baseline variants is maximized

     2 : find a maximum scoring set of variants in the query such that each
     variant can be matched by a subset of the baseline variants

     3 : (1 to 1 direct match) find a maximum scoring set of entry pairs
     such that each entry pair contains one query and one baseline variant
     that result in the same sequence. In this scheme, different scoring
     functions and matching mode have no difference.

### Help Information:

use `-h/--help` for detailed help message.

