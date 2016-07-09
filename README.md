# VarMatch
robust matching of small variant datasets using flexible scoring schemes

# Authors
- Chen Sun (The Pennsylvania State University)
- Paul Medvedev (The Pennsylvania State University)

# Release Date
### TBA
Any questions about VarMatch, please email to chensun at cse dot psu dot edu.

If you identify a bug in VarMatch, please either reported on 'github Issues' of VarMatch, or email directly to chensun at cse dot psu dot edu.



# Prerequisite
- GCC 4.7 or later for c++11 support
- Python 2.7 or later
- matplotlib*

> *matplotlib is only used for graphic visualization. you can use '-G' parameter to disable visualization function

> *matplotlib is not a prerequisite if either `-f`, `-G` or `-C` parameter is used

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
./varmatch -b baseline.vcf -q query.vcf -g ref.fa -o out -f
```
- `-b` baseline vcf file
- `-q` query vcf file
- `-g` genome fasta file
- `-o` output file prefix, default value is `out`
- `-f` fast mode*, equivalent to use parameters `-u 0 -m 0 -s 0 -C`

>*fast mode is suggested for ordinary analysis

### Detail Usage

```
./varmatch  -g <file> -b <file> -q <file> [-o <string>] [-t <int>] [-u <0|1>]
     [-m <0|1>] [-s <0|1|2|3>] [-h] [-G] [-C] [-f]
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


    `-G`, `--no_graph`
          disable graphic module
    `-C`, `--disable_curves`
          disable Precision-Recall curves, if use -G or --no_graph, then
          automatically disable these curves
    `-f`, `--fast_mode`
          In this mode, automatically disable graphic module and precision-
          recall curves, only performs one matching criterion.
           Fast mode is equivalent to use following parameters compulsively: -G
          -u 0 -m 0 -s 0


### Help Information:

use `-h/--help` for detailed help message.

