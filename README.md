# VarMatch
robust matching of small variant datasets using flexible scoring schemes

# Authors
- Chen Sun (The Pennsylvania State University)
- Paul Medvedev (The Pennsylvania State University)

# Release Date
### TBA


# Prerequisite
- GCC 4.7 or later for c++11 support
- Python 2.7 or later
- matplotlib*

> *matplotlib is only used for graphic visualization. you can use '-G' parameter to disable visualization function

> *matplotlib is not a prerequisite if either `-f`, `-G` or `-C` parameter is used

# Installation
### Quick Install Instruction:
You can build VarMatch from source. 
```
git clone https://github.com/medvedevgroup/varmatch.git
cd varmatch
make all
```

### Uninstall
`cd` to the directory of VarMatch
```
make clean
```

# Test Data Set

- Links to a test data set (~15M) : [https://github.com/medvedevgroup/varmatch/blob/master/test_data.txt](https://github.com/medvedevgroup/varmatch/blob/master/test_data.txt)

- Links to the data used for bencharmking in the paper: [https://github.com/medvedevgroup/varmatch/blob/master/data.txt](https://github.com/medvedevgroup/varmatch/blob/master/data.txt)

# Usage
### Quick Usage:

*compare two vcf files to match variants*

```
./varmatch -b baseline.vcf -q query.vcf -g ref.fa -o output -f
```
- `-b` baseline vcf file
- `-q` query vcf file
- `-g` genome fasta file
- `-o` output file prefix, default value is `out`
- `-f` fast mode*, equivalent to use parameters `-u 0 -m 0 -s 0 -C`

>*fast mode is suggested for ordinary analysis.

>VarMatch accept baseline and query in VCF file format (e.g. xx.vcf), it does not accept gz file (e.g. xx.vcf.gz) in current version.

>see [Results of VarMatch](#results-of-varmatch) section for intepretation of results in output directory. 

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

# Results of VarMatch

### varmatch (.match) file
You can find varmatch (.match) files in VarMatch output directory, filename is in the format of query`x`.`u`\_`m`\_`s`.match

- `x` is the id of queries.
- `u` is the value of parameter `-u`, `--score_unit`
- `m` is the value of parameter `-m`, `--match_mode`
- `s` is the value of parameter `-s`, `--score_scheme`

For instance, if you use one query VCF file and use `-f` parameter, there is query1.0_0_0.match in your output file.

varmatch file contains the information of matched VCF entries from baseline and query VCF file.

Lines in varmatch file started with `#` are comment lines. They contain general information of baseline and query VCF file, and also general information of varmatch file. The first 2 lines of .match file starts with `###`:

- `###VCF1` is the baseline VCF filename

- `###VCF2` is the query VCF filename

varmatch files contains at least 9 fields:

  1. `CHROM` field represents chromosome ID
  2. `POS` field represents genome position on reference genome sequence
  3. `REF` field represents reference allele
  4. `ALT` field represents alternative alleles, multiple alleles are separated by `/`
  5. `VCF1` field represents variants from baseline. If it is a direct match, this column is `.`. If it is not a direct match, this column contains variants separated by `;`. Each variant contains three information: reference genome position, reference allele, alternative alleles.
  6. `VCF2` field represents variants from query. If it is a direct match, this column is `.`.
  7. `PHASE1` field represents phasing information of variants from baseline. If it is a direct match, this column is `.`.
  8. `PHASE2` field represents phasing information of variants from query. If it is a direct match, this column is `.`.
  9. `SCORE` field represents the total score of variants from baseline and query, calculated based on given score unit, and score scheme.

The meaning of each line in varmatch file:

> variants in `VCF1` column is equivalent to variants in `VCF2` column. If applying them separately on `REF` sequence, which is a substring of reference genome sequence starts at position `POS` of chromosome `CHROM`, can get the same donor sequences in `ALT` column.
The phasing information of variants in `VCF1` and `VCF2` are separately in `PHASE1` and `PHASE2`. 

varmatch (.match) file format gives a standard representation of equivalent variants, especially for complex variants.

If you have any suggestions of improving .match file format, please [contact me](#contact).

### .stat file

It contains some statistical information.

# License

See [license.txt](https://github.com/medvedevgroup/varmatch/blob/master/license.txt)

# Contact

chensun@cse.psu.edu

You also can report bugs or suggest features using issue tracker at GitHub [https://github.com/medvedevgroup/varmatch](https://github.com/medvedevgroup/varmatch)

# Acknowledgements

If using VarMatch, please cite:
Sun, C., & Medvedev, P. (2016). VarMatch: robust matching of small variant datasets using flexible scoring schemes. Bioinformatics, 33(9), 1301-1308.

Corresponding BiBTex:
```
@article{sun2016varmatch,
  title={VarMatch: robust matching of small variant datasets using flexible scoring schemes},
  author={Sun, Chen and Medvedev, Paul},
  journal={Bioinformatics},
  volume={33},
  number={9},
  pages={1301--1308},
  year={2016},
  publisher={Oxford University Press}
}
```

This project has been supported in part by NSF awards DBI-1356529, CCF-1439057, IIS-1453527, and IIS-1421908.
