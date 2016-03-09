# VarMatch
matching equivalent genetic variants with different complex representations

# Prerequisite
- python 2.7 or later
- GCC 4.7 or later for c++11 support

# Installation
**Quick Install Instruction:**
You can build VarMatch from source. 
```
git clone --recursive https://github.com/medvedevgroup/varmatch.git
cd varmatch
make all
```
:warning: Use recursive clone(`--recursive`) as VarMatch dependes on submodules.

**Install Without Normalization:**
VarMatch use `vt normalize` to normalize variants, but normalization is not necessary. If you do not want normalization, use following command to build:
```
git clone https://github.com/medvedevgroup/varmatch.git
cd varmatch
make vm
```
And when using VarMatch, use `-n` parameter to indicate that you do not need normalization.

# Usage
:one: **Quick Usage:**

*illustrate common parameters for quick usage of VarMatch*

```
./varmatch -r one.vcf -q another.vcf -g chr1.fa -o ./output
```
- -r reference vcf file
- -q query vcf file
- -g genome fasta file
- -o output directory, default value is `./output`

:two: **Compare two vcf files with multi chromosomes**

*illustrate how to use --multi_genome parameter*

*This parameter is used if vcf file contains variants of multi chromosome*

```
./varmatch --multi_genome=chromosome_list.txt -r one.vcf -q another.vcf
```

- --multi_genome chromosome list file contains chromosome name and FASTA file absolute path, separated by `\t`.

An example of chromosome list file is as follows:

>1&nbsp;&nbsp;&nbsp;&nbsp;/home/varmatch/human/chr1.fa

>2&nbsp;&nbsp;&nbsp;&nbsp;/home/varmatch/human/chr2.fa

>17&nbsp;&nbsp;&nbsp;/home/varmatch/human/backup/chr17.fa

>X&nbsp;&nbsp;&nbsp;&nbsp;/home/varmatch/human/chrxx.fa

>Y&nbsp;&nbsp;&nbsp;&nbsp;/home/anotherpath/human/chrY/human.y.fa

:three: **Remove duplicates in one vcf file**

*illustrate how to use --remove_dup module*

*This function module can be used to remove duplicates in database or single vcf file*

```
./varmatch -g genome.fa --remove_dup database.vcf -o ./output
```

- --remove_dup vcf file name

:four: **Using multi thread**

```
./varmatch -r one.vcf -q another.vcf -g chr1.fa -o ./output -t 8
```

- -t thread number

use `-h/--help` for detailed help message.


:five: **Detailed Usage:**
```
TBA
```

# Output
**Standard Output:**
```
######### Matching Result ################

 ref total: Number of VCF entries in reference VCF file
 que total: Number of VCF entries in query VCF file
 ref matches: Number of matched VCF entries in reference VCF file
 que matches: Number of matched VCF entries in query VCF file
 ref mismatch: Number of mismatch in reference VCF file
 alt mismatch: Number of mismatch in query VCF file

```
# VarMatch Parallel Mode
VarMatch supports parallel computing with `-t ` parameter.

Parallel computing is based on C++11 standard.
