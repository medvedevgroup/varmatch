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
**Quick Usage:**
```
./varmatch -r one.vcf -q another.vcf -g chr1.fa
```
- -r reference vcf file
- -q query vcf file
- -g genome fasta file

**Compare two vcf files with multi chromosomes**
```
./varmatch --multi_genome=chromosome_list.txt -r one.vcf -q another.vcf
```

- --multi_genome chromosome list file contains chromosome name and FASTA file absolute path, separated by TAB
A example of chromosome list file is as follows:

>1   /home/varmatch/human/chr1.fa

>2   /home/varmatch/human/chr2.fa

>17  /home/varmatch/human/backup/chr17.fa

>X   /home/varmatch/human/chrxx.fa

>Y   /home/anotherpath/human/chrY/human.y.fa

**Remove duplicates in one vcf file**

```

```

use `-h/--help` for detailed help message.

**Detailed Usage:**
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
