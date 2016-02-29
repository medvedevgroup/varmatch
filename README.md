# VarMatch
matching equivalent genetic variants with different complex representations

# Prerequite
python 2.7 or later
GCC 4.7 or later for c++11 support

# Installation
You can build VarMatch from source. Use recursive clone for VarMatch submodules.
```
git clone --recursive https://github.com/medvedevgroup/varmatch.git
cd varmatch
make all
```

VarMatch use vt normalize to normalize variants, but normalization is not necessary. If you do not want normalization, use following command to build:
```
git clone https://github.com/medvedevgroup/varmatch.git
cd varmatch
make varmatch
```
And when using VarMatch, use `-n` parameter to indicate that you do not need normalization.

# Usage
**Quick Usage:**
```
varmatch -r one.vcf -q another.vcf -g chr1.fa
```
-r reference vcf file name
-q query vcf file name
-g genome fasta file name

use -h/--help for detailed help message.

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
Parallel computing is based on C++11 Concurrency standard.