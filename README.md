# VCF-Compare
VCF entries cross-validating tool
# Prerequite
python 3.0+
# Installation
Directly download the directory
# Usage
```
python ./vcfcompare.py [-h] -r REFERENCE -q QUERY -g GENOME [-p FALSE_POSITIVE] [-n FALSE_NEGATIVE] [-t TRUE_POSITIVE] [-o OUTPUT] [-d] [-c CHR] [-s STAT]
		-r/--reference: reference VCF file path
		-q/--query: query VCF file path
		-g/--genome: genome sequence file, based on which reference and query VCF are called
```
**Quick usage**
python ./vcfcompare.py -r one.vcf -q another.vcf -g chr1.fa

use -h/--help for detailed help message.

# Output
Three numbers separated by tab in standard output:
```
Matched_Variants_Number    Matched_Entries_Number_In_Reference    Matched_Entries_Number_In_Query
```


