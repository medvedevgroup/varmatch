# VCF-Compare
VCF entries cross-validating tool
# Prerequite
python 3.0+
# Installation
Directly download the directory
# Usage
**Quick Usage:**
```
python ./vcfcompare.py -r one.vcf -q another.vcf -g chr1.fa
```
use -h/--help for detailed help message.
Note VCF file need to be normalized first to get accurate match.

**Detailed Usage:**
```
usage: vcfcompare.py [-h] -r REFERENCE -q QUERY -g GENOME [-p FALSE_POSITIVE]
                     [-n FALSE_NEGATIVE] [-o OUTPUT] [-d] [-c CHR] [-s STAT]

optional arguments:
  -h, --help            show this help message and exit
  -r REFERENCE, --reference REFERENCE
                        reference vcf file path, usually larger than query vcf file
  -q QUERY, --query QUERY
                        query vcf file path
  -g GENOME, --genome GENOME
                        reference genome file path, fasta file format
  -p FALSE_POSITIVE, --false_positive FALSE_POSITIVE
                        false positive, i.e. mismatch vcf entries in query vcf file, default=false_positive.vcf
  -n FALSE_NEGATIVE, --false_negative FALSE_NEGATIVE
                        false negative, i.e. mismatch vcf entries in reference vcf file, default=false_negative.vcf
  -o OUTPUT, --output OUTPUT
                        output matched variants in stage 2 and 3, default=multi_match.out
  -d, --direct_search   if activate, only perform stage 1, default=not activate
  -c CHR, --chr CHR     chromosome name or id, used for parallel multi genome analysis
  -s STAT, --stat STAT  append statistics result into a file, useful for parallel multi genome analysis

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
