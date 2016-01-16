#pragma once // the same purpose as #include guards

#include <unordered_map>
#include "util.h"

typedef struct SNP {
	SNP(char snp_type_ = 'S', string ref_ = "", string alt_ = "") :
		snp_type(snp_type_), ref(ref_), alt(alt_){}
	char snp_type;
	string ref;
	string alt;
}SNP;

typedef unordered_map<int, SNP> SnpHash;
typedef unordered_map<int, string> VCFEntryHash;

class VCF
{
public:
	VCF();
	~VCF();

	SnpHash refpos_2_snp;
	VCFEntryHash refpos_2_vcf_entry;
	SnpHash querypos_2_snp;
	VCFEntryHash querypos_2_vcf_entry;
	void ReadVCF(string filename, SnpHash & pos_2_snp, VCFEntryHash & pos_2_vcf_entry);
};

