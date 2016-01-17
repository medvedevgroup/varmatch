#pragma once // the same purpose as #include guards

#define DEBUG
#include <unordered_map>
#include "util.h"
#include "threadguard.h"

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
private:
	void ReadVCF(string filename, SnpHash & pos_2_snp, VCFEntryHash & pos_2_vcf_entry);
public:
	VCF();
	~VCF();

	SnpHash refpos_2_snp;
	VCFEntryHash refpos_2_vcf_entry;
	SnpHash querypos_2_snp;
	VCFEntryHash querypos_2_vcf_entry;
	
	void ReadRefVCF(string filename);
	void ReadQueryVCF(string filename);
	void DirectSearch();
};

