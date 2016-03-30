#pragma once

#include "vcf.h"

typedef struct DiploidVariant {
	DiploidVariant(int pos_ = 0,
		char var_type_ = 'S',
		string ref_ = "",
		vector<string> alts_ = {"",""},
		string genotype_ = "0/0",
		bool heterozygous_ = false,
		bool multi_alts_ = false
		int flag_ = 1) :
		pos(pos_),
		var_type(var_type_),
		ref(ref_),
		alts(alt_),
		genotype(genotype_),
		bool heterozygous(heterozygous_),
		bool multi_alts(multi_alts_),
		flag(flag_){}

	int pos;
	char snp_type;
	string ref;
	vector<string> alts;
	string genotype;
	bool heterozygous;
	bool multi_alts;
	int flag;
}DiploidVariant;

// define outside of struct, idiomatic solution for lexicographical compare for structures
bool operator <(const DiploidVariant& x, const DiploidVariant& y);

bool operator ==(const DiploidVariant& x, const DiploidVariant& y);

class DiploidVCF : public VCF
{
private:

public:
	DiploidVCF(int thread_num_);
	~DiploidVCF();

	// for public access
	void Compare(string ref_vcf,
		string query_vcf,
		string genome_seq,
		bool direct_search,
		string output_prefix,
		bool match_genotype,
		bool normalization);

};
