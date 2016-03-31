#pragma once

#include "vcf.h"

typedef struct DiploidVariant {
	DiploidVariant(int pos_ = 0,
		vector<char> var_types_ = {'S','S'},
		string ref_ = "",
		vector<string> alts_ = {"",""},
		string genotype_ = "0/0",
		bool heterozygous_ = false,
		bool multi_alts_ = false,
		int flag_ = 0) :
		pos(pos_),
		var_types(var_types_),
		ref(ref_),
		alts(alt_),
		genotype(genotype_),
		bool heterozygous(heterozygous_),
		bool multi_alts(multi_alts_),
		flag(flag_){}

	int pos;
	vector<char> var_types;
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

typedef vector<unordered_map<int, DiploidVariant > > VariantHash;
typedef vector<map<int, DiploidVariant > > VariantMap;

class DiploidVCF : public VCF
{
private:
	// data structure for direct search
	VariantHash refpos_2_var;
	VariantHash querypos_2_var;
	// private
	void VCF::ReadRefVCF(string filename) {
		ReadDiploidVCF(filename, this->refpos_2_var);
	}

	// private
	void VCF::ReadQueryVCF(string filename) {
		ReadDiploidVCF(filename, this->querypos_2_var);
	}

	void DirectSearchInThread(unordered_map<int, DiploidVariant> & ref_snps, unordered_map<int, DiploidVariant> & query_snps, int thread_index);
	void DirectSearchMultiThread();

protected:
	bool CompareSequence(string s1, string s2);
	bool scoring_basepair;

	bool ReadDiploidVCF(string filename, VariantHash & pos_2_var);
	bool NormalizeDiploidVariant(DiploidVariant & var);
	map<int, vector<DiploidVariant> > cluster_vars_map;
	bool VariantMatch(vector<DiploidVariant> & variant_list);
	void MatchWithIndel(vector<DiploidVariant> & variant_list,
		const string subsequence,
		const int offset,
		int index,
		map<int, DiploidVariant> separate_pos_var[],
		map<int, int> choices[], // 4 vectors
		int pos_boundries[],
		map<int, int> max_matches[],  // 4 vectors
		int max_score);
	bool CompareSequence(string s1, string s2);
	int CheckPrefix(const string subsequence,
		const int offset,
		map<int, DiploidVariant> separate_pos_var[],
		map<int, int> choices[],
		int pos_boundries[]);

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
