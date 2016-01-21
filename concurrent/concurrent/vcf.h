#pragma once // the same purpose as #include guards

#include <unordered_map>
#include <algorithm>
#include <vector>
#include <map>
#include <iostream>
#include <tuple>
#include <chrono>
#include "util.h"
#include "threadguard.h"
using namespace std;


typedef struct SNP {
    SNP(int pos_ = 0, char snp_type_ = 'S', string ref_ = "", string alt_ = "") :
		pos(pos_), snp_type(snp_type_), ref(ref_), alt(alt_){}

	int pos;
	char snp_type;
	string ref;
	string alt;
}SNP;

// define outside of struct, idiomatic solution for lexicographical compare for structures
bool operator <(const SNP& x, const SNP& y);

bool operator ==(const SNP& x, const SNP& y);

typedef vector<unordered_map<int, vector<SNP> > > SnpHash;
typedef unordered_map<int, string> VCFEntryHash;

typedef vector<map<int, vector<SNP> > > SnpMap;

class VCF
{
private:
	int thread_num;
	vector<int> pos_boundries; // boundries for split multi hash table
	string genome_sequence; // genome sequence from fasta file
	bool boundries_decided; // before deciding boundries, can not read vcf file, because do not know how to split
    bool complex_search;

	void ReadVCF(string filename, SnpHash & pos_2_snps, VCFEntryHash & pos_2_vcf_entry);
	void DirectSearchInThread(unordered_map<int, vector<SNP> > & ref_snps, unordered_map<int, vector<SNP> > & query_snps);
	void ComplexSearchInThread(map<int, vector<SNP> > & ref_snps, map<int, vector<SNP> > & query_snps);
	bool CompareSnps(SNP r, SNP q);
	
	//template function can only be defined in head file
	template <typename T>
	vector<vector<T>> CreateCombinations(vector<T> dict, int k) {
		vector<vector<T>> result;
		int n = dict.size();
		vector<bool> v(n);
		fill(v.begin(), v.end() - n + k, true);

		do {
			vector<T> t;
			for (int i = 0; i < n; ++i) {
				if (v[i]) t.push_back(dict[i]);
			}
			result.push_back(t);
		} while (prev_permutation(v.begin(), v.end()));
		return result;
	}

	bool ComplexMatch(SNP s, vector<SNP> comb);
	bool GreedyComplexMatch(SNP r_snp, map<int, vector<SNP> > & query_snps, vector<SNP> & deleted_ref_snps, vector<SNP> & deleted_que_snps);
	bool ExponentialComplexMatch(SNP r_snp, map<int, vector<SNP> > & query_snps, vector<SNP> & deleted_ref_snps, vector<SNP> & deleted_que_snps);

	string ModifySequenceBySnp(string sequence, SNP s, int offset);
	string ModifySequenceBySnpList(string sequence, vector<SNP> s, int offset);
	void FindVariantsInRange(int start, int end, map<int, vector<SNP> > snp_map, vector<SNP> & candidate_query_list);
	unsigned int EditDistance(const std::string& s1, const std::string& s2);
public:
	VCF(int thread_num_ = 0);
	~VCF();

	SnpHash refpos_2_snp;
	VCFEntryHash refpos_2_vcf_entry;
	SnpHash querypos_2_snp;
	VCFEntryHash querypos_2_vcf_entry;

	SnpMap refpos_snp_map;
	SnpMap querypos_snp_map;
	
	void ReadRefVCF(string filename);
	void ReadQueryVCF(string filename);
	void ReadGenomeSequence(string filename);
	void DecideBoundries();
	
	void DirectSearchMultiThread();
	void ComplexSearchMultiThread();
	void ClusteringSnps();
	void ClusteringSearch();

	int GetRefSnpNumber();
	int GetQuerySnpNumber();
};

