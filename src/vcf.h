#pragma once // the same purpose as #include guards

#include <unordered_map>
#include <algorithm>
#include <vector>
#include <map>
#include <list>
#include <tuple>
#include <chrono>
#include <cstdint>
#include <limits>
#include <thread>
#include <mutex>
#include "util.h"

#include <iostream>
using namespace std;


typedef struct SNP {
    SNP(int pos_ = 0, char snp_type_ = 'S', string ref_ = "", string alt_ = "", int flag_=1) :
		pos(pos_), snp_type(snp_type_), ref(ref_), alt(alt_), flag(flag_){}

	int pos;
	char snp_type;
	string ref;
	string alt;
	int flag;
}SNP;

// define outside of struct, idiomatic solution for lexicographical compare for structures
bool operator <(const SNP& x, const SNP& y);

bool operator ==(const SNP& x, const SNP& y);

typedef vector<unordered_map<int, vector<SNP> > > SnpHash;
typedef vector<map<int, vector<SNP> > > SnpMap;

class VCF
{
private:
    int debug_f;

    bool complex_search;
	bool clustering_search;

    void ReadVCF(string filename, SnpHash & pos_2_snps);
	void DirectSearchInThread(unordered_map<int, vector<SNP> > & ref_snps,
							unordered_map<int, vector<SNP> > & query_snps);
	void ComplexSearchInThread(map<int, vector<SNP> > & ref_snps,
							map<int, vector<SNP> > & query_snps);
	
	//template function can only be defined in head file
	template <typename T>
	vector<vector<T>> CreateCombinationsWithTarget(vector<T> dict, int k, vector<int> changes, int target) {
	//vector<vector<T>> CreateCombinations(vector<T> dict, int k) {
		vector<vector<T>> result;
		int n = dict.size();
		vector<bool> v(n);
		fill(v.begin(), v.end() - n + k, true);

		do {
			vector<T> t;
            int sum = 0;
			for (int i = 0; i < n; ++i) {
				if (v[i]){
                    t.push_back(dict[i]);
                    sum += changes[i];
                }
			}
            if(sum == target){
			    result.push_back(t);
            }
		} while (prev_permutation(v.begin(), v.end()));
		return result;
	}

	bool ComplexMatch(SNP s, vector<SNP> comb);
	bool GreedyComplexMatch(SNP r_snp,
							map<int, vector<SNP> > & query_snps,
							vector<SNP> & deleted_ref_snps,
							vector<SNP> & deleted_que_snps);
	bool ExponentialComplexMatch(SNP r_snp,
								map<int, vector<SNP> > & query_snps,
								map<int, vector<SNP> >::iterator & qit_start,
								vector<SNP> & deleted_ref_snps,
								vector<SNP> & deleted_que_snps);

	void FindVariantsInRange_NlgN(int start,
								int end,
								map<int, vector<SNP> > snp_map,
								vector<SNP> & candidate_query_list,
								vector<int> & candidate_changes);
	void FindVariantsInRange_Linear(int start,
									int end,
									map<int, vector<SNP> > snp_map,
									map<int, vector<SNP> >::iterator & qit_start,
									vector<SNP> & candidate_query_list,
									vector<int> & candidate_changes);
	unsigned int EditDistance(const std::string& s1, const std::string& s2);

    void ClusteringSearchInThread(int start, int end, int thread_index);


	void ClusteringSnpsOldAlgorithm(int threshold = 400, int lower_bound = 10);
    
    //-------------------------following can be public--------------------------
    // but for a better OO design, made them private
    string ref_mismatch_filename;
    string que_mismatch_filename;

	// data structure for direct search
	SnpHash refpos_2_snp;
	SnpHash querypos_2_snp;

	// data structure for complex search
	SnpMap refpos_snp_map;
	SnpMap querypos_snp_map;

	void ReadRefVCF(string filename);
	void ReadQueryVCF(string filename);
	void DirectSearchMultiThread();
	void ComplexSearchMultiThread();
	virtual void ClusteringSnps();
	// default value better be in declaration, or definition, but never both
	void ClusteringSearchMultiThread();

	int GetRefSnpNumber();
	int GetQuerySnpNumber();

	//---------------------------above can be public:---------------------------

protected:
    vector<int> pos_boundries; // boundries for split multi hash table
    bool boundries_decided; // before deciding boundries, can not read vcf file, because do not know how to split
    // for inherit
	int thread_num;
    string chromosome_name;
	string genome_sequence; // genome sequence from fasta file
    const static int MAX_REPEAT_LEN = 1000;

    // data structure for clustering search
    vector<SNP> data_list;
    vector<int> cluster_list;
    map<int, vector<SNP> > cluster_snps_map;

    // storing complex match results
    std::mutex complex_match_mutex;
    vector<vector<string>> complex_match_records;

    // for output
    string output_stat_filename;
    string output_simple_filename;
    string output_complex_filename;
	
	bool CompareSnps(SNP r, SNP q);
    void DecideBoundries();
	string ModifySequenceBySnp(const string sequence, SNP s, int offset);
	string ModifySequenceBySnpList(const string sequence, vector<SNP> s, int offset);
	bool CheckTandemRepeat(string sequence, int unit_threshold);
	void ReadGenomeSequence(string filename);
	
    bool MatchSnpLists(vector<SNP> & ref_snp_list,
            vector<SNP> & query_snp_list,
            vector<SNP> & mixed_list,
            const string subsequence,
            int offset,
            int thread_index);
	
    template <typename D>
	vector<vector<D>> CreateCombinations(vector<D> dict, int k) {
		vector<vector<D>> result;
		int n = dict.size();
		vector<bool> v(n);
		fill(v.begin(), v.end() - n + k, true);
		do {
			vector<D> t;
			int sum = 0;
			for (int i = 0; i < n; ++i) {
				if (v[i]) {
					t.push_back(dict[i]);
				}
			}
			result.push_back(t);
		} while (prev_permutation(v.begin(), v.end()));
		return result;
	}
    bool CheckVariantOverlap(vector<SNP> snp_list);


public:
	VCF(int thread_num_ = 0);
	~VCF();

    // for public access
	void Compare(string ref_vcf,
            string query_vcf,
            string genome_seq,
            bool direct_search,
            string output_prefix);
};

