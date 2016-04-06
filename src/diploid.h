#pragma once

#include "vcf.h"
// data structure for direct search
typedef struct DiploidVariant {
    DiploidVariant(int pos_ = -1,
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
        alts(alts_),
        genotype(genotype_),
        heterozygous(heterozygous_),
        multi_alts(multi_alts_),
        flag(flag_){}

    int pos;
    vector<char> var_types;
    string ref;
    vector<string> alts;
    string genotype;
    bool heterozygous;
    bool multi_alts;
    int flag; //in DiploidVariant, flag = 0 is reference, flag = 1 is query
}DiploidVariant;

// define outside of struct, idiomatic solution for lexicographical compare for structures
bool operator <(const DiploidVariant& x, const DiploidVariant& y);

bool operator ==(const DiploidVariant& x, const DiploidVariant& y);

class DiploidVCF : public VCF
{
private:

    typedef vector<unordered_map<int, DiploidVariant > > VariantHash;
    typedef vector<map<int, DiploidVariant > > VariantMap;

	VariantHash refpos_2_var;
	VariantHash querypos_2_var;

	vector<int> complex_ref_match_num;
	vector<int> complex_que_match_num;

    void ReadRefVCF(string filename);
    void ReadQueryVCF(string filename);

	void DirectSearchInThread(unordered_map<int, DiploidVariant> & ref_snps,
                           unordered_map<int, DiploidVariant> & query_snps,
                           int thread_index);
	void DirectSearchMultiThread();
    void ClusteringVariants();
    bool ClusteringMatchInThread(int, int, int);
	void ClusteringMatchMultiThread();

protected:
	bool scoring_basepair;
	bool overlap_match;
	map<int, vector<DiploidVariant> > cluster_vars_map;

	void DecideBoundaries();

	bool ReadDiploidVCF(string filename, VariantHash & pos_2_var);
	bool NormalizeDiploidVariant(DiploidVariant & var);

	bool VariantMatch(vector<DiploidVariant> & variant_list, int thread_index);

    bool FindBestMatch(vector<DiploidVariant> & variant_list,
		const string subsequence,
		const int offset,
		int index,
		map<int, DiploidVariant> separate_pos_var[],
		vector<vector<int>> max_choices[],  // 4 vectors
		int & max_score,
		bool & max_heterozygosity,
		string max_paths[]); //only two

    int CheckPrefix(const string subsequence,
		const int offset,
		map<int, DiploidVariant> separate_pos_var[],
		map<int, int> choices[],
		string cur_paths[]);

	bool RecurrentVariantMatch(vector<DiploidVariant> & variant_list, int thread_index);
	void RecurrentMatchWithIndel(vector<DiploidVariant> & variant_list,
		const string subsequence,
		const int offset,
		int index,
		map<int, DiploidVariant> separate_pos_var[],
		map<int, int> choices[], // 4 vectors
		map<int, int> max_matches[],  // 4 vectors
		int & max_score,
		string max_paths[]);

	vector<vector<vector<int>>> Combine(vector<int> & positions,
                                     vector<bool> & multi_indicators,
                                     int k);

	void FindComb(vector<int> & positions,
		vector<bool> & multi_indicators,
		int start,
		int k,
		vector<vector<int> > & sol,
		vector<vector<vector<int>>> & all_sol);

    void ModifyRefMultiVar(const string & ref,
                           int offset,
                           map<int, DiploidVariant> & pos_var,
                           vector<vector<int>> pos_choice,
                           string & donor,
                           int & score);

    void DivisiveHierarchicalClustering(list<vector<DiploidVariant> > & snp_clusters);

    bool VariantMatchWithOverlap(vector<DiploidVariant> & variant_list, int thread_index);

    bool FindBestMatchWithOverlap(vector<DiploidVariant> & variant_list,
                                const string subsequence,
                                const int offset,
                                int index,
                                map<int, DiploidVariant> separate_pos_var[],
                                set<int> selected_positions[]);

public:
	DiploidVCF(int thread_num_);
	~DiploidVCF();

    int test();
    // for public access
	void Compare(string ref_vcf,
		string query_vcf,
		string genome_seq,
		bool direct_search,
		string output_prefix,
		bool match_genotype,
		bool normalization,
		bool scoring_basepair,
		bool overlap_match);

};
