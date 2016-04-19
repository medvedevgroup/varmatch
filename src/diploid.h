#pragma once

#include "vcf.h"
// data structure for direct search
typedef struct DiploidVariant {
    DiploidVariant(int pos_ = -1,
        string ref_ = "",
        vector<string> alts_ = {"",""},
        string genotype_ = "0/0",
        bool heterozygous_ = false,
        bool multi_alts_ = false,
        int flag_ = 0) :
        pos(pos_),
        ref(ref_),
        alts(alts_),
        genotype(genotype_),
        heterozygous(heterozygous_),
        multi_alts(multi_alts_),
        flag(flag_){}

    int pos;
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


class VariantSelection{
public:
    int score;
    int separate_score[2];
    int min_genome_pos; // min(donor_sequence[0], donor_sequence[2])
    bool haplotypes_consistent;
    int genome_position[2]; // current reference position that has been covered
    int donor_length[2];
    vector<int> pos_vectors[2]; // selected variants, not necessary now
    vector<int> phasing_vectors[2]; // phasing vector for corresponding variant, for D_0
    int cur_var;

    VariantSelection(){
        score = 0;
        cur_var = -1;
        min_genome_pos = -1;
        haplotypes_consistent = false;
        for(int i = 0; i < 2; i++){
            separate_score[i] = 0;
            genome_position[i] = -1;
            donor_length[i] = 0;
            pos_vectors[i] = vector<int>();
            phasing_vectors[i] = vector<int>();
        }
    }
};

class DiploidVCF : public VCF
{
private:

    int total_ref_complex;
	int total_que_complex;

    typedef vector<unordered_map<int, DiploidVariant > > VariantHash;
    typedef vector<map<int, DiploidVariant > > VariantMap;

	VariantHash refpos_2_var;
	VariantHash querypos_2_var;

	vector<int> complex_ref_match_num;
	vector<int> complex_que_match_num;

	vector<DiploidVariant> variant_list;
	vector<DiploidVariant> ref_variant_list;
	vector<DiploidVariant> que_variant_list;

    void ReadRefVCF(string filename);
    void ReadQueryVCF(string filename);

	void DirectSearchInThread(unordered_map<int, DiploidVariant> & ref_snps,
                           unordered_map<int, DiploidVariant> & query_snps,
                           int thread_index);
	void DirectSearchMultiThread();
    void ClusteringVariants();
    bool ClusteringMatchInThread(int, int, int);
	void ClusteringMatchMultiThread();

	ofstream offf;
	const time_t ctt = time(0);

protected:
	bool scoring_basepair;
	bool overlap_match;
	bool variant_check;
	map<int, vector<DiploidVariant> > cluster_vars_map;

	void DecideBoundaries();

	bool ReadDiploidVCF(string filename, vector<DiploidVariant> & x_variant_list, int flag);
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

    bool FindBestDiploidMatch(vector<DiploidVariant> & variant_list,
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

    vector<vector<vector<int>>> DiploidCombine(vector<int> & positions,
                                               vector<bool> & heter_indicators,
                                               vector<bool> & multi_indicators,
                                               int k);

    void FindDiploidComb(vector<int> & positions,
        vector<bool> & heter_indicators,
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
                                map<int, int> selected_positions[]);

    bool VariantMatchPathCreation(vector<DiploidVariant> & variant_list, int thread_index, int cluster_id);
    bool CollapseSelections(VariantSelection & selection,
                            list<VariantSelection> & variant_selections);

    int CheckDonorSequences(vector<DiploidVariant> separate_var_list[],
                                      VariantSelection & selection,
                                      const string & subsequence,
                                      int offset,
                                      string donor_sequences[]);

      bool AddVariantToSelection(list<VariantSelection> & variant_selections,
       VariantSelection selection,
       DiploidVariant variant,
       int haplotype,
       vector<DiploidVariant> separate_var_list[],
       const string & subsequence,
       int offset,
       VariantSelection & best_selection);

    void SortVariantList();
    void ReadGenome(string filename);
    void LinearClusteringVariants();
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
		bool overlap_match,
		bool variant_check);

};
