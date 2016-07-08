#pragma once

#include "vcf.h"

// data structure for direct search
class DiploidVariant {
public:
    DiploidVariant(int pos_ = -1,
        string ref_ = "",
        vector<string> alts_ = {"",""},
        bool heterozygous_ = false,
        bool multi_alts_ = false,
        int mdl_ = 0,
        int mil_ = 0,
        int flag_ = 0) :
        pos(pos_),
        ref(ref_),
        alts(alts_),
        heterozygous(heterozygous_),
        multi_alts(multi_alts_),
        mdl(mdl_),
        mil(mil_),
        flag(flag_){}

    int pos;
    string ref;
    vector<string> alts;
    bool heterozygous;
    bool multi_alts;
    int mdl;
    int mil;
    int flag; //in DiploidVariant, flag = 0 is reference, flag = 1 is query

//    int get_pos() const{return pos};
//    string get_ref() const{return ref};
//    vector<string> get_alts() const{return alts};
//    bool get_heterozygous() const{return heterozygous};
//    bool get_multi_alts() const{return multi_alts};

    bool operator <(const DiploidVariant& y) const {
        return pos < y.pos;
    }

    // this is based on the assumption that all sequence are in upper case
    bool operator ==(const DiploidVariant& y) {
        if (pos == y.pos && ref == y.ref) {
            if(heterozygous == y.heterozygous && multi_alts == y.multi_alts){
                if (multi_alts && heterozygous) {
                    int match_times = 0;
                    for (int i = 0; i < 2; i++) {
                        for (int j = 0; j < 2; j++) {
                            if (alts[i] == y.alts[j])
                                match_times++;
                        }
                    }
                    if (match_times >= 2)
                        return true;
                }
                else if(alts[0] == y.alts[0]){
                    return true;
                }
            }
        }
        return false;
    }

    bool DirectCompare(const DiploidVariant& y){
        if (pos == y.pos && ref == y.ref) {
            if (multi_alts && heterozygous && y.multi_alts && y.heterozygous) {
                int match_times = 0;
                for (int i = 0; i < 2; i++) {
                    for (int j = 0; j < 2; j++) {
                        if (alts[i] == y.alts[j])
                            match_times++;
                    }
                }
                if (match_times > 0)
                    return true;
            }
            else if(alts[0] == y.alts[0]){
                return true;
            }
        }
        return false;
    }

    bool CompareNoGenotype(const DiploidVariant & y){
        if(pos == y.pos && ref == y.ref){
            if(alts[0] == y.alts[0]) return true;
            if(multi_alts){
                if(alts[1] == y.alts[0]) return true;
                if(y.multi_alts && alts[1] == y.alts[1]){
                    return true;
                }
            }
            if(y.multi_alts && alts[0] == y.alts[1]){
                return true;
            }
        }
        return false;
    }

};

// define outside of struct, idiomatic solution for lexicographical compare for structures
//bool operator <(const DiploidVariant& x, const DiploidVariant& y);

//bool operator ==(const DiploidVariant& x, const DiploidVariant& y);

class VariantSelection{
public:
    int score;
    int separate_score[2];
    int min_genome_pos; // min(donor_sequence[0], donor_sequence[2])
    bool haplotypes_consistent;
    int genome_position[2]; // genome position that has been considered, exclusive
    int donor_length[2];
    string donor_sequences[4];
    vector<int> pos_vectors[2]; // selected variants, not necessary now
    vector<int> phasing_vectors[2]; // phasing vector for corresponding variant, for D_0
    int cur_var;
    bool overlap_detected;

    VariantSelection(){
        score = 0;
        cur_var = -1;
        min_genome_pos = -1;
        haplotypes_consistent = false;
        overlap_detected = false;
        for(int i = 0; i < 2; i++){
            separate_score[i] = 0;
            genome_position[i] = 0;
            donor_length[i] = 0;
            pos_vectors[i] = vector<int>();
            phasing_vectors[i] = vector<int>();
            donor_sequences[i] = "";
        }
        donor_sequences[2] = "";
        donor_sequences[3] = "";
    }

    bool operator< (const VariantSelection& rhs) const // sort by min_genome_position
    {
        return min_genome_pos < rhs.min_genome_pos;
    }
};

class DiploidVCF : public VCF
{
private:

    typedef vector<unordered_map<int, DiploidVariant > > VariantHash;
    typedef vector<map<int, DiploidVariant > > VariantMap;

	VariantHash refpos_2_var;
	VariantHash querypos_2_var;

	vector<DiploidVariant> variant_list;
	vector<DiploidVariant> ref_variant_list;
	vector<DiploidVariant> que_variant_list;

    int ReadRefVCF(string filename);
    int ReadQueryVCF(string filename);

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

    vector<int> complex_ref_match_num;
	vector<int> complex_que_match_num;
    int total_ref_complex;
	int total_que_complex;

	bool scoring_basepair;
	bool overlap_match;
	bool variant_check;
	map<int, vector<DiploidVariant> > cluster_vars_map;

	void DecideBoundaries();

	int ReadDiploidVCF(string filename, vector<DiploidVariant> & x_variant_list, int flag);
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
    bool CollapseSelections(VariantSelection selection,
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

    int NormalizeVariantSequence(int pos,
                                 string & parsimonious_ref,
                                 string & parsimonious_alt0,
                                 string & parsimonious_alt1);

    int ExtendingDonorSequences(vector<DiploidVariant> separate_var_list[],
                                      VariantSelection & selection,
                                      const string & subsequence,
                                      int offset,
                                      int flag);

    bool CollapsePrefixMatchSelection(VariantSelection selection,
                                    list<VariantSelection> & variant_selections);

    void ReverseLinearClusteringVariants();
    bool AcceleratedVariantMatchPathCreation(vector<DiploidVariant> & variant_list, int thread_index, int cluster_id);
    bool VariantMatchPathCreationByDonor(vector<DiploidVariant> & variant_list, int thread_index, int cluster_id);

    int CheckDonorSequencesWithOverlap(vector<DiploidVariant> separate_var_list[],
                                      VariantSelection & selection,
                                      const string & subsequence,
                                      int offset,
                                      string donor_sequences[]);

    void PrintVariant(DiploidVariant var);

public:
	DiploidVCF(int thread_num_);
	~DiploidVCF();

	const static int VAR_LEN = 100;

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
