#pragma once

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
#include <set>

#include "util.h"
#include "diploidvariant.h"
//#include "tbb/task_scheduler_init.h"
//#include "tbb/blocked_range.h"
//#include "tbb/parallel_for.h"
//#include "tbb/concurrent_vector.h"

typedef struct VariantIndicator{
    VariantIndicator(int chr_id_ = -1,
    int var_id_ = -1,
    bool refer_ = true) :
    chr_id(chr_id_),
    var_id(var_id_),
    refer(refer_){}

    char chr_id;
    int var_id;
    bool refer;
}VariantIndicator;

typedef struct Interval {
    int start;
    int end;
    Interval() : start(0), end(0) {}
    Interval(int s, int e) : start(s), end(e) {}
}Interval;

class SequencePath{
public:
    SequencePath(int n, int v)
    {
        reference_length = n;
        for(int i = 0; i < 4; i++){
            string_sequences[i].resize(n, ".");
            // default value is "."
            donor_sequences[i] = "";
        }
        current_genome_pos = -1;
        score = 0;
        removable = false;
        same_donor_len = false;
        current_equal_donor_pos[0] = -1;
        current_equal_donor_pos[1] = -1;
        reached_sync_num = 0;

        for(int i = 0; i < v; i++){
            choice_vector.push_back(-89);
        }
    }
    int reference_length;
    vector<string> string_sequences[4];
    map<int, pair<int, int>> choice_made[2]; // this can be used to indicate if choice is made and which choice
    // one choice is a pair: variant id, phasing index
    int current_genome_pos;
    string donor_sequences[4];
    int current_equal_donor_pos[2];
    int score;
    bool removable;
    bool same_donor_len;
    int reached_sync_num;
    vector<int> choice_vector;
};

class WholeGenome{
private:
    int chrom_num;
    int thread_num;
    string ref_vcf_filename;
    string que_vcf_filename;
    int baseline_variant_total_num;
    int query_variant_total_num;
    vector<string> baseline_variant_strings;
    vector<string> query_variant_strings;
    bool detail_results;

    //int thread_num; VCF->DiploidVariant->WholeGenome
protected:
    map<string, int> chrid_by_chrname;
    map<int, string> chrname_by_chrid;
    map<string, int> chrname_dict;
    map<int, string> genome_sequences;
    vector<DiploidVariant> ** ref_variant_by_chrid;
    vector<DiploidVariant> ** que_variant_by_chrid;
    vector<vector<VariantIndicator>> ** variant_cluster_by_chrid;
    // so here cluster is represented as vector<vector<VariantIndicator>>
    // and we create a list of pointers point to cluster
    // and we hold the point to that list

    vector<vector<VariantIndicator>> variants_by_cluster;

    vector<string> *** match_records_by_mode_by_thread;

    //vector<int> *** baseline_matches_by_mode_by_thread;
    //vector<int> *** query_matches_by_mode_by_thread;
    vector<int> *** baseline_total_match_num;
    vector<int> *** query_total_match_num;

    vector<int> *** baseline_total_edit_distance;
    vector<int> *** query_total_edit_distance;

    //map<float, int> *** tp_qual_num_by_mode_by_thread;
    //map<float, int> *** fp_qual_num_by_mode_by_thread;

    //map<float, int> query_total_qual_num;

    string output_prefix;
    string output_dir;
    // copy the above into this.

    int score_unit_indicator;
    int match_mode_indicator;
    int score_scheme_indicator;

    vector<int> score_unit_list;
    vector<int> match_mode_list;
    vector<int> score_scheme_list;
    vector<int> mode_index_list;

    vector<double> threshold_list;
    int threshold_num;

    vector<float> per_list;

    bool ReadWholeGenomeSequence(string filename);
    bool ReadGenomeSequenceList(string filename);
    int ReadWholeGenomeVariant(string filename, bool flag);
    bool ReadVariantFileList(string filename);
    int ReadReferenceVariants(string filename);
    int ReadQueryVariants(string filename);
    bool ParallelClustering(); // parallel by chr id
    bool ParallelMatching(); // parallel by task
    bool TBBMatching();

    void SingleThreadClustering(int chr_id);
    //bool MatchingSingleCluster(int cluster_index, int thread_index, int match_mode);

    //override
    bool ClusteringMatchInThread(int start, int end, int thread_index);
    void ClusteringMatchMultiThread();
    int NormalizeVariantSequence(int pos,
                             string & parsimonious_ref,
                             string & parsimonious_alt0,
                             string & parsimonious_alt1,
                             int chr_id);

    struct compInterval {
        bool operator()(const Interval &a, const Interval &b) const {
            return a.start<b.start;
        }
    };

    vector<Interval> merge(vector<Interval> &intervals) {
        sort(intervals.begin(),intervals.end(),compInterval());
        vector<Interval> results;
        for(int i=0; i<intervals.size(); i++) {
            if(results.empty() || results.back().end < intervals[i].start)  // no overlap
                results.push_back(intervals[i]);
            else   // overlap
                results.back().end = max(results.back().end, intervals[i].end);
        }
        return results;
    }

    int PathNeedDecision(SequencePath& sp, multimap<int, int> * choices_by_pos[], int pos);
    int PathExtendOneStep(SequencePath& sp,
                          multimap<int, int> * choices_by_pos[],
                          const string & reference_sequence,
                          vector<int> & sync_points,
                          int match_mode,
                          int & variant_need_decision);

    bool PathMakeDecision(SequencePath& sp,
                                 vector<DiploidVariant> & variant_list,
                                 multimap<int, int> * choices_by_pos[],
                                 list<SequencePath> & sequence_path_list,
                                 const string & reference_sequence,
                                 int score_unit,
                                 int match_mode,
                                 int score_scheme);

    bool VariantMakeDecision(SequencePath& sp,
                             vector<DiploidVariant> & variant_list,
                             list<SequencePath> & sequence_path_list,
                             const string & reference_sequence,
                             int score_unit,
                             int match_mode,
                             int score_scheme,
                             int variant_index);

    bool VariantMakeDecisionNoGenotype(SequencePath& sp,
                         vector<DiploidVariant> & variant_list,
                         list<SequencePath> & sequence_path_list,
                         const string & reference_sequence,
                         int score_unit,
                         int match_mode,
                         int score_scheme,
                         int variant_index);

    bool AppendChangedSp(SequencePath& sp,
                         vector<DiploidVariant> & variant_list,
                         list<SequencePath> & sequence_path_list,
                         const string & reference_sequence,
                         int score_unit,
                         int match_mode,
                         int score_scheme,
                         int variant_index,
                         int c);

    bool AppendChangedSpNoGenotype(SequencePath& sp,
                         vector<DiploidVariant> & variant_list,
                         list<SequencePath> & sequence_path_list,
                         const string & reference_sequence,
                         int score_unit,
                         int match_mode,
                         int score_scheme,
                         int variant_index,
                         int c);

    bool PathMakeDecisionBackup(SequencePath& sp,
                                 vector<DiploidVariant> & variant_list,
                                 multimap<int, int> * choices_by_pos[],
                                 list<SequencePath> & sequence_path_list,
                                 const string & reference_sequence,
                                 int score_unit,
                                 int match_mode,
                                 int score_scheme);

    bool MatchingSingleClusterBaseExtending(int cluster_index,
                                            int thread_index,
                                            vector<DiploidVariant> & variant_list,
                                            string & subsequence,
                                            int offset,
                                            multimap<int, int> * choices_by_pos[],
                                            vector<int> & sync_points,
                                            int chr_id,
                                            int score_unit,
                                            int match_mode,
                                            int score_scheme,
                                            int threshold_index);

    bool DonorLengthEqual(SequencePath & a, SequencePath & b);
    void ConvergePaths(list<SequencePath> & path_list);
    int CheckPathEqualProperty(SequencePath & sp, int match_mode);

    int ScoreEditDistance(DiploidVariant & dv, int allele_indicator);
    int EditDistance(const std::string& s1, const std::string& s2);
    bool PathMakeDecisionNoGenotype(SequencePath& sp,
                                 vector<DiploidVariant> & variant_list,
                                 multimap<int, int> * choices_by_pos[],
                                 list<SequencePath> & sequence_path_list,
                                 const string & reference_sequence,
                                 int score_unit,
                                 int match_mode,
                                 int score_scheme);

    void ConstructMatchRecord(SequencePath & best_path,
                               vector<DiploidVariant> & variant_list,
                               string & subsequence,
                               int offset,
                               int thread_index,
                               int chr_id,
                               int mode_index,
                               int threshold_index);

    void ConstructMatchRecordBackup(SequencePath & best_path,
                               vector<DiploidVariant> & variant_list,
                               string & subsequence,
                               int offset,
                               int thread_index,
                               int chr_id,
                               int mode_index,
                               int threshold_index);

    void ConstructMatchRecordNoGenotype(SequencePath & best_path,
                                       vector<DiploidVariant> & variant_list,
                                       string & subsequence,
                                       int offset,
                                       int thread_index,
                                       int chr_id,
                                       int mode_index,
                                       int threshold_index);

    void ConstructMatchRecordNoGenotypeBackup(SequencePath & best_path,
                                       vector<DiploidVariant> & variant_list,
                                       string & subsequence,
                                       int offset,
                                       int thread_index,
                                       int chr_id,
                                       int mode_index,
                                       int threshold_index);

    int CalculateScore(DiploidVariant & dv,
                       int choice,
                       int score_unit,
                       int match_mode,
                       int score_scheme);

    int GetIndexFromMatchScore(int score_unit, int match_mode, int score_scheme);
    bool ClearQuery();

    inline void ToUpper(string & s){
        transform(s.begin(), s.end(), s.begin(), ::toupper);
    }

    bool CheckTandemRepeat(string sequence, int unit_threshold);

    bool MatchVariantListInThread(int thread_index, 
        int threshold_index,
        int chr_id,
        vector<DiploidVariant> & variant_list,
        int cluster_id);

    void initialize_score_matrix(int **score, char **trackBack, int M, int N);
    int needleman_wunsch(string S1, string S2, string &R1, string &R2);
    void GenerateAltVector(string ref, string alt, vector<string> & alt_vector);

    int CalculateEditDistance(DiploidVariant & dv,
                                int choice,
                                int match_mode);

public:
    WholeGenome(int thread_num_,
                string output_dir_,
                bool pr_curves_);

    ~WholeGenome();

    void ReadRef(string genome_seq, 
      string ref_vcf);

    void Compare(string query_vcf,
        string output_prefix,
        bool detail_results,
        int score_unit_,
        int match_mode_,
        int score_scheme_);

    void DirectMatch(string ref_vcf,
                string query_vcf,
                int match_mode_,
                string output_prefix);

    int test(); // for direct test
    void PrintPath(SequencePath & sp);

    const static int MATCH_MODE_NUM = 16;
    const static int VAR_LEN = 100;
    const static int MAX_REPEAT_LEN = 1000;
    const static int ROC_SAMPLE_NUM = 5;
    const static int MEANING_CHOICE_BOUND = -10;
    const static int NOT_USE = -9;
    const static int EASY_MATCH_VAR_NUM = 5;
};
