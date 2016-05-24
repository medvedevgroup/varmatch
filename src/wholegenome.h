#pragma once

#include "vcf.h"
#include "diploid.h"
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

    int chr_id;
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
    SequencePath(int n)
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
};

class WholeGenome: public DiploidVCF{
private:
    int chrom_num;

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
    // copy the above into this.

    bool ReadWholeGenomeSequence(string filename);
    bool ReadGenomeSequenceList(string filename);
    int ReadWholeGenomeVariant(string filename, int flag);
    bool ReadVariantFileList(string filename);
    int ReadReferenceVariants(string filename);
    int ReadQueryVariants(string filename);
    bool ParallelClustering(); // parallel by chr id
    bool ParallelMatching(); // parallel by task
    bool TBBMatching();

    void SingleThreadClustering(int chr_id);
    bool MatchingSingleCluster(int cluster_index, int thread_index);

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

    bool PathNeedDecision(SequencePath& sp, multimap<int, int> * choices_by_pos[], int pos);
    int PathExtendOneStep(SequencePath& sp,
                          multimap<int, int> * choices_by_pos[],
                          const string & reference_sequence,
                          vector<int> & sync_points);

    bool PathMakeDecision(SequencePath& sp,
                                 vector<DiploidVariant> & variant_list,
                                 multimap<int, int> * choices_by_pos[],
                                 list<SequencePath> & sequence_path_list,
                                 const string & reference_sequence);

    bool MatchingSingleClusterBaseExtending(int cluster_index, int thread_index);
    bool DonorLengthEqual(SequencePath & a, SequencePath & b);
    void ConvergePaths(list<SequencePath> & path_list);
    int CheckPathEqualProperty(SequencePath & sp);


public:
    WholeGenome(int thread_num_);
    ~WholeGenome();
    void Compare(string ref_vcf,
        string query_vcf,
        string genome_seq,
        string output_prefix,
        bool match_genotype,
        bool normalization,
        bool score_basepair,
        bool variant_check);

    void DirectMatch(string ref_vcf,
                string query_vcf,
                bool match_genometype,
                bool normalization);

    int test(); // for direct test
    void PrintPath(SequencePath & sp);
};
