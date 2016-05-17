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
};
