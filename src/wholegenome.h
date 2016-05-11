#pragma once

#include "vcf.h"
#include "diploid.h"

typedef struct VariantIndicator{
    int chr_id;
    int var_id;
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
};
