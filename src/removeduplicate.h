#pragma once

#include "vcf.h"

typedef unordered_map<string, string> VCFEntryHash;

class RemoveDuplicate: public VCF
{
private:
    map<int, vector<SNP> > nondup_pos_snp_map;
    VCFEntryHash nondup_vcfentry_hash; // id is pos_ref_alt with uppercase

    int GetThreadIndex(int pos);
    int ReadVCFWithoutDup(string filename);
    void ClusteringRemoveDuplicateInThread(int start, int end, int thread_index);
    void ClusteringRemoveDuplicateMultiThread();
    void ClusteringSnps() override;
    void DivisiveHierarchicalClustering(list<vector<SNP>>& snp_clusters);
    void FindMatches(vector<SNP> snp_list, int thread_index);
    bool FindOneMatch(vector<SNP> snp_list, const string subsequence, int offset, int thread_index);

public:
    RemoveDuplicate(int thread_num_);
    ~RemoveDuplicate();

    void Deduplicate(string vcf_filename,
            string genome_filename,
            bool direct_search,
            string output_prefix);

};
