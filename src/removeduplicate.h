#include "vcf.h"

typedef vector< unordered_map<string, string>> VCFEntryHash;

class RemoveDuplicate: public VCF
{
private:
    SnpMap nondup_pos_snp_map_list;
    VCFEntryHash nondup_vcfentry_hash_list; // id is pos_ref_alt with uppercase

    int GetThreadIndex(int pos);
    int ReadVCFWithoutDup(string filename);
    void ClusteringRemoveDuplicateInThread();
    void ClusteringRemoveDuplicateMultiThread();

public:
    RemoveDuplicate(int thread_num_);
    ~RemoveDuplicate();

    Deduplicate(string vcf_filename
            string genome_filename,
            bool direct_search,
            string output_prefix);

}