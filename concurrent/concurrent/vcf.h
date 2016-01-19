#pragma once // the same purpose as #include guards

#define DEBUG
#include <unordered_map>
#include <algorithm>
#include <vector>
#include <map>
#include <iostream>
#include "util.h"
#include "threadguard.h"
using namespace std;


typedef struct SNP {
	SNP(char snp_type_ = 'S', string ref_ = "", string alt_ = "") :
		snp_type(snp_type_), ref(ref_), alt(alt_){}
	char snp_type;
	string ref;
	string alt;
}SNP;

typedef vector<unordered_map<int, vector<SNP>>> SnpHash;
typedef unordered_map<int, string> VCFEntryHash;

class VCF
{
private:
	void ReadVCF(string filename, SnpHash & pos_2_snps, VCFEntryHash & pos_2_vcf_entry);
	int thread_num;
	vector<int> pos_boundries;
public:
	VCF(int thread_num_ = 0);
	~VCF();

	SnpHash refpos_2_snp;
	VCFEntryHash refpos_2_vcf_entry;
	SnpHash querypos_2_snp;
	VCFEntryHash querypos_2_vcf_entry;

	vector<int> refpos_list;
	vector<int> querypos_list;
	
	void ReadRefVCF(string filename);
	void ReadQueryVCF(string filename);
	void ReadGenomeSequence(string filename, string & genome_sequence);
	//void DirectSearchInThread(int start_index, int end_index);
	void DirectSearch();
	void ComplexSearch();
	void ClusteringSnps();
	void ClusteringSearch();
};

