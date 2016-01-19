#pragma once // the same purpose as #include guards

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

typedef vector<unordered_map<int, vector<SNP> > > SnpHash;
typedef unordered_map<int, string> VCFEntryHash;

class VCF
{
private:
	int thread_num;
	vector<int> pos_boundries; // boundries for split multi hash table
	string genome_sequence; // genome sequence from fasta file
	bool boundries_decided; // before deciding boundries, can not read vcf file, because do not know how to split

	void ReadVCF(string filename, SnpHash & pos_2_snps, VCFEntryHash & pos_2_vcf_entry);
	void DirectSearchInThread(unordered_map<int, vector<SNP> > & ref_snps, unordered_map<int, vector<SNP> > & query_snps);
	bool CompareSnps(SNP r, SNP q);
public:
	VCF(int thread_num_ = 0);
	~VCF();

	SnpHash refpos_2_snp;
	VCFEntryHash refpos_2_vcf_entry;
	SnpHash querypos_2_snp;
	VCFEntryHash querypos_2_vcf_entry;
	
	void ReadRefVCF(string filename);
	void ReadQueryVCF(string filename);
	void ReadGenomeSequence(string filename);
	void DecideBoundries();
	
	void DirectSearchMultiThread();
	void ComplexSearch();
	void ClusteringSnps();
	void ClusteringSearch();

	int GetRefSnpNumber();
	int GetQuerySnpNumber();
};

