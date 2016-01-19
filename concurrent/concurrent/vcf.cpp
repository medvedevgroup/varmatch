#include "vcf.h"

VCF::VCF(int thread_num_)
{
	if (thread_num_ == 0) {
		thread_num = 1;
	}
	else {
		thread_num = min(thread_num_, (int)thread::hardware_concurrency());
	}
	dout << "Thread Number: " << thread_num << endl;
}


VCF::~VCF()
{
}



void VCF::ReadVCF(string filename, SnpHash & pos_2_snp, VCFEntryHash& pos_2_vcf_entry) {
	ifstream vcf_file;
	vcf_file.open(filename.c_str());
	if (!vcf_file.good()) {
		cout << "[Error] VCF::ReadVCF can not open vcf file" << endl;
	}

	while (!vcf_file.eof()) { // alternative way is vcf_file != NULL
		string line;
		getline(vcf_file, line, '\n');
		//dout << line << endl;
		if (line.length() <= 1) continue;
		if (line[0] == '#') continue;
		auto columns = split(line, '\t');
		auto pos = atoi(columns[1].c_str());
		auto ref = columns[3];
		auto alt = columns[4];
		auto quality = columns[6];

		if (alt.find(",") != string::npos) continue; // can not deal with multi alt yet
		//todo(Chen) deal with multi alt

		char snp_type = 'S'; 
		if (ref.length() > alt.length()) {
			snp_type = 'D';
		}
		else if (ref.length() < alt.length()) {
			snp_type = 'I';
		}

		//pos_2_snp[pos].push_back(SNP(snp_type, ref, alt));
		pos_2_vcf_entry[pos] = line;
		pos_list.push_back(pos);
	}
	vcf_file.close();
	return;
}

void VCF::ReadGenomeSequence(string filename, string & genome_sequence) {
	ifstream genome_file;
	genome_file.open(filename.c_str());
	if (!genome_file.good()) {
		cout << "[Error] VCF::ReadGenomeSequence can not open fasta file" << endl;
	}

	genome_sequence = "";

	while(!genome_file.eof()) {
		string line;
		getline(genome_file, line, '\n');
		if (line.length() <= 1) continue;
		if (line[0] == '>') continue;
		genome_sequence += line;
	}
	genome_file.close();
	return;
}

void VCF::ReadRefVCF(string filename) {
	ReadVCF(filename, refpos_2_snp, refpos_2_vcf_entry, refpos_list);
}

void VCF::ReadQueryVCF(string filename) {
	ReadVCF(filename, querypos_2_snp, querypos_2_vcf_entry, querypos_list);
}

//void VCF::DirectSearchInThread(int start_index, int end_index) {}

// directly match by position
void VCF::DirectSearch() {
	//sort(refpos_list.begin(), refpos_list.end());
	//sort(querypos_list.begin(), querypos_list.end());
}

// match by overlapping reference region
void VCF::ComplexSearch() {}

// clustering snps
void VCF::ClusteringSnps() {}

// match by cluster
void VCF::ClusteringSearch() {}
