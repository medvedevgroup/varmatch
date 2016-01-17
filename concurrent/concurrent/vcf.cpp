#include "vcf.h"

VCF::VCF()
{
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

		pos_2_snp[pos] = SNP(snp_type, ref, alt);
		pos_2_vcf_entry[pos] = line;
	}
	vcf_file.close();
	return;
}

void VCF::ReadRefVCF(string filename) {
	ReadVCF(filename, this->refpos_2_snp, this->refpos_2_vcf_entry);
}

void VCF::ReadQueryVCF(string filename) {
	ReadVCF(filename, this->querypos_2_snp, this->querypos_2_vcf_entry);
}

void VCF::DirectSearch() {

}
