#include "wholegenome.h"
#include "tbb/task_scheduler_init.h"
#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"

using namespace std;

WholeGenome::WholeGenome(int thread_num_):DiploidVCF(thread_num_){
    chrom_num = 24;
    dout << "WholeGenome() Thread Number: " << thread_num << endl;
    ref_variant_by_chrid = new vector<DiploidVariant>*[chrom_num];
    que_variant_by_chrid = new vector<DiploidVariant>*[chrom_num];
    variant_cluster_by_chrid = new vector<vector<VariantIndicator>> *[chrom_num];
	for (int j = 0; j < chrom_num; j++) {
		ref_variant_by_chrid[j] = new vector<DiploidVariant>;
		que_variant_by_chrid[j] = new vector<DiploidVariant>;
		variant_cluster_by_chrid[j] = new vector<vector<VariantIndicator>>;
	}

    // chr_id starts from 0
	for(int j = 1; j <= 22; j++){
        string chr_name = to_string(j);
        chrname_dict[chr_name] = j-1;
        chr_name = "chr"+chr_name;
        chrname_dict[chr_name] = j-1;
	}
	chrname_dict["X"] = 22;
	chrname_dict["chrX"] = 22;
	chrname_dict["Y"] = 23;
	chrname_dict["chrY"] = 23;

}

WholeGenome::~WholeGenome(){

    for(int j = 0; j < chrom_num; j++){
        ref_variant_by_chrid[j]->clear();
        que_variant_by_chrid[j]->clear();
        variant_cluster_by_chrid[j]->clear();
        delete ref_variant_by_chrid[j];
        delete que_variant_by_chrid[j];
        delete variant_cluster_by_chrid[j];
    }
    delete[] ref_variant_by_chrid;
    delete[] que_variant_by_chrid;
    delete[] variant_cluster_by_chrid;
}

bool WholeGenome::ReadWholeGenomeSequence(string filename){
    std::ifstream input(filename);
    if(!input.good()){
        std::cerr << "Error opening '"<<filename<<"'. Bailing out." << std::endl;
        return false;
    }

    std::string line, name, content;
    while( std::getline( input, line ).good() ){
        if( line.empty() || line[0] == '>' ){ // Identifier marker
            if( !name.empty() ){ // Print out what we read from the last entry
                //std::cout << name << " : " << content << std::endl;
                if(chrname_dict.find(name) == chrname_dict.end()){
                    cout << "[VarMatch] Error: detected chromosome name: " << name <<" does not exist in human genome." << endl;
                    return false;
                }
                int chr_id = chrname_dict[name];
                chrid_by_chrname[name] = chr_id;
                chrname_by_chrid[chr_id] = name;
                genome_sequences[chr_id] = content;
                name.clear();
            }
            if( !line.empty() ){
                name = split(line, ' ')[0].substr(1);
            }
            content.clear();
        } else if( !name.empty() ){
            if( line.find(' ') != std::string::npos ){ // Invalid sequence--no spaces allowed
                name.clear();
                content.clear();
            } else {
                content += line;
            }
        }
    }
    if( !name.empty() ){ // Print out what we read from the last entry
        //std::cout << name << " : " << content << std::endl;
        if(chrname_dict.find(name) == chrname_dict.end()){
            cout << "[VarMatch] Error: detected chromosome name: " << name <<" does not exist in human genome." << endl;
            return false;
        }
        int chr_id = chrname_dict[name];
        chrid_by_chrname[name] = chr_id;
        chrname_by_chrid[chr_id] = name;
        genome_sequences[chr_id] = content;
    }
    // test

    for(auto it = genome_sequences.begin(); it != genome_sequences.end(); ++it){
        cout << it->first << ":" << (it->second).length();
    }
    return true;
}

bool WholeGenome::ReadGenomeSequenceList(string filename){

}

int WholeGenome::ReadWholeGenomeVariant(string filename, int flag){
    int total_num = 0;
	ifstream vcf_file;
	vcf_file.open(filename.c_str());
	if (!vcf_file.good()) {
		cout << "[VarMatch] Error: can not open vcf file" << endl;
		return -1;
	}
	int genotype_index = -1;
	char genotype_separator = '/';
	string chr_name = ".";
	//int genome_sequence_length = genome_sequence.length();
	while (!vcf_file.eof()) { // alternative way is vcf_file != NULL
		string line;
		getline(vcf_file, line, '\n');
		// check ineligible lines
		//dout << line << endl;
		if ((int)line.length() <= 1) continue;
		//if (line.find_first_not_of(' ') == std::string::npos) continue;

		if (line[0] == '#') {
//			if (line[1] == '#') continue;
//			auto head_names = split(line, '\t');
//			if (head_names.size() < 10 && match_genotype) {
//				cout << "[VarMatch] Warning: not enough information in VCF file for genotype matching." << endl;
//				cout << "[VarMatch] \tVCF file name " << filename << endl;
//				cout << "[VarMatch] \tAutomatically turn off genotype matching module." << endl;
//				match_genotype = false;
//			}
			continue;
		}
		auto columns = split(line, '\t');
		if (columns.size() < 10) {
			if(match_genotype){
                cout << "[VarMatch] Warning: not enough information in VCF file for genotype matching." << endl;
                cout << "[VarMatch] \tAutomatically turn off genotype matching module " << filename << endl;
                match_genotype = false;
                continue;
            }
            if(columns.size() < 6){
                cout << "[VarMatch] Warning: not enough information in VCF file for variant matching." << endl;
                cout << "[VarMatch] skip current variant: " << line << endl;
                continue;
            }
		}
		if (chr_name == ".") chr_name = columns[0];
		auto pos = atoi(columns[1].c_str()) - 1; // 0-based coordinate

		auto ref = columns[3];
		auto alt_line = columns[4];
		auto quality = columns[5];

		ToUpper(ref);
		ToUpper(alt_line);

		bool is_heterozygous_variant = false;
		bool is_multi_alternatives = false;

		if (columns.size() >= 10) {
			if (genotype_index < 0) {
                auto formats = split(columns[8], ':');
                for (int i = 0; i < formats.size(); i++) {
                    if (formats[i] == "GT") {
                        genotype_index = i;
                        break;
                    }
                }
                if(genotype_index < 0){
                    cout << "[VarMatch] VCF entry does not contain genotype information" << endl;
                    continue;
                }
			}
			auto additionals = split(columns[9], ':');
            vector<string> genotype_columns = split(additionals[genotype_index], genotype_separator);

            if(genotype_columns.size() != 2){
                genotype_separator = '|';
                genotype_columns = split(additionals[genotype_index], genotype_separator);
            }

			// normalize format of genotype: sorted, separated by |
			if (genotype_columns.size() != 2) {
				cout << "[VarMatch] Warning Unrecognized Genotype: " << additionals[genotype_index] << endl;
				continue;
			}
			else {
				if (genotype_columns[0] != genotype_columns[1]) {
					is_heterozygous_variant = true;
				}
			}

            if (genotype_columns[1] == "0" && genotype_columns[0] == "0" && match_genotype) {
                continue;
            }
		}

		vector<string> alt_list;
		if (alt_line.find(",") != std::string::npos) {
			alt_list = split(alt_line, ',');
			is_multi_alternatives = true;
		}
		else {
			alt_list.push_back(alt_line);
		}

        int snp_ins = max(0, (int)alt_list[0].length() - (int)ref.length());
        int snp_del = max(0, (int)ref.length() - (int)alt_list[0].length());
        if(is_multi_alternatives){
            snp_ins = max(snp_ins, (int)alt_list[1].length() - (int)ref.length());
            snp_del = max(snp_del, (int)ref.length() - (int)alt_list[1].length());
        }

		DiploidVariant dv(pos, ref, alt_list, is_heterozygous_variant, is_multi_alternatives, snp_del, snp_ins, flag);
		if (normalization) {
			NormalizeDiploidVariant(dv);
		}
        if(chrid_by_chrname.find(chr_name) != chrid_by_chrname.end()){
            int chr_id = chrid_by_chrname[chr_name];
            if(flag == 0){
                ref_variant_by_chrid[chr_id]->push_back(dv);
            }else{
                que_variant_by_chrid[chr_id]->push_back(dv);
            }
        }else{
            cout << "[VarMatch] Error: detected chromosome name does not exist: " << chr_name << endl;
            return 0;
        }

        total_num++;
	}
	vcf_file.close();
	return total_num;
}

bool WholeGenome::ReadVariantFileList(string filename){

}

bool WholeGenome::ParallelClustering(){

}

bool WholeGenome::ParallelMatching()
{

}

int WholeGenome::ReadReferenceVariants(string filename){
    return ReadWholeGenomeVariant(filename, 0);
}

int WholeGenome::ReadQueryVariants(string filename){
    return ReadWholeGenomeVariant(filename, 1);
}

void WholeGenome::Compare(string ref_vcf,
	string query_vcf,
	string genome_seq,
	string output_prefix,
	bool match_genotype,
	bool normalization,
	bool score_basepair,
	bool variant_check)
{
	ref_vcf_filename = ref_vcf;
	que_vcf_filename = query_vcf;
	this->normalization = normalization;
	this->scoring_basepair = score_basepair;
	this->variant_check = variant_check;
	this->match_genotype = match_genotype;

	output_stat_filename = output_prefix + ".stat";
    output_complex_filename = output_prefix + ".match";

    ReadWholeGenomeSequence(genome_seq);
    int ref_variant_num = ReadReferenceVariants(ref_vcf);
    int que_variant_num = ReadQueryVariants(query_vcf);
    dout << ref_variant_num << "," << que_variant_num << endl;
}























