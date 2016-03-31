#include "diploid.h"


bool operator <(const DiploidVariant& x, const DiploidVariant& y) {
	return x.pos < y.pos;
}

bool operator ==(const DiploidVariant& x, const DiploidVariant& y) {
	if (x.pos == y.pos && x.heterozygous == y.heterozygous && x.multi_alts == y.multi_alts) {
		if (x.multi_alts && x.heterozygous) {
			int match_times = 0;
			for (int i = 0; i < 2; i++) {
				for (int j = 0; j < 2; j++) {
					if (x.alts[i] == y.alts[j])
						match_times++;
				}
			}
			if (match_times >= 2)
				return true;
		}
		else if(x.alts[0] == y.alts[0]){
			return true;
		}
	}
	return false;
}


DiploidVCF::DiploidVCF(int thread_num_)
{
	genome_sequence = "";
	boundries_decided = false;
	clustering_search = false;
    scoring_basepair = false;
	match_genotype = true;
    if (thread_num_ <= 0) {
		thread_num = 1;
	}
	else {
		thread_num = min(thread_num_, (int)thread::hardware_concurrency());
	}
	dout << "Thread Number: " << thread_num << endl;
    chromosome_name = ".";
}

DiploidVCF::~DiploidVCF()
{
}

// private
void DiploidVCF::ReadRefVCF(string filename) {
    ReadDiploidVCF(filename, this->refpos_2_var);
}

// private
void DiploidVCF::ReadQueryVCF(string filename) {
    ReadDiploidVCF(filename, this->querypos_2_var);
}
// protected
bool DiploidVCF::NormalizeDiploidVariant(DiploidVariant & var) {
	int pos = var.pos;
	string parsimonious_ref = var.ref;
	string parsimonious_alt0 = var.alts[0];
	string parsimonious_alt1 = var.alts[0];
	if (var.heterozygous && var.multi_alts)
		parsimonious_alt1 = var.alts[1];

	int left_index = pos;
	if (genome_sequence.size() == 0) return false;
	if (parsimonious_ref.size() == 1 && parsimonious_alt0.size() == 1 && parsimonious_alt1.size() == 1) return true;
	if (toupper(genome_sequence[left_index]) != toupper(parsimonious_ref[0])) {
		dout << "[Error] genome sequence, subsequence, offset does not match." << endl;
		return false;
	}
	bool change_in_allels = true;
	while (change_in_allels) {
		change_in_allels = false;
		if (toupper(parsimonious_ref.back()) == toupper(parsimonious_alt0.back()) && toupper(parsimonious_ref.back()) == toupper(parsimonious_alt1.back())) {
			if ((parsimonious_ref.size() > 1 && parsimonious_alt0.size() > 1 && parsimonious_alt1.size() > 1) || left_index > 0) {
				parsimonious_ref.pop_back();
				parsimonious_alt0.pop_back();
				parsimonious_alt1.pop_back();
				change_in_allels = true;
			}
			else {
				return false;
			}
		}
		if (parsimonious_ref.length() == 0 || parsimonious_alt0.length() == 0 || parsimonious_alt1.length() == 0) {
			left_index--;
			char left_char = genome_sequence[left_index];
			parsimonious_ref = left_char + parsimonious_ref;
			parsimonious_alt0 = left_char + parsimonious_alt0;
			parsimonious_alt1 = left_char + parsimonious_alt1;
		}
	}
	while (toupper(parsimonious_ref[0]) == toupper(parsimonious_alt0[0]) && toupper(parsimonious_ref[0]) == toupper(parsimonious_alt1[0]) && parsimonious_ref.size() > 1 && parsimonious_alt0.size() > 1 && parsimonious_alt1.size() > 1) {
		parsimonious_ref.erase(0, 1);
		parsimonious_alt0.erase(0, 1);
		parsimonious_alt1.erase(0, 1);
	}
	var.pos = left_index;
	var.ref = parsimonious_ref;
	var.alts[0] = parsimonious_alt0;
	if (var.heterozygous && var.multi_alts)
		var.alts[1] = parsimonious_alt1;
	return true;
}

// protected
bool DiploidVCF::ReadDiploidVCF(string filename, VariantHash & pos_2_var) {
	if (!boundries_decided) {
		dout << "[VarMatch] DiploidVCF::ReadDiploidVCF cannot read vcf file before read genome file" << endl;
		return false;
	}

	ifstream vcf_file;
	vcf_file.open(filename.c_str());
	if (!vcf_file.good()) {
		cout << "[VarMatch] Error: can not open vcf file" << endl;
		return false;
	}

	if (normalization) {
		dout << "normalize variants while read" << endl;
	}
	string previous_line;
	while (!vcf_file.eof()) { // alternative way is vcf_file != NULL
		string line;
		getline(vcf_file, line, '\n');
		// check ineligible lines
		//dout << line << endl;
		if ((int)line.length() <= 1) continue;
		if (line.find_first_not_of(' ') == std::string::npos) continue;

		if (line[0] == '#') {
			if (line[1] == '#') continue;
			auto head_names = split(line, '\t');
			if (match_genotype && head_names.size() < 10) {
				cout << "[VarMatch] Warning: not enough information in VCF file for genotype matching." << endl;
				cout << "[VarMatch] \tVCF file name " << filename << endl;
				cout << "[VarMatch] \tAutomatically turn off genotype matching module." << endl;
				match_genotype = false;
			}
			continue;
		}
		auto columns = split(line, '\t');
		if (match_genotype && columns.size() < 10) {
			cout << "[VarMatch] Warning: not enough information in VCF file for genotype matching." << endl;
			cout << "[VarMatch] \tskip current variant " << filename << endl;
			continue;
		}
		if (chromosome_name == ".") chromosome_name = columns[0];
		auto pos = atoi(columns[1].c_str()) - 1;
		auto ref = columns[3];
		auto alt_line = columns[4];
		auto quality = columns[6];

		if (ref == ".") ref = "";
		if (alt_line == ".") alt_line = "";
		//decide which thread to use
		int thread_index = 0;
		for (int i = 0; i < pos_boundries.size(); i++) {
			if (pos < pos_boundries[i]) {
				thread_index = i;
				break;
			}
		}

		int genotype_index = -1;
		string genotype = "0/0";
		vector<string> genotype_columns;

		bool is_heterozygous_variant = false;
		bool is_multi_alternatives = false;

		if (columns.size() >= 10) {
			auto formats = split(columns[8], ':');
			for (int i = 0; i < formats.size(); i++) {
				if (formats[i] == "GT") {
					genotype_index = i;
					break;
				}
			}
			if (genotype_index < 0) {
				cout << "[VarMatch] Warning: not enough information in VCF file for genotype matching." << endl;
				cout << "[VarMatch] \tskip current variant " << filename << endl;
				continue;
			}
			auto additionals = split(columns[9], ':');
			genotype = additionals[genotype_index];

			if (genotype.find("/") != std::string::npos) {
				genotype_columns = split(genotype, '/');
			}
			else if (genotype.find("|") != std::string::npos) {
				genotype_columns = split(genotype, '|');
			}
			else {
				cout << "[VarMatch] Error: Unrecognized Genotype: " << genotype << endl;
				continue;
			}

			// normalize format of genotype: sorted, separated by |
			if (genotype_columns.size() != 2) {
				cout << "[VarMatch] Warning Unrecognized Genotype: " << genotype << endl;
			}
			else {
				sort(genotype_columns.begin(), genotype_columns.end());
				genotype = genotype_columns[0] + "|" + genotype_columns[1];
				if (genotype_columns[0] != genotype_columns[1]) {
					is_heterozygous_variant = true;
				}
			}
		}
		if (match_genotype && genotype == "0|0") {
			continue;
		}
		vector<string> alt_list;
		if (alt_line.find(",") != std::string::npos) {
			alt_list = split(alt_line, ',');
			is_multi_alternatives = true;
		}
		else {
			alt_list.push_back(alt_line);
		}
		
		vector<char> var_types;
		for (auto alt_it = alt_list.begin(); alt_it != alt_list.end(); ++alt_it) {
			string alt = *alt_it;
			char snp_type = 'S';
			if ((int)ref.length() > (int)alt.length()) {
				snp_type = 'D';
			}
			else if ((int)ref.length() < (int)alt.length()) {
				snp_type = 'I';
			}
			var_types.push_back(snp_type);
		}

		DiploidVariant dv(pos, var_types, ref, alt_list, genotype, is_heterozygous_variant, is_multi_alternatives);
		if (normalization) {
			NormalizeDiploidVariant(dv);
		}

		if (pos_2_var[thread_index].find(pos) != pos_2_var[thread_index].end()) {
			dout << "[VarMatch] Warning: overlap variants detected in " << filename << endl;
		}
		pos_2_var[thread_index][pos] = dv;
	}
	vcf_file.close();
	return true;
}

//private
void DiploidVCF::DirectSearchInThread(unordered_map<int, DiploidVariant> & ref_snps, unordered_map<int, DiploidVariant> & query_snps, int thread_index) {
	// handle heterozygous variants
	auto rit = ref_snps.begin();
	auto rend = ref_snps.end();
	for (; rit != rend;) {
		auto r_pos = rit->first;
		DiploidVariant r_var = rit->second;
		auto qit = query_snps.find(r_pos);
		if (qit != query_snps.end()) {
			DiploidVariant q_var = qit->second;
			if (r_var == q_var) {
				string matching_result = chromosome_name + '\t' + to_string(r_var.pos + 1) + "\t" + r_var.ref + "\t";
				auto alt_string = r_var.alts[0];
				if (r_var.multi_alts)
					alt_string += "," + r_var.alts[1];
				matching_result += alt_string;
				direct_match_records[thread_index]->push_back(matching_result);
				rit = ref_snps.erase(rit);
				query_snps.erase(qit);
			}
			else {
				++rit;
			}
		}
		else {
			++rit;
		}
	}
}

// directly match by position
// private
void DiploidVCF::DirectSearchMultiThread() {

	direct_match_records = new vector<string>*[thread_num];
	for (int j = 0; j < thread_num; j++) {
		direct_match_records[j] = new vector<string>;
	}

	vector<thread> threads;
	//spawn threads
	unsigned i = 0;
	for (; i < thread_num - 1; i++) {
		threads.push_back(thread(&DiploidVCF::DirectSearchInThread, this, ref(refpos_2_var[i]), ref(querypos_2_var[i]), i));
	}
	// also you need to do a job in main thread
	// i equals to (thread_num - 1)
	if (i != thread_num - 1) {
		dout << "[Error] thread number not match" << endl;
	}
	DirectSearchInThread(refpos_2_var[i], querypos_2_var[i], i);

	// call join() on each thread in turn before this function?
	std::for_each(threads.begin(), threads.end(), std::mem_fn(&std::thread::join));

	threads.clear();

	ofstream output_simple_file;
	output_simple_file.open(output_simple_filename);
	output_simple_file << "##VCF1:" << ref_vcf_filename << endl;
	output_simple_file << "##VCF2:" << que_vcf_filename << endl;
	output_simple_file << "#CHROM\tPOS\tREF\tALT" << endl;
	for (int i = 0; i < thread_num; i++) {
		for (int j = 0; j < direct_match_records[i]->size(); j++) {
			output_simple_file << direct_match_records[i]->at(j) << endl;
		}
	}
	output_simple_file.close();
	for (int j = 0; j < thread_num; j++) {
		delete direct_match_records[j];
	}
	delete[] direct_match_records;
}

bool DiploidVCF::VariantMatch(vector<DiploidVariant> & variant_list, int thread_index) {
	sort(variant_list.begin(), variant_list.end());
	map<int, DiploidVariant> separate_pos_var[2];
	bool separate_contians_indel[2];
	// separate into ref and que
	int min_pos = genome_sequence.length() + 1;
	int max_pos = -1;
	for (int i = 0; i < variant_list.size(); i++) {
		int flag = variant_list[i].flag; // flag indicate if the variant is from ref set or query set
		int pos = variant_list[i].pos;
		separate_pos_var[flag][pos] = variant_list[i];
		auto ref_sequence = variant_list[i].ref;
		auto alt_sequences = variant_list[i].alts;

		min_pos = min(pos, min_pos);
		max_pos = max((int)(pos + ref_sequence.length()), max_pos);

		if (ref_sequence.length() != alt_sequences[0].length())
			separate_contians_indel[flag] = true;
		if (variant_list[i].multi_alts) {
			if (ref_sequence.length() != alt_sequences[1].length()) {
				separate_contians_indel[flag] = true;
			}
		}
	}

	min_pos = max(min_pos - 1, 0);
	max_pos = min(max_pos + 1, (int)genome_sequence.length());

	if (!separate_contians_indel[0] && !separate_contians_indel[1]) {
		// There is no way that there will be a match
		return false;
	}

	string subsequence = genome_sequence.substr(min_pos, max_pos-min_pos);
	int offset = min_pos;
	// 0 for ref, 1 for query, same as flag
    map<int, int> choices[4];
	for(int i = 0; i < 2; i++){
        for(int j = 0; j < 2; j++){
            for(auto it = separate_pos_var[i].begin(); it != separate_pos_var[i].end(); ++it){
                auto pos = it->first;
                choices[i*2+j][pos] = -1;
            }
        }
    }
	map<int, int> max_matches[4];
	int max_score = 0;
	MatchWithIndel(variant_list,
		subsequence,
		offset,
		0,
		separate_pos_var,
		choices,
		max_matches,
		max_score);
	if (max_score == 0) {
		return false;
	}
	
	// here can construct matching result
	//for(int i = 0; i < 2; i++){
    //    for(int j = 0; j < 2; j++){
    //        if(i == 0){
    //            dout << "ref: " << endl;
    //        }else{
    //            dout << "alt: " << endl;
    //        }
    //        auto c = max_matches[i*2+j];
    //        for(auto it = c.begin(); it !=c.end(); ++it){
    //            dout << it->first << "," << it->second << " ";
    //        }
   //         dout << endl;
    //    }
    //}

	//vector<int> complex_ref_match_num;
	//vector<int> complex_que_match_num;
	map<int, bool> separate_pos_matched[2];
	for (int i = 0; i < 2; i++) {
		for (auto it = separate_pos_var[i].begin(); it != separate_pos_var[i].end(); ++it) {
			separate_pos_matched[i][it->first] = false;
		}
	}
	for(int i = 0; i < 2; i++){
		// i= 0 ref, =1 alt
	    for(int j = 0; j < 2; j++){
	        auto c = max_matches[i*2+j];
	        for(auto it = c.begin(); it !=c.end(); ++it){
				if (it->second > 0) {
					separate_pos_matched[i][it->first] = true;
				}
	        }
	    }
	}
	for (auto it = separate_pos_matched[0].begin(); it != separate_pos_matched[0].end(); ++it) {
		if (it->second) {
			complex_ref_match_num[thread_index] ++;
		}
	}
	for (auto it = separate_pos_matched[1].begin(); it != separate_pos_matched[1].end(); ++it) {
		if (it->second) {
			complex_que_match_num[thread_index] ++;
		}
	}
	cout << complex_ref_match_num[thread_index] << "," << complex_que_match_num[thread_index] << endl;
    return true;
}

void DiploidVCF::MatchWithIndel(vector<DiploidVariant> & variant_list,
	const string subsequence,
	const int offset,
	int index,
	map<int, DiploidVariant> separate_pos_var [],
	map<int, int> choices [], // 4 vectors
	map<int, int> max_matches[],  // 4 vectors
	int & max_score) {

    int prefix_match = CheckPrefix(subsequence, offset, separate_pos_var, choices);
    if (prefix_match < 0) return;
	// if prefix_match == 0, just prefix match
	if (prefix_match > 0) { // sequence direct match
		int score = prefix_match;
		if (max_score < score) {
			//cout << "higher score: " << score << endl;
            max_score = score; 
			for (int i = 0; i < 4; i++) {
				max_matches[i] = choices[i];
			}
		}
	}
	if (index >= variant_list.size()) return;
	auto variant = variant_list[index];
	int flag = variant.flag;
	int pos = variant.pos;
	int choice_end = 1;
	if (variant.multi_alts) choice_end = 2;
	for (int choice = 0; choice <= choice_end; choice++) {
		if(pos == 0){
            dout << "error pos = 0 " << endl;
        }
        choices[flag * 2][pos] = choice;
		if (choice == 0) {
			choices[flag * 2 + 1][pos] = 0;
		}
		else if (choice == 1) { // include
			if (variant.multi_alts) { // if multi_alts, then the other alleles should be included
				choices[flag * 2 + 1][pos] = 2;
			}
			else if (variant.heterozygous) { // if heterozygous but not multi_alts, then reference should be included
				choices[flag * 2 + 1][pos] = 0;
			}
			else { // homozygous one
				choices[flag * 2 + 1][pos] = 1;
			}
		}
		else {
			choices[flag * 2 + 1][pos] = 1;
		}

		MatchWithIndel(variant_list,
			subsequence,
			offset,
			index + 1,
			separate_pos_var,
			choices,
			max_matches,
			max_score);

		choices[flag * 2][pos] = -1;
		choices[flag * 2 + 1][pos] = -1;
	}
}

inline bool DiploidVCF::CompareSequence(string s1, string s2) {
	transform(s1.begin(), s1.end(), s1.begin(), ::toupper);
	transform(s2.begin(), s2.end(), s2.begin(), ::toupper);
	return s1 == s2;
}

// check if prefix match or equal
int DiploidVCF::CheckPrefix(const string subsequence,
	const int offset,
	map<int, DiploidVariant> separate_pos_var[],
	map<int, int> choices[])
{
	string paths[4] = { "" }; // 0 and 1 are ref, 2 and 3 are query path
	// create 4 paths
	for (int i = 0; i < 2; i++) {
		// create 
		for (int j = 0; j < 2; j++) {
			int index = i*2 + j;
			map<int, int> pos_choice = choices[index];
			string path = "";
			int start_pos = 0;
			auto it = pos_choice.begin();
            for (; it != pos_choice.end(); ++it) {
				int pos = it->first;
				int choice = it->second;
				auto variant = separate_pos_var[i][pos];
				string ref = variant.ref;
				auto alts = variant.alts;
				int offset_pos = pos - offset;
				if (offset_pos < start_pos) {
					return -1;
				}
				else if (offset_pos > start_pos) {
					path += subsequence.substr(start_pos, offset_pos - start_pos);
				}
				if(choice < 0)
                    break;
				if (choice == 0) {
					path += ref;
				}
				else if (choice == 1) {
					path += alts[0];
				}
				else {
					path += alts[1];
				}
                start_pos = max(start_pos, offset_pos + (int)ref.length());
			}
            if(it == pos_choice.end()){
                if(start_pos < subsequence.length()){
                    path += subsequence.substr(start_pos, subsequence.length()-start_pos);
                }
            }
			paths[index] = path;
		}
	}
    
	// check prefix match
	int const comb[2][4] = {
		{1,3,2,4},
		{1,4,2,3}
	};
	
	bool prefix_match = false;
	bool direct_match = false;
	for (int i = 0; i < 2; i++) {
		bool check_prefix_match[2] = { false };
		bool check_direct_match[2] = { false };
		for (int k = 0; k < 2; k++) {
			string s1 = paths[comb[i][k * 2]-1];
			string s2 = paths[comb[i][k * 2 + 1]-1];
            int min_len = min(s1.length(), s2.length());
			string s1_sub = s1.substr(0, min_len);
			string s2_sub = s2.substr(0, min_len);
			check_prefix_match[k] = CompareSequence(s1_sub, s2_sub);
			check_direct_match[k] = CompareSequence(s1, s2);
		}
		if (check_prefix_match[0] && check_prefix_match[1])
			prefix_match = true;
		if (check_direct_match[0] && check_direct_match[1])
			direct_match = true;
	}
	if (direct_match) {
        //for(int i = 0; i < 4; i++){
        //    dout << paths[i] << endl;
        //}
        //dout << endl;

		int score = 0;
		for (int i = 0; i < 2; i++) {
			auto pos_var = separate_pos_var[i];
			for (auto it = pos_var.begin(); it != pos_var.end(); ++it) {
				if(choices[i*2][it->first] <= 0 && choices[i*2+1][it->first] <= 0){
                    continue;
                }
                if (scoring_basepair) {
					score += it->second.ref.length();
				}
				else {
					score += 1;
				}
			}
		}
		return score;
	}
	if (prefix_match) return 0;
    return -1;
}

int DiploidVCF::test() {
	genome_sequence = "GTCAGCCGG";
	DiploidVariant d1(1, vector<char> ({'S', 'S'}), "T", vector<string> ({"A", "C"}), "1/2", true, true, 0);
	DiploidVariant d2(4, vector<char> ({'S', 'S'}), "G", vector<string> ({"C", ""}), "0/1", true, false, 0);
	DiploidVariant d3(5, vector<char> ({'S', 'S'}), "C", vector<string> ({"T", ""}), "0/1", true, false, 0); // this is false negative
	DiploidVariant d4(6, vector<char> ({'S', 'S'}), "C", vector<string> ({"G", ""}), "0/1", true, false, 0);
	DiploidVariant d5(7, vector<char> ({'S', 'S'}), "G", vector<string> ({"A", ""}), "0/1", true, false, 0);
	DiploidVariant d6(1, vector<char> ({'S', 'S'}), "T", vector<string> ({"A", "C"}), "1/2", true, true, 1);
	DiploidVariant d7(3, vector<char> ({'D', 'S'}), "AG", vector<string> ({"A", ""}), "0/1", true, false, 1);
	DiploidVariant d8(7, vector<char> ({'I', 'S'}), "G", vector<string> ({"GA", ""}), "0/1", true, false, 1);
    
    complex_ref_match_num.push_back(0);
    complex_que_match_num.push_back(0);

	vector<DiploidVariant> var_list = { d1,d2,d3,d4,d5,d6,d7,d8 };
	cout << VariantMatch(var_list, 0) << endl;
	return 0;
}


void DiploidVCF::ClusteringVariants() {

	vector<DiploidVariant> variant_list;
	// in DiploidVariant, flag = 0 is reference, flag = 1 is query
	for (int i = 0; i < refpos_2_var.size(); i++) {
		auto & m = refpos_2_var[i];
		for (auto it = m.begin(); it != m.end(); ++it) {
			auto v = it->second;
			if (v.flag != 0) {
				v.flag = 0;
			}
			variant_list.push_back(v);
		}
	}

	for (int i = 0; i < querypos_2_var.size(); i++) {
		auto & m = querypos_2_var[i];
		for (auto it = m.begin(); it != m.end(); ++it) {
			auto v = it->second;
			v.flag = 1;
			variant_list.push_back(v);
		}
	}

	if (variant_list.size() == 0)
		return;

	sort(variant_list.begin(), variant_list.end());

	int cluster_index = 0;
	int ins_len[2] = { 0 };
	int del_len[2] = { 0 };
	int c_start = 0;
	int c_end = 0;

	for (int i = 0; i < variant_list.size(); i++) {
		auto snp = variant_list[i];
		// check if need to separator clusters
		if (i > 0) {
			c_end = snp.pos;
			if (c_end - c_start >= 2) {
				string separator = genome_sequence.substr(c_start, c_end - c_start);
				int max_change = max(ins_len[0] + del_len[1], ins_len[1] + del_len[0]);
				if ((int)(separator.length()) > 2 * max_change &&
					((int)(separator.length()) > MAX_REPEAT_LEN || !CheckTandemRepeat(separator, max_change)))
				{
					cluster_index++;
					ins_len[0] = 0;
					del_len[0] = 0;
					ins_len[1] = 0;
					del_len[1] = 0;
					c_start = 0; // re-assign c_start
				}
			}
		}
		if (c_start < snp.pos + (int)(snp.ref.length())) c_start = snp.pos + (int)(snp.ref.length());

		// assign snp to cluster
		cluster_vars_map[cluster_index].push_back(snp);
		int ref_length = (int)(snp.ref.length());

		int flag = snp.flag;
		int alt0_length = snp.alts[0].length();
		int diff0_length = abs(ref_length - alt0_length);
		if (snp.multi_alts) {
			int alt1_length = snp.alts[1].length();
			int diff1_length = abs(ref_length - alt1_length);
			if (snp.var_types[0] == snp.var_types[1]) {
				int diff_length = max(diff0_length, diff1_length);
				if (snp.var_types[0] == 'I') {
					ins_len[flag] += diff_length;
				}
				else if (snp.var_types[0] == 'D') {
					del_len[flag] += diff_length;
				}
			}
			else {
				if (snp.var_types[0] == 'I') {
					ins_len[flag] += diff0_length;
				}
				else if (snp.var_types[0] == 'D') {
					del_len[flag] += diff0_length;
				}
				if (snp.var_types[1] == 'I') {
					ins_len[flag] += diff1_length;
				}
				else if (snp.var_types[1] == 'D') {
					del_len[flag] += diff1_length;
				}
			}
		}
		else {
			if (snp.var_types[0] == 'I') {
				ins_len[flag] += diff0_length;
			}
			else if (snp.var_types[0] == 'D') {
				del_len[flag] += diff0_length;
			}
		}
	}
}

// private
bool DiploidVCF::ClusteringMatchInThread(int start, int end, int thread_index) {
	for (int cluster_id = start; cluster_id < end; cluster_id++) {
		if (cluster_vars_map.find(cluster_id) != cluster_vars_map.end()) {
			auto & var_list = cluster_vars_map[cluster_id];
			if (var_list.size() > 20) {
				cout << "cluster to large" << endl;
				continue;
			}
			VariantMatch(var_list, thread_index);
		}
	}
}

// private
void DiploidVCF::ClusteringMatchMultiThread() {
	clustering_search = true;
	int start = cluster_vars_map.begin()->first;
	int cluster_number = cluster_vars_map.size();
	int cluster_end_boundary = start + cluster_number;
	int cluster_step = cluster_number / thread_num;
	if (cluster_step * thread_num < cluster_number) cluster_step++;
	int end = start + cluster_step;
	//initialize vector size, each allocating will have a lock
	complex_match_records = new vector<string>*[thread_num];
	for (int j = 0; j < thread_num; j++) {
		complex_match_records[j] = new vector<string>;
		complex_ref_match_num.push_back(0);
		complex_que_match_num.push_back(0);
	}
	
	vector<thread> threads;
	//spawn threads
	unsigned i = 0;
	for (; i < thread_num - 1; i++) {
		int variant_number = 0;
		for (int cluster_id = start; cluster_id < end; cluster_id++) {
			if (cluster_vars_map.find(cluster_id) != cluster_vars_map.end()) {
				variant_number += cluster_vars_map[cluster_id].size();
			}
		}
		threads.push_back(thread(&DiploidVCF::ClusteringMatchInThread, this, start, end, i));
		start = end;
		end = start + cluster_step;
	}
	// also you need to do a job in main thread
	// i equals to (thread_num - 1)
	if (i != thread_num - 1) {
		dout << "[Error] thread number not match" << endl;
	}
	if (start >= cluster_vars_map.size()) {
		dout << "[Error] index out of map range" << endl;
	}
	else {
		int variant_number = 0;
		for (int cluster_id = start; cluster_id < end; cluster_id++) {
			if (cluster_vars_map.find(cluster_id) != cluster_vars_map.end()) {
				variant_number += cluster_vars_map[cluster_id].size();
			}
		}
		ClusteringMatchInThread(start, end, i);
	}

	// call join() on each thread in turn before this function?
	std::for_each(threads.begin(), threads.end(), std::mem_fn(&std::thread::join));

	ofstream output_complex_file;
	output_complex_file.open(output_complex_filename);
	output_complex_file << "##VCF1:" << ref_vcf_filename << endl;
	output_complex_file << "##VCF2:" << que_vcf_filename << endl;
	output_complex_file << "#CHROM\tPOS\tREF\tALT\tVCF1\tVCF2" << endl;
	for (int i = 0; i < thread_num; i++) {
		for (int j = 0; j < complex_match_records[i]->size(); j++) {
			if (complex_match_records[i]->at(j).find_first_not_of(' ') != std::string::npos) {
				output_complex_file << complex_match_records[i]->at(j);
			}
		}
	}
	output_complex_file.close();

	for (int j = 0; j < thread_num; j++) {
		delete complex_match_records[j];
	}
	delete[] complex_match_records;
}

// for public access
void DiploidVCF::Compare(string ref_vcf,
	string query_vcf,
	string genome_seq,
	bool direct_search,
	string output_prefix,
	bool match_genotype,
	bool normalization,
	bool score_basepair) {



}

