// author: Chen Sun, chensun@cse.psu.edu
#include "diploid.h"

// inline function protected
// code reviewed by Channing
inline bool CompareSequence(string s1, string s2) {
	transform(s1.begin(), s1.end(), s1.begin(), ::toupper);
	transform(s2.begin(), s2.end(), s2.begin(), ::toupper);
	return s1 == s2;
}

bool operator <(const DiploidVariant& x, const DiploidVariant& y) {
	return x.pos < y.pos;
}

bool operator ==(const DiploidVariant& x, const DiploidVariant& y) {
	if (x.pos == y.pos && CompareSequence(x.ref, y.ref) && x.heterozygous == y.heterozygous && x.multi_alts == y.multi_alts) {
		if (x.multi_alts && x.heterozygous) {
			int match_times = 0;
			for (int i = 0; i < 2; i++) {
				for (int j = 0; j < 2; j++) {
					if (CompareSequence(x.alts[i], y.alts[j]))
						match_times++;
				}
			}
			if (match_times >= 2)
				return true;
		}
		else if(CompareSequence(x.alts[0], y.alts[0])){
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
// [todo] unit test normalization
// normalization modifies vt normalize algorithm
// code reviewed by Channing 4/2/2016
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
			if ((parsimonious_ref.size() > 1 && parsimonious_alt0.size() > 1 && parsimonious_alt1.size() > 1) || left_index > 0) { // when left_index == 0, can not make further changes
				parsimonious_ref.pop_back();
				parsimonious_alt0.pop_back();
				parsimonious_alt1.pop_back();
				change_in_allels = true;
			}
            // else do not make further changes
		}
		if (parsimonious_ref.length() == 0 || parsimonious_alt0.length() == 0 || parsimonious_alt1.length() == 0) {
			left_index--;
			char left_char = genome_sequence[left_index];
			parsimonious_ref = left_char + parsimonious_ref;
			parsimonious_alt0 = left_char + parsimonious_alt0;
			parsimonious_alt1 = left_char + parsimonious_alt1;
		}
	}
	while (toupper(parsimonious_ref[0]) == toupper(parsimonious_alt0[0]) &&
            toupper(parsimonious_ref[0]) == toupper(parsimonious_alt1[0]) &&
            parsimonious_ref.size() > 1 &&
            parsimonious_alt0.size() > 1 &&
            parsimonious_alt1.size() > 1)
    {
		parsimonious_ref.erase(0, 1);
		parsimonious_alt0.erase(0, 1);
		parsimonious_alt1.erase(0, 1);
        left_index ++; // left_index indicates variant position, if truncate the leftmost, then
	}
	var.pos = left_index;
	var.ref = parsimonious_ref;
	var.alts[0] = parsimonious_alt0;
	if (var.heterozygous && var.multi_alts)
		var.alts[1] = parsimonious_alt1;
	return true;
}

// protected
// code reviewed by Channing and Succulent on 4/2/2016
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
		auto pos = atoi(columns[1].c_str()) - 1; // 0-based coordinate
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
			dout << "[VarMatch] Warning: different variants sharing the same positon detected in " << filename << endl;
		}
		//if(dv.pos == -1) cout << "@:" << line << endl;
		pos_2_var[thread_index][pos] = dv;
	}
	vcf_file.close();
	return true;
}

// protected override
// code reviewed by Channing and Succulent on 4/2/2016
void DiploidVCF::DecideBoundaries() {
	int genome_size = genome_sequence.size();

    if(genome_size == 0){
        dout << "[VarMatch] Warning: no genome sequence detected when decide boundries. " << endl;
    }

	int distance = genome_size / thread_num;
	for (int i = 0; i < thread_num - 1; i++) {
		pos_boundries.push_back((i + 1)*distance);
	}
	pos_boundries.push_back(genome_size);

	for (int i = 0; i < thread_num; i++) {
		refpos_2_var.push_back(unordered_map<int, DiploidVariant>());
		querypos_2_var.push_back(unordered_map<int, DiploidVariant>());
	}

	boundries_decided = true;
}

//private
void DiploidVCF::DirectSearchInThread(unordered_map<int, DiploidVariant> & ref_snps,
                                      unordered_map<int, DiploidVariant> & query_snps,
                                      int thread_index) {
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

bool DiploidVCF::RecurrentVariantMatch(vector<DiploidVariant> & variant_list, int thread_index) {
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
	if (separate_pos_var[0].size() == 0 || separate_pos_var[1].size() == 0) {
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
	string max_paths[2];
	int max_score = 0;
	RecurrentMatchWithIndel(variant_list,
		subsequence,
		offset,
		0,
		separate_pos_var,
		choices,
		max_matches,
		max_score,
		max_paths);
	if (max_score == 0) {
		return false;
	}

	// matched, print out matches
	bool multiple_match = true;
	if (CompareSequence(max_paths[1], subsequence) || CompareSequence(max_paths[1], max_paths[0])) {
		multiple_match = false;
	}
	string alt_record = max_paths[0];
	if (multiple_match)
		alt_record += "/" + max_paths[1];
	string match_record = chromosome_name + "\t" + to_string(offset) + "\t" + subsequence + "\t" + alt_record;
	string vcf_record[2] = { "" };
	string phase_record[4] = { "" };

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
	for (int i = 0; i < 2; i++) {
		auto final_iter = separate_pos_matched[i].end();
		--final_iter;
		for (auto it = separate_pos_matched[i].begin(); it != separate_pos_matched[i].end(); ++it) {
			if (it->second) {
				int pos = it->first;
				DiploidVariant variant = separate_pos_var[i][pos];
				string alt1_string = variant.alts[0];
				if (variant.multi_alts) {
					alt1_string = variant.alts[1];
				}
				else if(! variant.heterozygous) {
					alt1_string = variant.ref;
				}
				string variant_record = to_string(pos) + "," + variant.ref + "," + variant.alts[0];
				if (multiple_match)
					variant_record += "/" + alt1_string;
				vcf_record[i] += variant_record;
				//cout << pos << ":" << max_matches[i*2+1][pos] << endl;
				phase_record[i * 2] += to_string(max_matches[i * 2][pos]);
				phase_record[i * 2 + 1] += to_string(max_matches[i * 2 + 1][pos]);
				if (it != final_iter) {
					vcf_record[i] += ";";
					phase_record[i * 2] += ",";
					phase_record[i * 2 + 1] += ",";
				}
			}
		}
	}
	match_record += "\t" + vcf_record[0] + "\t" + vcf_record[1];
	if (multiple_match) {
		match_record += "\t" + phase_record[0] + "/" + phase_record[1] + "\t" + phase_record[2] + "/" + phase_record[3];
	}
	else {
		match_record += "\t.\t.";
	}
	match_record += "\t" + to_string(max_score) + "\n";
	cout << match_record ;

	for (int i = 0; i < 2; i++)
	{
		if (i == 0) {
			cout << "ref: ";
		}
		else {
			cout << "alt: ";
		}
		cout << separate_pos_var[i].size() << endl;
		for (auto it = separate_pos_var[i].begin(); it != separate_pos_var[i].end(); ++it) {
			auto v = it->second;
			cout << v.pos << "," << v.ref << "," << v.alts[0];
			if (v.multi_alts) {
				cout << v.alts[1];
			}
			cout << ";";
		}
		cout << endl;
	}
	cout << endl;

	complex_match_records[thread_index]->push_back(match_record);
    return true;
}

void DiploidVCF::RecurrentMatchWithIndel(vector<DiploidVariant> & variant_list,
	const string subsequence,
	const int offset,
	int index,
	map<int, DiploidVariant> separate_pos_var [],
	map<int, int> choices [], // 4 vectors
	map<int, int> max_matches[],  // 4 vectors
	int & max_score,
	string max_paths[]) {

	string cur_paths[2];
    int prefix_match = CheckPrefix(subsequence, offset, separate_pos_var, choices, cur_paths);
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
			for (int i = 0; i < 2; i++) {
				max_paths[i] = cur_paths[i];
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

		RecurrentMatchWithIndel(variant_list,
			subsequence,
			offset,
			index + 1,
			separate_pos_var,
			choices,
			max_matches,
			max_score,
			max_paths);

		choices[flag * 2][pos] = -1;
		choices[flag * 2 + 1][pos] = -1;
	}
}

// check if prefix match or equal
int DiploidVCF::CheckPrefix(const string subsequence,
	const int offset,
	map<int, DiploidVariant> separate_pos_var[],
	map<int, int> choices[],
	string cur_paths[])
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
					//return -1;
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
        for(int i = 0; i < 4; i++){
            dout << paths[i] << endl;
        }
        dout << endl;

		int score = 0;
		for (int i = 0; i < 2; i++) {
			cur_paths[i] = paths[i];
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

// code reviewed by Channing 4/3/2016
vector<vector<vector<int>>> DiploidVCF::Combine(vector<int> & positions, vector<bool> & multi_indicators, int k) {
	vector<vector<int>> sol;
	vector<vector<vector<int>>> all_sol;
	if (k == 0 || k > positions.size()) {
		return all_sol;
	}
	FindComb(positions,
		multi_indicators,
		0,
		k,
		sol,
		all_sol);
	return all_sol;
}

// code review by Channing 4/3/2016
// [TODO] unit test
void DiploidVCF::FindComb(vector<int> & positions,
	vector<bool> & multi_indicators,
	int start,
	int k,
	vector<vector<int> > & sol,
	vector<vector<vector<int>>> & all_sol)
{
	if (k == 0) {
		all_sol.push_back(sol);
		return;
	}
	int n = positions.size();
	for (int i = start; i <= n - k; i++) {
		sol.push_back(vector<int>({ positions[i], 0 }));
		FindComb(positions, multi_indicators, i + 1, k - 1, sol, all_sol);
		sol.pop_back();
		if (multi_indicators[i]) { // try second allele
			sol.push_back(vector<int>({ positions[i], 1 }));
			FindComb(positions, multi_indicators, i + 1, k - 1, sol, all_sol);
			sol.pop_back();
		}
	}
}

// code reviewed by Chen on 4/4/2016
bool DiploidVCF::VariantMatch(vector<DiploidVariant> & variant_list, int thread_index) {
    if(variant_list.size() <= 1) return false;
	sort(variant_list.begin(), variant_list.end());
	map<int, DiploidVariant> separate_pos_var[2];
	bool separate_contians_indel[2];
	// separate into ref and que
	int min_pos = genome_sequence.length() + 1;
	int max_pos = -1;
	for (int i = 0; i < variant_list.size(); i++) {
		int flag = variant_list[i].flag; // flag indicate if the variant is from ref set(0) or query set(1)
		int pos = variant_list[i].pos;
		separate_pos_var[flag][pos] = variant_list[i];
		auto ref_sequence = variant_list[i].ref;
		auto alt_sequences = variant_list[i].alts;

		min_pos = min(pos, min_pos);
		max_pos = max((int)(pos + ref_sequence.length()), max_pos);
        //dout << pos << "," << ref_sequence << "," << alt_sequences[0] << "," << flag << endl;
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
		//dout << "[VarMatch] Warning : skip this cluster as no indel detected" << endl;
		return false;
	}
	if (separate_pos_var[0].size() == 0 || separate_pos_var[1].size() == 0) {
		return false;
	}

	string subsequence = genome_sequence.substr(min_pos, max_pos - min_pos);
	int offset = min_pos;
	vector<vector<int>> max_choices[4]; // -1 for ref, 0 for alts[0], 1 for alts[1] (only applied to multi_alts)
	string max_paths[2];
	int max_score = 0;
	bool max_heterozygosity = false;
    FindBestMatch(variant_list,
                  subsequence,
                  offset,
                  0,
                  separate_pos_var,
                  max_choices,
                  max_score,
                  max_heterozygosity,
                  max_paths);

	if (max_score == 0) {
		return false;
	}

	// matched, print out matches
	bool multiple_match = max_heterozygosity;

    vector<string> alt_list;
    vector<char> var_types;
    string temp_gt = "1/1";
    alt_list.push_back(max_paths[0]);
    if(multiple_match)
        alt_list.push_back(max_paths[1]);
	DiploidVariant dv(offset, var_types, subsequence, alt_list, temp_gt, true, multiple_match);
	NormalizeDiploidVariant(dv);

	string alt_record = dv.alts[0];
	if (multiple_match)
		alt_record += "/" + dv.alts[1];
	string match_record = chromosome_name + "\t" + to_string(dv.pos+1) + "\t" + dv.ref + "\t" + alt_record;
	string vcf_record[2] = { "" };
	string phase_record[4] = { "" };

	for (auto it = max_choices[0].begin(); it != max_choices[0].end(); ++it) {
        complex_ref_match_num[thread_index] ++;
	}
	for (auto it = max_choices[2].begin(); it != max_choices[2].end(); ++it) {
        complex_que_match_num[thread_index] ++;
	}
	for (int i = 0; i < 2; i++) {
		auto final_iter = max_choices[i*2].size()-1;
		for (int k = 0; k < max_choices[i*2].size(); k++) {
            int pos = max_choices[i*2][k][0];
            DiploidVariant variant = separate_pos_var[i][pos];
            string alt1_string = variant.alts[0];
            if (variant.multi_alts) {
                alt1_string = variant.alts[1];
            }
            else if (variant.heterozygous) {
                alt1_string = variant.ref;
            }
            string variant_record = to_string(pos+1) + "," + variant.ref + "," + variant.alts[0];
            if (multiple_match)
                variant_record += "/" + alt1_string;
            vcf_record[i] += variant_record;
            //cout << pos << ":" << max_matches[i*2+1][pos] << endl;
            if(multiple_match){
                phase_record[i * 2] += to_string(max_choices[i*2][k][1]+1);
                phase_record[i * 2 + 1] += to_string(max_choices[i * 2 + 1][k][1]+1);
            }
            if (k != final_iter) {
                vcf_record[i] += ";";
                if(multiple_match){
                    phase_record[i * 2] += ",";
                    phase_record[i * 2 + 1] += ",";
                }
            }
		}
	}
	match_record += "\t" + vcf_record[0] + "\t" + vcf_record[1];
	if (multiple_match) {
		match_record += "\t" + phase_record[0] + "/" + phase_record[1] + "\t" + phase_record[2] + "/" + phase_record[3];
	}
	else {
		match_record += "\t.\t.";
	}
	match_record += "\t" + to_string(max_score) + "\n";
	// dout << match_record;

	// for (int i = 0; i < 2; i++)
	// {
	// 	if (i == 0) {
	// 		dout << "ref: ";
	// 	}
	// 	else {
	// 		dout << "alt: ";
	// 	}
	// 	dout << separate_pos_var[i].size() << endl;
	// 	for (auto it = separate_pos_var[i].begin(); it != separate_pos_var[i].end(); ++it) {
	// 		auto v = it->second;
	// 		dout << v.pos << "," << v.ref << "," << v.alts[0];
	// 		if (v.multi_alts) {
	// 			dout << v.alts[1];
	// 		}
	// 		dout << ";";
	// 	}
	// 	dout << endl;
	// }
	// dout << endl;

	complex_match_records[thread_index]->push_back(match_record);
	return true;
}

// code reviewed by Chen on 4/3/2016
bool DiploidVCF::FindBestMatch(vector<DiploidVariant> & variant_list,
	const string subsequence,
	const int offset,
	int index,
	map<int, DiploidVariant> separate_pos_var[],
	vector<vector<int>> max_choices[],  // 4 vectors
	int & max_score,
	bool & max_heterozygosity,
	string max_paths[])
{
	vector<int> positions[2]; // 0 from ref, 1 from query
	vector<bool> indicators[2]; // 0 from ref, 1 from query, indicate if multi_alts(true) or not(false)
	for (int i = 0; i < 2; i++) {
		for (auto it = separate_pos_var[i].begin(); it != separate_pos_var[i].end(); ++it) {
			auto v = it->second;
			positions[i].push_back(v.pos);
			indicators[i].push_back(v.multi_alts);
		}
	}
	// construct ref combinations in hash table, key is donor sequence
    unordered_map<string, vector<vector<int>> > seq_choice_ref;
    unordered_map<string, int> seq_score_ref; // corresponding score, if same key, store the one with highest score
	for (int i = 1; i <= positions[0].size(); i++) { // i : how many variants are chosen
		vector<vector<vector<int>>> ref_choice_list = Combine(positions[0], indicators[0], i);

		for (auto rit = ref_choice_list.begin(); rit != ref_choice_list.end(); ++rit) { // iterate all combinations with i variants
            // each combination is a vector of pairs(position, alt_index), alt_index is 0 or 1 (if multi_alts)
            string donor;
            int score;
            ModifyRefMultiVar(subsequence, offset, separate_pos_var[0], *rit, donor, score);
            if(CompareSequence(donor, subsequence)) continue;
            if(seq_choice_ref.find(donor) != seq_choice_ref.end() && seq_score_ref[donor] > score){
                continue;
            }else{
                // either overwrite or insert new
                seq_choice_ref[donor] = *rit;
                seq_score_ref[donor] = score;
            }
            //dout << "ref-donor: " << donor << endl;
		}
	}
	// now all combinations are stored in hash table seq_choice_ref
	// search query
    for(int i = 1; i <= positions[1].size(); i++){
            // iterate all combinations with i variants
        vector<vector<vector<int>>> que_choice_list = Combine(positions[1], indicators[1], i);
        for (auto qit = que_choice_list.begin(); qit != que_choice_list.end(); ++qit){
            string donor;
            int score;
            ModifyRefMultiVar(subsequence, offset, separate_pos_var[1], *qit, donor, score);
            if(CompareSequence(donor, subsequence)) continue;
            //dout << "que-donor: " << donor << endl;
            if(seq_choice_ref.find(donor) != seq_choice_ref.end()){
                // first check if there is heterozygous alleles
                int total_score = seq_score_ref[donor] + score;
                if(total_score <= max_score) continue;
                bool local_heter = false;
                bool local_multi = false;
                vector<vector<int>> ref_var_choices = seq_choice_ref[donor];
                vector<vector<int>> que_var_choices = *qit;
                vector<vector<int>> ref_other_choices;
                vector<vector<int>> que_other_choices;
                // check and construct heterozygous alleles
                for(int ri = 0; ri < ref_var_choices.size(); ri++){
                    int ref_pos = ref_var_choices[ri][0];
                    DiploidVariant ref_variant = separate_pos_var[0][ref_pos];
                    if (ref_variant.multi_alts){
                        local_multi = true;
                        ref_other_choices.push_back(vector<int>({ref_pos, 1 - ref_var_choices[ri][1]}));
                    }else if(ref_variant.heterozygous){
                        local_heter = true;
                        ref_other_choices.push_back(vector<int>({ref_pos,-1}));
                    }else{
                        ref_other_choices.push_back(vector<int>({ref_pos, ref_var_choices[ri][1]}));
                    }
                }
                // if not find heter, continue checking
                for(int qi = 0; qi < que_var_choices.size(); qi++){
                    int que_pos = que_var_choices[qi][0];
                    DiploidVariant que_variant = separate_pos_var[1][que_pos];
                    if(que_variant.multi_alts){
                        local_multi = true;
                        que_other_choices.push_back(vector<int>({que_pos, 1- que_var_choices[qi][1]}));
                    }else if (que_variant.heterozygous){
                        local_heter = true;
                        que_other_choices.push_back(vector<int>({que_pos, -1}));
                    }else{
                        que_other_choices.push_back(vector<int>({que_pos, que_var_choices[qi][1]}));
                    }
                }


                if(local_multi){
                    // also check the other chromosome matches
                    int temp_score;
                    string ref_other_donor;
                    ModifyRefMultiVar(subsequence, offset, separate_pos_var[0], ref_other_choices, ref_other_donor, temp_score);
                    string que_other_donor;
                    ModifyRefMultiVar(subsequence, offset, separate_pos_var[1], que_other_choices, que_other_donor, temp_score);
                    if(CompareSequence(ref_other_donor, que_other_donor)){
                        max_choices[0] = ref_var_choices;
                        max_choices[1] = ref_other_choices;
                        max_choices[2] = que_var_choices;
                        max_choices[3] = que_other_choices;
                        max_paths[0] = donor;
                        max_paths[1] = ref_other_donor;
                        max_score = total_score;
                        max_heterozygosity = true;
                    }
                }else if(local_heter){
                    // also check the other chromosome matches
                    int temp_score;
                    string ref_other_donor;
                    ModifyRefMultiVar(subsequence, offset, separate_pos_var[0], ref_other_choices, ref_other_donor, temp_score);
                    string que_other_donor;
                    ModifyRefMultiVar(subsequence, offset, separate_pos_var[1], que_other_choices, que_other_donor, temp_score);
                    if(CompareSequence(ref_other_donor, que_other_donor)){
                        max_choices[0] = ref_var_choices;
                        max_choices[2] = que_var_choices;
                        max_paths[0] = donor;
                        max_score = total_score;
                        max_heterozygosity = false;
                    }
                }else{
                    max_choices[0] = ref_var_choices;
                    max_choices[2] = que_var_choices;
                    max_paths[0] = donor;
                    max_score = total_score;
                    max_heterozygosity = false;
                    //delay construct optimal solution at the very end.
                }
            }
        }
    }
    if(max_score > 0) return true;
    return false;
}

// code reviewed by Channing on 4/3/2016
void DiploidVCF::ModifyRefMultiVar(const string & genome,
                                   int offset,
                                   map<int, DiploidVariant> & pos_var,
                                   vector<vector<int>> pos_choice,
                                   string & donor,
                                   int & score) {
    donor = genome;
    score = 0; // if return before end of function, score = 0
    int local_score = 0;
    transform(donor.begin(), donor.end(), donor.begin(), ::toupper);
    int start_pos = 0;
    std::sort(pos_choice.begin(), pos_choice.end(),
          [](const std::vector<int>& a, const std::vector<int>& b) {
            return a[0]>b[0];}); // sorted by position in reverse order
    for(int i = 0; i < pos_choice.size(); i++){
        assert(pos_choice[i].size() == 2);
        int pos = pos_choice[i][0];
        int alt_index = pos_choice[i][1];
        int offset_pos = pos - offset;
        DiploidVariant variant = pos_var[pos];
        if(alt_index > 0 && !variant.multi_alts){
            dout << "[VarMatch] Warning: modify reference genome with allele not exist" << endl;
            return;
        }
        int offset_end = offset_pos + (int) variant.ref.length();
        string alt = "";
        if(alt_index >= 0){
            alt = variant.alts[alt_index];
        }else{
            alt = variant.ref;
        }
        int donor_length = donor.length();
		if(offset_pos > donor_length || offset_end > donor_length){
            dout << "[VarMatch] Warning: overlapping variants detected." << endl; // the most reason is overlapping variants
            return;
		}
		donor = donor.substr(0, offset_pos) + alt + donor.substr(offset_end, donor_length - offset_end);
		if(scoring_basepair){
            local_score += variant.ref.size();
		}else{
            local_score++;
		}
    }
    transform(donor.begin(), donor.end(), donor.begin(), ::toupper);
    //only assign score here, if fail to change reference, score will be 0
    score = local_score;
    return;
}

//int DiploidVCF::test() {
//	genome_sequence = "GTCAGCCGG";
//	DiploidVariant d1(1, vector<char> ({'S', 'S'}), "T", vector<string> ({"A", "C"}), "1/2", true, true, 0);
//	DiploidVariant d2(4, vector<char> ({'S', 'S'}), "G", vector<string> ({"C", ""}), "0/1", true, false, 0);
//	DiploidVariant d3(5, vector<char> ({'S', 'S'}), "C", vector<string> ({"T", ""}), "0/1", true, false, 0); // this is false negative
//	DiploidVariant d4(6, vector<char> ({'S', 'S'}), "C", vector<string> ({"G", ""}), "0/1", true, false, 0);
//	DiploidVariant d5(7, vector<char> ({'S', 'S'}), "G", vector<string> ({"A", ""}), "0/1", true, false, 0);
//	DiploidVariant d6(1, vector<char> ({'S', 'S'}), "T", vector<string> ({"A", "C"}), "1/2", true, true, 1);
//	DiploidVariant d7(3, vector<char> ({'D', 'S'}), "AG", vector<string> ({"A", ""}), "0/1", true, false, 1);
//	DiploidVariant d8(7, vector<char> ({'I', 'S'}), "G", vector<string> ({"GA", ""}), "0/1", true, false, 1);
//
//    complex_ref_match_num.push_back(0);
//    complex_que_match_num.push_back(0);
//    complex_match_records = new vector<string>*[1];
//    complex_match_records[0] = new vector<string>;
//	//vector<DiploidVariant> var_list = { d2,d3,d4,d5,d7,d8 };
//	vector<DiploidVariant> var_list = { d1,d2,d3,d4,d5,d6,d7,d8 };
//	cout << VariantMatch(var_list, 0) << endl;
//	return 0;
//}

int DiploidVCF::test() {
	genome_sequence = "AATtgtgtg";

	DiploidVariant d1(1, vector<char>({ 'D', 'S' }), "AT", vector<string>({ "A", "AA" }), "1/2", true, true, 0);
	DiploidVariant d2(1, vector<char>({ 'D', 'S' }), "AT", vector<string>({ "A", "AT" }), "0/1", true, true, 1);
	DiploidVariant d3(7, vector<char>({ 'I', 'S' }), "T", vector<string>({ "AA", "" }), "0/1", true, false, 1);

	complex_ref_match_num.push_back(0);
	complex_que_match_num.push_back(0);
    complex_match_records = new vector<string>*[1];
    complex_match_records[0]= new vector<string>;
	//vector<DiploidVariant> var_list = { d2,d3,d4,d5,d7,d8 };
	vector<DiploidVariant> var_list = { d1,d2,d3 };
	cout << VariantMatch(var_list, 0) << endl;
	return 0;
}

// [TODO] unit test clustering variant position
// code reviewed by Chen on 4/4/2016
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
            //if(v.pos == -1) cout << "@@@@@@@@@@@@@" << endl;
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
		c_start = max(c_start, snp.pos + (int)snp.ref.length() );

		// assign snp to cluster
		//if(snp.pos == -1) cout << "@@@@@@@@@@@@@" << endl;
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

void DiploidVCF::DivisiveHierarchicalClustering(list<vector<DiploidVariant>> & snp_clusters){
    // I use list of vectors instead of vector of vectors, to take advantage of member func of list
    if(snp_clusters.size() == 0) return;
    bool flag = true;
    list<bool> potential_list;
    for(int i = 0; i < snp_clusters.size(); i++){
        potential_list.push_back(true);
    }
    while(flag){
        flag = false;
        int list_size = snp_clusters.size();
        for(int i = 0; i < list_size; i++){
            auto front_cluster = snp_clusters.front();
            auto front_posential = potential_list.front(); // record if this can be separated
            // at the very beginning, all clusters can be separated
            // all newly separated sub-clusters can be separated
            // if one cluster marked not separated, then it can never be separated

            snp_clusters.pop_front();
            potential_list.pop_front();
            if(front_cluster.size() == 1){
                snp_clusters.push_back(front_cluster);
                potential_list.push_back(false);
                continue;
            }

            if(! front_posential){
                snp_clusters.push_back(front_cluster);
                potential_list.push_back(front_posential);
                continue;
            }
            int max_start = -1;
            int max_end = -1;
            int max_length = -1;
            int start = front_cluster[0].pos + (int)front_cluster[0].ref.length();
            // find the largest gap, see if we can separate from that gap
            for(int k = 0; k < front_cluster.size(); k++){
                auto snp = front_cluster[k];
                auto snp_pos = snp.pos;
                if(max_length < snp_pos - start){
                    max_length = snp_pos - start;
                    max_start = start;
                    max_end = snp_pos;
                }
            }

            if(max_length <= 0){
                snp_clusters.push_back(front_cluster);
                potential_list.push_back(false);
                continue;
            }
            int left_ins = 0;
            int left_del = 0;
            int right_ins = 0;
            int right_del = 0;
            vector<DiploidVariant> left_snp_list;
            vector<DiploidVariant> right_snp_list;
            string separator = genome_sequence.substr(max_start, max_end-max_start);
            for(int k = 0; k < front_cluster.size(); k++){
                DiploidVariant snp = front_cluster[k];
                int snp_diff = abs((int)snp.ref.length() - (int)snp.alts[0].length());
                if(snp.multi_alts && snp.var_types[0] == snp.var_types[1]){
                    snp_diff = max(snp_diff, abs((int)snp.ref.length() - (int)snp.alts[1].length()));
                }
                if(snp.pos <= max_start){
                    if(snp.var_types[0] == 'I'){
                        left_ins += snp_diff;
                    }else if(snp.var_types[0] == 'D'){
                        left_del += snp_diff;
                    }
                    left_snp_list.push_back(snp);
                }else{
                    if(snp.var_types[0] == 'I'){
                        right_ins += snp_diff;
                    }else if(snp.var_types[0] == 'D'){
                        right_del += snp_diff;
                    }
                    right_snp_list.push_back(snp);
                }
                if(snp.multi_alts && snp.var_types[0] != snp.var_types[1]){
                    snp_diff =  abs((int)snp.ref.length() - (int)snp.alts[1].length());
                    if(snp.pos <= max_start){
                        if(snp.var_types[1] == 'I'){
                            left_ins += snp_diff;
                        }else if(snp.var_types[1] == 'D'){
                            left_del += snp_diff;
                        }
                        //left_snp_list.push_back(snp);
                    }else{
                        if(snp.var_types[1] == 'I'){
                            right_ins += snp_diff;
                        }else if(snp.var_types[1] == 'D'){
                            right_del += snp_diff;
                        }
                        //right_snp_list.push_back(snp);
                    }
                }
            }
            //check
            if(left_snp_list.size() == 0 || right_snp_list.size() == 0){
                snp_clusters.push_back(front_cluster);
                potential_list.push_back(false);
                continue;
            }

            vector<int> change_list = {left_ins, left_del, right_ins, right_del};
            int max_change = 0;
            for(int k = 0; k < change_list.size(); k++){
                if (max_change < change_list[k]) max_change = change_list[k];
            }
            if ((int)separator.length() > 2 * max_change &&
                    ((int)separator.length() > MAX_REPEAT_LEN || !CheckTandemRepeat(separator, max_change)))
            {
                flag = true;
                snp_clusters.push_back(left_snp_list);
                potential_list.push_back(true);
                snp_clusters.push_back(right_snp_list);
                potential_list.push_back(true);
            }else{
                snp_clusters.push_back(front_cluster);
                potential_list.push_back(false);
                continue;
            }
        }
    }
    return;
}

// private
// code reviewed
bool DiploidVCF::ClusteringMatchInThread(int start, int end, int thread_index) {
    // end exclusive
	for (int cluster_id = start; cluster_id < end; cluster_id++) {
		if (cluster_vars_map.find(cluster_id) != cluster_vars_map.end()) {
			auto & var_list = cluster_vars_map[cluster_id];
			if (var_list.size() <= 1) continue;
			if (var_list.size() > 20) {
				list<vector<DiploidVariant>> cluster_list;
				cluster_list.push_back(var_list);
				DivisiveHierarchicalClustering(cluster_list);
				for(auto lt = cluster_list.begin(); lt != cluster_list.end(); ++lt){
                    auto cur_cluster = * lt;
                    if(cur_cluster.size() > 20){
                        dout << "[VarMatch] Warning: Large size cluster skipped" << endl;
                        continue;
                    }
                    VariantMatch(cur_cluster, thread_index);
				}
			}else{
                VariantMatch(var_list, thread_index);
			}
		}
	}
}

// private
void DiploidVCF::ClusteringMatchMultiThread() {
	clustering_search = true;
	int start = cluster_vars_map.begin()->first; // start cluster id
	int cluster_number = cluster_vars_map.size(); // cluster number
	int cluster_end_boundary = start + cluster_number; // end cluster id, exclusive
	int cluster_step = cluster_number / thread_num; // assign clusters to threads
	if (cluster_step * thread_num < cluster_number) cluster_step++;
	int end = start + cluster_step;
	//initialize vector size
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
	output_complex_file << "#CHROM\tPOS\tREF\tALT\tVCF1\tVCF2\tPHASE1\tPHASE2\tSCORE" << endl;
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

	int total_ref_complex = 0;
	int total_que_complex = 0;
	for (int i = 0; i < complex_ref_match_num.size(); i++)
		total_ref_complex += complex_ref_match_num[i];
	for (int i = 0; i < complex_que_match_num.size(); i++)
		total_que_complex += complex_que_match_num[i];

	cout << "complex match: " << total_ref_complex << "," << total_que_complex << endl;
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

	ref_vcf_filename = ref_vcf;
	que_vcf_filename = query_vcf;
	this->match_genotype = match_genotype;
	this->normalization = normalization;
	this->scoring_basepair = score_basepair;
	output_stat_filename = output_prefix + ".stat";
	output_simple_filename = output_prefix + ".simple";
	output_complex_filename = output_prefix + ".complex";

	//------------read genome sequence and decide boundary according to thread number
	dsptime();
	dout << " Read genome sequence file... " << endl;
	ReadGenomeSequence(genome_seq);
	dsptime();
	dout << " Finish reading genome sequence file." << endl;
	//------------read ref and query vcf file
	dsptime();
	dout << " Read reference vcf file... " << endl;
	ReadRefVCF(ref_vcf);
	dsptime();
	dout << " Read query vcf file... " << endl;
	ReadQueryVCF(query_vcf);
	dsptime();
	dout << " Finish reading all vcf file." << endl;

	//------------check vcf entry number before matching
	//int ref_total_indel_num, que_total_indel_num;
	int ref_total_num = refpos_2_var[0].size();
	int que_total_num = querypos_2_var[0].size();
	dout << "total num: " << ref_total_num << "," << que_total_num << endl;

	//------------direct search
	dsptime();
	dout << " Direct search ... " << endl;
	DirectSearchMultiThread();
	dsptime();
	dout << " Finish direct search." << endl;
	//int ref_direct_left_indel_num, que_direct_left_indel_num;
	int ref_direct_left_num = refpos_2_var[0].size();
	int que_direct_left_num = querypos_2_var[0].size();
	int ref_direct_match_num = ref_total_num - ref_direct_left_num;
	int que_direct_match_num = que_total_num - que_direct_left_num;
	cout << "direct match: " << ref_direct_match_num << "," << que_direct_match_num << endl;

	if (direct_search) {
		return;
	}

	//-------------clustering search
	dsptime();
	dout << " Clustering snps ... " << endl;
	ClusteringVariants();
	dsptime();
	dout << " Finish clustering." << endl;
	dsptime();
	dout << " Clustering search ... " << endl;
	ClusteringMatchMultiThread();
	dsptime();
	dout << " Finish clustering search." << endl;
	return;
}

