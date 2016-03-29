#include "vcf.h"


bool operator <(const SNP& x, const SNP& y) {
	return x.pos < y.pos;
}

bool operator ==(const SNP& x, const SNP& y) {
	if (x.pos == y.pos && x.snp_type == y.snp_type && x.alt == y.alt && x.genotype == y.genotype) {
		return true;
	}
	return false;
}

VCF::VCF(int thread_num_)
{
    debug_f = 0;
	genome_sequence = "";
	boundries_decided = false;
	clustering_search = false;
	match_genotype = false;
    if (thread_num_ <= 0) {
		thread_num = 1;
	}
	else {
		thread_num = min(thread_num_, (int)thread::hardware_concurrency());
	}
	dout << "Thread Number: " << thread_num << endl;
    chromosome_name = ".";
}

VCF::~VCF()
{
}

bool VCF::NormalizeSnp(int pos, string ref, string alt, string & parsimonious_ref, string & parsimonious_alt) {
	parsimonious_ref = ref;
	parsimonious_alt = alt;
	int left_index = pos;
	if (genome_sequence.size() == 0) return false;
	//if (parsimonious_ref.size() == 1 || parsimonious_alt.size() == 1) return true;
	if (toupper(genome_sequence[left_index]) != toupper(parsimonious_ref[0])) {
		dout << "[Error] genome sequence, subsequence, offset does not match." << endl;
		return false;
	}
	bool change_in_allels = true;
	while (change_in_allels) {
		change_in_allels = false;
		if (toupper(parsimonious_ref.back()) == toupper(parsimonious_alt.back())) {
			if((parsimonious_ref.size() > 1 && parsimonious_alt.back() > 1) || left_index > 0){
                parsimonious_ref.pop_back();
				parsimonious_alt.pop_back();
				change_in_allels = true;
			}
			else {
				return false;
			}
		}
		if (parsimonious_ref.length() == 0 || parsimonious_alt.length() == 0) {
			left_index--;
			char left_char = genome_sequence[left_index];
			parsimonious_ref = left_char + parsimonious_ref;
			parsimonious_alt = left_char + parsimonious_alt;
		}
	}
	while (toupper(parsimonious_ref[0]) == toupper(parsimonious_alt[0]) && parsimonious_ref.size() > 1 && parsimonious_alt.size() > 1) {
		parsimonious_ref.erase(0, 1);
		parsimonious_alt.erase(0, 1);
	}
    //if(parsimonious_ref != ref){
    //    cout << ref << "," << alt << "," << parsimonious_ref << "," << parsimonious_alt << endl;
    //}
	return true;
}

void VCF::ReadVCF(string filename, SnpHash & pos_2_snp) {
	if (!boundries_decided) {
		dout << "[Error] VCF::ReadVCF cannot read vcf file before read genome file" << endl;
		return;
	}

	ifstream vcf_file;
	vcf_file.open(filename.c_str());
	if (!vcf_file.good()) {
		cout << "[VarMatch] Error: can not open vcf file" << endl;
		return;
	}

    if(normalization){
        dout << "normalize while read" << endl;
    }
    string previous_line;
	while (!vcf_file.eof()) { // alternative way is vcf_file != NULL
		string line;
		getline(vcf_file, line, '\n');
		// check ineligible lines
        //dout << line << endl;
		if ((int)line.length() <= 1) continue;
		if (line.find_first_not_of(' ') == std::string::npos) continue;

		if (line[0] == '#'){
			if(line[1] == '#') continue;
			auto head_names = split(line, '\t');
			if(match_genotype && head_names.size() < 10){
				cout << "[VarMatch] Warning: not enough information in VCF file for genotype matching." << endl;
				cout << "[VarMatch] \tVCF file name " << filename << endl;
				cout << "[VarMatch] \tAutomatically turn off genotype matching module." << endl;
				match_genotype = false;
			}
			continue;
		}
		auto columns = split(line, '\t');
		if(match_genotype && columns.size() < 10){
			cout << "[VarMatch] Warning: not enough information in VCF file for genotype matching." << endl;
			cout << "[VarMatch] \tskip current variant " << filename << endl;
			continue;
		}
		if(chromosome_name == ".") chromosome_name = columns[0];
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
		string genotype = "1/1";
		vector<string> genotype_columns;

		if (match_genotype){
			auto formats = split(columns[8], ':');
			for(int i = 0; i < formats.size(); i++){
				if(formats[i] == "GT"){
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
			
			if(genotype.find("/") != std::string::npos){
				genotype_columns = split(genotype, '/');
			}else if(genotype.find("|") != std::string::npos){
				genotype_columns = split(genotype, '|');
			}else{
				cout << "[VarMatch] Error: Unrecognized Genotype: " << genotype << endl;
				continue;
			}
			// normalize format of genotype: sorted, separated by |
			if(genotype_columns.size() != 2){
				cout << "[VarMatch] Warning Unrecognized Genotype: " << genotype << endl;
                //cout << genotype_columns.size() << endl;
                //dout << line << endl;
			}else{
				sort(genotype_columns.begin(), genotype_columns.end());
				genotype = genotype_columns[0]+"|"+genotype_columns[1];
			}
		}

		vector<string> alt_list;
		if (alt_line.find(",") != std::string::npos){
			alt_list = split(alt_line, ',');
		}else{
			alt_list.push_back(alt_line);
		}
		if(!match_genotype){
			for(auto alt_it = alt_list.begin(); alt_it != alt_list.end(); ++alt_it){
				string alt = *alt_it;
				char snp_type = 'S';
				if ((int)ref.length() > (int)alt.length()) {
					snp_type = 'D';
				}
				else if ((int)ref.length() < (int)alt.length()) {
					snp_type = 'I';
				}
				if (normalization) {
					string norm_ref, norm_alt;
					NormalizeSnp(pos, ref, alt, norm_ref, norm_alt);
					pos_2_snp[thread_index][pos].push_back(SNP(pos, snp_type, norm_ref, norm_alt, genotype));
				}
				else {
					pos_2_snp[thread_index][pos].push_back(SNP(pos, snp_type, ref, alt, genotype));
				}
			}
		}else{
			//append variants according to genotype
            if(genotype == "0|0") continue;
            vector<char> snp_type_list;
            for(int i = 0; i < alt_list.size(); i++){
            	string alt = alt_list[i];
				char snp_type = 'S';
				if ((int)ref.length() > (int)alt.length()) {
					snp_type = 'D';
				}
				else if ((int)ref.length() < (int)alt.length()) {
					snp_type = 'I';
				}
				snp_type_list.push_back(snp_type);
            }

			if (alt_list.size() > 1) {
				haplotype_matching_check[thread_index][pos] = 0;
			}
            
            int genotype_val = atoi(genotype_columns[0].c_str()) - 1;

            if(genotype_val >= 0){ // if genotype == -1, it is reference which does not need to be added
				if (genotype_val >= alt_list.size()) {
					cout << "[VarMatch] Warning: Unrecognized Genotype. " << genotype_val << endl;
				}
				if (normalization) {
					string norm_ref, norm_alt;
					NormalizeSnp(pos, ref, alt_list[genotype_val], norm_ref, norm_alt);
					pos_2_snp[thread_index][pos].push_back(SNP(pos, snp_type_list[genotype_val], norm_ref, norm_alt, genotype));
				}
				else {
					pos_2_snp[thread_index][pos].push_back(SNP(pos, snp_type_list[genotype_val], ref, alt_list[genotype_val], genotype));
				}
            }
            if(genotype_columns[0] != genotype_columns[1]){
            	// add another alt, one genotype corresponding to one alt
            	genotype_val = atoi(genotype_columns[1].c_str()) - 1;

            	if(genotype_val >= 0){
					if (genotype_val >= alt_list.size()) {
						cout << "[VarMatch] Warning: Unrecognized Genotype[2]. " << genotype_val << endl;
					}
					if (normalization) {
						string norm_ref, norm_alt;
						NormalizeSnp(pos, ref, alt_list[genotype_val], norm_ref, norm_alt);
						pos_2_snp[thread_index][pos].push_back(SNP(pos, snp_type_list[genotype_val], norm_ref, norm_alt, genotype));
					}
					else {
						pos_2_snp[thread_index][pos].push_back(SNP(pos, snp_type_list[genotype_val], ref, alt_list[genotype_val], genotype));
					}
            	}
            }
		}
	}
	vcf_file.close();
	return;
}

void VCF::ReadGenomeSequence(string filename) {
	ifstream genome_file;
	genome_file.open(filename.c_str());
	if (!genome_file.good()) {
		cout << "[VarMatch] can not open FASTA file: ";
		cout << filename << endl;
		return;
	}

	genome_sequence = "";

	while(!genome_file.eof()) {
		string line;
		getline(genome_file, line, '\n');
		if ((int)line.length() <= 1) continue;
		if (line[0] == '>') continue;
		genome_sequence += line;
	}
	genome_file.close();
	// boundries can get after knowing genome sequence.
	DecideBoundries();
	return;
}

void VCF::DecideBoundries() {
	int genome_size = genome_sequence.size();

	int distance = genome_size / thread_num;
	for (int i = 0; i < thread_num - 1; i++) {
		pos_boundries.push_back((i + 1)*distance);
	}
	pos_boundries.push_back(genome_size);

	// initialize two for copy
	unordered_map<int, vector<SNP> > ref_h;
	unordered_map<int, vector<SNP> > que_h;
	map<int, vector<SNP> > ref_m;
	map<int, vector<SNP> > que_m;

	for (int i = 0; i < thread_num; i++) {
		refpos_2_snp.push_back(ref_h);
		querypos_2_snp.push_back(que_h);
		refpos_snp_map.push_back(ref_m);
		querypos_snp_map.push_back(que_m);
	}

	boundries_decided = true;

}

void VCF::ReadRefVCF(string filename) {
	ReadVCF(filename, this->refpos_2_snp);
}

void VCF::ReadQueryVCF(string filename) {
	ReadVCF(filename, this->querypos_2_snp);
}

bool VCF::CompareSnps(SNP r, SNP q) {
	if(r.pos != q.pos) return false;

	// directly match genotype
	if(match_genotype && r.genotype != q.genotype) return false;
    auto ref_ref = r.ref;
	transform(ref_ref.begin(), ref_ref.end(), ref_ref.begin(), ::toupper);
	auto ref_alt = r.alt;
	transform(ref_alt.begin(), ref_alt.end(), ref_alt.begin(), ::toupper);
	auto que_ref = q.ref;
	transform(que_ref.begin(), que_ref.end(), que_ref.begin(), ::toupper);
	auto que_alt = q.alt;
	transform(que_alt.begin(), que_alt.end(), que_alt.begin(), ::toupper);
	if (ref_ref == que_ref && ref_alt == que_alt) return true;
	return false;
}

void VCF::DirectSearchInThread(unordered_map<int, vector<SNP> > & ref_snps, unordered_map<int, vector<SNP> > & query_snps, int thread_index) {
	// handle heterozygous variants
    auto rit = ref_snps.begin();
	auto rend = ref_snps.end();
	for (; rit != rend;) {
		auto r_pos = rit->first;
		auto & r_snps = rit->second;
		auto qit = query_snps.find(r_pos);
		if (qit != query_snps.end()) {
			auto & q_snps = qit->second;

			for (auto r_snp_it = r_snps.begin(); r_snp_it != r_snps.end(); ) {
				bool matched_r_snp = false;
				for (auto q_snp_it = q_snps.begin(); q_snp_it != q_snps.end(); ) {
					if (CompareSnps(*r_snp_it, *q_snp_it)) {
                        auto temp_snp = *r_snp_it;
                        string matching_result = chromosome_name + '\t' + to_string(temp_snp.pos+1) + "\t" + temp_snp.ref + "\t" + temp_snp.alt;
                        direct_match_records[thread_index]->push_back(matching_result);
						matched_r_snp = true;
						q_snps.erase(q_snp_it);
						break;
					}
					else {
						++q_snp_it;
					}
				}
				if (matched_r_snp) {
					r_snps.erase(r_snp_it);
				}
				else {
					++r_snp_it;
				}
			}
			if (r_snps.size() == 0) {
				rit = ref_snps.erase(rit);
			}
			else {
				++rit;
			}
			if (q_snps.size() == 0) {
				query_snps.erase(qit);
			}
		}else{
            ++rit;
        }
	}
}

// directly match by position
void VCF::DirectSearchMultiThread() {

    direct_match_records = new vector<string>* [thread_num];
    for(int j = 0; j < thread_num; j++){
        direct_match_records[j] = new vector<string>;
    }

	vector<thread> threads;
	//spawn threads
	unsigned i = 0;
	for (; i < thread_num - 1; i++) {
		threads.push_back( thread(&VCF::DirectSearchInThread, this, ref(refpos_2_snp[i]), ref(querypos_2_snp[i]), i));
	}
	// also you need to do a job in main thread
	// i equals to (thread_num - 1)
	if (i != thread_num - 1) {
		dout << "[Error] thread number not match" << endl;
	}
	DirectSearchInThread(refpos_2_snp[i], querypos_2_snp[i],i);

	// call join() on each thread in turn before this function?
	std::for_each(threads.begin(), threads.end(), std::mem_fn(&std::thread::join));

    threads.clear();

    ofstream output_simple_file;
    output_simple_file.open(output_simple_filename);
    output_simple_file << "##VCF1:" << ref_vcf_filename << endl;
    output_simple_file << "##VCF2:" << que_vcf_filename << endl;
    output_simple_file << "#CHROM\tPOS\tREF\tALT" << endl;
    for(int i = 0; i < thread_num; i++){
        for (int j = 0; j < direct_match_records[i]->size(); j++){
            output_simple_file << direct_match_records[i]->at(j) << endl;
        }
    }
    output_simple_file.close();
    for(int j = 0; j < thread_num; j++){
        delete direct_match_records[j];
    }
    delete [] direct_match_records;

}

string VCF::ModifySequenceBySnp(const string sequence, SNP s, int offset) {
	string result = "";
	int snp_pos = s.pos - offset;
	int snp_end = snp_pos + (int)s.ref.length();
	if(snp_end > (int)sequence.length()){
        dout << "[Error] snp end greater than sequence length" << endl;
    }
	result += sequence.substr(0, snp_pos);
	result += s.alt;
	result += sequence.substr(snp_end, sequence.length() - snp_end);
	transform(result.begin(), result.end(), result.begin(), ::toupper);
	return result;
}
string VCF::ModifySequenceBySnpList(const string sequence, vector<SNP> s, int offset) {
	string result = sequence;
	int start_pos = 0;
    if(s.size() == 1){
        return ModifySequenceBySnp(sequence, s[0], offset);
    }
    sort(s.begin(), s.end());
	for (int i = s.size()-1; i >= 0; i--) {
		int snp_pos = s[i].pos - offset;
		int snp_end = snp_pos + (int)s[i].ref.length();
		string snp_alt = s[i].alt;
        int result_length = (int)result.length();
        if(snp_pos > result_length || snp_end > result_length){
            result = sequence;
            transform(result.begin(), result.end(), result.begin(), ::toupper);
            return result;
        }
        result = result.substr(0, snp_pos) + s[i].alt + result.substr(snp_end, result_length-snp_end);
    }
	transform(result.begin(), result.end(), result.begin(), ::toupper);
	return result;
}

bool VCF::CheckVariantOverlap(vector<SNP> snp_list){
    if (snp_list.size() <= 1) return false;
    int previous_ends = -1;
    for(int i = 0; i < snp_list.size(); i++){
        if(snp_list[i].pos < previous_ends) return true;
        if( previous_ends < snp_list[i].pos + (int)snp_list[i].ref.length()){
            previous_ends = snp_list[i].pos + (int)snp_list[i].ref.length();
        }
    }
    return false;
}

void f(){
    this_thread::sleep_for(chrono::seconds(2));
    cout << "Hello World" << endl;
}

bool VCF::CheckTandemRepeat(string sequence, int unit_threshold) {
	int sequence_length = (int)sequence.length();
    if(sequence_length == 1) return true;
	transform(sequence.begin(), sequence.end(), sequence.begin(), ::toupper);
    int end_index = sequence_length / 2 + 1;
	bool final_checking = false;
    int repeat_threshold = min(end_index-1, unit_threshold);
	for (int repeat_length = 1; repeat_length <= end_index; repeat_length++) {
		bool is_tandem_repeat = true;
        int repeat_time = 1;
		string repeat_region = sequence.substr(0, repeat_length);
		int start_position = repeat_length;
		while (start_position < sequence_length) {
			if (start_position + repeat_length > sequence_length)
				break;
			string matching_region = sequence.substr(start_position, repeat_length);
			if (matching_region != repeat_region) {
				is_tandem_repeat = false;
				break;
			}
			start_position += repeat_length;
            repeat_time ++;
		}
		if (is_tandem_repeat && repeat_time > 1) {
            final_checking = true;
			break;
		}
    }
	return final_checking;
}

/*
	clustering snps
	algorithm description, please refer to paper method
*/
void VCF::ClusteringSnps() {
    // handle heterozygous snps
	for (int i = 0; i < refpos_2_snp.size(); i++) {
		auto & m = refpos_2_snp[i];
		for (auto it = m.begin(); it != m.end(); ++it) {
			auto & v = it->second;
			for (int k = 0; k < v.size(); k++) {
				if (v[k].flag != 1) {
					v[k].flag = 1;
				}
				data_list.push_back(v[k]);
			}
		}
	}
	for (int i = 0; i < querypos_2_snp.size(); i++) {
		auto & m = querypos_2_snp[i];
		for (auto it = m.begin(); it != m.end(); ++it) {
			auto & v = it->second;
			for (int k = 0; k < v.size(); k++) {
				v[k].flag = -1;
				data_list.push_back(v[k]);
			}
		}
	}

	if (data_list.size() == 0)
		return;

	sort(data_list.begin(), data_list.end());

	int cluster_index = 0;
	int ins_ref = 0;
	int del_ref = 0;
	int ins_que = 0;
	int del_que = 0;
	int c_start = 0;
    int c_end = 0;

    for (int i = 0; i < data_list.size(); i++) {
        auto snp = data_list[i];
		// check if need to separator clusters
        if (i > 0) {
			c_end = snp.pos;
            if(c_end-c_start >= 2){
                string separator = genome_sequence.substr(c_start, c_end - c_start);
                int max_change = max(ins_ref + del_que, ins_que + del_ref);
                if ((int)(separator.length()) > 2 * max_change &&
                    ((int)(separator.length()) > MAX_REPEAT_LEN || !CheckTandemRepeat(separator, max_change)))
                {
                    cluster_index++;
                    ins_ref = 0;
                    del_ref = 0;
                    ins_que = 0;
                    del_que = 0;
                    c_start = 0; // re-assign c_start
                }
            }
		}
        if(c_start < snp.pos + (int)(snp.ref.length())) c_start = snp.pos + (int)(snp.ref.length());

        // assign snp to cluster
		cluster_snps_map[cluster_index].push_back(snp);
		int ref_length = (int)(snp.ref.length());
		int alt_length = (int)(snp.alt.length());
		int diff_length = abs(ref_length - alt_length);
		if (snp.flag == 1) {
			if (snp.snp_type == 'I') {
				ins_ref += diff_length;
			}
			else if (snp.snp_type == 'D') {
				del_ref += diff_length;
			}
		}
		else {
			if (snp.snp_type == 'I') {
				ins_que += diff_length;
			}
			else if (snp.snp_type == 'D') {
				del_que += diff_length;
			}
		}
	}
}
bool VCF::MatchSnpListsByPathConstruction(vector<SNP> & ref_snp_list,
	vector<SNP> & query_snp_list,
	vector<SNP> & mixed_list,
	const string subsequence,
	int offset,
	int thread_index)
{

}


bool VCF::MatchSnpLists(vector<SNP> & ref_snp_list,
                        vector<SNP> & query_snp_list,
                        vector<SNP> & mixed_list,
                        const string subsequence,
                        int offset,
                        int thread_index)
{
    // handle heterozygous snps
	map<string, vector<SNP> > ref_choice_snps;
	sort(mixed_list.begin(), mixed_list.end());

	for (int i = ref_snp_list.size(); i >= 1; i--) {
		vector<vector<SNP> > combinations = CreateCombinations(ref_snp_list, i);
        for (int k = 0; k < combinations.size(); k++) {
			auto c = combinations[k];
			if(CheckVariantOverlap(c)) continue;
            string ref_sequence = ModifySequenceBySnpList(subsequence, c, offset);
            ref_choice_snps[ref_sequence] = c;
		}
	}
	for (int i = query_snp_list.size(); i >= 1; i--) {
		vector<vector<SNP> > combinations = CreateCombinations(query_snp_list, i);
		for (int k = 0; k < combinations.size(); k++) {
			auto c = combinations[k];
			if(CheckVariantOverlap(c)) continue;
            string que_sequence = ModifySequenceBySnpList(subsequence, c, offset);
			if (ref_choice_snps.find(que_sequence) != ref_choice_snps.end()) {
				// delete all matched
                auto r = ref_choice_snps[que_sequence];
				sort(r.begin(), r.end());
                string matching_result = "";
                matching_result += chromosome_name;

                string parsimonious_ref = subsequence;
                string parsimonious_alt = que_sequence;
                if(parsimonious_ref == parsimonious_alt){
                    dout << "[Error] in variant, ref == alt";
                }
                int min_parsimonious_len = min(parsimonious_ref.size(), parsimonious_alt.size());
				// normalize
				int left_index = offset;
				if (toupper(genome_sequence[left_index]) != toupper(subsequence[0])) {
					dout << "[Error] genome sequence, subsequence, offset does not match." << endl;
				}
				bool change_in_allels = true;
				while (change_in_allels) {
					change_in_allels = false;
					if (toupper(parsimonious_ref.back()) == toupper(parsimonious_alt.back())) {
						if ((parsimonious_ref.size() > 1 && parsimonious_alt.back() > 1) || left_index > 0) {
							parsimonious_ref.pop_back();
							parsimonious_alt.pop_back();
							change_in_allels = true;
						}
						else {
							return false;
						}
					}
					if (parsimonious_ref.length() == 0 || parsimonious_alt.length() == 0) {
						left_index--;
						char left_char = genome_sequence[left_index];
						parsimonious_ref = left_char + parsimonious_ref;
						parsimonious_alt = left_char + parsimonious_alt;
					}
				}
				while (toupper(parsimonious_ref[0]) == toupper(parsimonious_alt[0]) && parsimonious_ref.size() > 1 && parsimonious_alt.size() > 1) {
					parsimonious_ref.erase(0, 1);
					parsimonious_alt.erase(0, 1);
				}

                matching_result += "\t" + parsimonious_ref + "\t" + parsimonious_alt;

                string ref_matching_variants = "";

				for (int m = 0; m < r.size(); m++) {
					SNP r_snp = r[m];
					for (auto n = mixed_list.begin(); n != mixed_list.end(); n++) {
						SNP m_snp = *n;
						if (m_snp.pos == r_snp.pos &&
                            m_snp.ref == r_snp.ref &&
                            m_snp.alt == r_snp.alt &&
                            m_snp.flag == r_snp.flag)
                        {
							mixed_list.erase(n);
                            break;
                        }
					}
                    for (auto n = ref_snp_list.begin(); n != ref_snp_list.end(); n++){
                        SNP m_snp = *n;
						if (m_snp.pos == r_snp.pos &&
                            m_snp.ref == r_snp.ref &&
                            m_snp.alt == r_snp.alt &&
                            m_snp.flag == r_snp.flag)
                        {
                        	// 1-based
                            ref_matching_variants += to_string(m_snp.pos+1) + "," + m_snp.ref + "," + m_snp.alt + ";";
							ref_snp_list.erase(n);
                            break;
                        }
                    }
				}
                matching_result += "\t" + ref_matching_variants;
                string que_matching_variants = "";
				sort(c.begin(), c.end());
				for (int m = 0; m < c.size(); m++) {
					SNP q_snp = c[m];
					for (auto n = mixed_list.begin(); n != mixed_list.end(); n++) {
						SNP m_snp = *n;
						if (m_snp.pos == q_snp.pos &&
                            m_snp.ref == q_snp.ref &&
                            m_snp.alt == q_snp.alt &&
                            m_snp.flag == q_snp.flag)
                        {
							mixed_list.erase(n);
                            break;
                        }
					}
                    for (auto n = query_snp_list.begin(); n != query_snp_list.end(); n++){
                        SNP m_snp = *n;
						if (m_snp.pos == q_snp.pos &&
                            m_snp.ref == q_snp.ref &&
                            m_snp.alt == q_snp.alt &&
                            m_snp.flag == q_snp.flag)
                        {
                        	// 1-based
                            que_matching_variants += to_string(m_snp.pos+1) + "," + m_snp.ref + "," + m_snp.alt + ";";
							query_snp_list.erase(n);
                            break;
                        }
                    }
				}
                matching_result += "\t" + que_matching_variants + "\n";
                complex_match_records[thread_index]->push_back(matching_result);
                return true;
			}
		}
	}
	return false;
}

void VCF::ClusteringSearchInThread(int start, int end, int thread_index) {
    for (int cluster_id = start; cluster_id < end; cluster_id++) {
		if (cluster_snps_map.find(cluster_id) != cluster_snps_map.end()) {
            auto & snp_list = cluster_snps_map[cluster_id];
            vector<SNP> candidate_ref_snps;
			vector<SNP> candidate_que_snps;
			vector<SNP> candidate_snps;
            int min_pos = std::numeric_limits<int>::max();
			int max_pos = 0;
			for (int i = 0; i < snp_list.size(); i++) {
				auto s = snp_list[i];
				if (s.flag == 1) {
					candidate_ref_snps.push_back(s);
				}
				else if(s.flag == -1) {
					candidate_que_snps.push_back(s);
				}
                candidate_snps.push_back(s);
				if (min_pos > s.pos) min_pos = s.pos;
				if (max_pos < s.pos + (int)(s.ref.length())) max_pos = s.pos + (int)(s.ref.length());
			}

			min_pos = max(0, min_pos - 1);
			max_pos = min(max_pos + 1, (int)genome_sequence.length());
            string subsequence = genome_sequence.substr(min_pos, max_pos-min_pos);

            if (candidate_ref_snps.size() == 0 || candidate_que_snps.size() == 0) continue;
			if (candidate_ref_snps.size() <= 1 && candidate_que_snps.size() <= 1) continue;
            if(candidate_ref_snps.size() > 10 || candidate_que_snps.size() > 10){
                vector<SNP> cluster_ref_snps;
                vector<SNP> cluster_que_snps;
                int ins_ref = 0;
                int del_ref = 0;
                int ins_que = 0;
                int del_que = 0;
                int c_start = std::numeric_limits<int>::max();
                int c_end = std::numeric_limits<int>::max();
                for(int i = 0; i < candidate_snps.size(); i++){
                    candidate_snps[i].pos += (int)candidate_snps[i].ref.length();
                }

                sort(candidate_snps.begin(), candidate_snps.end());

                for (int i = candidate_snps.size()-1; i >= 0; i--) {
                    auto snp = candidate_snps[i];
                    // check if need to separator clusters
                    if (i < candidate_snps.size() - 1) {
                        int c_start = snp.pos;
                        if(c_start < c_end){
                            string separator = genome_sequence.substr(c_start, c_end - c_start);
                            int max_change = max(ins_ref + del_que, ins_que + del_ref);
                            if ((int)separator.length() > 2 * max_change && !CheckTandemRepeat(separator, max_change))
                            {
                                while(cluster_ref_snps.size() > 0 &&
                                        cluster_que_snps.size() > 0 &&
                                        MatchSnpLists(cluster_ref_snps, cluster_que_snps, snp_list, subsequence, min_pos, thread_index));
                                cluster_ref_snps.clear();
                                cluster_que_snps.clear();
                                ins_ref = 0;
                                del_ref = 0;
                                ins_que = 0;
                                del_que = 0;
                            }
                        }
                    }

                    if(c_end > snp.pos- (int)snp.ref.length()) c_end = snp.pos - (int)snp.ref.length();
                    // assign snp to cluster
                    snp.pos -= (int)snp.ref.length();
                    if(snp.flag == 1){
                        cluster_ref_snps.push_back(snp);
                    }else{
                        cluster_que_snps.push_back(snp);
                    }
                    int ref_length = (int)snp.ref.length();
                    int alt_length = (int)snp.alt.length();
                    int diff_length = abs(ref_length - alt_length);
                    if (snp.flag == 1) {
                        if (snp.snp_type == 'I') {
                            ins_ref += diff_length;
                        }
                        else if (snp.snp_type == 'D') {
                            del_ref += diff_length;
                        }
                    }
                    else {
                        if (snp.snp_type == 'I') {
                            ins_que += diff_length;
                        }
                        else if (snp.snp_type == 'D') {
                            del_que += diff_length;
                        }
                    }
                }

                //if separating cluster does not work, try heuristic, if still not work, discard this cluster
                if(cluster_ref_snps.size() > 20 || cluster_que_snps.size() > 20){
                    // final check by variant length, if not applicable, skip it and give a warning.
                    if (cluster_ref_snps.size() > cluster_que_snps.size()){

                        int ref_sum_del_len = 0;
                        int ref_sum_ins_len = 0;
                        for(int j = 0; j < cluster_ref_snps.size(); j++){
                            int len_change = cluster_ref_snps[j].ref.size() -  cluster_ref_snps[j].alt.size();
                            if (len_change > 0){
                                ref_sum_del_len += len_change;
                            }else if(len_change < 0){
                                ref_sum_ins_len -= len_change;
                            }
                        }
                        bool skip_flag = false;
                        for(int j = 0; j < cluster_que_snps.size(); j++){
                            int len_change = cluster_que_snps[j].ref.size() - cluster_que_snps[j].alt.size();
                            if(len_change > 0){
                                if (ref_sum_del_len < len_change){
                                    skip_flag = true;
                                    break;
                                }
                            }else if(len_change < 0){
                                if (ref_sum_ins_len < len_change * -1){
                                    skip_flag = true;
                                    break;
                                }
                            }
                        }
                        if (skip_flag) continue;
                    }else{
                        int que_sum_del_len = 0;
                        int que_sum_ins_len = 0;
                        for(int j = 0; j < cluster_que_snps.size(); j++){
                            int len_change = cluster_que_snps[j].ref.size() -  cluster_que_snps[j].alt.size();
                            if (len_change > 0){
                                que_sum_del_len += len_change;
                            }else if(len_change < 0){
                                que_sum_ins_len -= len_change;
                            }
                        }
                        bool skip_flag = false;
                        for(int j = 0; j < cluster_ref_snps.size(); j++){
                            int len_change = cluster_ref_snps[j].ref.size() - cluster_ref_snps[j].alt.size();
                            if(len_change > 0){
                                if (que_sum_del_len < len_change){
                                    skip_flag = true;
                                    break;
                                }
                            }else if(len_change < 0){
                                if (que_sum_ins_len < len_change * -1){
                                    skip_flag = true;
                                    break;
                                }
                            }
                        }
                        if(skip_flag) continue;

                    }
                    cout << "[Warning] large cluster found, skip it." << endl;
                    continue;
                }

                while(cluster_ref_snps.size() > 0 &&
                        cluster_que_snps.size() > 0 &&
                        MatchSnpLists(cluster_ref_snps, cluster_que_snps, snp_list, subsequence, min_pos, thread_index));

            }
            else
            {
                while(candidate_ref_snps.size() > 0 &&
                        candidate_que_snps.size() > 0 &&
                        MatchSnpLists(candidate_ref_snps, candidate_que_snps, snp_list, subsequence, min_pos, thread_index));
            }
		}
		else {
			break;
		}
	}
}

// match by cluster
void VCF::ClusteringSearchMultiThread() {
	clustering_search = true;
	int start = cluster_snps_map.begin()->first;
	int cluster_number = cluster_snps_map.size();
	int cluster_end_boundary = start + cluster_number;
	int cluster_step = cluster_number / thread_num;
	if (cluster_step * thread_num < cluster_number) cluster_step++;
	int end = start + cluster_step;
    //initialize vector size, each allocating will have a lock
    complex_match_records = new vector<string>* [thread_num];
    for(int j = 0; j < thread_num; j++){
        complex_match_records[j] = new vector<string>;
    }

	vector<thread> threads;
	//spawn threads
	unsigned i = 0;
	for (; i < thread_num - 1; i++) {
        int variant_number = 0;
        for (int cluster_id = start; cluster_id < end; cluster_id++) {
		    if (cluster_snps_map.find(cluster_id) != cluster_snps_map.end()) {
                variant_number += cluster_snps_map[cluster_id].size();
            }
        }
		threads.push_back(thread(&VCF::ClusteringSearchInThread, this, start, end, i));
		start = end;
		end = start + cluster_step;
	}
	// also you need to do a job in main thread
	// i equals to (thread_num - 1)
	if (i != thread_num - 1) {
		dout << "[Error] thread number not match" << endl;
	}
	if (start >= cluster_snps_map.size()) {
		dout << "[Error] index out of map range" << endl;
	}
	else {
        int variant_number = 0;
        for (int cluster_id = start; cluster_id < end; cluster_id++) {
		    if (cluster_snps_map.find(cluster_id) != cluster_snps_map.end()) {
                variant_number += cluster_snps_map[cluster_id].size();
            }
        }
		ClusteringSearchInThread(start, end, i);
	}

	// call join() on each thread in turn before this function?
	std::for_each(threads.begin(), threads.end(), std::mem_fn(&std::thread::join));

    ofstream output_complex_file;
    output_complex_file.open(output_complex_filename);
    output_complex_file << "##VCF1:" << ref_vcf_filename << endl;
    output_complex_file << "##VCF2:" << que_vcf_filename << endl;
    output_complex_file << "#CHROM\tPOS\tREF\tALT\tVCF1\tVCF2" << endl;
    for(int i = 0; i < thread_num; i++){
        for (int j = 0; j < complex_match_records[i]->size(); j++){
            if(complex_match_records[i]->at(j).find_first_not_of(' ') != std::string::npos){
                output_complex_file << complex_match_records[i]->at(j);
            }
        }
    }
    output_complex_file.close();

    for(int j = 0; j < thread_num; j++){
        delete complex_match_records[j];
    }
    delete [] complex_match_records;
}

int VCF::GetRefSnpNumber(int & indel_num) {
	int result = 0;
	indel_num = 0;
	if (clustering_search) {
		for (auto it = cluster_snps_map.begin(); it != cluster_snps_map.end(); it++) {
			auto v = it->second;
			for (int i = 0; i < v.size(); i++) {
				if (v[i].flag == 1) {
					result++;
					if (v[i].ref.length() != v[i].alt.length())
						indel_num++;
				}
			}
		}
	}else{
	    for (int i = 0; i < refpos_2_snp.size(); i++) {
			for (auto it = refpos_2_snp[i].begin(); it != refpos_2_snp[i].end(); ++it) {
				auto v = it->second;
				result += v.size();
				for (int k = 0; k < v.size(); k++) {
					if (v[k].ref.length() != v[k].alt.length())
						indel_num++;
				}
			}
	    }
    }
	return result;
}

int VCF::GetQuerySnpNumber(int & indel_num) {
	int result = 0;
	indel_num = 0;
	if (clustering_search) {
		for (auto it = cluster_snps_map.begin(); it != cluster_snps_map.end(); it++) {
			auto v = it->second;
			for (int i = 0; i < v.size(); i++) {
				if (v[i].flag == -1) {
					result++;
					if (v[i].ref.length() != v[i].alt.length())
						indel_num++;
				}
			}
		}
	}else{
	    for (int i = 0; i < querypos_2_snp.size(); i++) {
			for (auto it = querypos_2_snp[i].begin(); it != querypos_2_snp[i].end(); ++it) {
				auto v = it->second;
				result += v.size();
				for (int k = 0; k < v.size(); k++) {
					if (v[k].ref.length() != v[k].alt.length())
						indel_num++;
				}
			}
	    }
    }
	return result;
}

void VCF::Compare(string ref_vcf,
        string query_vcf,
        string genome_seq,
        bool direct_search,
        string output_prefix,
        bool match_genotype,
		bool normalization){

    ref_vcf_filename = ref_vcf;
    que_vcf_filename = query_vcf;
    this->match_genotype = match_genotype;
	this->normalization = normalization;
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
    int ref_total_indel_num, que_total_indel_num;
	int ref_total_num = GetRefSnpNumber(ref_total_indel_num);
    int que_total_num = GetQuerySnpNumber(que_total_indel_num);
    dout << " referece vcf entry number [total, indel]: " << ref_total_num << "," << ref_total_indel_num << endl;
	dout << " query vcf entry number: [total, indel] " << que_total_num << "," << que_total_indel_num << endl;


	//------------direct search
	dsptime();
	dout << " Direct search ... " << endl;
	DirectSearchMultiThread();
	dsptime();
	dout << " Finish direct search." << endl;
    int ref_direct_left_indel_num, que_direct_left_indel_num;
    int ref_direct_left_num = GetRefSnpNumber(ref_direct_left_indel_num);
    int que_direct_left_num = GetQuerySnpNumber(que_direct_left_indel_num);
    int ref_direct_match_num = ref_total_num - ref_direct_left_num;
    int que_direct_match_num = que_total_num - que_direct_left_num;
    int ref_direct_match_indel_num = ref_total_indel_num - ref_direct_left_indel_num;
    int que_direct_match_indel_num = que_total_indel_num - que_direct_left_indel_num;
	dout << " referece vcf entry direct match number [total, indel]: " << ref_direct_match_num << "," << ref_direct_match_indel_num << endl;
	dout << " query vcf entry direct match number [total, indel]: " << que_direct_match_num  << "," << que_direct_match_indel_num << endl;

	if (direct_search){
	    dout << " referece vcf entry mismatch number [total, indel]: " << ref_direct_left_num << "," << ref_direct_left_indel_num << endl;
	    dout << " query vcf entry mismatch number [total, indel]: " << que_direct_left_num << "," << que_direct_left_indel_num << endl;
        ofstream output_stat_file;
        output_stat_file.open(output_stat_filename);
        output_stat_file << ref_total_num << endl;
        output_stat_file << que_total_num << endl;
        output_stat_file << ref_direct_match_num << endl;
        output_stat_file << que_direct_match_num << endl;
        output_stat_file << ref_direct_left_num << endl;
        output_stat_file << que_direct_left_num << endl;
        //=====================================================
        output_stat_file << ref_total_indel_num << endl;
		output_stat_file << que_total_indel_num << endl;
		output_stat_file << ref_direct_match_indel_num << endl;
		output_stat_file << que_direct_match_indel_num << endl;
		output_stat_file << ref_direct_left_indel_num << endl;
		output_stat_file << que_direct_left_indel_num << endl;
        output_stat_file.close();

        return;
    }

	//-------------clustering search
	dsptime();
	dout << " Clustering snps ... " << endl;
	ClusteringSnps();
	dsptime();
	dout << " Finish clustering." << endl;
	dsptime();
	dout << " Clustering search ... " << endl;
	ClusteringSearchMultiThread();
	dsptime();
	dout << " Finish clustering search." << endl;
    int ref_cluster_left_indel_num, que_cluster_left_indel_num;
	int ref_cluster_left_num = GetRefSnpNumber(ref_cluster_left_indel_num);
    int que_cluster_left_num = GetQuerySnpNumber(que_cluster_left_indel_num);
    int ref_cluster_match_num = ref_direct_left_num - ref_cluster_left_num;
    int que_cluster_match_num = que_direct_left_num - que_cluster_left_num;
    int ref_cluster_match_indel_num = ref_direct_left_indel_num - ref_cluster_left_indel_num;
    int que_cluster_match_indel_num = que_direct_left_indel_num - que_cluster_left_indel_num;

    dout << " referece vcf entry cluster match number [total, indel]: " << ref_cluster_match_num << "," << ref_cluster_match_indel_num << endl;
	dout << " query vcf entry cluster match number [total, indel]: " << que_cluster_match_num << "," << que_cluster_match_indel_num << endl;

	dout << " referece vcf entry mismatch number [total, indel]: " << ref_cluster_left_num << "," << ref_cluster_left_indel_num << endl;
	dout << " query vcf entry mismatch number [total, indel]: " << que_cluster_left_num  << "," << que_cluster_left_indel_num << endl;

    //write stat file
    ofstream output_stat_file;
    output_stat_file.open(output_stat_filename);
    output_stat_file << ref_total_num << endl;
    output_stat_file << que_total_num << endl;
    output_stat_file << ref_direct_match_num << endl;
    output_stat_file << que_direct_match_num << endl;
    output_stat_file << ref_cluster_match_num << endl;
    output_stat_file << que_cluster_match_num << endl;
    output_stat_file << ref_cluster_left_num << endl;
    output_stat_file << que_cluster_left_num << endl;
	//=====================================================
	output_stat_file << ref_total_indel_num << endl;
	output_stat_file << que_total_indel_num << endl;
	output_stat_file << ref_direct_match_indel_num << endl;
	output_stat_file << que_direct_match_indel_num << endl;
	output_stat_file << ref_cluster_match_indel_num << endl;
	output_stat_file << que_cluster_match_indel_num << endl;
	output_stat_file << ref_direct_left_indel_num << endl;
	output_stat_file << que_direct_left_indel_num << endl;
    output_stat_file.close();

    return;
}
