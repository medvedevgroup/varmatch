#include "vcf.h"


bool operator <(const SNP& x, const SNP& y) {
	return x.pos < y.pos;
}

bool operator ==(const SNP& x, const SNP& y) {
	if (x.pos == y.pos && x.snp_type == y.snp_type && x.alt == y.alt) {
		return true;
	}
	return false;
}

VCF::VCF(int thread_num_)
{
	genome_sequence = "";
	boundries_decided = false;
    complex_search = false;
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
	if (!boundries_decided) {
		cout << "[Error] VCF::ReadVCF cannot read vcf file before read genome file" << endl;
		return;
	}

	ifstream vcf_file;
	vcf_file.open(filename.c_str());
	if (!vcf_file.good()) {
		cout << "[Error] VCF::ReadVCF can not open vcf file" << endl;
		return;
	}

	while (!vcf_file.eof()) { // alternative way is vcf_file != NULL
		string line;
		getline(vcf_file, line, '\n');
		//dout << line << endl;
		if (line.length() <= 1) continue;
		if (line[0] == '#') continue;
		auto columns = split(line, '\t');
		auto pos = atoi(columns[1].c_str()) - 1;
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

		//decide which thread to use
		int index = 0;
		for (int i = 0; i < pos_boundries.size(); i++) {
			if (pos < pos_boundries[i]) {
				index = i;
				break;
			}
		}

		pos_2_snp[index][pos].push_back(SNP(pos, snp_type, ref, alt));
		pos_2_vcf_entry[pos] = line;
	}
	vcf_file.close();
	return;
}

void VCF::ReadGenomeSequence(string filename) {
	ifstream genome_file;
	genome_file.open(filename.c_str());
	if (!genome_file.good()) {
		cout << "[Error] VCF::ReadGenomeSequence can not open fasta file" << endl;
		return;
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
	ReadVCF(filename, refpos_2_snp, refpos_2_vcf_entry);
}

void VCF::ReadQueryVCF(string filename) {
	ReadVCF(filename, querypos_2_snp, querypos_2_vcf_entry);
}

bool VCF::CompareSnps(SNP r, SNP q) {
	if (r.snp_type == q.snp_type && r.alt == q.alt) return true;
	return false;
}

void VCF::DirectSearchInThread(unordered_map<int, vector<SNP> > & ref_snps, unordered_map<int, vector<SNP> > & query_snps) {
	auto rit = ref_snps.begin();
	auto rend = ref_snps.end();
	for (; rit != rend;) {
		auto r_pos = rit->first;
		auto & r_snps = rit->second;
		auto qit = query_snps.find(r_pos);
		if (qit != query_snps.end()) {
			
			auto & q_snps = qit->second;

			if (r_snps.size() != 1 || q_snps.size() != 1) {
				cout << "[Error] snp vector size not right" << endl;
			}

			vector<vector<SNP>::iterator> r_deleted_snps;
			vector<vector<SNP>::iterator> q_deleted_snps;
			for (auto r_snp_it = r_snps.begin(); r_snp_it != r_snps.end(); ++r_snp_it) {
				for (auto q_snp_it = q_snps.begin(); q_snp_it != q_snps.end(); ++q_snp_it) {
					if (CompareSnps(*r_snp_it, *q_snp_it)) {
						r_deleted_snps.push_back(r_snp_it);
						q_deleted_snps.push_back(q_snp_it);
					}
				}
			}
			for (int i = 0; i < r_deleted_snps.size(); i++) {
				r_snps.erase(r_deleted_snps[i]);
			}
			if (r_snps.size() == 0) {
				rit = ref_snps.erase(rit);
			}
			else {
				++rit;
			}
			for (int i = 0; i < q_deleted_snps.size(); i++) {
				q_snps.erase(q_deleted_snps[i]);
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
	vector<thread> threads;
	//spawn threads
	unsigned i = 0;
	for (; i < thread_num - 1; i++) {
		threads.push_back( thread(&VCF::DirectSearchInThread, this, ref(refpos_2_snp[i]), ref(querypos_2_snp[i])) );
	}
	// also you need to do a job in main thread
	// i equals to (thread_num - 1)
	if (i != thread_num - 1) {
		dout << "[Error] thread number not match" << endl;
	}
	DirectSearchInThread(refpos_2_snp[i], querypos_2_snp[i]);

	// call join() on each thread in turn before this function?
	std::for_each(threads.begin(), threads.end(), std::mem_fn(&std::thread::join));

    threads.clear();
}

string VCF::ModifySequenceBySnp(string sequence, SNP s, int offset) {
	// [todo] unit test
	string result = "";
	int snp_pos = s.pos - offset;
	int snp_end = snp_pos + s.ref.length();
	if(snp_end > sequence.length()){
        dout << "[Error] snp end greater than sequence length" << endl;
    }
    //if(snp_end > sequence.length()){
    //    cout << "snp end greater than sequence length" << endl;
    //    cout << snp_end << "\t" << sequence.length() << endl;
    //}
	result += sequence.substr(0, snp_pos);
	result += s.alt;
	result += sequence.substr(snp_end, sequence.length() - snp_end);
	transform(result.begin(), result.end(), result.begin(), ::toupper);
	return result;
}
string VCF::ModifySequenceBySnpList(string sequence, vector<SNP> s, int offset) {
	// [todo] unit test
	string result = "";
	int start_pos = 0;
	for (int i = 0; i < s.size(); i++) {
		int snp_pos = s[i].pos - offset;
		int snp_end = snp_pos + s[i].ref.length();
		string snp_alt = s[i].alt;
		result += sequence.substr(start_pos, snp_pos - start_pos);
		result += snp_alt;
		start_pos = snp_end;
	}
	if (start_pos < sequence.length()) {
		result += sequence.substr(start_pos, sequence.length() - start_pos);
	}
	transform(result.begin(), result.end(), result.begin(), ::toupper);
	return result;
}

bool VCF::ComplexMatch(SNP s, vector<SNP> comb) {
	//size of comb >= 1
    assert(comb.size() >= 1);
	sort(comb.begin(), comb.end());
	int ref_left = s.pos;
	int ref_right = ref_left + s.ref.length();

	int comb_size = comb.size();
	int que_left = comb[0].pos;
	int que_right = comb[comb_size - 1].pos + comb[comb_size - 1].ref.length();

	int genome_left = min(ref_left, que_left);
	int genome_right = max(ref_right, que_right);

	genome_left = max(0, genome_left - 10);
	genome_right = min(genome_right + 10, (int)genome_sequence.length());

	string subsequence = genome_sequence.substr(genome_left, genome_right - genome_left);
	return ModifySequenceBySnp(subsequence, s, genome_left) == ModifySequenceBySnpList(subsequence, comb, genome_left);
}

bool VCF::ExponentialComplexMatch(SNP r_snp, map<int, vector<SNP> > & query_snps, vector<SNP> & deleted_ref_snps, vector<SNP> & deleted_que_snps) {

	int ref_start_pos = r_snp.pos;
	auto ref_ref = r_snp.ref;
	auto ref_alt = r_snp.alt;
    int ref_change = ref_alt.length() - ref_ref.length();
	int ref_end_pos = ref_start_pos + ref_ref.size();
	vector<SNP> candidate_query_list;
    vector<int> candidate_changes;
	FindVariantsInRange(ref_start_pos, ref_end_pos, query_snps, candidate_query_list, candidate_changes);

	int candidate_size = candidate_query_list.size();
	if (candidate_size == 0) return false;
	//check all combinations, from largest, one single match is enough
	bool flag = false;
	for (int k = candidate_query_list.size(); k >= 1; --k) {
		vector<vector<SNP>> combinations = CreateCombinations(candidate_query_list, k, candidate_changes, ref_change);
		//vector<vector<SNP>> combinations = CreateCombinations(candidate_query_list, k);
		bool matched = false;
		// check combinations with k elements
		for (auto cit = combinations.begin(); cit != combinations.end(); ++cit) {
			auto comb = *cit;
			if (ComplexMatch(r_snp, *cit)) {
				matched = true;

				// delete corresponding snps
				deleted_ref_snps.push_back(r_snp);
				for (auto combit = comb.begin(); combit != comb.end(); ++combit) {
					deleted_que_snps.push_back(*combit);
				}
				break;
			}
		}
		if (matched) {
			flag = true;
			break;
		}
	}
	return flag;
}

unsigned int VCF::EditDistance(const std::string& s1, const std::string& s2)
{
	const std::size_t len1 = s1.size(), len2 = s2.size();
	std::vector<unsigned int> col(len2 + 1), prevCol(len2 + 1);

	for (unsigned int i = 0; i < prevCol.size(); i++)
		prevCol[i] = i;
	for (unsigned int i = 0; i < len1; i++) {
		col[0] = i + 1;
		for (unsigned int j = 0; j < len2; j++)
			// note that std::min({arg1, arg2, arg3}) works only in C++11,
			// for C++98 use std::min(std::min(arg1, arg2), arg3)
			col[j + 1] = std::min({ prevCol[1 + j] + 1, col[j] + 1, prevCol[j] + (s1[i] == s2[j] ? 0 : 1) });
		col.swap(prevCol);
	}
	return prevCol[len2];
}

bool VCF::GreedyComplexMatch(SNP r_snp, map<int, vector<SNP> > & query_snps, vector<SNP> & deleted_ref_snps, vector<SNP> & deleted_que_snps) {
    dout << "==" << endl;
	int ref_start_pos = r_snp.pos;
	auto ref_ref = r_snp.ref;
	auto ref_alt = r_snp.alt;
	int ref_end_pos = ref_start_pos + ref_ref.size();

	auto itlow = query_snps.lower_bound(ref_start_pos);
	auto itup = query_snps.upper_bound(ref_end_pos);
	// theoretically, since all snps are independent, we do not need the following two if-statement
	// but this is engineering
	if (itlow != query_snps.begin()) {
		itlow--;
	}
	if (itup != query_snps.end()) {
		itup++;
	}

	// we have candidate query snps stored separatedly in two data structures
	// one according to positions
	map<int, vector<SNP> > candidate_query_map;
	// one purely stored
	vector<SNP> comb;

    int que_left = std::numeric_limits<int>::max();
    int que_right = 0;
	for (auto it = itlow; it != itup; ++it) {
		auto v = it->second;
        bool flag = false;
		for (int i = 0; i < v.size(); i++) {
			int snp_start = v[i].pos;
			int snp_end = snp_start + v[i].ref.length();
			//if(ref_start_pos <= snp_start && ref_end_pos >= snp_start){
            if (min(ref_end_pos, snp_end) - max(ref_start_pos, snp_start) > 0) {
				comb.push_back(v[i]);
                if (que_left > snp_start ) que_left = snp_start;
                if (que_right < snp_end) que_right = snp_end;
                flag = true;
                dout << "candidates:" << ref_start_pos << "," << ref_end_pos << ": ";
                dout << snp_start << "," << snp_end << endl;
                dout << v[i].pos << "\t" << v[i].ref << "\t" << v[i].alt << endl;
			}
		}
        if(flag){
		    candidate_query_map[it->first] = v;
        }
	}

	if (comb.size() < 1) return false;
	sort(comb.begin(), comb.end());

	int ref_left = r_snp.pos;
	int ref_right = ref_left + r_snp.ref.length();

	int comb_size = comb.size();
	//int que_left = comb[0].pos;
	//int que_right = comb[comb_size - 1].pos + comb[comb_size - 1].ref.length();

	int genome_left = min(ref_left, que_left);
	int genome_right = max(ref_right, que_right);

	genome_left = max(0, genome_left - 1);
	genome_right = min(genome_right + 1, (int)genome_sequence.length());

	string subsequence = genome_sequence.substr(genome_left, genome_right - genome_left);
	string ref_subseq = ModifySequenceBySnp(subsequence, r_snp, genome_left);
	string que_subseq = subsequence;
	int edit_distance = EditDistance(ref_subseq, que_subseq);
    int len_distance = abs(ref_subseq.length() - que_subseq.length());
	int que_offset = genome_left;

	vector<SNP> deleted_snps;
    dout << "genome sequence:" << que_subseq << endl;
    dout << r_snp.pos << "\t" << r_snp.ref << "\t" << r_snp.alt << endl;
    dout << "ref:" << ref_subseq << endl;
	for (auto it = candidate_query_map.begin(); it != candidate_query_map.end(); ++it) {
		auto v = it->second;
		SNP min_s;
		int min_distance = std::numeric_limits<int>::max();
        int min_len_distance = min_distance;
		string min_subseq;
		for (int i = 0; i < v.size(); i++) {
			auto s = v[i];
            dout << que_subseq << "\t" << s.pos << "\t" << que_offset << endl;
            dout << s.pos << "\t" << s.ref << "\t" << s.alt << endl;
			auto subseq = ModifySequenceBySnp(que_subseq, s, que_offset);
			int distance = EditDistance(ref_subseq, subseq);
            int len_dis = abs(ref_subseq.size()- subseq.size());
			if (distance < min_distance && len_dis <= min_len_distance) {
				min_distance = distance;
				min_s = s;
				min_subseq = subseq;
                min_len_distance = len_dis;
			}
		}
		if (min_distance < edit_distance && min_len_distance <= len_distance) {
			//dout << "distance: " << edit_distance << "," << min_distance << endl;
            edit_distance = min_distance;
			que_subseq = min_subseq;
			que_offset += (min_s.ref.length() - min_s.alt.length());
			deleted_snps.push_back(min_s);
		}
	}

	if (edit_distance == 0) {
		deleted_ref_snps.push_back(r_snp);
		for (int i = 0; i < deleted_snps.size(); i++)
			deleted_que_snps.push_back(deleted_snps[i]);

		return true;
	}
	else {
		return false;
	}

	return true;
}

void VCF::FindVariantsInRange(int start, int end, map<int, vector<SNP> > snp_map, vector<SNP> & candidate_query_list, vector<int>& candidate_changes) {
	auto itlow = snp_map.lower_bound(start);
	auto itup = snp_map.upper_bound(end);
	// theoretically, since all snps are independent, we do not need the following two if-statement
	// but this is engineering
	if (itlow != snp_map.begin()) {
		itlow--;
	}
	if (itup != snp_map.end()) {
		itup++;
	}

	for (auto it = itlow; it != itup; ++it) {
		auto v = it->second;
		for (int i = 0; i < v.size(); i++) {
			int snp_start = v[i].pos;
			int snp_end = snp_start + v[i].ref.length();
            int change = v[i].alt.length() - v[i].ref.length();
			//if(end > snp_start && start >= snp_start){
            if (min(end, snp_end) - max(start, snp_start) > 0) {
				candidate_query_list.push_back(v[i]);
                candidate_changes.push_back(change);
			}
		}
	}
}

void VCF::ComplexSearchInThread(map<int, vector<SNP> > & ref_snps, map<int, vector<SNP> > & query_snps) {
	// linear algorithm
	vector<SNP> deleted_ref_snps;
	vector<SNP> deleted_que_snps;
	// for each position in ref, i.e. a vector
	//dout << ref_snps.size() << endl;
    //dout << query_snps.size() << endl;
    //int i = 0;
    for (auto rit = ref_snps.begin(); rit != ref_snps.end(); ++rit) {
        //i ++;
        //dout << i << endl;
		int ref_start_pos = rit->first;
		auto ref_snp_list = rit->second;
		auto qit_start = query_snps.begin();

		// for each snp in the vector
		for (int i = 0; i < ref_snp_list.size(); i++) {
			auto ref_snp = ref_snp_list[i];

			ExponentialComplexMatch(ref_snp_list[i], query_snps, deleted_ref_snps, deleted_que_snps);

			//GreedyComplexMatch(ref_snp_list[i], query_snps, deleted_ref_snps, deleted_que_snps);
		}
	}

	// delete all snps, first check position, then delete by value matching
	// [todo] this actually should be one separated function

	for (auto it = deleted_ref_snps.begin(); it != deleted_ref_snps.end(); ++it) {
		auto snp = *it;
		auto pos = snp.pos;
		auto & v = ref_snps[pos];
		auto vit = v.begin();
		while (vit != v.end()) {
			if (snp == *vit) {
				vit = v.erase(vit);
				break;
			}
			else {
				++vit;
			}
		}
		if (v.size() == 0) {
			ref_snps.erase(ref_snps.find(pos));
		}
	}

	for (auto it = deleted_que_snps.begin(); it != deleted_que_snps.end(); ++it) {
		auto snp = *it;
		auto pos = snp.pos;
		auto & v = query_snps[pos];
		auto vit = v.begin();
		while (vit != v.end()) {
			if (snp == *vit) {
				vit = v.erase(vit);
                break;
			}
			else {
				++vit;
			}
		}
		if (v.size() == 0) {
			query_snps.erase(query_snps.find(pos));
		}
	}
}

void f(){
    this_thread::sleep_for(chrono::seconds(2));
    cout << "Hello World" << endl;
}
// match by overlapping reference region
void VCF::ComplexSearchMultiThread() {
	if (GetRefSnpNumber() == 0 || GetQuerySnpNumber() == 0) return;
    complex_search = true;

	// transfer data from hash to map
	for (int i = 0; i < refpos_2_snp.size(); i++) {
		auto & pos_snp_hash = refpos_2_snp[i];
		for (auto rit = pos_snp_hash.begin(); rit != pos_snp_hash.end(); ++rit) {
			//refpos_snp_map[i][rit->first] = rit->second;
            auto p = rit->first;
            auto v = rit->second;
            for(int j  = 0; j < v.size(); j++){
                refpos_snp_map[i][p].push_back(v[j]);
            }
		}
	}
	refpos_2_snp.clear();

	for (int i = 0; i < querypos_2_snp.size(); i++) {
		auto & pos_snp_hash = querypos_2_snp[i];
		for (auto qit = pos_snp_hash.begin(); qit != pos_snp_hash.end(); ++qit) {
			//querypos_snp_map[i][qit->first] = qit->second;
            auto p = qit->first;
            auto v = qit->second;
            for(int j = 0; j < v.size(); j++){
                querypos_snp_map[i][p].push_back(v[j]);
            }
		}
	}
	querypos_2_snp.clear();


	vector<thread> threads;
	//spawn threads
	unsigned i = 0;
	for (; i < thread_num - 1; i++) {
        if(i >= refpos_snp_map.size() || i >= querypos_snp_map.size()){
            dout << "[Error] index out of map range" << endl;
            continue;
        }
        if(refpos_snp_map[i].size() == 0 || querypos_snp_map[i].size() == 0){
            continue;
        }
		//threads.push_back(thread(f));
        //dout << "create new thread" << endl;
        threads.push_back(thread(&VCF::ComplexSearchInThread, this, ref(refpos_snp_map[i]), ref(querypos_snp_map[i])));
	}
	// also you need to do a job in main thread
	// i equals to (thread_num - 1)
	if (i != thread_num - 1) {
		dout << "[Error] thread number not match" << endl;
	}
    if(i >= refpos_snp_map.size() || i >= querypos_snp_map.size()){
        dout << "[Error] index out of map range" << endl;
    }else if(refpos_snp_map[i].size() != 0 && querypos_snp_map[i].size() != 0){
	    ComplexSearchInThread(refpos_snp_map[i], querypos_snp_map[i]);
    }

	// call join() on each thread in turn before this function?
	std::for_each(threads.begin(), threads.end(), std::mem_fn(&std::thread::join));

    //[todo] check boundary
    
    //----------------------------change direction-------------------------------
    threads.clear();
    i = 0;
	for (; i < thread_num - 1; i++) {
        if(i >= refpos_snp_map.size() || i >= querypos_snp_map.size()){
            dout << "[Error] index out of map range" << endl;
            continue;
        }
        if(refpos_snp_map[i].size() == 0 || querypos_snp_map[i].size() == 0){
            continue;
        }
		//threads.push_back(thread(f));
        //dout << "create new thread" << endl;
        threads.push_back(thread(&VCF::ComplexSearchInThread, this, ref(querypos_snp_map[i]), ref(refpos_snp_map[i])));
	}
	// also you need to do a job in main thread
	// i equals to (thread_num - 1)
	if (i != thread_num - 1) {
		dout << "[Error] thread number not match" << endl;
	}
    if(i >= refpos_snp_map.size() || i >= querypos_snp_map.size()){
        dout << "[Error] index out of map range" << endl;
    }else if(refpos_snp_map[i].size() != 0 && querypos_snp_map[i].size() != 0){
	    ComplexSearchInThread(querypos_snp_map[i], refpos_snp_map[i]);
    }

	// call join() on each thread in turn before this function?
	std::for_each(threads.begin(), threads.end(), std::mem_fn(&std::thread::join));
    //[todo] check boundries
    //-----------------------------------------------------------------------------

}

// clustering snps
void VCF::ClusteringSnps() {}

// match by cluster
void VCF::ClusteringSearch() {}

int VCF::GetRefSnpNumber() {
	int result = 0;
    if(complex_search){
        for (int i = 0; i < refpos_snp_map.size(); i++){
            result += refpos_snp_map[i].size();
        }
    }else{
	    for (int i = 0; i < refpos_2_snp.size(); i++) {
		    result += refpos_2_snp[i].size();
	    }
    }
	return result;
}

int VCF::GetQuerySnpNumber() {
	int result = 0;
    if(complex_search){
        for(int i = 0; i < querypos_snp_map.size(); i++){
            result += querypos_snp_map[i].size();
        }
    }else{
	    for (int i = 0; i < querypos_2_snp.size(); i++) {
		    result += querypos_2_snp[i].size();
	    }   
    }
	return result;
}
