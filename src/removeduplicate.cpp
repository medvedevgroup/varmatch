#include "removeduplicate.h"

RemoveDuplicate::RemoveDuplicate(int thread_num_):VCF(thread_num_){}
RemoveDuplicate::~RemoveDuplicate(){}

int RemoveDuplicate::GetThreadIndex(int pos){
    for(int i = 0; i < pos_boundries.size(); i++){
        if(pos < pos_boundries[i]){
            return i;
        }
    }
}

int RemoveDuplicate::ReadVCFWithoutDup(string filename){
    if(!boundries_decided){
        cout << "[Error: RemoveDuplicate] ReadVCFWithoutDup can not read vcf file before read genome file" << endl;
        return -1;
    }

    ifstream vcf_file;
    vcf_file.open(filename.c_str());
    if (!vcf_file.good()) {
        cout << "[Error] RemoveDuplicate::ReadVCFWithoutDup can not open vcf file" << endl;
        return -1;
    }
    int var_num = 0;
    int nodup_var_num = 0;
    while(!vcf_file.eof()){
        string line;
        getline(vcf_file, line, '\n');
        if (line.length() <= 1) continue;
        if (line[0] == '#') continue;
        auto columns = split(line, '\t');
        if(chromosome_name == ".") chromosome_name = columns[0];
        auto pos = atoi(columns[1].c_str()) - 1;
        string ref = columns[3];
        string alt = columns[4];
        string quality = columns[6];

        vector<string> alt_list;
        if(alt.find(",") != string::npos){
            continue;
            // deal with multi alt
            alt_list = split(alt, ',');
        }else{
            alt_list.push_back(alt);
        }

        //int thread_index = GetThreadIndex(pos);

        char snp_type;
        for (auto it = alt_list.begin(); it != alt_list.end(); ++it){
            snp_type = 'S';
            string a = *it;
            if(ref.length() > alt.length()){
                snp_type = 'D';
            }else if(ref.length() < alt.length()){
                snp_type = 'I';
            }
            var_num ++;
            string varid = to_string(pos) + "_" + ref + "_" + a;
            transform(varid.begin(), varid.end(), varid.begin(), ::toupper);
            if(nondup_vcfentry_hash.find(varid) != nondup_vcfentry_hash.end()){
                nodup_var_num ++;
                nondup_vcfentry_hash[varid] = line;
                nondup_pos_snp_map[pos].push_back(SNP(pos, snp_type, ref, a));
            }
        }
    }

    vcf_file.close();
    return var_num;
}

void RemoveDuplicate::ClusteringSnps() {
    //int num = 0;
    for (auto it = nondup_pos_snp_map.begin(); it != nondup_pos_snp_map.end(); ++it) {
        auto & v = it->second;
        for (int k = 0; k < v.size(); k++) {
            data_list.push_back(v[k]);
        }
    }
    if (data_list.size() == 0)
        return;
    sort(data_list.begin(), data_list.end());

    int cluster_index = 0;
    int ins_total = 0;
    int del_total = 0;
    int c_start = 0;
    int c_end = 0;
    
    for (int i = 0; i < data_list.size(); i++) {
        auto snp = data_list[i];
        // check if need to separator clusters
        if (i > 0) {
            c_end = snp.pos;
            if(c_end-c_start >= 2){
                string separator = genome_sequence.substr(c_start, c_end - c_start);
                int max_change = max(ins_total, del_total);
                if (separator.length() > 2 * max_change &&
                    (separator.length() > MAX_REPEAT_LEN || !CheckTandemRepeat(separator, max_change))) 
                {
                    cluster_index++;
                    ins_total = 0;
                    del_total = 0;
                    c_start = 0; // re-assign c_start
                }
            }
        }

        if(c_start < snp.pos + snp.ref.length()) c_start = snp.pos + snp.ref.length();
        // assign snp to cluster
        cluster_snps_map[cluster_index].push_back(snp);
        int ref_length = snp.ref.length();
        int alt_length = snp.alt.length();
        int diff_length = abs(ref_length - alt_length);

        if (snp.snp_type == 'I') {
            ins_total += diff_length;
        }
        else if (snp.snp_type == 'D') {
            del_total += diff_length;
        }
    }
}

void RemoveDuplicate::DivisiveHierarchicalClustering(list<vector<SNP> > & snp_clusters){
    // 
    if(snp_clusters.size() == 0) return;
    bool flag = true;
    list<bool> potential_list;
    for(int i = 0; i < potential_list.size(); i++){
        potential_list.push_back(true);
    }
    while(flag){
        flag = false;
        list_size = snp_clusters.size();
        for(int i = 0; i < list_size; i++){
            auto front_cluster = snp_clusters.front();
            auto front_posential = potential_list.front();
            snp_clusters.pop_front();
            potential_list.pop_front();
            if(! front_posential){
                snp_clusters.push_back(front_cluster);
                potential_list.push_back(front_posential);
                continue;
            }

            int max_start = -1;
            int max_end = -1;
            int max_length = -1;
            int start = front_cluster[0].pos + front_cluster[0].ref.length();
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
            vector<SNP> left_snp_list;
            vector<SNP> right_snp_list;
            string separator = genome_sequence.substr(max_start, max_end-max_start);
            for(int k = 0; k < front_cluster.size(); k++){
                auto snp = front_cluster[k];
                int snp_diff = abs(snp.ref.length() - snp.alt.length());
                if(snp.pos <= max_start){
                    if(snp.type == 'I'){
                        left_ins += snp_diff;
                    }else if(snp.type == 'D'){
                        left_del += snp_diff;
                    }
                    left_snp_list.push_back(snp);
                }else{
                    if(snp.type == 'I'){
                        right_ins += snp_diff;
                    }else if(snp.type == 'D'){
                        right_del += snp_diff;
                    }
                    right_snp_list.push_back(snp);
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
            if (separator.length() > 2 * max_change &&
                    (separator.length() > MAX_REPEAT_LEN || !CheckTandemRepeat(separator, max_change)))
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

bool RemoveDuplicate::FindOneMatch(vector<SNP> snp_list, 
    const string subsequence,
    int offset,
    int thread_index)
{
    if(snp_list.size() <= 1) return false;
    unordered_map<string, vector<SNP>> donor_snps;
    for(int i = 1; i < snp_list.size(); i++){
        vector<vector<SNP> > combinations = CreateCombinations(snp_list, i);
        for(int k = 0; k < combinations.size(); k++){
            vector<SNP> comb = combinations[k];
            string alt_sequence = ModifySequenceBySnpList(subsequence, comb, offset);
            if(donor_snps.find(temp_sequence) != donor_snps.end()){
                string matching_result = "";
                matching_result += chromosome_name;
                string parsimonious_ref = subsequence;
                string parsimonious_alt = alt_sequence;
                if(parsimonious_ref == parsimonious_alt){
                    dout << "[Error:RemoveDuplicate::FindOneMatch] in variant, ref == alt";
                }
                int min_parsimonious_len = min(parsimonious_ref.size(), parsimonious_alt.size());
                int chop_left = 0;
                int chop_right = 0;
                for(int i = 0; i < min_parsimonious_len; i++){
                    if(toupper(parsimonious_ref[i]) == toupper(parsimonious_alt[i])){
                        chop_left ++;
                    }else{
                        break;
                    }
                }
                for(int i = min_parsimonious_len-1; i >= 0; i--){
                    if(toupper(parsimonious_ref[i]) == toupper(parsimonious_alt[i])){
                        chop_right ++;
                    }else{
                        break;
                    }
                }
                // 1-based
                matching_result += "\t" + to_string(chop_left + offset + 1);

                parsimonious_ref = parsimonious_ref.substr(chop_left, parsimonious_ref.length() - chop_left - chop_right);
                parsimonious_alt = parsimonious_alt.substr(chop_left, parsimonious_alt.length() - chop_left - chop_right);
                matching_result += "\t" + parsimonious_ref + "\t" + parsimonious_alt;

                string set_matching_string = "";
                for(int m = 0; m < comb.size; m++){
                    auto m_snp = comb[m];
                    set_matching_string += to_string(m_snp.pos+1) + "," + m_snp.ref + "," + m_snp.alt + ";";
                }
                matching_result += "\t"+set_matching_string;

                set_matching_string = "";
                for(int m = 0; m < donor_snps[alt_sequence].size; m++){
                    auto m_snp = donor_snps[alt_sequence][m];
                    set_matching_string += to_string(m_snp.pos+1) + "," + m_snp.ref + "," + m_snp.alt + ";";
                }
                matching_result += "\t"+set_matching_string + "\n";

                if(thread_num == 1){
                    std::lock_guard<std::mutex> guard(complex_match_mutex);
                    complex_match_records[thread_index].push_back(matching_result);
                }else{
                    complex_match_records[thread_index].push_back(matching_result);
                }

                return true;
            }
        }
    }
    return false;
}

void RemoveDuplicate::FindMatches(snp_list, thread_index){
    if(snp_list.size() <= 1) return;
    int min_pos = 0;
    int max_pos = 0;
    sort(snp_list.begin(), snp_list.end());
    min_pos = snp_list[0].pos;
    for(int i = 0; i < snp_list.size(); i++){
        int temp_pos = snp_list[i].pos + snp_list[i].ref.length();
        if(max_pos < temp_pos) max_pos = temp_pos;
    }
    min_pos = max(0, min_pos - 1);
    max_pos = min(max_pos + 1, (int)genome_sequence.length());
    string subsequence = genome_sequence.substr(min_pos, max_pos-min_pos);
    while(snp_list.size() > 1 && CheckMatchCombination(snp_list, subsequence, min_pos, thread_index));
}

void RemoveDuplicate::ClusteringRemoveDuplicateInThread(int start, int end, int thread_index){
    for (int cluster_id = start; cluster_id < end; cluster_id++) {
        if (cluster_snps_map.find(cluster_id) == cluster_snps_map.end()) continue;
        auto & snp_list = cluster_snps_map[cluster_id];
        if(snp_list.size() > 20){
            //DivisiveHierarchicalClustering
            list<vector<SNP> > snp_clusters = {snp_list};
            DivisiveHierarchicalClustering(snp_clusters);
            for(int i = 0; i < snp_clusters.size(); i++){
                //FindMatches for snp_clusters[i]
                FindMatches(snp_clusters[i], thread_index);
            }
        }else{
            //FindMatches
            FindMatches(snp_list, thread_index);
        }
    }
}

void RemoveDuplicate::ClusteringRemoveDuplicateMultiThread(){
    int start = cluster_snps_map.begin()->first;
    int cluster_number = cluster_snps_map.size();
    int cluster_end_boundary = start + cluster_number;
    int cluster_step = cluster_number / thread_num;
    if (cluster_step * thread_num < cluster_number) cluster_step++;
    int end = start + cluster_step;
    
    //initialize vector size, each allocating will have a lock
    vector<string> temp_vector;
    for(int j = 0; j < thread_num; j++){
        complex_match_records.push_back(temp_vector);
    }

    vector<thread> threads;
    //spawn threads
    unsigned i = 0;
    for (; i < thread_num - 1; i++) {
        //threads.push_back(thread(f));
        //dout << "create new thread" << endl;
        int variant_number = 0;
        for (int cluster_id = start; cluster_id < end; cluster_id++) {
            if (cluster_snps_map.find(cluster_id) != cluster_snps_map.end()) {
                variant_number += cluster_snps_map[cluster_id].size();
            }
        }
        threads.push_back(thread(&RemoveDuplicate::ClusteringRemoveDuplicateInThread, this, start, end, i));
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
        ClusteringRemoveDuplicateInThread(start, end, i);
    }

    // call join() on each thread in turn before this function?
    std::for_each(threads.begin(), threads.end(), std::mem_fn(&std::thread::join));
}

void RemoveDuplicate::Deduplicate(string vcf_filename
            string genome_filename,
            bool direct_search,
            string output_prefix)
{
    output_stat_filename = output_prefix + ".stat";
    output_simple_filename = output_prefix + ".simple";
    output_complex_filename = output_prefix + ".complex";

        //------------read genome sequence and decide boundary according to thread number
    dsptime();
    dout << " Read genome sequence file... " << endl;
    ReadGenomeSequence(genome_filename);
    dsptime();
    dout << " Finish reading genome sequence file." << endl;

    dsptime();
    dout << " Read vcf file and remove simple duplications... " << endl;
    ReadRefVCF(vcf_filename);
    dsptime();

    //-------------clustering search
    dsptime();
    dout << " Clustering snps ... " << endl;
    ClusteringSnps();
    dsptime();
    dout << " Finish clustering." << endl;
    dsptime();
    dout << " Detect complex duplications..." << endl;
    ClusteringRemoveDuplicateMultiThread();
    dsptime();
    dout << " Output complex duplications..." << endl;
    
    ofstream output_complex_file;
    output_complex_file.open(output_complex_filename);
    output_complex_file << "#CHR\tPOS\tREF\tALT\tSet1\tSet2" << endl;
    for(int i = 0; i < complex_match_records.size(); i++){
        for (int j = 0; j < complex_match_records[i].size(); j++){
            if(complex_match_records[i][j].find_first_not_of(' ') != std::string::npos){
                output_complex_file << complex_match_records[i][j];
            }
        }
    }
    output_complex_file.close();
    complex_match_records.clear();
    
    return;
}