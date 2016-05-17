#include "wholegenome.h"

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
    int real_chrom_num = 0;
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
                real_chrom_num++;
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
        real_chrom_num++;
    }
    // test

    chrom_num = real_chrom_num;
    dout << "detected chromosome num: " << chrom_num << endl;
//    for(auto it = genome_sequences.begin(); it != genome_sequences.end(); ++it){
//        cout << it->first << ":" << (it->second).length();
//    }
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
	//int genome_sequence_length = genome_sequence.length();
	while (!vcf_file.eof()) { // alternative way is vcf_file != NULL
		string line;
		getline(vcf_file, line, '\n');
		// check ineligible lines
		//dout << line << endl;
		if ((int)line.length() <= 1) continue;
		//if (line.find_first_not_of(' ') == std::string::npos) continue;

		if (line[0] == '#') {
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
		string chr_name = columns[0];
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

        if(snp_ins > VAR_LEN || snp_del > VAR_LEN){
            //dout << "[VarMatch] skip large INDEL with length > 50 bp" << endl;
            continue;
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
            int chr_id = chrname_dict[chr_name];
            if(flag == 0){
                ref_variant_by_chrid[chr_id]->push_back(dv);
            }else{
                que_variant_by_chrid[chr_id]->push_back(dv);
            }
        }

        total_num++;
	}
	vcf_file.close();
	return total_num;
}

bool WholeGenome::ReadVariantFileList(string filename){

}

bool WholeGenome::ParallelClustering(){
    // parallel by chr
    int parallel_steps = chrom_num / thread_num;
    if(parallel_steps*thread_num < chrom_num) parallel_steps += 1;
    int chr_id = 0;
    for(int i = 0; i < parallel_steps; i++){
        vector<thread> threads;
        for(int j = 0; j < thread_num-1 && chr_id < chrom_num-1; j++){
            if(chrname_by_chrid.find(chr_id) != chrname_by_chrid.end()){
                threads.push_back(thread(&WholeGenome::SingleThreadClustering, this, chr_id));
            }
            chr_id ++;
        }
        if(chr_id < chrom_num){
            if(chrname_by_chrid.find(chr_id) != chrname_by_chrid.end()){
                SingleThreadClustering(chr_id);
            }
            chr_id ++;
        }
        std::for_each(threads.begin(), threads.end(), std::mem_fn(&std::thread::join));
        threads.clear();
    }


    for(int i = 0; i < chrom_num; i++){
        variants_by_cluster.insert(variants_by_cluster.end(), variant_cluster_by_chrid[i]->begin(), variant_cluster_by_chrid[i]->end());
    }

    // test output
    dout << endl;
    map<int, int> size_num;
    map<int, int> size_chrid;
    for(int i = 0; i < chrom_num; i++){
        dout << i << ": " << variant_cluster_by_chrid[i]->size() << endl;
        for(int j = 0; j < variant_cluster_by_chrid[i]->size(); j++){
            int temp_size = variant_cluster_by_chrid[i]->at(j).size();
            if(size_num.find(temp_size) != size_num.end()){
                size_num[temp_size] ++;
            }else{
                size_num[temp_size] = 1;
            }
            if(size_chrid.find(temp_size) == size_chrid.end()){
                size_chrid[temp_size] = i;
            }
        }
    }

    cout << endl;
    for(auto it = size_num.begin(); it != size_num.end(); ++it){
        dout << it->first << ": " << it->second << endl;
    }

//    cout << endl;
//    cout << "size and location:" << endl;
//    for(auto it = size_chrid.begin(); it != size_chrid.end(); ++it){
//        dout << it->first << ": " << it->second << endl;
//    }

    return true;
}

bool WholeGenome::ParallelMatching(){

}

bool WholeGenome::TBBMatching()
{

}

bool WholeGenome::ClusteringMatchInThread(int start, int end, int thread_index) {

	for (int cluster_id = start; cluster_id < end; cluster_id++) {
        if(cluster_id >= variants_by_cluster.size()) break;
        MatchingSingleCluster(cluster_id, thread_index);
	}
	return true;
}

// private
void WholeGenome::ClusteringMatchMultiThread() {
	int start = 0;
	int cluster_number = variants_by_cluster.size(); // cluster number
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
		threads.push_back(thread(&WholeGenome::ClusteringMatchInThread, this, start, end, i));
		start = end;
		end = start + cluster_step;
	}
	// also you need to do a job in main thread
	// i equals to (thread_num - 1)
	if (i != thread_num - 1) {
		dout << "[Error] thread number not match" << endl;
	}
	if (start >= variants_by_cluster.size()) {
		dout << "[Error] index out of map range" << endl;
	}
	else {
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
				output_complex_file << chromosome_name << "\t" << complex_match_records[i]->at(j);
			}
		}
	}
	output_complex_file.close();

	for (int j = 0; j < thread_num; j++) {
		delete complex_match_records[j];
	}
	delete[] complex_match_records;

	total_ref_complex = 0;
	total_que_complex = 0;
	for (int i = 0; i < complex_ref_match_num.size(); i++)
		total_ref_complex += complex_ref_match_num[i];
	for (int i = 0; i < complex_que_match_num.size(); i++)
		total_que_complex += complex_que_match_num[i];

	cout << "complex match: " << total_ref_complex << "," << total_que_complex << endl;
}


bool WholeGenome::MatchingSingleCluster(int cluster_index, int thread_index){

    vector<VariantIndicator> vi_list = variants_by_cluster[cluster_index];
    if(vi_list.size() <= 1) return false;
    // create variant_list from vi_list;
    vector<DiploidVariant> variant_list;
    int chr_id = -1;
    for(int i = 0; i < vi_list.size(); i++){
        VariantIndicator vi = vi_list[i];
        chr_id = vi.chr_id;
        int var_id = vi.var_id;
        if(vi.refer){
            variant_list.push_back(ref_variant_by_chrid[chr_id]->at(var_id));
        }else{
            variant_list.push_back(que_variant_by_chrid[chr_id]->at(var_id));
        }
    }
    if(chr_id == -1 || chr_id >= chrom_num){
        cout << "[VarMatch] Error in matching single cluster" << endl;
        return false;
    }


    // after this, all are copied from diploid.cpp
    // Warning: remember to change all genome_sequence to genome_sequences[chr_id]
    // Please!! use search!!


    //*****************************************************************************************************************
    //*****************************************************************************************************************
    //*****************************************************************************************************************
    //*****************************************************************************************************************
    //*****************************************************************************************************************


    sort(variant_list.begin(), variant_list.end()); // here we need to sort
    vector<DiploidVariant> separate_var_list[2];
	// separate into ref and que
	int total_mil = 0;
	int total_mdl = 0;
	int min_pos = genome_sequences[chr_id].length() + 1;
	int max_pos = -1;
	for (int i = 0; i < variant_list.size(); i++) {
		int flag = variant_list[i].flag; // flag indicate if the variant is from ref set(0) or query set(1)
		int pos = variant_list[i].pos;
		separate_var_list[flag].push_back(variant_list[i]);
		total_mil += variant_list[i].mil;
		total_mdl += variant_list[i].mdl;
		auto ref_sequence = variant_list[i].ref;
		auto alt_sequences = variant_list[i].alts;
		min_pos = min(pos, min_pos);
		max_pos = max((int)(pos + ref_sequence.length()), max_pos);
	}
	min_pos = max(min_pos - 1, 0);
	max_pos = min(max_pos + 1, (int)genome_sequences[chr_id].length()); //exclusive
	if (separate_var_list[0].size() == 0 || separate_var_list[1].size() == 0) {
		return false;
	}
	if (separate_var_list[0].size() == 1 && separate_var_list[1].size() == 1){
        // try direct match to save time
        if(separate_var_list[0][0] == separate_var_list[1][0]){
            complex_ref_match_num[thread_index]++;
            complex_que_match_num[thread_index]++;

            DiploidVariant tv = separate_var_list[0][0];
            string match_record = to_string(tv.pos+1) + "\t" + tv.ref + "\t" + tv.alts[0];
            if(tv.multi_alts) match_record += "/" + tv.alts[1];
            match_record += "\t.\t.\t.\t.\t.\n";
            complex_match_records[thread_index]->push_back(match_record);
            // output match result
            return true;
        }
        // if not match, still can match by changing genome
	}else if(separate_var_list[0].size() == 1 || separate_var_list[1].size() == 1){
        int flag = 0;
        if(separate_var_list[1].size() == 1) flag = 1;
        int r_flag = 1-flag;
        if(separate_var_list[r_flag].size() > 4){
            int total_r_mdl = 0;
            int total_r_mil = 0;

            for(int k = 0; k < separate_var_list[r_flag].size(); k++){
                DiploidVariant var = separate_var_list[r_flag][k];
                int var_mdl = var.mdl;
                int var_mil = var.mil;
                int ref_length = var.ref.length();
                total_r_mdl += var_mdl;
                total_r_mil += var_mil;
            }

            if(max(separate_var_list[flag][0].mdl, separate_var_list[flag][0].mil) > max(total_r_mdl, total_r_mil)) return false;
        }
	}

	// remove singular variant
    vector<bool> appliable_flag[2];
    int total_change = total_mil+total_mdl;
    for(int i = 0; i < 2; i++){
        for(int k = 0; k < separate_var_list[i].size(); k++){
            DiploidVariant cur_var = separate_var_list[i][k];
            int max_change = max(cur_var.mil, cur_var.mdl);
            if(max_change > total_change-max_change){
                appliable_flag[i].push_back(false);
            }else{
                appliable_flag[i].push_back(true);
            }
        }
    }
	string subsequence = genome_sequences[chr_id].substr(min_pos, max_pos - min_pos);

	ToUpper(subsequence); // subsequence only contains upper char
	int offset = min_pos;
	int subsequence_length = max_pos - min_pos;
	list<VariantSelection> variant_selections; // sorted by last matched donor length
	VariantSelection best_selection;
	VariantSelection dummy;

    bool overlap_detected = false;

    for(int i = 0; i < 2; i++){
        int largest_pos = 0;
        for(int k = 0; k < separate_var_list[i].size(); k++){
            auto var = separate_var_list[i][k];
            if(var.pos <= largest_pos){
                overlap_detected = true;
                break;
            }
            largest_pos = max(largest_pos, (int)(var.pos+var.ref.length()));
        }
        if(overlap_detected) break;
    }
    dummy.overlap_detected = overlap_detected;

    variant_selections.push_back(dummy);

    map<string, int> score_by_consistent_donor; // donor should be sorted

    while(variant_selections.size() != 0){
        VariantSelection current_selection = variant_selections.front();
        variant_selections.pop_front();

        bool get_ref_var = true;
        int ref_var_taken = current_selection.phasing_vectors[0].size();
        int que_var_taken = current_selection.phasing_vectors[1].size();
        if(ref_var_taken >= separate_var_list[0].size()){
            get_ref_var = false;
        }else if(que_var_taken < separate_var_list[1].size()){
              if(current_selection.genome_position[0] > current_selection.genome_position[1]){
                get_ref_var = false;
              }else if( current_selection.genome_position[0] == current_selection.genome_position[1]){
                if(min(current_selection.donor_length[0], current_selection.donor_length[1]) > min(current_selection.donor_length[2], current_selection.donor_length[3])){
                    get_ref_var = false;
                }
              }
        }

        DiploidVariant current_variant;
        bool can_take_variant = true;
        if(get_ref_var){
            can_take_variant = appliable_flag[0][ref_var_taken];
            current_variant = separate_var_list[0][ref_var_taken];
        }else{
            can_take_variant = appliable_flag[1][que_var_taken];
            current_variant = separate_var_list[1][que_var_taken];
        }

        int current_flag = current_variant.flag;

//            cout << "current selection" << endl;
//            PrintSelection(current_selection);
//            cout << "add variant";
//            PrintVariant(current_variant);

        bool added = false;
        // make choose decision before not choose decision, save del times
        if(can_take_variant){
            added = AddVariantToSelection(variant_selections,
                                current_selection,
                                current_variant,
                                0,
                                separate_var_list,
                                subsequence,
                                offset,
                                best_selection);
    //            cout << "added state : " << added << endl;
    //            PrintSelectionsList(variant_selections);

            if(current_variant.heterozygous){
                added = AddVariantToSelection(variant_selections,
                                    current_selection,
                                    current_variant,
                                    1,
                                    separate_var_list,
                                    subsequence,
                                    offset,
                                    best_selection);
    //                cout << "added state : " << added << endl;
    //                PrintSelectionsList(variant_selections);
            }
        }

       added= AddVariantToSelection(variant_selections,
                            current_selection,
                            current_variant,
                            -1,
                            separate_var_list,
                            subsequence,
                            offset,
                            best_selection);
//            cout << "added state : " << added << endl;
//            PrintSelectionsList(variant_selections);

    }
//    dout << best_selection.score << endl;
    if (best_selection.score <= 0) return false;
//    cout << "best selection: " << endl;
//    PrintSelection(best_selection);

    complex_ref_match_num[thread_index] += best_selection.separate_score[0];
    complex_que_match_num[thread_index] += best_selection.separate_score[1];

    bool multiple_match = true;



    if(best_selection.donor_sequences[0] == best_selection.donor_sequences[1]) multiple_match = true;


    string parsimonious_ref = subsequence;
    string parsimonious_alt0 = best_selection.donor_sequences[0];
    string parsimonious_alt1 = best_selection.donor_sequences[1];

    int parsimonious_pos = NormalizeVariantSequence(offset,
                             parsimonious_ref,
                             parsimonious_alt0,
                             parsimonious_alt1,
                             chr_id);

    string match_record = to_string(parsimonious_pos+1) + "\t" + parsimonious_ref + "\t" + parsimonious_alt0;
    if(multiple_match) match_record += "/" + parsimonious_alt1;

    string vcf_record[2];
    string phasing_record[2];

	for (int i = 0; i < 2; i++) {
		auto final_iter = separate_var_list[i].size()-1;
		vector<int> phasing_vector = best_selection.phasing_vectors[i];
		for (int k = 0; k < separate_var_list[i].size(); k++) {
            int phasing = phasing_vector[k];
            if(phasing == -1) continue;
            DiploidVariant variant = separate_var_list[i][k];
            string alt_string = variant.alts[0];
            if(variant.multi_alts){
                alt_string += "/" + variant.alts[1];
            }
            string phasing_string = "";
            if(phasing == 0){
                phasing_string += "1";
                if(variant.heterozygous){
                    if(variant.multi_alts){
                        phasing_string += "|2";
                    }else{
                        phasing_string += "|0";
                    }
                }else{
                    phasing_string += "|1";
                }
            }else if(phasing == 1){
                if(variant.multi_alts){
                    phasing_string += "2|1";
                }else{
                    phasing_string += "0|1";
                }
            }
            string variant_record = to_string(variant.pos+1) + "," + variant.ref + "," + alt_string;
            vcf_record[i] += variant_record;
            phasing_record[i] += phasing_string;
            if (k != final_iter) {
                vcf_record[i] += ";";
                phasing_record[i] += ";";
            }
		}
	}
	match_record += "\t" + vcf_record[0] + "\t" + vcf_record[1];
    match_record += "\t" + phasing_record[0] + "\t" + phasing_record[1];
	match_record += "\t" + to_string(best_selection.score) + "\n";

	complex_match_records[thread_index]->push_back(match_record);
    // add matching result
    return true;
}

int WholeGenome::NormalizeVariantSequence(int pos, string & parsimonious_ref, string & parsimonious_alt0, string & parsimonious_alt1, int chr_id) {

	int left_index = pos;
	if (genome_sequences[chr_id].size() == 0) return -1;
	if (parsimonious_ref.size() == 1 && parsimonious_alt0.size() == 1 && parsimonious_alt1.size() == 1) return true;

	bool change_in_allels = true;
	while (change_in_allels) {
		change_in_allels = false;
		if (parsimonious_ref.back() == parsimonious_alt0.back() && parsimonious_ref.back() == parsimonious_alt1.back() ) {
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
			char left_char = toupper(genome_sequences[chr_id][left_index]);
			parsimonious_ref = left_char + parsimonious_ref;
			parsimonious_alt0 = left_char + parsimonious_alt0;
			parsimonious_alt1 = left_char + parsimonious_alt1;
		}
	}
	while (parsimonious_ref[0] == parsimonious_alt0[0] &&
            parsimonious_ref[0] == parsimonious_alt1[0] &&
            parsimonious_ref.size() > 1 &&
            parsimonious_alt0.size() > 1 &&
            parsimonious_alt1.size() > 1)
    {
		parsimonious_ref.erase(0, 1);
		parsimonious_alt0.erase(0, 1);
		parsimonious_alt1.erase(0, 1);
        left_index ++; // left_index indicates variant position, if truncate the leftmost, then
	}
	return left_index;
}

void WholeGenome::SingleThreadClustering(int chr_id) {
	int ins_len[2] = { 0 };
	int del_len[2] = { 0 };
	int c_start = 0;
	int c_end = 0;
    sort(ref_variant_by_chrid[chr_id]->begin(), ref_variant_by_chrid[chr_id]->end());
    sort(que_variant_by_chrid[chr_id]->begin(), que_variant_by_chrid[chr_id]->end());
    int ref_size = ref_variant_by_chrid[chr_id]->size();
    int que_size = que_variant_by_chrid[chr_id]->size();
    //dout << chr_id << "," << ref_size << "," << que_size << endl;

    int ref_index = 0;
    int que_index = 0;
    bool not_first = false;
    DiploidVariant snp;
    vector<VariantIndicator> vi_list;
    while (ref_index < ref_size || que_index < que_size) {
		bool take_que = true;
		if(ref_index < ref_size && que_index < que_size){
            if(ref_variant_by_chrid[chr_id]->at(ref_index).pos < que_variant_by_chrid[chr_id]->at(que_index).pos){
                take_que = false;
            }
		}else if(ref_index < ref_size){
            take_que = false;
		}
        int var_index;
		if(take_que){

            snp = que_variant_by_chrid[chr_id]->at(que_index);
            //cout << "q |" << que_index << "," << snp.pos << endl;
            var_index = que_index;
            que_index++;
		}else{
            snp = ref_variant_by_chrid[chr_id]->at(ref_index);
            //cout << "r |" << ref_index << "," << snp.pos << endl;
            var_index = ref_index;
            ref_index++;
		}
		// check if need to separator clusters
		if (not_first) {
			c_end = snp.pos;
			if (c_end - c_start >= 2) {
                int separator_length = c_end - c_start;
				string separator = genome_sequences[chr_id].substr(c_start, separator_length);
				int max_change = max(ins_len[0] + del_len[1], ins_len[1] + del_len[0]);
				bool separate_cluster = false;
				if(max_change == 0){
                    separate_cluster = true;
				}
				else if (separator_length > 2 * max_change &&
					(separator_length > MAX_REPEAT_LEN || !CheckTandemRepeat(separator, max_change)))
				{
				    separate_cluster = true;
				}

				if(separate_cluster){
                    variant_cluster_by_chrid[chr_id]->push_back(vi_list);
                    vi_list.clear();
					ins_len[0] = 0;
					del_len[0] = 0;
					ins_len[1] = 0;
					del_len[1] = 0;
					c_start = 0; // re-assign c_start
				}
			}
		}
		c_start = max(c_start, snp.pos + (int)snp.ref.length() );
        VariantIndicator current_variant_indicator(chr_id, var_index, !take_que);
        vi_list.push_back(current_variant_indicator);
		//cluster_vars_map[cluster_index].push_back(snp);
		if(!not_first) not_first = true;
		int ref_length = (int)(snp.ref.length());
		int flag = snp.flag;
//        DiploidVariant snp = front_cluster[k];
//        int rq = snp.flag;
        ins_len[flag] += snp.mil;
        del_len[flag] += snp.mdl;
	}
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

    ParallelClustering();
    ClusteringMatchMultiThread();
}

void WholeGenome::DirectMatch(string ref_vcf,
    string query_vcf,
    bool match_genometype,
    bool normalization)
{

    this->normalization = normalization;
    this->match_genotype = match_genometype;
    int ref_variant_num = ReadReferenceVariants(ref_vcf);
    int que_variant_num = ReadQueryVariants(query_vcf);
    dout << ref_variant_num << "," << que_variant_num << endl;
    int match_num = 0;
    for(int i = 0; i < chrom_num; i++){
        if(ref_variant_by_chrid[i]->size() == 0 || que_variant_by_chrid[i]->size() == 0)
            continue;
        map<int, DiploidVariant> ref_variant_by_pos;
        for(int j = 0; j < ref_variant_by_chrid[i]->size(); j++){
            DiploidVariant var = ref_variant_by_chrid[i]->at(j);
            int pos = var.pos;
            ref_variant_by_pos[pos] = var;
        }

        for(int j = 0; j < que_variant_by_chrid[i]->size(); j++){
            DiploidVariant var = que_variant_by_chrid[i]->at(j);
            int pos = var.pos;
            if(ref_variant_by_pos.find(pos) == ref_variant_by_pos.end())
                continue;
            DiploidVariant ref_var = ref_variant_by_pos[pos];
            if (var == ref_var){
                match_num ++;
            }
        }
    }
    dout << "matched variants: " << match_num << endl;
}























