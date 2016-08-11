#include "mergevariants.h"

using namespace std;

// constructor
WholeGenome::WholeGenome(int thread_num_,
    string output_dir_,
    bool pr_curves){

    thread_num = thread_num_;
    chrom_num = 24;

    output_dir = output_dir_;

    accumulated_query_number = 0;

    global_baseline_id = 0;

    //thread_num = thread_num_;
    //dout << "WholeGenome() Thread Number: " << thread_num << endl;

    ref_variant_by_chrid = new vector<DiploidVariant>*[chrom_num];
    //variant_matching_condition = new unordered_map<string, long long> * [chrom_num];

	for (int j = 0; j < chrom_num; j++) {
		ref_variant_by_chrid[j] = new vector<DiploidVariant>;
        //variant_matching_condition[j] = new unordered_map<string, long long>;
	}

    que_variant_by_chrid = new vector<DiploidVariant>*[chrom_num];
    for (int j = 0; j < chrom_num; j++) {
        que_variant_by_chrid[j] = new vector<DiploidVariant>;
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

    if(pr_curves){
        per_list = {0.0, 0.1, 0.2, 0.3, 0.9};
    }else{
        per_list = {0.0};
    }

}

inline int WholeGenome::GetIndexFromMatchScore(int score_unit, int match_mode, int score_scheme){
    int result = 0;
    result |= score_unit & 1;
    result <<= 1;
    result |= match_mode & 1;
    result <<= 2;
    result |= score_scheme & 3;
    return result;
}

// distructor
WholeGenome::~WholeGenome(){

    for(int j = 0; j < chrom_num; j++){
        ref_variant_by_chrid[j]->clear();
        delete ref_variant_by_chrid[j];
        que_variant_by_chrid[j]->clear();
        delete que_variant_by_chrid[j];
        //variant_matching_condition[j]->clear();
        //delete variant_matching_condition[j];
    }
    delete[] ref_variant_by_chrid;
    delete[] que_variant_by_chrid;
    //delete[] variant_matching_condition;
}

bool WholeGenome::ReadWholeGenomeSequence(string filename){
    std::ifstream input(filename);
    if(!input.good()){
        std::cerr << "Error opening '"<<filename<<"'. Bailing out." << std::endl;
        return false;
    }

    std::string line, name, content;
    int real_chrom_num = 0;
    int chr_id = 0;
    int current_id = -1;
    while( std::getline( input, line ).good() ){
        if( line.empty() || line[0] == '>' ){ // Identifier marker
            if( !name.empty() ){ // Print out what we read from the last entry
                //std::cout << name << " : " << content << std::endl;
                if(chrname_dict.find(name) == chrname_dict.end()){
                    cout << "[VarMatch] Error: detected chromosome name: " << name <<" does not exist in human genome." << endl;
                    return false;
                }
                //int chr_id = chrname_dict[name];
                if(chrid_by_chrname.find(name) == chrid_by_chrname.end()){
                    chrid_by_chrname[name] = chr_id;
                    chr_id++;
                }
                current_id = chrid_by_chrname[name];
                chrname_by_chrid[current_id] = name;
                genome_sequences[current_id] = content;
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
        if(chrid_by_chrname.find(name) == chrid_by_chrname.end()){
            chrid_by_chrname[name] = chr_id;
            chr_id++;
        }
        current_id = chrid_by_chrname[name];
        chrname_by_chrid[current_id] = name;
        genome_sequences[current_id] = content;
        real_chrom_num++;
    }
    // test

    chrom_num = real_chrom_num;
    //dout << "detected chromosome num: " << chrom_num << endl;
//    for(auto it = genome_sequences.begin(); it != genome_sequences.end(); ++it){
//        cout << it->first << ":" << (it->second).length();
//    }
    return true;
}

bool WholeGenome::ReadGenomeSequenceList(string filename){

}

int WholeGenome::ReadWholeGenomeVariant(string filename, bool flag){
    int total_num = 0;
    int long_num = 0;
    double QUAL_LOWER_BOUND = 0.1;

	ifstream vcf_file;
	vcf_file.open(filename.c_str());
	if (!vcf_file.good()) {
		cout << "[VarMatch] Error: can not open vcf file" << endl;
		return -1;
	}

    vector<float> quality_list;

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
			if(match_mode_indicator != 1){
                cout << "[VarMatch] Warning: not enough information in VCF file for genotype matching." << endl;
                cout << "[VarMatch] \tAutomatically turn off genotype matching module " << filename << endl;
                match_mode_indicator = 1;
                //continue;
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
        double quality;
        if(columns[5] == "."){
            quality = 1.0;
        }else{
            quality = stod(columns[5]);
        }

        if(flag){
            quality_list.push_back(quality);
        }

		ToUpper(ref);
		ToUpper(alt_line);

		bool is_heterozygous_variant = false;
		bool is_multi_alternatives = false;
        bool is_zero_one_var = false;
        
        vector<string> genotype_columns;
		
        if (match_mode_indicator != 1) { // match mode indicator is -1 or 0
			if (genotype_index < 0) {
                // change genotype index
                auto formats = split(columns[8], ':');
                for (int i = 0; i < formats.size(); i++) {
                    if (formats[i] == "GT") {
                        genotype_index = i;
                        break;
                    }
                }
                // if GT not found
                if(genotype_index < 0){
                    if(match_mode_indicator != 1 && match_mode_indicator != 1){
                        cout << "[VarMatch] Warning: VCF entry does not contain genotype information." << endl;
                        cout << "[VarMatch] \tAutomatically turn off genotype matching mode. " << endl;
                        match_mode_indicator = 1;
                    }
                }
			}

            
            if(match_mode_indicator != 1){

    			auto additionals = split(columns[9], ':');
                genotype_columns = split(additionals[genotype_index], genotype_separator);

                if(genotype_columns.size() != 2){
                    
                    if(genotype_separator == '/'){
                        genotype_separator = '|';
                    }else{
                        genotype_separator = '/';
                    }
                    genotype_columns = split(additionals[genotype_index], genotype_separator);
                }

    			// normalize format of genotype: sorted, separated by |
    			if (genotype_columns.size() != 2) {
    				cout << "[VarMatch] Warning: Unrecognized Genotype: " << additionals[genotype_index] << endl;
                    cout << "[VarMatch] \tAutomatically turn off genotype matching mode." << endl;
                    match_mode_indicator = 1;
    			}
    			else {
    				if (genotype_columns[0] != genotype_columns[1]) {
    					is_heterozygous_variant = true;
    				}
                    if (genotype_columns[1] == "0" && genotype_columns[0] == "0") {
                        //cout << "Skip Variants when both genotype is refernce allele: " << line << endl;   
                        continue;
                    }
                    if(genotype_columns[0] == "0" || genotype_columns[1] == "0"){
                        is_zero_one_var = true;
                    }
    			}
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

        if(alt_list.size() > 2){
            if(match_mode_indicator != 1){
            vector<string> temp_alt_list = alt_list;
            alt_list.clear();
            for(int i = 0; i < 2; i++){
                int alt_indicator = stoi(genotype_columns[i]);
                if(alt_indicator == 0) continue;
                alt_list.push_back(temp_alt_list[alt_indicator-1]);
            }
            }else{
                vector<string> temp_alt_list = alt_list;
                alt_list.clear();
                alt_list.push_back(temp_alt_list[0]);
                alt_list.push_back(temp_alt_list[1]);
            }
        }

        int snp_ins = max(0, (int)alt_list[0].length() - (int)ref.length());
        int snp_del = max(0, (int)ref.length() - (int)alt_list[0].length());
        if(is_multi_alternatives){
            snp_ins = max(snp_ins, (int)alt_list[1].length() - (int)ref.length());
            snp_del = max(snp_del, (int)ref.length() - (int)alt_list[1].length());
        }

        if(snp_ins > VAR_LEN || snp_del > VAR_LEN){
            //dout << "[VarMatch] skip large INDEL with length > " << VAR_LEN << "| "<< line <<endl;
            long_num ++;
            continue;
        }

		DiploidVariant dv(pos, ref, alt_list, is_heterozygous_variant, is_multi_alternatives, snp_del, snp_ins, flag, quality, is_zero_one_var);
		//if (normalization) {
			//NormalizeDiploidVariant(dv);
		//}
        if(chrid_by_chrname.find(chr_name) != chrid_by_chrname.end()){
            int chr_id = chrid_by_chrname[chr_name];
            if(flag == false){
                dv.matching_condition = 1;
                dv.accumulated_vcf_num = 1;
                global_baseline_id ++;
                dv.unique = global_baseline_id;
                ref_variant_by_chrid[chr_id]->push_back(dv);

                //baseline_variant_strings.push_back(line);
            }else{
                dv.matching_condition = 0;
                dv.accumulated_vcf_num = 1;
                que_variant_by_chrid[chr_id]->push_back(dv);
                query_variant_strings.push_back(line);
            }
        }else{
            //cout << "[VarMatch] skip current variant as no corresponding reference genome sequence found." << endl;
            cout << "[VarMatch] Warning: variant skipped, no reference genome sequence found for chromosome name: " << chr_name << "." << endl;
            continue;
            int chr_id = chrname_dict[chr_name];
            if(flag == false){
                ref_variant_by_chrid[chr_id]->push_back(dv);
                //baseline_variant_strings.push_back(line);
            }else{
                que_variant_by_chrid[chr_id]->push_back(dv);
                query_variant_strings.push_back(line);
            }
        }

        total_num++;
	}
	vcf_file.close();

    if(flag){
        sort(quality_list.begin(), quality_list.end());
        auto qual_lower_it = lower_bound(quality_list.begin(), quality_list.end(), QUAL_LOWER_BOUND);
        int qual_lower_index = qual_lower_it - quality_list.begin();
        int rest_size = quality_list.size() - qual_lower_index;

        vector<float> temp_percentage_list;
        temp_percentage_list.push_back(0.0);
        threshold_list.push_back(0.0);
        
        for(int i = 1; i < per_list.size(); i++){
            int additional_index = (int)(rest_size * per_list[i]);
            int real_index = qual_lower_index + additional_index;
            if(real_index >= quality_list.size()) real_index = quality_list.size() - 1;
            double quality = quality_list[real_index];
            threshold_list.push_back(quality);

            auto quality_lowit = lower_bound(quality_list.begin(), quality_list.end(), quality);
            int quality_low_index = quality_lowit - quality_list.begin();
            // following program will retain variants >= quality threshold
            
            int quality_size = quality_low_index + 1; // counting number, +/- 1 does not matter
            if(quality_size > quality_list.size()) quality_size = quality_list.size();
            double percentage = (double)quality_size/ quality_list.size();
            temp_percentage_list.push_back(percentage);
        }
        threshold_num = threshold_list.size();
        // revice percentage
        per_list = temp_percentage_list;
    }
    //cout << flag << "," << total_num << "," << long_num << endl;
	return total_num;
}

bool WholeGenome::ReadVariantFileList(string filename){

}

int WholeGenome::ScoreEditDistance(DiploidVariant & dv, int allele_indicator){
    return EditDistance(dv.ref, dv.alts[allele_indicator]);
}

inline int WholeGenome::EditDistance(const std::string& s1, const std::string& s2)
{
	const std::size_t len1 = s1.size(), len2 = s2.size();
	std::vector<unsigned int> col(len2+1), prevCol(len2+1);

	for (unsigned int i = 0; i < prevCol.size(); i++)
		prevCol[i] = i;
	for (unsigned int i = 0; i < len1; i++) {
		col[0] = i+1;
		for (unsigned int j = 0; j < len2; j++)
                        // note that std::min({arg1, arg2, arg3}) works only in C++11,
                        // for C++98 use std::min(std::min(arg1, arg2), arg3)
			col[j+1] = std::min({ prevCol[1 + j] + 1, col[j] + 1, prevCol[j] + (s1[i]==s2[j] ? 0 : 1) });
		col.swap(prevCol);
	}
	return prevCol[len2];
}

// Needleman Wunsch Initialization
inline void WholeGenome::initialize_score_matrix(int **score, char **trackBack, int M, int N)
{
    for (int i = 0; i < M+1; i++)
    {
        score[0][i] = i * -1;
        trackBack[0][i] = '-';
    }

    for (int i = 0; i < N+1; i++)
    {
        score[i][0] = i * -1;
        trackBack[i][0] = '|';
    }

    trackBack[0][0] = '*';
}

int WholeGenome::needleman_wunsch(string S1, string S2, string &R1, string &R2)
{
    int M = S1.length();
    int N = S2.length();
    /*
    0MMMMMMMMMMMMMMMM
    N
    N
    N
    N
    N
    N
    so the matrix is N*M
    */
    int **score = new int *[N+1];
    for (int i = 0; i <= N; i++)
    {
        score[i] = new int [M+1];
    }

    char **trackBack = new char *[N+1];
    // * for match, - for ->, | for moving downward
    for (int i = 0; i <= N; i++)
    {
        trackBack[i] = new char [M+1];
    }
    R1 = "";
    R2 = "";
    initialize_score_matrix(score, trackBack, M, N);

    for (int i = 1; i <=N; i++)
    {
        for (int k = 1; k <= M; k++)
        {
            char S1_k = S1[k-1];
            char S2_i = S2[i-1];
            int matchingCost = score[i-1][k-1];
            if(S1_k != S2_i) matchingCost--;
            int rightCost = score[i][k-1] - 1;
            int downCost = score[i-1][k] - 1;
            if (matchingCost > rightCost && matchingCost > downCost)
            {
                score[i][k] = matchingCost;
                trackBack[i][k] = '*';
            }else if(rightCost >= downCost)
            {
                score[i][k] = rightCost;
                trackBack[i][k] = '-';
            }else
            {
                score[i][k] = downCost;
                trackBack[i][k] = '|';
            }
        }
    }

    //trackBack
    int n = N;
    int m = M;
    while(n > 0 || m > 0)
    {
        if (trackBack[n][m] == '*')
        {
            R1 += S1[m-1];
            R2 += S2[n-1];
            n--;
            m--;
        }else if(trackBack[n][m] == '-')
        {
            R1 += S1[m-1];
            R2 += '-';
            m--;
        }else if(trackBack[n][m] == '|')
        {
            R1 += '-';
            R2 += S2[n-1];
            n--;
        }
    }
    reverse(R1.begin(), R1.end());
    reverse(R2.begin(), R2.end());

    int result = score[N][M];

    for (int i = 0; i <= N; i++)
    {
        delete score[i];
        delete trackBack[i];
    }

    delete score;
    delete trackBack;

    return result;
}

void WholeGenome::GenerateAltVector(string ref, string alt, vector<string> & alt_vector){
    if(ref.size() == 0) return;
    string ref_match = "";
    string alt_match = "";
    needleman_wunsch(ref, alt, ref_match, alt_match);
    int current_ref_index = -1;
    for(int i = 0; i < ref.size(); i++){
        alt_vector.push_back("");
    }

    for(int i = 0; i < ref_match.size(); i++){
        if(ref_match[i] == '-'){
            if(current_ref_index < 0){
                alt_vector[0].push_back(alt_match[i]);
            }else{
                alt_vector[current_ref_index].push_back(alt_match[i]);
            }
        }else if(alt_match[i] == '-'){
            // pass
            current_ref_index ++;
        }else{
            current_ref_index ++;
            if(current_ref_index >= ref.size()){
                alt_vector[ref.size()-1].push_back(alt_match[i]);
            }
            alt_vector[current_ref_index].push_back(alt_match[i]);
        }
    }
    return;
}

bool WholeGenome::ParallelClustering(){
    // parallel by chr
    variant_cluster_by_chrid = new vector<vector<VariantIndicator>> *[chrom_num];
    for (int j = 0; j < chrom_num; j++) {
        variant_cluster_by_chrid[j] = new vector<vector<VariantIndicator>>;
    }

    int parallel_steps = chrom_num / thread_num;
    if(parallel_steps*thread_num < chrom_num) parallel_steps += 1;
    int chr_id = 0;
    for(int i = 0; i < parallel_steps; i++){
        vector<thread> threads;
        for(int j = 0; j < thread_num-1 && chr_id < chrom_num-1; j++){
            if(chrname_by_chrid.find(chr_id) != chrname_by_chrid.end()){
                if(ref_variant_by_chrid[chr_id]->size() > 0 && que_variant_by_chrid[chr_id]->size() > 0){
                    threads.push_back(thread(&WholeGenome::SingleThreadClustering, this, chr_id));
                }
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
        if(variant_cluster_by_chrid[i]->size() > 0){
            variants_by_cluster.insert(variants_by_cluster.end(), variant_cluster_by_chrid[i]->begin(), variant_cluster_by_chrid[i]->end());
        }
    }

    // test output
    //dout << endl;
    map<int, int> size_num;
    map<int, int> size_chrid;
    for(int i = 0; i < chrom_num; i++){
        //dout << i << ": " << variant_cluster_by_chrid[i]->size() << endl;
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

    //cout << endl;
    //for(auto it = size_num.begin(); it != size_num.end(); ++it){
    //    dout << it->first << ": " << it->second << endl;
    //}

//    cout << endl;
//    cout << "size and location:" << endl;
//    for(auto it = size_chrid.begin(); it != size_chrid.end(); ++it){
//        dout << it->first << ": " << it->second << endl;
//    }
        // clean at the end of function

    for(int j = 0; j < chrom_num; j++){
        variant_cluster_by_chrid[j]->clear();
        delete variant_cluster_by_chrid[j];
    }
    delete[] variant_cluster_by_chrid;

    return true;
}

bool WholeGenome::ParallelMatching(){

}

bool WholeGenome::TBBMatching()
{

}


bool WholeGenome::CheckTandemRepeat(string sequence, int unit_threshold) {
    int sequence_length = (int)sequence.length();
    //cout << sequence_length << "," << unit_threshold << endl;
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

// preprocess
bool WholeGenome::MatchVariantListInThread(int thread_index, 
    int threshold_index,
    int chr_id,
    vector<DiploidVariant> & variant_list,
    unordered_map<int, VariantIndicator> & baseline_unique_indicator,
    int cluster_id){
    //===================================================
    //cout << cluster_id << endl;
    sort(variant_list.begin(), variant_list.end());
    // decide reference sequence
    vector<DiploidVariant> separate_var_list[2];
    vector<Interval> intervals;
    // separate into ref and que
    int total_mil = 0;
    int total_mdl = 0;
    int min_pos = genome_sequences[chr_id].length() + 1;
    int max_pos = -1;
    for (int i = 0; i < variant_list.size(); i++) {
        int flag = 0;
        if (variant_list[i].flag) flag = 1; // flag indicate if the variant is from ref set(0) or query set(1)
        int pos = variant_list[i].pos;
        separate_var_list[flag].push_back(variant_list[i]);
        total_mil += variant_list[i].mil;
        total_mdl += variant_list[i].mdl;
        auto ref_sequence = variant_list[i].ref;
        auto alt_sequences = variant_list[i].alts;
        min_pos = min(pos, min_pos);
        max_pos = max((int)(pos + ref_sequence.length()), max_pos);

        int end_pos = pos + ref_sequence.length() - 1; // included end position!!
        intervals.push_back(Interval(pos, end_pos));
    }
    min_pos = max(min_pos - 1, 0);
    max_pos = min(max_pos + 1, (int)genome_sequences[chr_id].length()); //exclusive

    if (separate_var_list[0].size() == 0 || separate_var_list[1].size() == 0) {
        //dout << separate_var_list[0].size() << ", " << separate_var_list[1].size() << endl;
        // no match exist, need to update reference
        if(threshold_index == 0){
            for(auto baseline_var : separate_var_list[0]){
                int baseline_unique = baseline_var.unique;
                auto vi = baseline_unique_indicator[baseline_unique];
                ref_variant_by_chrid[vi.chr_id]->at(vi.var_id).matching_condition <<= 1;
                ref_variant_by_chrid[vi.chr_id]->at(vi.var_id).accumulated_vcf_num += 1;
            }

            for(auto query_var: separate_var_list[1]){
                // add to the baseline latter on
                query_var.unique = chr_id;
                mismatch_query_by_mode_by_thread[thread_index][0]->push_back(query_var);
            }
        }
        return false;
    }
    if (separate_var_list[0].size() == 1 && separate_var_list[1].size() == 1){
        // try direct match to save time
        if(separate_var_list[0][0] == separate_var_list[1][0]){

            DiploidVariant tv = separate_var_list[0][0];
            string match_record = chrname_by_chrid[chr_id] + "\t" + to_string(tv.pos+1) + "\t" + tv.ref + "\t" + tv.alts[0];
            if(tv.multi_alts) match_record += "/" + tv.alts[1];
            match_record += "\t.\t.\t.\t.\t.\n";
            // here we need to push back for all mode_index
            //complex_match_records[thread_index]->push_back(match_record);

            int edit_distance = CalculateEditDistance(tv, 0, 0);

            if(threshold_index == 0){
                for(auto baseline_var : separate_var_list[0]){
                    int baseline_unique = baseline_var.unique;
                    auto vi = baseline_unique_indicator[baseline_unique];
                    ref_variant_by_chrid[vi.chr_id]->at(vi.var_id).matching_condition <<= 1;
                    ref_variant_by_chrid[vi.chr_id]->at(vi.var_id).matching_condition |= 1;
                    ref_variant_by_chrid[vi.chr_id]->at(vi.var_id).accumulated_vcf_num += 1;
                }
            }

            for(int mi = 0; mi < mode_index_list.size(); mi ++){
                int mode_i = mode_index_list[mi];
                //if(mi == 0){
                

                // this line should be recovered
                match_records_by_mode_by_thread[thread_index][mode_i]->push_back(match_record);
                
                //}else{
                //    match_records_by_mode_by_thread[thread_index][mode_i]->push_back("$"+to_string(match_records_by_mode_by_thread[thread_index][0]->size()));
                    // use dollor to represent that it is the same
                //}
                baseline_total_match_num[thread_index][threshold_index]->at(mode_i)++;
                query_total_match_num[thread_index][threshold_index]->at(mode_i)++;

                baseline_total_edit_distance[thread_index][threshold_index]->at(mode_i) +=  edit_distance;
                query_total_edit_distance[thread_index][threshold_index]->at(mode_i) += edit_distance;
                //calculate the edit distance
            }
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

            if(max(separate_var_list[flag][0].mdl, separate_var_list[flag][0].mil) > max(total_r_mdl, total_r_mil)){ 
                // no match, need to update the reference
                if(threshold_index == 0){
                    for(auto baseline_var : separate_var_list[0]){
                        int baseline_unique = baseline_var.unique;
                        auto vi = baseline_unique_indicator[baseline_unique];
                        ref_variant_by_chrid[vi.chr_id]->at(vi.var_id).matching_condition <<= 1;
                        ref_variant_by_chrid[vi.chr_id]->at(vi.var_id).accumulated_vcf_num += 1;
                    }

                    for(auto query_var: separate_var_list[1]){
                        // add to the baseline latter on
                        query_var.unique = chr_id;
                        mismatch_query_by_mode_by_thread[thread_index][0]->push_back(query_var);
                    }
                }
                return false;
            }
        }
    }

    separate_var_list[0].clear();
    separate_var_list[1].clear();
    // remove singular variant
    // [todo] try removing this filter to see running time changes
    vector<bool> appliable_flag;
    int total_change = total_mil+total_mdl;

    if(variant_list.size() > EASY_MATCH_VAR_NUM){
        for(int k = 0; k < variant_list.size(); k++){
            DiploidVariant cur_var = variant_list[k];
            int max_change = max(cur_var.mil, cur_var.mdl);
            if(max_change > total_change-max_change){
                appliable_flag.push_back(false);
                //dout << "this variant is removed" << endl;
            }else{
                appliable_flag.push_back(true);
            }
        }
    }else{
        for(int k = 0; k < variant_list.size(); k++){
            appliable_flag.push_back(true);
        }
    }

    string subsequence = genome_sequences[chr_id].substr(min_pos, max_pos - min_pos);

    ToUpper(subsequence); // subsequence only contains upper char
    int offset = min_pos;
    int subsequence_length = max_pos - min_pos;

    // have subsequence in hand
    //generate decision point
    multimap<int, int> * choices_by_pos[2];
    // choice by pos is to also equal to var by pos
    for(int i = 0; i < 2; i++){
        choices_by_pos[i] = new multimap<int, int>();
    }

    for(int index = 0; index < variant_list.size(); index++){
        if(!appliable_flag[index]) continue;
        // remove decision point if not applicable
        int pos = variant_list[index].pos - offset;
        int flag = 0;
        if(variant_list[index].flag) flag = 1;
        choices_by_pos[flag]->insert(pair<int, int>(pos, index));
        //dout << pos << index << endl;
    }

    vector<Interval> mergered_intervals = merge(intervals);
//    unordered_map<int, bool> sync_points;
//    for(int i = 0; i < mergered_intervals.size(); i++){
//        sync_points[mergered_intervals[i].end-offset] = true;
//    }
    vector<int> sync_points;
    for(int i = 0; i < mergered_intervals.size(); i++){
        sync_points.push_back(mergered_intervals[i].end-offset);
    }

    if(sync_points.back() < subsequence.size() - 1){
        sync_points.push_back(subsequence.size()-1);
    }

    int score_unit;
    int match_mode;
    int score_scheme;

    for(int i = 0; i < score_unit_list.size(); i++){
        score_unit = score_unit_list[i];
        for(int j = 0; j < match_mode_list.size(); j++){
            match_mode = match_mode_list[j];
            for(int k = 0; k < score_scheme_list.size(); k++){
                score_scheme = score_scheme_list[k];

                // change only when threshold_index == 0 && score_unit == 0 && match_mode == 0 && score_scheme == 0;
                bool method2 = MatchingSingleClusterBaseExtending(
                                            cluster_id,
                                            thread_index,
                                            variant_list,
                                            subsequence,
                                            offset,
                                            choices_by_pos,
                                            sync_points,
                                            chr_id,
                                            score_unit,
                                            match_mode,
                                            score_scheme,
                                            threshold_index,
                                            baseline_unique_indicator);
            }
        }
    }

    for(int i = 0; i < 2; i++){
        delete choices_by_pos[i];
    }
    //delete choices_by_pos;

    return true;
}

// transfer indicator to variant 
bool WholeGenome::ClusteringMatchInThread(int start, int end, int thread_index) {

	for (int cluster_id = start; cluster_id < end; cluster_id++) {
        if(cluster_id >= variants_by_cluster.size()) break;
        //dout << cluster_id << endl;
        //bool method1 = MatchingSingleCluster(cluster_id, thread_index);
        vector<VariantIndicator> vi_list = variants_by_cluster[cluster_id];
        if(vi_list.size() <= 1) {
            // probabolity need to update reference
            if(vi_list.size() == 0) continue;
            auto vi = vi_list[0];
            if(vi.refer){
                ref_variant_by_chrid[vi.chr_id]->at(vi.var_id).matching_condition <<= 1;
                ref_variant_by_chrid[vi.chr_id]->at(vi.var_id).accumulated_vcf_num += 1;
            }else{
                // add to baseline later on
                DiploidVariant query_var = que_variant_by_chrid[vi.chr_id]->at(vi.var_id);
                query_var.unique = vi.chr_id;
                mismatch_query_by_mode_by_thread[thread_index][0]->push_back(query_var);
            }
            continue;
        }
        // create variant_list from vi_list;

        

        for(int t = 0; t < threshold_num; t++){

            double quality_threshold = threshold_list[t];
            unordered_map<int, VariantIndicator> baseline_unique_indicator;

            vector<DiploidVariant> variant_list;
            int chr_id = -1;
            for(int i = 0; i < vi_list.size(); i++){
                VariantIndicator vi = vi_list[i];
                chr_id = vi.chr_id;
                int var_id = vi.var_id;
                DiploidVariant var;
                if(vi.refer){
                    var = ref_variant_by_chrid[chr_id]->at(var_id);
                    if(t == 0){ // only fill it when t = 0, leave it clear otherwise
                        int unique = var.unique;
                        baseline_unique_indicator[unique] = vi;
                    }
                }else{
                    var = que_variant_by_chrid[chr_id]->at(var_id);
                }
                if(var.qual < quality_threshold)
                    continue;
                variant_list.push_back(var);
            }
            if(chr_id == -1 || chr_id >= chrom_num){
                cout << "[VarMatch] Error in matching single cluster" << endl;
                continue;
            }

            // we do not consider quality threshold in current stage
            // so only when t = 0 do we keep our matching record

            MatchVariantListInThread(thread_index, 
                                    t,
                                    chr_id,
                                    variant_list,
                                    baseline_unique_indicator,
                                    cluster_id);

        }

        //if(method1 != method2){
        //    cout << "not same result for cluster :" << cluster_id << ": " << method1 << "," << method2 << endl;
        //}

	}
	return true;
}


// to reduce memory usage of paths, move all functions about SequencePath out into WholeGenome with a parameter SequencePath
int WholeGenome::PathNeedDecision(SequencePath& sp, multimap<int, int> * choices_by_pos[], int pos){
    for(int i = 0; i < 2; i++){
        if(choices_by_pos[i]->find(pos) != choices_by_pos[i]->end()){
            pair<multimap<int, int>::iterator, multimap<int, int>::iterator> var_range;
            var_range = choices_by_pos[i]->equal_range(pos);

            for(auto it = var_range.first; it != var_range.second; ++it){
                int var_index = (*it).second;
                if(sp.choice_vector[var_index] <= MEANING_CHOICE_BOUND) return var_index;
            }
            // you need to make choices now
            // if(sp.choice_made[i].find(pos) == sp.choice_made[i].end()){
            //     // no choice made at current pos
            //     return true;
            // }
        }
    }
    return -1;
}

// if match_mode == 1, i.e. variant match mode, only check one sequence
// otherwise, check two sequences
int WholeGenome::CheckPathEqualProperty(SequencePath & sp, int match_mode)
{

    if(match_mode == 0){
        //bool equal_sequences = false;
        // same ref position, same donor length, same donor sequence, keep
        if(sp.donor_sequences[0].length() == sp.donor_sequences[2].length() &&
           sp.donor_sequences[1].length() == sp.donor_sequences[3].length()){
            if(sp.donor_sequences[0] == sp.donor_sequences[2] && sp.donor_sequences[1] == sp.donor_sequences[3]){
                sp.same_donor_len = true;
                sp.current_equal_donor_pos[0] = sp.donor_sequences[0].length()-1;
                sp.current_equal_donor_pos[1] = sp.donor_sequences[1].length()-1;
                return 0;
            }else{
                //dout << "delete this path at pos: " << sp.current_genome_pos << " for not equal donor sequence";
                //PrintPath(sp);
                return -1;
            }
        }else{
            sp.same_donor_len = false;
            int min_donor_identical_len[2];
            for(int i = 0; i < 2; i++){
                // compare each strain
                min_donor_identical_len[i] = min(sp.donor_sequences[0+i].length(), sp.donor_sequences[2+i].length());
                for(int k = sp.current_equal_donor_pos[i]+1; k < min_donor_identical_len[i]; k++){
                    if(sp.donor_sequences[0+i][k] != sp.donor_sequences[2+i][k]){
                        return -1;
                    }
                }
                sp.current_equal_donor_pos[i] = min_donor_identical_len[i]-1;
            }
            return 0;
        }
    }else{
        if(sp.donor_sequences[0].length() == sp.donor_sequences[2].length()){
            if(sp.donor_sequences[0] == sp.donor_sequences[2]){
                sp.same_donor_len = true;
                sp.current_equal_donor_pos[0] = sp.donor_sequences[0].length()-1;
                //sp.current_equal_donor_pos[1] = sp.donor_sequences[1].length()-1;
                return 0;
            }else{
                //dout << "delete this path at pos: " << sp.current_genome_pos << " for not equal donor sequence";
                //PrintPath(sp);
                return -1;
            }
        }else{
            sp.same_donor_len = false;
            int min_donor_identical_len[2];
            //for(int i = 0; i < 2; i++)
            int i = 0;
            {
                // compare each strain
                min_donor_identical_len[i] = min(sp.donor_sequences[0+i].length(), sp.donor_sequences[2+i].length());
                for(int k = sp.current_equal_donor_pos[i]+1; k < min_donor_identical_len[i]; k++){
                    if(sp.donor_sequences[0+i][k] != sp.donor_sequences[2+i][k]){
                        return -1;
                    }
                }
                sp.current_equal_donor_pos[i] = min_donor_identical_len[i]-1;
            }
            return 0;
        }
    }
}

// one step is not one nt, but to the next sync point
// i.e. one step, one sync point
int WholeGenome::PathExtendOneStep(SequencePath& sp,
                                   multimap<int, int> * choices_by_pos[],
                                   const string & reference_sequence,
                                   vector<int> & sync_points,
                                   int match_mode,
                                   int & variant_need_decision){
    //-1 operation fail, path deleted
    //0 operation succeed
    //1 operation fail, need to make decision first, then extend
    //2 path reached end, need to check if good

    if(sp.reached_sync_num >= sync_points.size()) return -1;

    int start_pos = sp.current_genome_pos + 1;
    int end_pos = sync_points[sp.reached_sync_num]; // the next sync point, end pos included

    for(int next_genome_pos = start_pos; next_genome_pos <= end_pos; next_genome_pos++){

        // before make decision, we need to check if the equal property still holds
        int variant_need_decision_ = PathNeedDecision(sp, choices_by_pos, next_genome_pos);
        if(variant_need_decision_ >= 0){

            // check equal property
            int statu = CheckPathEqualProperty(sp, match_mode);
            if(statu == -1) return -1;
            variant_need_decision = variant_need_decision_;
            return 1; // need decision on next position
        }

        // else extend one nt
        for(int i = 0; i < 4; i++){

            if(match_mode == 1){
                if(i%2 != 0) continue;
            }

            if(sp.string_sequences[i][next_genome_pos] == "."){
                sp.donor_sequences[i] += reference_sequence[next_genome_pos];
            }else{
                sp.donor_sequences[i] += sp.string_sequences[i][next_genome_pos];
            }
        }
        sp.current_genome_pos = next_genome_pos;
    }

    // reaches the end of end_pos
    sp.reached_sync_num ++;

    if(sp.reached_sync_num >= sync_points.size()){
        // last sync point is the end of ref genome sequence
        if(sp.donor_sequences[0] == sp.donor_sequences[2] &&
           sp.donor_sequences[1] == sp.donor_sequences[3]){
            return 2;
       }else{
            //dout << "delete this path at pos: " << sp.current_genome_pos << " for reach end but not equal";
            //PrintPath(sp);
            return -1;
       }
    }
    return CheckPathEqualProperty(sp, match_mode);
    // first try to converge, then extend

}

int WholeGenome::CalculateScore(DiploidVariant & dv,
                                int choice,
                                int score_unit,
                                int match_mode,
                                int score_scheme){
    int score = 0;
    if(choice <= NOT_USE) return score;
    if(score_unit == 0){
        score = 1;
    }else if(score_unit == 1){
        if(match_mode == 0){
            if(choice == -1){
                score += ScoreEditDistance(dv, 0);
            }else if(choice == -2){
                score += ScoreEditDistance(dv, 1);
            }else if(choice == 0){
                score += ScoreEditDistance(dv, 0);
                if(dv.multi_alts && !dv.zero_one_var){
                    score += ScoreEditDistance(dv, 1);
                }
            }else{
                score += ScoreEditDistance(dv, 0);
                score += ScoreEditDistance(dv, 1);
            }
        }else{
            score += ScoreEditDistance(dv, choice);
        }
    }

    if(score_scheme == 0){
        return score;
    }else if(score_scheme == 1 || score_scheme == 2){
        if(dv.flag == false && score_scheme == 1){
            return score;
        }else if(dv.flag && score_scheme == 2){
            return score;
        }else{
            return 0;
        }
    }
}


// this is the special function to calculate edit distance
int WholeGenome::CalculateEditDistance(DiploidVariant & dv,
                                int choice,
                                int match_mode){
    int score = 0;
    if(choice <= NOT_USE) return score;

    if(match_mode == 0){
        if(choice == -1){
            score += ScoreEditDistance(dv, 0);
        }else if(choice == -2){
            score += ScoreEditDistance(dv, 1);
        }else if(choice == 0){
            score += ScoreEditDistance(dv, 0);
            if(dv.multi_alts && !dv.zero_one_var){
                score += ScoreEditDistance(dv, 1);
            }
        }else{
            score += ScoreEditDistance(dv, 0);
            score += ScoreEditDistance(dv, 1);
        }
    }else{
        score += 2 * ScoreEditDistance(dv, choice);
    }

    return score;
}


// function no longer used, move to VariantMakeDecisionNoGenotype
// no genotype means you can maintain only one strand
// for simplicity, also work on original SequencePath data structure
// when making decision, only decide one path
// when extending, only extend one path
// when comparing, only compare one path
bool WholeGenome::PathMakeDecisionNoGenotype(SequencePath& sp,
                                 vector<DiploidVariant> & variant_list,
                                 multimap<int, int> * choices_by_pos[],
                                 list<SequencePath> & sequence_path_list,
                                 const string & reference_sequence,
                                 int score_unit,
                                 int match_mode,
                                 int score_scheme)
{
    int pos = sp.current_genome_pos+1;
    vector<pair<int, int>> candidate_choices[2];
    for(int i = 0; i < 2; i++){
        // because if it's (-1,-1), it will do nothing, so it's ok to have this one...
        candidate_choices[i].push_back(pair<int, int>(-1, -1));
        // to maintain existance
        // in this position, make choice of not use any variants, no matter if there is variant

        pair<multimap<int, int>::iterator, multimap<int, int>::iterator> var_range;
        var_range = choices_by_pos[i]->equal_range(pos);

        for(auto it = var_range.first; it != var_range.second; ++it){
            int var_index = (*it).second;
            DiploidVariant var = variant_list[var_index];
            // check if current var influence
            string ref = var.ref; //even we do not know the offset, we know ref start from pos of reference_sequence
            string alts[2];
            alts[0] = var.alts[0];
            alts[1] = alts[0];
            if(var.multi_alts){ //here do not have to change anything
                alts[1] = var.alts[1];
            }

            // not just purely consider if a vqriant can be applied, but if a choice
            bool choice_applicable = true;
            for(int k = 0; k < ref.length(); k++){
            // for each ref char
                int y = 0;
                // for each strain
                if(sp.string_sequences[i*2+y][k+pos] != "."){
                    // decision in this area has already been made
                    if(k >= alts[y].length()){
                        choice_applicable = false;
                        break;
                    }else{
                        if(ref[k] != alts[y][k]){
                            choice_applicable = false;
                            break;
                        }
                    }
                }
            }

            if(choice_applicable){
                candidate_choices[i].push_back(pair<int, int>(var_index, 0));
            }

            if(var.multi_alts){ // here do not have to change anything

                //if heterozygous, then there is another choice, check if it is applicable
                string temp = alts[0];
                alts[0] = alts[1];
                alts[1] = temp;

                choice_applicable = true;
                for(int k = 0; k < ref.length(); k++){
                // for each ref char
                    //for(int y = 0; y < 2; y++)
                    int y = 0;
                    // for each strain
                    if(sp.string_sequences[i*2+y][k+pos] != "."){
                        // decision in this area has already been made
                        if(k >= alts[y].length()){
                            // should be a deletion
                            choice_applicable = false;
                            break;
                        }else{
                            // should be equal at current position
                            // can be an insertion, as long as current position is the same
                            if(ref[k] != alts[y][k]){
                                choice_applicable = false;
                                break;
                            }
                        }
                    }
                }

                if(choice_applicable){
                    candidate_choices[i].push_back(pair<int, int>(var_index, 1));
                }
            }
        }
    }

    //dout << candidate_choices[0].size() << "," << candidate_choices[1].size() << endl;

    for(int i = 0; i < candidate_choices[0].size(); i++){
        for(int j = 0; j < candidate_choices[1].size(); j++){
            // iterate all choices
            SequencePath path = sp;
            pair<int, int> var_choice[2];
            var_choice[0] = candidate_choices[0][i];
            var_choice[1] = candidate_choices[1][j];
            for(int x = 0; x < 2; x++){
                // iterate truth and predict
                int var_index = var_choice[x].first;
                if(var_index != -1){
                    DiploidVariant var = variant_list[var_index];
                    // if(var.flag != x){
                    //     dout << "Error" << endl;
                    // }
                    string ref = var.ref;
                    string alts[2];
                    int c = var_choice[x].second;
                    alts[0] = var.alts[c];
                    path.score += CalculateScore(var,
                                                 c,
                                                 score_unit,
                                                 match_mode,
                                                 score_scheme);

                    ToUpper(ref);
                    ToUpper(alts[0]);
                    int y = 0;

                    int k = 0;
                    for(; k < ref.length()-1; k++){
                        if(k < alts[y].length()){
                            if(ref[k] != alts[y][k]){
                                path.string_sequences[x*2+y][pos+k] = alts[y].substr(k,1);
                            }
                            // else change nothing
                        }else{
                            path.string_sequences[x*2+y][pos+k] = "";
                        }
                    }
                    // hence k == ref.length()-1, the last position
                    if(k < alts[y].length()){
                        string alt_part = alts[y].substr(k, alts[y].length()-k);
                        if(alt_part.length() > 1){
                            if(alt_part[0] == ref[k]){
                                if(path.string_sequences[x*2+y][pos+k] == "."){
                                    path.string_sequences[x*2+y][pos+k] = alt_part;
                                }else{
                                    path.string_sequences[x*2+y][pos+k] += alt_part.substr(1, alt_part.size() - 1);
                                }
                            }else{
                                path.string_sequences[x*2+y][pos+k] = alt_part;
                            }
                        }else{
                            if(ref[k] != alts[y][k]){
                                path.string_sequences[x*2+y][pos+k] = alt_part;
                            }
                        }
                    }else{
                        path.string_sequences[x*2+y][pos+k] = "";
                    }

                }
                path.choice_made[x][pos] = var_choice[x];
            }
            sequence_path_list.push_back(path);
        }
    }

    //expected number of inserted paths are 2,3,4,6,x...
    return true;
}

bool WholeGenome::AppendChangedSp(SequencePath& sp,
                         vector<DiploidVariant> & variant_list,
                         list<SequencePath> & sequence_path_list,
                         const string & reference_sequence,
                         int score_unit,
                         int match_mode,
                         int score_scheme,
                         int variant_index,
                         int c)
{
    int pos = sp.current_genome_pos+1;

    SequencePath path = sp;

    if(c == NOT_USE){
        path.choice_vector[variant_index] = c;
        sequence_path_list.push_back(path);
        return true;
    }

    pair<int, int> var_choice[2];
    int x = 0;
    int var_index = variant_index;
    DiploidVariant var = variant_list[var_index];
    if(var.flag) x = 1;
    string ref = var.ref;
    string alts[2];

    if(c == -1){
        alts[0] = ref;
        alts[1] = var.alts[0];
    }else if(c == -2){
        alts[0] = ref;
        alts[1] = var.alts[1];
    }else if(c >= 0){
        // c == 0 or 1
        alts[0] = var.alts[c];
        alts[1] = alts[0];

        if(var.multi_alts && !var.zero_one_var){
            // choose 1 or 0
            alts[1] = var.alts[1- c];
        }else{
            // c is 0, choose 0 or -1
            if(var.heterozygous) alts[1] = ref;
        }
    }else{
        dout << "Unrecognized choice" << endl;
    }
    path.score += CalculateScore(var,
                                 c,
                                 score_unit,
                                 match_mode,
                                 score_scheme);
    ToUpper(ref);
    ToUpper(alts[0]);
    ToUpper(alts[1]);
    for(int y = 0; y < 2; y++){
        // iterate two alts
        string alt = alts[y];
        if(alt == ref) continue;
        vector<string> alt_vector;
        GenerateAltVector(ref, alt, alt_vector);

        int k = 0;
        for(; k < ref.length()-1; k++){

            if(alt_vector[k].size() != 1 || ref[k] != alt_vector[k][0]){
                path.string_sequences[x*2+y][pos+k] = alt_vector[k];
            }
            // else changes nothing

        }
        // hence k == ref.length()-1, the last position
        assert(k == ref.length()-1);
        string alt_part = alt_vector[k];
        if(alt_part.length() > 0){
            if(alt_part.length() > 1){
                if(alt_part[0] == ref[k]){
                    if(path.string_sequences[x*2+y][pos+k] == "."){
                        path.string_sequences[x*2+y][pos+k] = alt_part;
                    }else{
                        path.string_sequences[x*2+y][pos+k] += alt_part.substr(1, alt_part.size() - 1);
                    }
                }else{
                    path.string_sequences[x*2+y][pos+k] = alt_part;
                }
            }else{
                if(ref[k] != alt_vector[k][0]){
                    path.string_sequences[x*2+y][pos+k] = alt_part;
                }
            }
        }else{
            path.string_sequences[x*2+y][pos+k] = "";
        }
    }

    // choice made
    path.choice_vector[variant_index] = c;
    //dout << "after decision at variant " << variant_index << endl;
    //PrintPath(path);
    sequence_path_list.push_back(path);
    return true;
}


// Question: when you make decision, do you also need to align?
// Answer: No, as it makes no difference, so currently you can skip alignment
bool WholeGenome::VariantMakeDecision(SequencePath& sp,
                         vector<DiploidVariant> & variant_list,
                         list<SequencePath> & sequence_path_list,
                         const string & reference_sequence,
                         int score_unit,
                         int match_mode,
                         int score_scheme,
                         int variant_index)
{

    int pos = sp.current_genome_pos+1;

    int var_index = variant_index;
    DiploidVariant var = variant_list[var_index];

    // also this variant may not be used
    AppendChangedSp(sp,
                variant_list,
                sequence_path_list,
                reference_sequence,
                score_unit,
                match_mode,
                score_scheme,
                var_index,
                NOT_USE);

    int i = 0;
    if(var.flag) i = 1;
    //PrintVariant(var);

    // check if current var influence
    string ref = var.ref; //even we do not know the offset, we know ref start from pos of reference_sequence
    string alts[2];
    alts[0] = var.alts[0];
    alts[1] = alts[0];
    if(var.multi_alts && !var.zero_one_var){
        alts[1] = var.alts[1];
    }else if(var.heterozygous){
        alts[1] = ref;
    }

    // not just purely consider if a vqriant can be applied, but if a choice
    int skiped_y = -1;
    if(alts[1] == ref) skiped_y = 1;

    bool choice_applicable = true;
    for(int k = 0; k < ref.length(); k++){
    // for each ref char
        for(int y = 0; y < 2; y++){
            // for each strain
            if(y == skiped_y) continue;
            if(sp.string_sequences[i*2+y][k+pos] != "."){
                // decision in this area has already been made
                if(k >= alts[y].length()){
                    choice_applicable = false;
                    break;
                }else{
                    if(ref[k] != alts[y][k]){
                        choice_applicable = false;
                        break;
                    }
                }
            }
        }
        if(!choice_applicable) break;
    }

    if(choice_applicable){
        //candidate_choices[i].push_back(pair<int, int>(var_index, 0));
        AppendChangedSp(sp,
                        variant_list,
                        sequence_path_list,
                        reference_sequence,
                        score_unit,
                        match_mode,
                        score_scheme,
                        var_index,
                        0);
    }

    if(var.heterozygous){

        //if heterozygous, then there is another choice, check if it is applicable

        string temp = alts[0];
        alts[0] = alts[1];
        alts[1] = temp;

        skiped_y = -1;
        if(alts[0] == ref) skiped_y = 0;

        choice_applicable = true;
        for(int k = 0; k < ref.length(); k++){
        // for each ref char
            for(int y = 0; y < 2; y++){
                // for each strain
                if(skiped_y == y) continue;
                if(sp.string_sequences[i*2+y][k+pos] != "."){
                    // decision in this area has already been made
                    if(k >= alts[y].length()){
                        // should be a deletion
                        choice_applicable = false;
                        break;
                    }else{
                        // should be equal at current position
                        // can be an insertion, as long as current position is the same
                        if(ref[k] != alts[y][k]){
                            choice_applicable = false;
                            break;
                        }
                    }
                }
            }
            if(!choice_applicable) break;
        }

        if(choice_applicable){
            if(var.multi_alts && !var.zero_one_var){
                //candidate_choices[i].push_back(pair<int, int>(var_index, 1));
                AppendChangedSp(sp,
                        variant_list,
                        sequence_path_list,
                        reference_sequence,
                        score_unit,
                        match_mode,
                        score_scheme,
                        var_index,
                        1);
            }else{
                //candidate_choices[i].push_back(pair<int, int>(var_index, -1));
                AppendChangedSp(sp,
                        variant_list,
                        sequence_path_list,
                        reference_sequence,
                        score_unit,
                        match_mode,
                        score_scheme,
                        var_index,
                        -1);
            }
        }
    }

    if(var.multi_alts && var.zero_one_var){
        // here contains another two combinations  alt1/ref and ref/alt1
        alts[0] = var.alts[1];
        alts[1] = ref;

        choice_applicable = true;
        int y = 0;
        for(int k = 0; k < ref.length(); k++){
        // for each ref char
            {
                // for each strain
                if(sp.string_sequences[i*2+y][k+pos] != "."){
                    // decision in this area has already been made
                    if(k >= alts[y].length()){
                        // should be a deletion
                        choice_applicable = false;
                        break;
                    }else{
                        // should be equal at current position
                        // can be an insertion, as long as current position is the same
                        if(ref[k] != alts[y][k]){
                            choice_applicable = false;
                            break;
                        }
                    }
                }
            }
            if(!choice_applicable) break;
        }

        if(choice_applicable){
            //candidate_choices[i].push_back(pair<int, int>(var_index, 1));
            AppendChangedSp(sp,
                variant_list,
                sequence_path_list,
                reference_sequence,
                score_unit,
                match_mode,
                score_scheme,
                var_index,
                1);
        }


        alts[0] = ref;
        alts[1] = var.alts[1];

        choice_applicable = true;
        y = 1;
        for(int k = 0; k < ref.length(); k++){
        // for each ref char
            {
                // for each strain
                if(sp.string_sequences[i*2+y][k+pos] != "."){
                    // decision in this area has already been made
                    if(k >= alts[y].length()){
                        // should be a deletion
                        choice_applicable = false;
                        break;
                    }else{
                        // should be equal at current position
                        // can be an insertion, as long as current position is the same
                        if(ref[k] != alts[y][k]){
                            choice_applicable = false;
                            break;
                        }
                    }
                }
            }
            if(!choice_applicable) break;
        }

        if(choice_applicable){
            //candidate_choices[i].push_back(pair<int, int>(var_index, -2));
            AppendChangedSp(sp,
                variant_list,
                sequence_path_list,
                reference_sequence,
                score_unit,
                match_mode,
                score_scheme,
                var_index,
                -2);
        }
    }
        
}

// no genotype means you only need to maintain one strand
// for simplicity, also work on original SequencePath data structure
// when making decision, only decide one path
// when extending, only extend one path
// when comparing, only compare one path
bool WholeGenome::VariantMakeDecisionNoGenotype(SequencePath& sp,
                         vector<DiploidVariant> & variant_list,
                         list<SequencePath> & sequence_path_list,
                         const string & reference_sequence,
                         int score_unit,
                         int match_mode,
                         int score_scheme,
                         int variant_index)
{

    int pos = sp.current_genome_pos+1;

    int var_index = variant_index;
    DiploidVariant var = variant_list[var_index];

    // also this variant may not be used
    AppendChangedSpNoGenotype(sp,
                variant_list,
                sequence_path_list,
                reference_sequence,
                score_unit,
                match_mode,
                score_scheme,
                var_index,
                NOT_USE);

    int i = 0;
    if(var.flag) i = 1;
    //PrintVariant(var);

    // check if current var influence
    string ref = var.ref; //even we do not know the offset, we know ref start from pos of reference_sequence
    string alts[2];
    alts[0] = var.alts[0];
    alts[1] = alts[0];
    if(var.multi_alts && !var.zero_one_var){
        alts[1] = var.alts[1];
    }else if(var.heterozygous){
        alts[1] = ref;
    }

    bool choice_applicable = true;
    for(int k = 0; k < ref.length(); k++){
    // for each ref char
        int y = 0;
        {
            if(sp.string_sequences[i*2+y][k+pos] != "."){
                // decision in this area has already been made
                if(k >= alts[y].length()){
                    choice_applicable = false;
                    break;
                }else{
                    if(ref[k] != alts[y][k]){
                        choice_applicable = false;
                        break;
                    }
                }
            }
        }
        if(!choice_applicable) break;
    }

    if(choice_applicable){
        //candidate_choices[i].push_back(pair<int, int>(var_index, 0));
        AppendChangedSpNoGenotype(sp,
                        variant_list,
                        sequence_path_list,
                        reference_sequence,
                        score_unit,
                        match_mode,
                        score_scheme,
                        var_index,
                        0);
    }

    // if variants is 0/1, then it does not make sense to apply reference, as it is the same as not_use
    // if variants is 0/1 but contains multi alts, then should try another alt
    // if variants is 1/2 , then should try another alt
    // if variants is 1/1 or 2/2 then should not try another alt
    // but here we do not care the phasing
    // so as long as variant has multi_alts, use another alt

    if(var.multi_alts){

        //if it contains multi alts, then there is another choice, check if it is applicable

        string temp = alts[0];
        alts[0] = alts[1];
        alts[1] = temp;

        choice_applicable = true;
        for(int k = 0; k < ref.length(); k++){
        // for each ref char
            int y = 0;
            {
                if(sp.string_sequences[i*2+y][k+pos] != "."){
                    // decision in this area has already been made
                    if(k >= alts[y].length()){
                        // should be a deletion
                        choice_applicable = false;
                        break;
                    }else{
                        // should be equal at current position
                        // can be an insertion, as long as current position is the same
                        if(ref[k] != alts[y][k]){
                            choice_applicable = false;
                            break;
                        }
                    }
                }
            }
            if(!choice_applicable) break;
        }

        if(choice_applicable){
            AppendChangedSp(sp,
                    variant_list,
                    sequence_path_list,
                    reference_sequence,
                    score_unit,
                    match_mode,
                    score_scheme,
                    var_index,
                    1);
        }
    }   
}


bool WholeGenome::AppendChangedSpNoGenotype(SequencePath& sp,
                         vector<DiploidVariant> & variant_list,
                         list<SequencePath> & sequence_path_list,
                         const string & reference_sequence,
                         int score_unit,
                         int match_mode,
                         int score_scheme,
                         int variant_index,
                         int c)
{
    int pos = sp.current_genome_pos+1;
    SequencePath path = sp;
    if(c == NOT_USE){
        path.choice_vector[variant_index] = c;
        sequence_path_list.push_back(path);
        return true;
    }

    pair<int, int> var_choice[2];
    int x = 0;
    int var_index = variant_index;
    DiploidVariant var = variant_list[var_index];
    if(var.flag) x = 1;
    string ref = var.ref;
    string alts[2];

    if(c == 0 || c == 1){
        // c == 0 or 1
        alts[0] = var.alts[c];
    }else{
        dout << "Unrecognized choice" << endl;
    }
    path.score += CalculateScore(var,
                                 c,
                                 score_unit,
                                 match_mode,
                                 score_scheme);
    ToUpper(ref);
    ToUpper(alts[0]);
    int y = 0;

    string alt = alts[y];
    vector<string> alt_vector;
    GenerateAltVector(ref, alt, alt_vector);
    int k = 0;
    for(; k < ref.length()-1; k++){
        if(alt_vector[k].size() != 1 || ref[k] != alt_vector[k][0]){
            path.string_sequences[x*2+y][pos+k] = alt_vector[k];
        }
        // else changes nothing
    }
    // hence k == ref.length()-1, the last position
    assert(k == ref.length()-1);
    string alt_part = alt_vector[k];
    if(alt_part.length() > 0){
        if(alt_part.length() > 1){
            if(alt_part[0] == ref[k]){
                if(path.string_sequences[x*2+y][pos+k] == "."){
                    path.string_sequences[x*2+y][pos+k] = alt_part;
                }else{
                    path.string_sequences[x*2+y][pos+k] += alt_part.substr(1, alt_part.size() - 1);
                }
            }else{
                path.string_sequences[x*2+y][pos+k] = alt_part;
            }
        }else{
            if(ref[k] != alt_vector[k][0]){
                path.string_sequences[x*2+y][pos+k] = alt_part;
            }
        }
    }else{
        path.string_sequences[x*2+y][pos+k] = "";
    }
    // choice made
    path.choice_vector[variant_index] = c;
    //dout << "after decision at variant " << variant_index << endl;
    //PrintPath(path);
    sequence_path_list.push_back(path);
    return true;
}

// this function is no longer used, because you can not make decison for one position at once,
// there might be multiple variants in one position,
// so a better way to do this is to make decision for one variant at a time
// previously I just want to save some time, but ignore the multiple variant condition
bool WholeGenome::PathMakeDecision(SequencePath& sp,
                                 vector<DiploidVariant> & variant_list,
                                 multimap<int, int> * choices_by_pos[],
                                 list<SequencePath> & sequence_path_list,
                                 const string & reference_sequence,
                                 int score_unit,
                                 int match_mode,
                                 int score_scheme)
{
    int pos = sp.current_genome_pos+1;

    vector<pair<int, int>> candidate_choices[2];
    for(int i = 0; i < 2; i++){

        // because if it's (-1,-1), it will do nothing, so it's ok to have this one...
        candidate_choices[i].push_back(pair<int, int>(-1, -1));
        // in this position, make choice of not use any variants, no matter if there is variant

        pair<multimap<int, int>::iterator, multimap<int, int>::iterator> var_range;
        var_range = choices_by_pos[i]->equal_range(pos);

        for(auto it = var_range.first; it != var_range.second; ++it){
            int var_index = (*it).second;
            DiploidVariant var = variant_list[var_index];
            //PrintVariant(var);

            // check if current var influence
            string ref = var.ref; //even we do not know the offset, we know ref start from pos of reference_sequence
            string alts[2];
            alts[0] = var.alts[0];
            alts[1] = alts[0];
            if(var.multi_alts && !var.zero_one_var){
                alts[1] = var.alts[1];
            }else if(var.heterozygous){
                alts[1] = ref;
            }

            // not just purely consider if a vqriant can be applied, but if a choice
            int skiped_y = -1;
            if(alts[1] == ref) skiped_y = 1;

            bool choice_applicable = true;
            for(int k = 0; k < ref.length(); k++){
            // for each ref char
                for(int y = 0; y < 2; y++){
                    // for each strain
                    if(y == skiped_y) continue;
                    if(sp.string_sequences[i*2+y][k+pos] != "."){
                        // decision in this area has already been made
                        if(k >= alts[y].length()){
                            choice_applicable = false;
                            break;
                        }else{
                            if(ref[k] != alts[y][k]){
                                choice_applicable = false;
                                break;
                            }
                        }
                    }
                }
                if(!choice_applicable) break;
            }

            if(choice_applicable){
                candidate_choices[i].push_back(pair<int, int>(var_index, 0));
            }

            if(var.heterozygous){

                //if heterozygous, then there is another choice, check if it is applicable

                string temp = alts[0];
                alts[0] = alts[1];
                alts[1] = temp;

                skiped_y = -1;
                if(alts[0] == ref) skiped_y = 0;

                choice_applicable = true;
                for(int k = 0; k < ref.length(); k++){
                // for each ref char
                    for(int y = 0; y < 2; y++){
                        // for each strain
                        if(skiped_y == y) continue;
                        if(sp.string_sequences[i*2+y][k+pos] != "."){
                            // decision in this area has already been made
                            if(k >= alts[y].length()){
                                // should be a deletion
                                choice_applicable = false;
                                break;
                            }else{
                                // should be equal at current position
                                // can be an insertion, as long as current position is the same
                                if(ref[k] != alts[y][k]){
                                    choice_applicable = false;
                                    break;
                                }
                            }
                        }
                    }
                    if(!choice_applicable) break;
                }

                if(choice_applicable){
                    if(var.multi_alts && !var.zero_one_var){
                        candidate_choices[i].push_back(pair<int, int>(var_index, 1));
                    }else{
                        candidate_choices[i].push_back(pair<int, int>(var_index, -1));
                    }
                }
            }

            if(var.multi_alts && var.zero_one_var){
                // here contains another two combinations  alt1/ref and ref/alt1
                alts[0] = var.alts[1];
                alts[1] = ref;

                choice_applicable = true;
                int y = 0;
                for(int k = 0; k < ref.length(); k++){
                // for each ref char
                    {
                        // for each strain
                        if(sp.string_sequences[i*2+y][k+pos] != "."){
                            // decision in this area has already been made
                            if(k >= alts[y].length()){
                                // should be a deletion
                                choice_applicable = false;
                                break;
                            }else{
                                // should be equal at current position
                                // can be an insertion, as long as current position is the same
                                if(ref[k] != alts[y][k]){
                                    choice_applicable = false;
                                    break;
                                }
                            }
                        }
                    }
                    if(!choice_applicable) break;
                }

                if(choice_applicable){
                    candidate_choices[i].push_back(pair<int, int>(var_index, 1));
                }


                alts[0] = ref;
                alts[1] = var.alts[1];

                choice_applicable = true;
                y = 1;
                for(int k = 0; k < ref.length(); k++){
                // for each ref char
                    {
                        // for each strain
                        if(sp.string_sequences[i*2+y][k+pos] != "."){
                            // decision in this area has already been made
                            if(k >= alts[y].length()){
                                // should be a deletion
                                choice_applicable = false;
                                break;
                            }else{
                                // should be equal at current position
                                // can be an insertion, as long as current position is the same
                                if(ref[k] != alts[y][k]){
                                    choice_applicable = false;
                                    break;
                                }
                            }
                        }
                    }
                    if(!choice_applicable) break;
                }

                if(choice_applicable){
                    candidate_choices[i].push_back(pair<int, int>(var_index, -2));
                }
            }
        }
    }

    //dout << candidate_choices[0].size() << "," << candidate_choices[1].size() << endl;

    for(int i = 0; i < candidate_choices[0].size(); i++){
        for(int j = 0; j < candidate_choices[1].size(); j++){
            // iterate all choices
            SequencePath path = sp;
            pair<int, int> var_choice[2];
            var_choice[0] = candidate_choices[0][i];
            var_choice[1] = candidate_choices[1][j];
            for(int x = 0; x < 2; x++){
                // iterate truth and predict
                int var_index = var_choice[x].first;
                if(var_index != -1){
//                    string temp_sequence = reference_sequence.substr(pos, 1);
//                    path.string_sequences[x*2][pos] = temp_sequence;
//                    path.string_sequences[x*2+1][pos] = temp_sequence;
//                }else{
                    // set score


                    DiploidVariant var = variant_list[var_index];
                    // if(var.flag != x){
                    //     dout << "Error" << endl;
                    // }
                    string ref = var.ref;
                    string alts[2];

                    int c = var_choice[x].second;
                    if(c == -1){
                        alts[0] = ref;
                        alts[1] = var.alts[0];
                    }else if(c == -2){
                        alts[0] = ref;
                        alts[1] = var.alts[1];
                    }else{
                        // c == 0 or 1
                        alts[0] = var.alts[c];
                        alts[1] = alts[0];

                        if(var.multi_alts && !var.zero_one_var){
                            // choose 1 or 0
                            alts[1] = var.alts[1- c];
                        }else{
                            // c is 0, choose 0 or -1
                            if(var.heterozygous) alts[1] = ref;
                        }
                    }

                    path.score += CalculateScore(var,
                                                 c,
                                                 score_unit,
                                                 match_mode,
                                                 score_scheme);

                    ToUpper(ref);
                    ToUpper(alts[0]);
                    ToUpper(alts[1]);
                    for(int y = 0; y < 2; y++){
                        // iterate two alts
                        string alt = alts[y];
                        vector<string> alt_vector;
                        GenerateAltVector(ref, alt, alt_vector);

                        int k = 0;
                        for(; k < ref.length()-1; k++){

                            if(alt_vector[k].size() != 1 || ref[k] != alt_vector[k][0]){
                                path.string_sequences[x*2+y][pos+k] = alt_vector[k];
                            }
                            // else changes nothing

                        }
                        // hence k == ref.length()-1, the last position
                        assert(k == ref.length()-1);
                        string alt_part = alt_vector[k];
                        if(alt_part.length() > 0){
                            if(alt_part.length() > 1){
                                if(alt_part[0] == ref[k]){
                                    if(path.string_sequences[x*2+y][pos+k] == "."){
                                        path.string_sequences[x*2+y][pos+k] = alt_part;
                                    }else{
                                        path.string_sequences[x*2+y][pos+k] += alt_part.substr(1, alt_part.size() - 1);
                                    }
                                }else{
                                    path.string_sequences[x*2+y][pos+k] = alt_part;
                                }
                            }else{
                                if(ref[k] != alt_vector[k][0]){
                                    path.string_sequences[x*2+y][pos+k] = alt_part;
                                }
                            }
                        }else{
                            path.string_sequences[x*2+y][pos+k] = "";
                        }
                    }
                }
                path.choice_made[x][pos] = var_choice[x];
            }
            // choice made
            dout << "after decision at pos " << pos << endl;
            PrintPath(path);
            sequence_path_list.push_back(path);
        }
    }

    //expected number of inserted paths are 2,3,4,6,x...
    return true;
}

bool WholeGenome::PathMakeDecisionBackup(SequencePath& sp,
                                 vector<DiploidVariant> & variant_list,
                                 multimap<int, int> * choices_by_pos[],
                                 list<SequencePath> & sequence_path_list,
                                 const string & reference_sequence,
                                 int score_unit,
                                 int match_mode,
                                 int score_scheme)
{
    //expected number of inserted paths are 2,3,4,6,x...
    return true;
}

void WholeGenome::PrintPath(SequencePath & sp){
    cout << "- Sequence Path:" << endl;
    cout << "@ String Sequences:" << endl;
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < sp.string_sequences[i].size(); j++){
            cout << sp.string_sequences[i][j] << " ";
        }
        cout << endl;
    }
    cout << "@ Donor Sequences:" << endl;
    for(int i = 0; i < 4; i++){
        cout << sp.donor_sequences[i] << endl;
    }
    cout << "@ Removable: " << sp.removable << endl;

    for(int i = 0; i < sp.choice_vector.size(); i++){
        cout << sp.choice_vector[i] << ",";
    }
    cout << endl;
}

// next: while until current path list is empty
// if extend, add to next path list
// if need decision, make decision, append to current list
// if reach end, compare with best path
bool WholeGenome::MatchingSingleClusterBaseExtending(int cluster_index,
                                                    int thread_index,
                                                    vector<DiploidVariant> & variant_list,
                                                    string & subsequence,
                                                    int offset,
                                                    multimap<int, int> * choices_by_pos[],
                                                    vector<int> & sync_points,
                                                    int chr_id,
                                                    int score_unit,
                                                    int match_mode,
                                                    int score_scheme,
                                                    int threshold_index,
                                                    unordered_map<int, VariantIndicator> & baseline_unique_indicator){
    //--------------for unit test------------------------------
    //dout << variant_list.size() << endl;

    //int chr_id = 0;
    //-------------end unit test-------------------------------

    //cout << chr_id << "," << cluster_index << "," << variant_list.size() << endl;

    // so a legal sync_points vector contains at least two
    // first is the end of variant, there should be at least one variant
    // second is the end of subsequence, there should be at least one nt not influenced by a variant

    list<SequencePath> current_path_list;
    list<SequencePath> next_path_list;
    SequencePath sp(subsequence.length(), variant_list.size());
    SequencePath best_path = sp;
    current_path_list.push_back(sp);
    int considered_path_num = 0;
    int mode_index = GetIndexFromMatchScore(score_unit, match_mode, score_scheme);
    while(current_path_list.size() != 0){
        bool reach_sync_point = true;
        // extend path before reaches sync points
        //cout << "\t" << current_path_list.size() << "," << considered_path_num << endl;
        
        while(current_path_list.size() != 0){
            SequencePath path = current_path_list.front();
            current_path_list.pop_front();
            considered_path_num ++;
            if (considered_path_num >= MAX_SEARCH_SPACE){
                //cout << "direct matching large variant cluster" << endl;
                //cout << variant_list.size() << "," << considered_path_num << endl;
                return MatchingSingleClusterDirectly(cluster_index, 
                    thread_index, 
                    mode_index, 
                    threshold_index, 
                    chr_id,
                    variant_list);
            }
            //dout << path.current_genome_pos << ":" << current_path_list.size() << endl;
            //PrintPath(path);
            int variant_need_decision = -1;
            int is_extend = PathExtendOneStep(path, choices_by_pos, subsequence, sync_points, match_mode, variant_need_decision);
            //cout << variant_need_decision << endl;
            //PrintPath(path);
            if(is_extend == -1){
                // discard path
                continue;
            }
            else if(is_extend == 0){
                next_path_list.push_back(path);
                // here the path is supposed to reach the next sync point
                // so it goes into next path list, and decrease the number of current path list
            }else if(is_extend == 1){
                if(match_mode == 0){
                    // PathMakeDecision(path,
                    //                  variant_list,
                    //                  choices_by_pos,
                    //                  current_path_list,
                    //                  subsequence,
                    //                  score_unit,
                    //                  match_mode,
                    //                  score_scheme);

                    VariantMakeDecision(path,
                                         variant_list,
                                         current_path_list,
                                         subsequence,
                                         score_unit,
                                         match_mode,
                                         score_scheme,
                                         variant_need_decision);
                }else{
                    // PathMakeDecisionNoGenotype(path,
                    //                            variant_list,
                    //                            choices_by_pos,
                    //                            current_path_list,
                    //                            subsequence,
                    //                            score_unit,
                    //                            match_mode,
                    //                            score_scheme);

                    VariantMakeDecisionNoGenotype(path,
                                                 variant_list,
                                                 current_path_list,
                                                 subsequence,
                                                 score_unit,
                                                 match_mode,
                                                 score_scheme,
                                                 variant_need_decision);
                }
            }else if(is_extend == 2){
                if(path.score > best_path.score){
                    best_path = path; // only when you reach the very end can you be considered as best path
                    //PrintPath(best_path);
                }
            }
        }
        current_path_list = next_path_list;
        next_path_list.clear();
        if(current_path_list.size() > 0){
            //int current_genome_pos = current_path_list.front().current_genome_pos;
            // after revise, we do not need this check
            //if(sync_points.find(current_genome_pos) != sync_points.end()){
                //dout << "converge paths at position: " << current_genome_pos << endl;
                //dout << "before converge: " << current_path_list.size() << endl;
                ConvergePaths(current_path_list);
                //dout << "after converge: " << current_path_list.size() << endl;
            //}
        }
    }
    current_path_list.clear();
    next_path_list.clear();
    // print best_path
    if(best_path.score <= 0) {
        // no match exist, need to update reference
        for(auto variant : variant_list){
            if(variant.flag){
                // query, add to baseline later
                variant.unique = chr_id;
                mismatch_query_by_mode_by_thread[thread_index][0]->push_back(variant);
            }else{
                // baseline, no match found
                int baseline_unique = variant.unique;
                auto vi = baseline_unique_indicator[baseline_unique];
                ref_variant_by_chrid[vi.chr_id]->at(vi.var_id).matching_condition <<= 1;
                ref_variant_by_chrid[vi.chr_id]->at(vi.var_id).accumulated_vcf_num += 1;
            }
        }
        return false;
    }
    // if(considered_path_num > max_path_num_by_mode_by_thread[thread_index][mode_index]){
    //     max_path_num_by_mode_by_thread[thread_index][mode_index] = considered_path_num;
    // }
    //dout << "new method: " << best_path.score << endl;

    //PrintPath(best_path);

    //==========================output ======================

    //return true;
    if(match_mode == 0){
        ConstructMatchRecord(best_path,
                             variant_list,
                             subsequence,
                             offset,
                             thread_index,
                             chr_id,
                             mode_index,
                             threshold_index,
                             baseline_unique_indicator);
    }else{
        // there is not need to track baseline usage for no genotype in current stage, because mode index will never be zero
        ConstructMatchRecordNoGenotype(best_path,
                                       variant_list,
                                       subsequence,
                                       offset,
                                       thread_index,
                                       chr_id,
                                       mode_index,
                                       threshold_index,
                                       baseline_unique_indicator);
    }
    return true;
}

int GetMatchmodeFromModeIndex(int mode_index){
    int result = mode_index;
    result >>= 2;
    result &= 1;
    return result;
}

bool WholeGenome::MatchingSingleClusterDirectly(int cluster_index, 
    int thread_index, 
    int mode_index, 
    int threshold_index, 
    int chr_id,
    vector<DiploidVariant> & variant_list){
    //cout << "direct match : " << variant_list.size() << endl;
    map<int, vector<DiploidVariant>> separate_var_list_map[2];
    for (int i = 0; i < variant_list.size(); i++) {
        int flag = 0;
        if (variant_list[i].flag) flag = 1; // flag indicate if the variant is from ref set(0) or query set(1)
        int pos = variant_list[i].pos;

        if(separate_var_list_map[flag].find(pos) == separate_var_list_map[flag].end()){
            separate_var_list_map[flag][pos] = vector<DiploidVariant>({variant_list[i]});
        }else{
            separate_var_list_map[flag][pos].push_back(variant_list[i]);
        }
    }

    //int match_mode = GetMatchmodeFromModeIndex(mode_index);

    for(auto baseline_it = separate_var_list_map[0].begin(); baseline_it != separate_var_list_map[0].end(); ++baseline_it){
        int pos = baseline_it->first;
        if(separate_var_list_map[1].find(pos) != separate_var_list_map[1].end()){
            vector<DiploidVariant> baseline_var_list = baseline_it->second;
            vector<DiploidVariant> query_var_list = separate_var_list_map[1][pos];
            for(auto baseline_var : baseline_var_list){
                for(auto query_var : query_var_list){
                    if(baseline_var == query_var){
                        DiploidVariant tv = baseline_var;
                        string match_record = chrname_by_chrid[chr_id] + "\t" + to_string(tv.pos+1) + "\t" + tv.ref + "\t" + tv.alts[0];
                        if(tv.multi_alts) match_record += "/" + tv.alts[1];
                        match_record += "\t.\t.\t.\t.\t.\n";
                        int edit_distance = CalculateEditDistance(tv, 0, 0);
                        match_records_by_mode_by_thread[thread_index][mode_index]->push_back(match_record);
                        baseline_total_match_num[thread_index][threshold_index]->at(mode_index)++;
                        query_total_match_num[thread_index][threshold_index]->at(mode_index)++;

                        baseline_total_edit_distance[thread_index][threshold_index]->at(mode_index) += edit_distance;
                        query_total_edit_distance[thread_index][threshold_index]->at(mode_index) += edit_distance;
                        break;
                    }
                }
            }
        }
    } 
    return true;
}



void WholeGenome::ConstructMatchRecord(SequencePath & best_path,
                                       vector<DiploidVariant> & variant_list,
                                       string & subsequence,
                                       int offset,
                                       int thread_index,
                                       int chr_id,
                                       int mode_index,
                                       int threshold_index,
                                       unordered_map<int, VariantIndicator> & baseline_unique_indicator){
    int truth_num = 0;
    int predict_num = 0;
    int truth_edit_distance = 0;
    int predict_edit_distance = 0;

    bool need_match_record = false;

    if (threshold_index == 0) need_match_record = true;

    bool multiple_match = false;

    if(best_path.donor_sequences[0] != best_path.donor_sequences[1]) multiple_match = true;

    string parsimonious_ref = subsequence;
    string parsimonious_alt0 = best_path.donor_sequences[0];
    string parsimonious_alt1 = best_path.donor_sequences[1];

    int parsimonious_pos = offset;
//    NormalizeVariantSequence(offset,
//                             parsimonious_ref,
//                             parsimonious_alt0,
//                             parsimonious_alt1,
//                             chr_id);

    string match_record = chrname_by_chrid[chr_id] + "\t" + to_string(parsimonious_pos+1) + "\t" + parsimonious_ref + "\t" + parsimonious_alt0;
    if(multiple_match) match_record += "/" + parsimonious_alt1;

    string vcf_record[2];
    string phasing_record[2];

    for (int i = 0; i < 2; i++) {
        for (int var_index = 0; var_index < variant_list.size(); var_index++) {
            DiploidVariant variant = variant_list[var_index];
            if(variant.flag != i) continue;
            //The exact wording from the C++ standard is (§4.7/4): "If the source type is bool,
            // the value false is converted to zero and the value true is converted to one."
            int phasing = best_path.choice_vector[var_index];
            if(phasing <= NOT_USE) { 

                
                if(threshold_index == 0 && mode_index == 0){
                    if(i == 0){
                        // put a zero there
                        int baseline_unique = variant.unique;
                        auto vi = baseline_unique_indicator[baseline_unique];
                        ref_variant_by_chrid[vi.chr_id]->at(vi.var_id).matching_condition <<= 1;
                        ref_variant_by_chrid[vi.chr_id]->at(vi.var_id).accumulated_vcf_num += 1;

                    }else{
                        // add to baseline later on
                        variant.unique = chr_id;
                        mismatch_query_by_mode_by_thread[thread_index][0]->push_back(variant);
                    }
                }
                continue; //this variant is not used
            }else{
                if(threshold_index == 0 && mode_index == 0){
                    if(i == 0){
                        // put a one there
                        int baseline_unique = variant.unique;
                        auto vi = baseline_unique_indicator[baseline_unique];
                        ref_variant_by_chrid[vi.chr_id]->at(vi.var_id).matching_condition <<= 1;
                        ref_variant_by_chrid[vi.chr_id]->at(vi.var_id).matching_condition |= 1;
                        ref_variant_by_chrid[vi.chr_id]->at(vi.var_id).accumulated_vcf_num += 1;
                    }else{
                        // do nothing
                    }
                }
            }
            int edit_distance = CalculateEditDistance(variant, phasing, 0);
            if(i == 0){
                truth_num ++;
                truth_edit_distance += edit_distance;
            }else{
                predict_num ++;
                predict_edit_distance += edit_distance;
            }

            if(need_match_record){
                string alt_string = variant.alts[0];
                if(variant.multi_alts){
                    alt_string += "/" + variant.alts[1];
                }
                string phasing_string = "";
                if(phasing == 0){
                    phasing_string += "1";
                    if(variant.heterozygous){
                        if(variant.multi_alts && !variant.zero_one_var){
                            phasing_string += "|2";
                        }else{
                            phasing_string += "|0";
                        }
                    }else{
                        phasing_string += "|1";
                    }
                }else if(phasing == 1){
                    if(variant.multi_alts && !variant.zero_one_var){
                        phasing_string += "2|1";
                    }else if(variant.multi_alts && variant.zero_one_var){
                        phasing_string += "2|0";
                    }else{
                        phasing_string += "0|1";
                    }
                }else if(phasing == -1){
                    phasing_string += "0|1";
                }else if(phasing == -2){
                    phasing_string += "0|2";
                }
                string variant_record = to_string(variant.pos+1) + "," + variant.ref + "," + alt_string;
                vcf_record[i] += variant_record;
                phasing_record[i] += phasing_string;
                vcf_record[i] += ";";
                phasing_record[i] += ";";
            }
        }
        if(need_match_record){
            vcf_record[i] = vcf_record[i].substr(0, vcf_record[i].size()-1);
            phasing_record[i] = phasing_record[i].substr(0, phasing_record[i].size()-1);
        }

    }

    if(need_match_record){
        match_record += "\t" + vcf_record[0] + "\t" + vcf_record[1];
        match_record += "\t" + phasing_record[0] + "\t" + phasing_record[1];
        match_record += "\t" + to_string(best_path.score) + "\n";

        //complex_match_records[thread_index]->push_back(match_record);
        

        // this line should be recovered
        match_records_by_mode_by_thread[thread_index][mode_index]->push_back(match_record);
    


    }

    baseline_total_match_num[thread_index][threshold_index]->at(mode_index) += truth_num;
    query_total_match_num[thread_index][threshold_index]->at(mode_index) += predict_num;

    baseline_total_edit_distance[thread_index][threshold_index]->at(mode_index) += truth_edit_distance;
    query_total_edit_distance[thread_index][threshold_index]->at(mode_index) += predict_edit_distance;
}

void WholeGenome::ConstructMatchRecordNoGenotype(SequencePath & best_path,
                                                 vector<DiploidVariant> & variant_list,
                                                 string & subsequence,
                                                 int offset,
                                                 int thread_index,
                                                 int chr_id,
                                                 int mode_index,
                                                 int threshold_index,
                                                 unordered_map<int, VariantIndicator> & baseline_unique_indicator){
    int truth_num = 0;
    int predict_num = 0;

    int truth_edit_distance = 0;
    int predict_edit_distance = 0;

    bool need_match_record = false;
    if(threshold_index == 0) need_match_record = true;

    bool multiple_match = false;
    string parsimonious_ref = subsequence;
    string parsimonious_alt0 = best_path.donor_sequences[0];
    string parsimonious_alt1 = best_path.donor_sequences[0];

    int parsimonious_pos = offset;

//    NormalizeVariantSequence(offset,
//                             parsimonious_ref,
//                             parsimonious_alt0,
//                             parsimonious_alt1,
//                             chr_id);

    string match_record = chrname_by_chrid[chr_id] + "\t" + to_string(parsimonious_pos+1) + "\t" + parsimonious_ref + "\t" + parsimonious_alt0;
    //if(multiple_match) match_record += "/" + parsimonious_alt1;

    string vcf_record[2];
    string phasing_record[2];

    for (int i = 0; i < 2; i++) {
        for (int var_index = 0; var_index < variant_list.size(); var_index++) {
            DiploidVariant variant = variant_list[var_index];
            if(variant.flag != i) continue;
            //The exact wording from the C++ standard is (§4.7/4): "If the source type is bool,
            // the value false is converted to zero and the value true is converted to one."
            int phasing = best_path.choice_vector[var_index];
            if(phasing <= NOT_USE) continue;
            int edit_distance = CalculateEditDistance(variant, phasing, 1);
            if(i == 0){
                truth_num ++;
                truth_edit_distance += edit_distance;
            }else{
                predict_num ++;
                predict_edit_distance += edit_distance;
            }

            if(need_match_record){
                string alt_string = variant.alts[0];
                if(variant.multi_alts){
                    alt_string += "/" + variant.alts[1];
                }
                string phasing_string = "";
                if(phasing == 0){
                    phasing_string += "1|1";
                }else if(phasing == 1){
                    phasing_string += "2|2";
                }
                string variant_record = to_string(variant.pos+1) + "," + variant.ref + "," + alt_string;
                vcf_record[i] += variant_record;
                phasing_record[i] += phasing_string;
                vcf_record[i] += ";";
                phasing_record[i] += ";";
            }
        }
        if(need_match_record){
            vcf_record[i] = vcf_record[i].substr(0, vcf_record[i].size()-1);
            phasing_record[i] = phasing_record[i].substr(0, phasing_record[i].size()-1);
        }

    }

    if(need_match_record){
       match_record += "\t" + vcf_record[0] + "\t" + vcf_record[1];
        match_record += "\t" + phasing_record[0] + "\t" + phasing_record[1];
       match_record += "\t" + to_string(best_path.score) + "\n";
       //complex_match_records[thread_index]->push_back(match_record);
       // this line should be recovered
       match_records_by_mode_by_thread[thread_index][mode_index]->push_back(match_record);
    }

    baseline_total_match_num[thread_index][threshold_index]->at(mode_index) += truth_num;
    query_total_match_num[thread_index][threshold_index]->at(mode_index) += predict_num;

    baseline_total_edit_distance[thread_index][threshold_index]->at(mode_index) += truth_edit_distance;
    query_total_edit_distance[thread_index][threshold_index]->at(mode_index) += predict_edit_distance;
}

// function no longer used, backup old method
void WholeGenome::ConstructMatchRecordBackup(SequencePath & best_path,
                                       vector<DiploidVariant> & variant_list,
                                       string & subsequence,
                                       int offset,
                                       int thread_index,
                                       int chr_id,
                                       int mode_index,
                                       int threshold_index){
    int truth_num = 0;
    int predict_num = 0;

    bool need_match_record = false;

    if (threshold_index == 0) need_match_record = true;

    bool multiple_match = false;

    if(best_path.donor_sequences[0] != best_path.donor_sequences[1]) multiple_match = true;

    string parsimonious_ref = subsequence;
    string parsimonious_alt0 = best_path.donor_sequences[0];
    string parsimonious_alt1 = best_path.donor_sequences[1];

    int parsimonious_pos = offset;
//    NormalizeVariantSequence(offset,
//                             parsimonious_ref,
//                             parsimonious_alt0,
//                             parsimonious_alt1,
//                             chr_id);

    string match_record = chrname_by_chrid[chr_id] + "\t" + to_string(parsimonious_pos+1) + "\t" + parsimonious_ref + "\t" + parsimonious_alt0;
    if(multiple_match) match_record += "/" + parsimonious_alt1;

    string vcf_record[2];
    string phasing_record[2];

	for (int i = 0; i < 2; i++) {
		for (auto it = best_path.choice_made[i].begin(); it != best_path.choice_made[i].end(); ++it) {
            pair<int, int> selection = it->second;
            int phasing = selection.second;
            if(selection.first == -1) continue;
            DiploidVariant variant = variant_list[selection.first];
            if(!variant.flag){
                truth_num++;
            }else{
                predict_num++;
            }

            if(need_match_record){
                string alt_string = variant.alts[0];
                if(variant.multi_alts){
                    alt_string += "/" + variant.alts[1];
                }
                string phasing_string = "";
                if(phasing == 0){
                    phasing_string += "1";
                    if(variant.heterozygous){
                        if(variant.multi_alts && !variant.zero_one_var){
                            phasing_string += "|2";
                        }else{
                            phasing_string += "|0";
                        }
                    }else{
                        phasing_string += "|1";
                    }
                }else if(phasing == 1){
                    if(variant.multi_alts && !variant.zero_one_var){
                        phasing_string += "2|1";
                    }else if(variant.multi_alts && variant.zero_one_var){
                        phasing_string += "2|0";
                    }else{
                        phasing_string += "0|1";
                    }
                }else if(phasing == -1){
                    phasing_string += "0|1";
                }else if(phasing == -2){
                    phasing_string += "0|2";
                }
                string variant_record = to_string(variant.pos+1) + "," + variant.ref + "," + alt_string;
                vcf_record[i] += variant_record;
                phasing_record[i] += phasing_string;
                vcf_record[i] += ";";
                phasing_record[i] += ";";
            }
		}
        if(need_match_record){
            vcf_record[i] = vcf_record[i].substr(0, vcf_record[i].size()-1);
            phasing_record[i] = phasing_record[i].substr(0, phasing_record[i].size()-1);
        }

	}

    if(need_match_record){
    	match_record += "\t" + vcf_record[0] + "\t" + vcf_record[1];
        match_record += "\t" + phasing_record[0] + "\t" + phasing_record[1];
    	match_record += "\t" + to_string(best_path.score) + "\n";

    	//complex_match_records[thread_index]->push_back(match_record);
    	
        // this line should be recovered


        match_records_by_mode_by_thread[thread_index][mode_index]->push_back(match_record);
    

    }

    baseline_total_match_num[thread_index][threshold_index]->at(mode_index) += truth_num;
    query_total_match_num[thread_index][threshold_index]->at(mode_index) += predict_num;
}

void WholeGenome::ConstructMatchRecordNoGenotypeBackup(SequencePath & best_path,
                                                 vector<DiploidVariant> & variant_list,
                                                 string & subsequence,
                                                 int offset,
                                                 int thread_index,
                                                 int chr_id,
                                                 int mode_index,
                                                 int threshold_index){
    int truth_num = 0;
    int predict_num = 0;

    bool need_match_record = false;
    if(threshold_index == 0) need_match_record = true;

    bool multiple_match = false;
    string parsimonious_ref = subsequence;
    string parsimonious_alt0 = best_path.donor_sequences[0];
    string parsimonious_alt1 = best_path.donor_sequences[0];

    int parsimonious_pos = offset;

//    NormalizeVariantSequence(offset,
//                             parsimonious_ref,
//                             parsimonious_alt0,
//                             parsimonious_alt1,
//                             chr_id);

    string match_record = chrname_by_chrid[chr_id] + "\t" + to_string(parsimonious_pos+1) + "\t" + parsimonious_ref + "\t" + parsimonious_alt0;
    //if(multiple_match) match_record += "/" + parsimonious_alt1;

    string vcf_record[2];
    string phasing_record[2];

	for (int i = 0; i < 2; i++) {
		for (auto it = best_path.choice_made[i].begin(); it != best_path.choice_made[i].end(); ++it) {
            pair<int, int> selection = it->second;
            int phasing = selection.second;
            if(selection.first == -1) continue;
            if (phasing == -1) continue;
            DiploidVariant variant = variant_list[selection.first];
            if(!variant.flag){
                truth_num++;
            }else{
                predict_num++;
            }

            if(need_match_record){
                string alt_string = variant.alts[0];
                if(variant.multi_alts){
                    alt_string += "/" + variant.alts[1];
                }
                string phasing_string = "";
                if(phasing == 0){
                    phasing_string += "1|1";
                }else if(phasing == 1){
                    phasing_string += "2|2";
                }
                string variant_record = to_string(variant.pos+1) + "," + variant.ref + "," + alt_string;
                vcf_record[i] += variant_record;
                phasing_record[i] += phasing_string;
                vcf_record[i] += ";";
                phasing_record[i] += ";";
            }
		}
        if(need_match_record){
            vcf_record[i] = vcf_record[i].substr(0, vcf_record[i].size()-1);
            phasing_record[i] = phasing_record[i].substr(0, phasing_record[i].size()-1);
        }

	}

    if(need_match_record){
	   match_record += "\t" + vcf_record[0] + "\t" + vcf_record[1];
        match_record += "\t" + phasing_record[0] + "\t" + phasing_record[1];
	   match_record += "\t" + to_string(best_path.score) + "\n";

	   //complex_match_records[thread_index]->push_back(match_record);
	

       // this line should be recovered
       match_records_by_mode_by_thread[thread_index][mode_index]->push_back(match_record);
    

    }

    baseline_total_match_num[thread_index][threshold_index]->at(mode_index) += truth_num;
    query_total_match_num[thread_index][threshold_index]->at(mode_index) += predict_num;
}

bool WholeGenome::DonorLengthEqual(SequencePath & a, SequencePath & b){
    bool truth_same = false;
    bool query_same = false;

    if(a.donor_sequences[0].length() == b.donor_sequences[0].length() &&
       a.donor_sequences[1].length() == b.donor_sequences[1].length()){
        truth_same = true;
    }
    else if(a.donor_sequences[0].length() == b.donor_sequences[1].length() &&
            a.donor_sequences[1].length() == b.donor_sequences[0].length()){
                truth_same = true;
            }


    if(a.donor_sequences[2].length() == b.donor_sequences[2].length() &&
       a.donor_sequences[3].length() == b.donor_sequences[3].length()){
        query_same = true;
    }
    else if(a.donor_sequences[2].length() == b.donor_sequences[3].length() &&
            a.donor_sequences[3].length() == b.donor_sequences[2].length()){
                query_same = true;
            }

    if(truth_same && query_same) return true;
    return false;
}

bool IsRemovable(SequencePath & s){ return s.removable;}

void WholeGenome::ConvergePaths(list<SequencePath> & path_list){
    //dout << "===========start converge===================" << endl;
    int path_num = path_list.size();
    if(path_num <= 1) return;
    for(list<SequencePath>::iterator i = path_list.begin(); i!= path_list.end(); ++i){
        SequencePath  ref_path = *i;
        if(ref_path.removable) continue;
        if(!ref_path.same_donor_len) continue;
        list<SequencePath>::iterator j = i;
        ++j;
        for(; j != path_list.end(); ++j){
            SequencePath que_path = *j;
            if(que_path.removable) continue;
            if(!que_path.same_donor_len) continue;
            //dout << "Comparing following paths: " << endl;
            //PrintPath(ref_path);
            //PrintPath(que_path);
            if(DonorLengthEqual(ref_path, que_path)){
                if(ref_path.score >= que_path.score){
                    (*j).removable = true;
                    //dout << "delete path: " << endl;
                    //PrintPath((*j));
                }else{
                    (*i).removable = true;
                    //dout << "delete path: " << endl;
                    //PrintPath((*i));
                    break;
                }
            }
            //dout << "-    -     -   -   -   -   -  - - -" << endl;
        }
    }

    path_list.remove_if(IsRemovable);
}

int WholeGenome::test() {
	genome_sequences[0] = "GTCAGCCGG";
	DiploidVariant d1(1, "T", vector<string> ({"A", "C"}), true, true, 0,0,0);
	DiploidVariant d2(4, "G", vector<string> ({"C", ""}), true, false, 0,0,0);
	DiploidVariant d3(5, "C", vector<string> ({"T", ""}), true, false, 0,0,0); // this is false negative
	DiploidVariant d4(6, "C", vector<string> ({"G", ""}), true, false, 0,0,0);
	DiploidVariant d5(7, "G", vector<string> ({"A", ""}), true, false, 0,0,0);
	DiploidVariant d6(1, "T", vector<string> ({"A", "C"}), true, true, 0,0,1);
	DiploidVariant d7(3, "AG", vector<string> ({"A", ""}), true, false, 1,0,1);
	DiploidVariant d8(7, "G", vector<string> ({"GA", ""}), true, false, 0,1,1);

    //complex_match_records = new vector<string>*[1];
    //complex_match_records[0] = new vector<string>;
	//vector<DiploidVariant> var_list = { d2,d3,d4,d5,d7,d8 };
	vector<DiploidVariant> var_list = { d1,d2,d3,d4,d5,d6,d7,d8 };
	//cout << MatchingSingleClusterBaseExtending(var_list, 0) << endl;
	//cout << complex_match_records[0]->at(0) << endl;
	return 0;
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
	//complex_match_records = new vector<string>*[thread_num];
	match_records_by_mode_by_thread = new vector<string>**[thread_num];
    mismatch_query_by_mode_by_thread = new vector<DiploidVariant>**[thread_num];

    //query_matches_by_mode_by_thread = new vector<int> ** [thread_num];

	for(int i = 0; i < thread_num; i++){
        match_records_by_mode_by_thread[i] = new vector<string>*[MATCH_MODE_NUM];
        mismatch_query_by_mode_by_thread[i] = new vector<DiploidVariant>*[MATCH_MODE_NUM];
        //max_path_num_by_mode_by_thread.push_back(vector<int>());
        for(int j = 0; j < MATCH_MODE_NUM; j++){
            match_records_by_mode_by_thread[i][j] = new vector<string>;
            mismatch_query_by_mode_by_thread[i][j] = new vector<DiploidVariant>;
            //max_path_num_by_mode_by_thread[i].push_back(0);
        }
	}

    baseline_total_match_num = new vector<int>** [thread_num];
    query_total_match_num = new vector<int> ** [thread_num];

    baseline_total_edit_distance = new vector<int>** [thread_num];
    query_total_edit_distance = new vector<int>** [thread_num];

    for(int i = 0; i < thread_num; i++){
        
        baseline_total_match_num[i] = new vector<int>* [ROC_SAMPLE_NUM];
        query_total_match_num[i] = new vector<int>* [ROC_SAMPLE_NUM];

        baseline_total_edit_distance[i] = new vector<int> * [ROC_SAMPLE_NUM];
        query_total_edit_distance[i] = new vector<int>* [ROC_SAMPLE_NUM];

        for(int j = 0; j < ROC_SAMPLE_NUM; j++){
            baseline_total_match_num[i][j] = new vector<int>;
            baseline_total_match_num[i][j]->resize(MATCH_MODE_NUM, 0);
            query_total_match_num[i][j] = new vector<int>;
            query_total_match_num[i][j]->resize(MATCH_MODE_NUM, 0);

            baseline_total_edit_distance[i][j] = new vector<int>;
            baseline_total_edit_distance[i][j]->resize(MATCH_MODE_NUM, 0);
            query_total_edit_distance[i][j] = new vector<int>;
            query_total_edit_distance[i][j]->resize(MATCH_MODE_NUM, 0);
        }
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

    //output all results
    cout << "writing results..." << endl;
    ofstream output_stat_file;
    output_stat_file.open(output_dir + "/" + output_prefix+".stat");

    cout << "=========VarMatch Result Stat.=======" << endl;
    string stat_head_string = "#score_unit\tmatch_mode\tscore_unit\tqual_threshold\tbaseline_match_num\tquery_match_num\tquery_total_num\tbaseline_total_ED\tquery_total_ED";
    cout << stat_head_string << endl;
    output_stat_file << "##Baseline:" << baseline_variant_total_num << endl;
    output_stat_file << "##Query:"<< query_variant_total_num << endl;
    output_stat_file << stat_head_string << endl;

    int score_unit;
    int match_mode;
    int score_scheme;

    for(int x = 0; x < score_unit_list.size(); x++){
        score_unit = score_unit_list[x];
        for(int y = 0; y < match_mode_list.size(); y++){
            match_mode = match_mode_list[y];
            for(int z = 0; z < score_scheme_list.size(); z++){
                score_scheme = score_scheme_list[z];
                int mode_index = GetIndexFromMatchScore(score_unit, match_mode, score_scheme);
                int total_ref_complex = 0;
                int total_que_complex = 0;

                string threshold_string = "";
                string baseline_match_num_string = "";
                string query_match_num_string = "";
                string query_total_num_string = "";

                string baseline_edit_distance_string = "";
                string query_edit_distance_string = "";

                for(int t = 0; t < threshold_num; t++){
                    
                    threshold_string += to_string(threshold_list[t]);

                    int baseline_match_num_by_threshold_by_mode = 0;
                    int query_match_num_by_threshold_by_mode = 0;

                    int baseline_edit_distance_by_threshold_by_mode = 0;
                    int query_edit_distance_by_threshold_by_mode = 0;

                    for(int i = 0; i < thread_num; i++){
                        baseline_match_num_by_threshold_by_mode += baseline_total_match_num[i][t]->at(mode_index);
                        query_match_num_by_threshold_by_mode += query_total_match_num[i][t]->at(mode_index);

                        baseline_edit_distance_by_threshold_by_mode += baseline_total_edit_distance[i][t]->at(mode_index);
                        query_edit_distance_by_threshold_by_mode += query_total_edit_distance[i][t]->at(mode_index);
                    }

                    baseline_match_num_string += to_string(baseline_match_num_by_threshold_by_mode);
                    query_match_num_string += to_string(query_match_num_by_threshold_by_mode);
                    query_total_num_string += to_string((int)(query_variant_total_num * (1-per_list[t])) );

                    baseline_edit_distance_string += to_string(baseline_edit_distance_by_threshold_by_mode);
                    query_edit_distance_string += to_string(query_edit_distance_by_threshold_by_mode);

                    if(t < threshold_num-1){
                        threshold_string += ",";
                        baseline_match_num_string += ",";
                        query_match_num_string += ",";
                        query_total_num_string += ",";
                        baseline_edit_distance_string += ",";
                        query_edit_distance_string += ",";
                    }
                                  
                }

                string total_match_num_string = to_string(score_unit) + "\t" +
                                                to_string(match_mode) + "\t" + 
                                                to_string(score_scheme) + "\t" +
                                                threshold_string + "\t" +
                                                baseline_match_num_string + "\t" + 
                                                query_match_num_string + "\t" + 
                                                query_total_num_string;// + "\t" + to_string(mode_index);
                cout << total_match_num_string << "\t" << baseline_edit_distance_string << "\t" << query_edit_distance_string << endl;;
                output_stat_file << total_match_num_string << endl;
            }
        }
    }
    output_stat_file.close();

    int bench_mode_index = GetIndexFromMatchScore(0, 0, 0);

    for(int x = 0; x < score_unit_list.size(); x++){
        score_unit = score_unit_list[x];
        for(int y = 0; y < match_mode_list.size(); y++){
            match_mode = match_mode_list[y];
            for(int z = 0; z < score_scheme_list.size(); z++){
                score_scheme = score_scheme_list[z];
                int mode_index = GetIndexFromMatchScore(score_unit, match_mode, score_scheme);
                string filename_index = to_string(score_unit) + "_" + to_string(match_mode) + "_" + to_string(score_scheme);
                
                ofstream output_complex_file;
                output_complex_file.open(output_dir + "/" + output_prefix+"."+filename_index+".match");

                output_complex_file << "##VCF1:" << ref_vcf_filename << endl;
                output_complex_file << "##VCF2:" << que_vcf_filename << endl;
                output_complex_file << "#CHROM\tPOS\tREF\tALT\tVCF1\tVCF2\tPHASE1\tPHASE2\tSCORE" << endl;

                for(int i = 0; i < thread_num; i++){
                    for(int k = 0; k < match_records_by_mode_by_thread[i][mode_index]->size(); k++){
                        if (match_records_by_mode_by_thread[i][mode_index]->at(k).find_first_not_of(' ') != std::string::npos) {
                            //if(match_records_by_mode_by_thread[i][mode_index]->at(k)[0] == '$'){
                                //int bench_mode_index = stoi(match_records_by_mode_by_thread[i][mode_index]->at(k).erase(0,1));
                                //output_complex_file << match_records_by_mode_by_thread[i][0]->at(k);
                            //}else{
                                output_complex_file << match_records_by_mode_by_thread[i][mode_index]->at(k);
                            //}
                        }
                    }
                }
                output_complex_file.close();
            }
        }
    }

    // int max_path_num = 0;
    // for(int i = 0; i < thread_num; i++){
    //     for(int j = 0; j < MATCH_MODE_NUM; j++){
    //         max_path_num = max(max_path_num, max_path_num_by_mode_by_thread[i][j]);
    //         //max_path_num += max_path_num_by_mode_by_thread[i][j];
    //     }
    // }

    //cout << "max path number: " << max_path_num << endl;

    //==== add mismatch_query_by_mode_by_thread into ref_variant_by_chrid by chr_id
    for(int thread_index = 0; thread_index < thread_num; thread_index++){
        for(int k = 0; k < mismatch_query_by_mode_by_thread[thread_index][0]->size(); k++){
            auto query_variant = mismatch_query_by_mode_by_thread[thread_index][0]->at(k);
            int chr_id = query_variant.unique;
            global_baseline_id ++;
            query_variant.unique = global_baseline_id;
            query_variant.flag = 0;
            //query_variant.matching_condition <<= 1;
            query_variant.matching_condition = 1;
            query_variant.accumulated_vcf_num ++;
            ref_variant_by_chrid[chr_id]->push_back(query_variant);
        }
    }

    // clear all matching records
	for(int i = 0; i < thread_num; i++){
        for(int j = 0; j < MATCH_MODE_NUM; j++){
            delete match_records_by_mode_by_thread[i][j];
            delete mismatch_query_by_mode_by_thread[i][j];

        }
        for(int j = 0; j < ROC_SAMPLE_NUM; j++){
            delete baseline_total_match_num[i][j];
            delete query_total_match_num[i][j];

            delete baseline_total_edit_distance[i][j];
            delete query_total_edit_distance[i][j];
        }
        delete[] match_records_by_mode_by_thread[i];
        delete[] baseline_total_match_num[i];
        delete[] query_total_match_num[i];
        
        delete[] baseline_total_edit_distance[i];
        delete[] query_total_edit_distance[i];
        delete[] mismatch_query_by_mode_by_thread[i];
	}
	delete[] match_records_by_mode_by_thread;
    delete[] baseline_total_match_num;
    delete[] query_total_match_num;

    delete[] baseline_total_edit_distance;
    delete[] query_total_edit_distance;
    delete[] mismatch_query_by_mode_by_thread;
    //max_path_num_by_mode_by_thread.clear();

}

//[TODO] unit test
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
		int flag = 0;
        if(snp.flag) flag = 1;
//        DiploidVariant snp = front_cluster[k];
//        int rq = snp.flag;
        ins_len[flag] += snp.mil;
        del_len[flag] += snp.mdl;
	}
    if(vi_list.size() > 0){
        variant_cluster_by_chrid[chr_id]->push_back(vi_list);
    }
}

int WholeGenome::ReadReferenceVariants(string filename){
    return ReadWholeGenomeVariant(filename, false);
}

int WholeGenome::ReadQueryVariants(string filename){
    return ReadWholeGenomeVariant(filename, true);
}

void WholeGenome::ReadRef(string genome_seq, string ref_vcf){

    ReadWholeGenomeSequence(genome_seq);
    baseline_variant_total_num = ReadReferenceVariants(ref_vcf);
    ref_vcf_filename = ref_vcf;

}

void WholeGenome::ReadDirectRef(string genome_seq, string ref_vcf){

    for(int i = 0; i < 22; i++){
        chrid_by_chrname[to_string(i+1)] = i;
        chrname_by_chrid[i] = to_string(i+1);
    }
    chrid_by_chrname["X"] = 22;
    chrname_by_chrid[22] = "X";
    chrid_by_chrname["Y"] = 23;
    chrname_by_chrid[23]="Y";
    baseline_variant_total_num = ReadReferenceVariants(ref_vcf);
    ref_vcf_filename = ref_vcf;

}

void WholeGenome::OutputMatchingMatrix(){
    ofstream output_matrix_file;
    output_matrix_file.open(output_dir + "/match.matrix");
    output_matrix_file << "##VCF1:" + ref_vcf_filename << endl;
    string title = "#CHROM\tPOS\tREF\tALT\tVCF1";
    for(int i = 0; i < query_vcf_filename_list.size(); i++){
        output_matrix_file << "##VCF" + to_string(i+2) + ":" + query_vcf_filename_list[i] << endl;
        title += "\tVCF"+to_string(i+2);
    }
    output_matrix_file << title << endl;

    int vcf_file_number = query_vcf_filename_list.size() + 1; // including the reference

    for(auto it = chrname_by_chrid.begin(); it != chrname_by_chrid.end(); ++it){
        int chr_id = it->first;
        string chrname = it->second;
        sort(ref_variant_by_chrid[chr_id]->begin(), ref_variant_by_chrid[chr_id]->end());
        for(int i = 0; i < ref_variant_by_chrid[chr_id]->size(); i++){
            DiploidVariant variant = ref_variant_by_chrid[chr_id]->at(i);
            string line = chrname + "\t" + to_string(variant.pos+1) + "\t" + variant.ref + "\t" + variant.alts[0];
            if(variant.multi_alts) line += ","+variant.alts[1];
            vector<string> match_indicator_list;
            long long matching_condition = variant.matching_condition;
            for(int k = 0; k < vcf_file_number; k++){
                int match_indicator = matching_condition & 1;
                match_indicator_list.push_back(to_string(match_indicator));
                matching_condition >>= 1;
            }
            // reverse match_indicator_list
            for(int k = match_indicator_list.size()-1; k >= 0; k--){
                line += "\t"+match_indicator_list[k];
            }

            output_matrix_file << line << endl;
        }
    }
    output_matrix_file.close();
}

void WholeGenome::Compare(string query_vcf,
	string output_prefix,
    bool detail_results,
    int score_unit_,
    int match_mode_,
    int score_scheme_)
{
    // initialize query variant data structure


	que_vcf_filename = query_vcf;
    query_vcf_filename_list.push_back(que_vcf_filename);

    this->output_prefix = output_prefix;
    this->detail_results = detail_results;

    score_unit_indicator = score_unit_;
    match_mode_indicator = match_mode_;
    score_scheme_indicator = score_scheme_;

    if(score_scheme_indicator == 3){
        DirectMatch(ref_vcf_filename, query_vcf, match_mode_, output_prefix);
        return;
    }

    query_variant_total_num = ReadQueryVariants(query_vcf);

    if(score_unit_indicator == -1){
        score_unit_list.push_back(0);
        score_unit_list.push_back(1);
    }else{
        score_unit_list.push_back(score_unit_indicator);
    }

    if(match_mode_indicator == -1){
        match_mode_list.push_back(0);
        match_mode_list.push_back(1);
    }else{
        match_mode_list.push_back(match_mode_indicator);
    }

    if(score_scheme_indicator == -1){
        score_scheme_list.push_back(0);
        score_scheme_list.push_back(1);
        score_scheme_list.push_back(2);
    }else{
        score_scheme_list.push_back(score_scheme_indicator);
    }

    for(int i = 0; i < score_unit_list.size(); i++){
        for(int j = 0; j < match_mode_list.size(); j++){
            for(int k = 0; k < score_scheme_list.size(); k++){
                int mode_index = GetIndexFromMatchScore(score_scheme_list[i], match_mode_list[j], score_scheme_list[k]);
                mode_index_list.push_back(mode_index);  // so that I can directly know how many mode, do not need to calculate all the time
            }
        }
    }

    cout << "Baseline VCF: " << ref_vcf_filename << endl;
    cout << "Query VCF: " << query_vcf << endl;
    cout << "========VCF Stat.==========" << endl;
    cout << "Total Number of VCF Entries: " << endl;
    cout << "Baseline: " << baseline_variant_total_num << "; Query: " << query_variant_total_num << endl;

    cout << "parallel clustering..." << endl;
    ParallelClustering();
    
    cout << "matching variants..." << endl;
    ClusteringMatchMultiThread();

    // most clustering results are cleared inside ParallelClustering function except the following one
    // which is needed for matching
    variants_by_cluster.clear();
    // clean at the end of function
    for(int j = 0; j < chrom_num; j++){
        que_variant_by_chrid[j]->clear();
        //delete que_variant_by_chrid[j];
    }
    //delete[] que_variant_by_chrid;

    query_variant_strings.clear();
    query_variant_total_num = 0;
    threshold_list.clear();
    threshold_num = 0;
    // The following three matching results are cleared inside ClusteringMatchMultiThread function
    // match_records_by_mode_by_thread;
    // baseline_total_match_num;
    // query_total_match_num;

    score_unit_list.clear();
    match_mode_list.clear();
    score_scheme_list.clear();
    mode_index_list.clear();

    return;
}

void WholeGenome::DirectMatch(string ref_vcf, string query_vcf, int match_mode_, string output_prefix)
{
    //dout << "direct match" << endl;
    match_mode_indicator = match_mode_;
    //int ref_variant_num = ReadReferenceVariants(ref_vcf);
    int que_variant_num = ReadQueryVariants(query_vcf);
    //dout << que_variant_num << endl;
    int match_num = 0;
    ofstream output_stat_file;
    output_stat_file.open(output_dir + "/" + output_prefix+".direct");
    for(int i = 0; i < chrom_num; i++){
        if(ref_variant_by_chrid[i]->size() == 0 || que_variant_by_chrid[i]->size() == 0)
            continue;
        //[TODO] not the right way to do it, at least need multimap
        multimap<int, int> ref_variant_by_pos;
        for(int j = 0; j < ref_variant_by_chrid[i]->size(); j++){
            DiploidVariant var = ref_variant_by_chrid[i]->at(j);
            int pos = var.pos;
            ref_variant_by_pos.insert(pair<int, int>(pos, j));
        }

        for(int j = 0; j < que_variant_by_chrid[i]->size(); j++){
            DiploidVariant var = que_variant_by_chrid[i]->at(j);
            int pos = var.pos;
            if(ref_variant_by_pos.find(pos) == ref_variant_by_pos.end())
                continue;

            pair<multimap<int, int>::iterator, multimap<int, int>::iterator> var_range;
            var_range = ref_variant_by_pos.equal_range(pos);

            for(auto it = var_range.first; it != var_range.second; ++it){
                int ref_index = (*it).second;
                DiploidVariant ref_var = ref_variant_by_chrid[i]->at(ref_index);
                if (match_mode_indicator != 1 && var == ref_var){
                    match_num ++;
                    string matched_variant = chrname_by_chrid[i] + "\t" + to_string(ref_var.pos) + "\t" + ref_var.ref + "\t";
                    output_stat_file << matched_variant << endl;
                    break;
                }else if(match_mode_indicator == 1 && var.CompareNoGenotype(ref_var)){
                    match_num ++;
                    break;
                }
            }
        }
    }
    output_stat_file.close();
    dout << "matched variants: " << match_num << endl;
}
