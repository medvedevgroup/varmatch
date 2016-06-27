#include <map>
#include <unordered_map>
#include <vector>
#include <iostream>
#include <algorithm>
#include <tclap/CmdLine.h>
#include "util.h"

using namespace std;

// all intervals are 0 based coordinate

typedef struct Interval {
    int start;
    int end;
    Interval() : start(0), end(0) {}
    Interval(int s, int e) : start(s), end(e) {}
}Interval;

struct compInterval {
    bool operator()(const Interval &a, const Interval &b) const {
        return a.start<b.start;
    }
};

typedef struct Args {
	string baseline_filename;
    vector<string> vcf_filenames;
}Args;

bool TclapParser(Args & args, int argc, char** argv){
	string version = "0.9";

	try {
		std::string desc = "Please cite our paper if you are using this program in your research. \n";
		TCLAP::CmdLine cmd(desc, ' ', version);

		TCLAP::ValueArg<std::string> arg_baseline_filename("b", "baseline", "VCF file", true, "", "file");
		TCLAP::MultiArg<std::string> arg_vcf_filenames("v", "vcf_files", "VCF file list", true, "file list");

        cmd.add(arg_vcf_filenames);
        cmd.add(arg_baseline_filename);

		cmd.parse(argc, argv);

		args.baseline_filename = arg_baseline_filename.getValue();
		args.vcf_filenames = arg_vcf_filenames.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		std::cerr << "error: " << e.error() << " for arg " << e.argId() << "\n";
		abort();
	}
	return true;
}

vector<Interval> merge(vector<Interval> &intervals) {
    sort(intervals.begin(),intervals.end(),compInterval());
    vector<Interval> results;
    for(int i=0; i<intervals.size(); i++) {
        if(results.empty() || results.back().end < intervals[i].start)  // no overlap
            results.push_back(intervals[i]);
        else   // overlap
            results.back().end = max(results.back().end, intervals[i].end);
    }
    return results;
}

int ReadWholeGenomeVariant(string filename,
	vector<vector<Interval>> & interval_list_list,
	vector<multimap<int, string>> & variant_hash_list,
	map<string, int> & chrname_index)
{
    int total_num = 0;
	ifstream vcf_file;
	vcf_file.open(filename.c_str());
	if (!vcf_file.good()) {
		cout << "[VarMatch] Error: can not open vcf file" << endl;
		return -1;
	}

	int genotype_index = -1;
	char genotype_separator = '/';
	int chr_num = 0;
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
		string chr_name = columns[0];
		
		if(chrname_index.find(chr_name) == chrname_index.end()){
			chrname_index[chr_name] = chr_num;
			chr_num++;
			interval_list_list.push_back(vector<Interval>());
			variant_hash_list.push_back(multimap<int, string>());
		}
		int chr_index = chrname_index[chr_name];

		auto pos = atoi(columns[1].c_str()) - 1; // 0-based coordinate
		auto ref = columns[3];
		int end_pos = pos + ref.size();
		auto alt_line = columns[4];

		bool is_heterozygous_variant = false;
		bool is_multi_alternatives = false;

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

        if(snp_ins > 0 || snp_del > 0){
            // this is an indel
        	int interval_start = pos - 10;
        	int interval_end = pos + 10;
        	interval_list_list[chr_index].push_back(Interval(interval_start, interval_end));
        }
        variant_hash_list[chr_index].insert(make_pair(pos, line)); 
        total_num++;
	}
	vcf_file.close();
	return total_num;
}

int FilterComplexVariant(int argc, char* argv[]){
	
	string input_filename = string(argv[1]);
	string output_filename = input_filename + ".cv.vcf";
	vector<vector<Interval>> interval_list_list;
	vector<multimap<int, string>> variant_hash_list;
	map<string, int> chrname_index;

	ReadWholeGenomeVariant(input_filename, 
		interval_list_list, 
		variant_hash_list,
		chrname_index);

	vector<vector<Interval>> merged_intervals_list;
	for(auto interval_list : interval_list_list){
		merged_intervals_list.push_back(merge(interval_list));
	}

	vector<map<int, int>> end_start_list;

	for(auto merged_intervals : merged_intervals_list){
		map<int, int> end_start;
		for(auto i : merged_intervals){
			end_start[i.end] = i.start;
		}
		end_start_list.push_back(end_start);
	}

	vector<unordered_map<int, string>> candidate_variant_hash_list;
	vector<vector<int>> candidate_variant_pos_list;
	cout << "filtering candidate variants..." << endl;
	for(int k = 0; k < variant_hash_list.size(); k++){
		cout << "filtering candidate variants on chromosome " << k << endl;
		auto variant_hash = variant_hash_list[k];
		auto end_start = end_start_list[k];
		//cout << variant_hash.size() << "," << end_start.size() << endl;

		unordered_map<int, string> candidate_variant_hash;
		vector<int> candidate_variant_pos;
		
		for(auto it = variant_hash.begin(); it != variant_hash.end(); ++it){
			int varp = it->first;
			
			auto lowit = end_start.lower_bound(varp);
			if(lowit == end_start.end()) continue;
			int interval_start = lowit->second;
			int interval_end = lowit->first;
			if(varp >= interval_start && varp < interval_end){
				// candidate variant
				candidate_variant_hash[varp] = it->second;
				candidate_variant_pos.push_back(varp);
			}
		}
		candidate_variant_hash_list.push_back(candidate_variant_hash);
		candidate_variant_pos_list.push_back(candidate_variant_pos);
	}
	cout << "filtered all candidate variants." << endl;

	ofstream cv_file;
    cv_file.open(output_filename);

    cout << "filtering complex variants..." << endl;
    for(int k = 0; k < variant_hash_list.size(); k++){
    	auto candidate_variant_pos = candidate_variant_pos_list[k];
    	auto candidate_variant_hash = candidate_variant_hash_list[k];
    	cout << "filtering complex variants on chromosome " << k << endl;
		for(int i = 0; i < candidate_variant_pos.size(); i++){
			int cur_pos = candidate_variant_pos[i];
			if(i > 0){
				int pre_pos = candidate_variant_pos[i-1];
				if(cur_pos - pre_pos <= 10){
					cv_file << candidate_variant_hash[cur_pos] << endl;
					continue;
				}
			}

			if(i < candidate_variant_pos.size() - 1){
				int next_pos = candidate_variant_pos[i+1];
				if(next_pos - cur_pos <= 10){
					cv_file << candidate_variant_hash[cur_pos] << endl;
					continue;
				}
			}
		}
	}
	cout << "finished" << endl;
	cv_file.close();
}



int FiltetCandidateVariant(int argc, char* argv[]){
    Args args;
    TclapParser(args, argc, argv);

    return 0;
}

int main(int argc, char* argv[]){
	return FilterComplexVariant(argc, argv);
}