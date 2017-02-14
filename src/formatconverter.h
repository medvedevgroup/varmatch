
#include <iostream>
#include <thread>
#include <tclap/CmdLine.h>
#include <unordered_map>

#include "util.h"


using namespace std;

/*
How member functions and member variables are organized in a class:
	1. public > protect > private
	2. constructor and destructor go first
	3. all member variables are in front of member functions
	4. all class member functions are sorted lexicographically (Except constructor and destructor)
*/
class FormatConverter{
	// convert .match file to intermediate VCF file
	// should be an match database shared by both baseline and query
	// then separately convert, accept two files: VCF and missed VCF
	
	// the whole converting is based on the following assumption:
	//	using chr+pos+ref+alts can uniquely identify a variant in a VCF
	//	no two VCF entries share the same combinations
public:
	FormatConverter(string match_filename){
		buildMatchingDatabase(match_filename);
	}

	bool convertVCF(const string& vcf_filename, const string& output_filename, string miss_filename = "", bool is_query = false){
		unordered_map<string, bool> miss_database;
		if(miss_filename != "")
			readMissVCF(miss_filename, miss_database);

		return convertVCF2Intermediate(vcf_filename, output_filename, miss_database, is_query);
	}


private:

	// database: matched, but only appeared in query
	unordered_map<string, bool> que_match_database;
	// database: matched, but only appeared in baseline
	unordered_map<string, bool> ref_match_database;
	// database to store matched file
	// we have the assumption that all variants has at most two alts for simplicity
	// key is chr_pos_ref_alt, alt sorted lexicographically
	unordered_map<string, bool> shared_database;

	bool add2Database(string line, string chr, unordered_map<string, bool>& database)
	//----
	//	
	//----
	// Code Reviewed 2/13/2017
	{
		vector<string> variant_columns = split(line, ';');
		for(auto variant: variant_columns){
			// corner case: empty string
			if(variant.size() < 3) continue;
			
			vector<string> cols = split(variant, ',');
			vector<string> alt_columns = split(cols[2], '/');
			sort(alt_columns.begin(), alt_columns.end());
			string unique_id = chr + "_" + cols[0] + "_" + cols[1] + "_" + alt_columns[0];
			if(alt_columns.size() > 1) unique_id += "," + alt_columns[1];
			database[unique_id] = true;
		}
	}

	bool addVcfInfo(const vector<string>& columns, string BK_val, string BD_val, string& line)
	//----
	// 
	//----
	{
		if(columns.size() < 10){
			//cout << columns.size() << endl;
			return false;
		}
		string info = columns[8];
		info += ":BK:BD";
		string value = columns[9];
		value += BK_val + ":" + BD_val;
		line = columns[0];
		for(size_t i = 1; i < columns.size()-2; i++) line += "\t" + columns[i];
		line += "\t" + info + "\t" + value + "\n";

		return true;
	}

	bool buildMatchingDatabase(const string & match_filename)
	//----
	// buildMatchingDatabase
	//	read .match file and write into the following three database
	//	1. que_match_database
	//  2. ref_match_database
	//  3. shared_database	
	//----
	// Paramters:
	//	const string& match_filename		filename of .match file	
	//----
	// Partial Effective:
	//	modify three unordered_map database
	//----
	// Code Reviewed 2/13/2017
	{
		ifstream vcf_file;
		vcf_file.open(match_filename.c_str());
		if (!vcf_file.good()) {
			cout << "[FormatConverter] Error: can not open .match file" << endl;
			return false;
		}

		while (!vcf_file.eof()) { // alternative way is vcf_file != NULL
			// read a line
			string line;
			getline(vcf_file, line, '\n');
			
			// check ineligible lines
			if ((int)line.length() <= 1) continue;
			//if (line.find_first_not_of(' ') == std::string::npos) continue;

			if (line[0] == '#') {
				continue;
			}

			vector<string> columns = split(line, '\t');
			auto chr = columns[0];

			auto vcf1 = columns[4];
			auto vcf2 = columns[5];
			
			if(vcf1 == ".")
			// direct match, go to shared_database 
			{
				auto pos = columns[1];
				auto ref = columns[2];
				auto alt_columns = split(columns[3], '/');
				sort(alt_columns.begin(), alt_columns.end());
				string unique_id = chr + "_" + pos + "_" + ref + "_" + alt_columns[0];
				if(alt_columns.size() > 1) unique_id += "," + alt_columns[1];
				shared_database[unique_id] = true;
			}else
			// go to que_match_database and ref_match_database separately
			{
				add2Database(vcf1, chr, ref_match_database);
				add2Database(vcf2, chr, que_match_database);
			}

		}
		vcf_file.close();

		//cout << "build matching database" << endl;
		return true;
	}

	bool convertVCF2Intermediate(const string& vcf_filename, const string& output_filename, unordered_map<string, bool>& miss_database, bool is_query)
	// read a VCF entry and directly converted to intermediate entry
	{
		// check if miss_database is empty
		bool miss_database_NOT_empty = ! miss_database.empty();
		
		// open input file and put all into memory
		// 	otherwise, I/O will be a bottleneck
		vector<string> dataset;

		ifstream vcf_file;
		vcf_file.open(vcf_filename.c_str());
		if (!vcf_file.good()) {
			cout << "[FormatConverter] Error: can not open missed VCF file" << endl;
			return false;
		}

		while (!vcf_file.eof()) { // alternative way is vcf_file != NULL
			// read a line
			string line;
			getline(vcf_file, line, '\n');	
			// check ineligible lines
			if (line.length() <= 1) continue;
			// put all into dataset, including comments
			dataset.emplace_back(line);
		}
		// close input file
		vcf_file.close();

		// open output file
		ofstream output_file;
		output_file.open(output_filename.c_str());
		if(! output_file.good()){
			cout << "[FormatConverter] Error: can not open output file" << endl;
			return false;
		}
		for(auto line: dataset){
			if(line[0] == '#'){
				if(line[1] == '#'){
					output_file << line << endl;
				}else{
					output_file << "##FORMAT=<ID=BK,Number=1,Type=String,Description=\"Sub-type for decision (match/mismatch type)\">" << endl;
					output_file << "##FORMAT=<ID=BD,Number=1,Type=String,Description=\"Decision for call (TP/FP/FN/N)\">" << endl;
					output_file << line << endl;
				}
				continue;
			}
			// parse line, create two version of unique id
			auto columns = split(line, '\t');
			auto alt_columns = split(columns[4], ',');
			sort(alt_columns.begin(), alt_columns.end());
			// short id, without alts

			// long id, with alts sorted lexicographically
			auto long_id = columns[0] + "_" + columns[1] + "_" + columns[3] + "_" + alt_columns[0];
			if(alt_columns.size() > 1) long_id += "," + alt_columns[1];

			string BK_val = ".";
			string BD_val = "N";
			string new_line = "";
			// all use long id
			// if miss_database not empty, each time should check if it exist
			if(miss_database_NOT_empty){
				if(miss_database.find(long_id) != miss_database.end()){
					bool add_info = addVcfInfo(columns, ".", "N", new_line);
					output_file << (add_info ? new_line : line);
					continue;
				}
			}
			//cout << long_id << endl;
			// each time should also check if exist in shared_database 
			//for(auto it = shared_database.begin(); it != shared_database.end(); ++it){
			//	cout << it->first << endl;
			//}
			if(shared_database.find(long_id) != miss_database.end()){
				bool add_info = addVcfInfo(columns, "gm", "TP", new_line);
				output_file << (add_info ? new_line : line);
				continue;
			}
			// check ref/que_match_database
			if(is_query){
				// check que_match_database
				if(que_match_database.find(long_id) != que_match_database.end()){
					bool add_info = addVcfInfo(columns, ".", "FP", new_line);
					output_file << (add_info ? new_line : line);
					continue;
				}
			}else{
				// check ref_match_database
				if(ref_match_database.find(long_id) != ref_match_database.end()){
					bool add_info = addVcfInfo(columns, ".", "FN", new_line);
					output_file << (add_info ? new_line : line);
					continue;
				}
			}

			bool add_info = addVcfInfo(columns, ".", "N", new_line);
			output_file << (add_info ? new_line : line);
		}

		// close output file
		output_file.close();
	}

	bool readMissVCF(const string & miss_filename, unordered_map<string, bool>& db)
	//----
	// readMissVCF
	// 	return a database of not accessed VCF entries
	// 	either not pass Quality threshold or length threshold
	//----
	// Parameter:
	//	const string&	miss_filename 		VCF entries not accessed
	// 	unordered_map<string, bool>& db 	database of not accessed VCF entries
	//----
	// Partial Effect:
	//	database of missed VCF, key is chr_pos_ref_alts, alts are separated by "," and sorted lexicographically
	//----
	// Code reviewed 2/13/2017
	{
		ifstream vcf_file;
		vcf_file.open(miss_filename.c_str());
		if (!vcf_file.good()) {
			cout << "[FormatConverter] Error: can not open missed VCF file" << endl;
			return false;
		}

		while (!vcf_file.eof()) { // alternative way is vcf_file != NULL
			// read a line
			string line;
			getline(vcf_file, line, '\n');
			
			// check ineligible lines
			if ((int)line.length() <= 1) continue;
			//if (line.find_first_not_of(' ') == std::string::npos) continue;

			if (line[0] == '#') {
				continue;
			}

			vector<string> columns = split(line, '\t');
			auto chr = columns[0];
			auto pos = columns[1];
			auto ref = columns[3];
			
			// create a unique string for alts
			vector<string> alt_columns = split(columns[4], ',');
			sort(alt_columns.begin(), alt_columns.end());
			string alt_unique = (alt_columns.empty() ? "" : alt_columns[0]);
			if(alt_columns.size() > 1) alt_unique += "," + alt_columns[1];

			auto unique_id = chr + "_" + pos + "_" + ref + "_" + alt_unique;

			// add to database
			db[unique_id] = true;
		}
		vcf_file.close();

		return true;
	}

};
