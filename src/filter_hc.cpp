#include <tclap/CmdLine.h>
#include <map>
#include <unordered_map>
#include <iostream>
#include "util.h"

using namespace std;

typedef struct Args {
	string bed_filename;
    vector<string> vcf_filenames;
}Args;

bool TclapParser(Args & args, int argc, char** argv){
	string version = "0.9";

	try {
		std::string desc = "Please cite our paper if you are using this program in your research. \n";
		TCLAP::CmdLine cmd(desc, ' ', version);

		TCLAP::ValueArg<std::string> arg_bed_filename("b", "bedfile", "bedfile", true, "", "file");
		TCLAP::MultiArg<std::string> arg_vcf_filenames("v", "vcf_files", "VCF file list", true, "file list");

        cmd.add(arg_vcf_filenames);
        cmd.add(arg_bed_filename);

		cmd.parse(argc, argv);

		args.bed_filename = arg_bed_filename.getValue();
		args.vcf_filenames = arg_vcf_filenames.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		std::cerr << "error: " << e.error() << " for arg " << e.argId() << "\n";
		abort();
	}
	return true;
}

void ReadBedfile(string bed_filename, 
	map<string, int> & chrname_2_index, 
	vector<map<int, int>> & chr_end_start){

    int chr_num = 0;

    ifstream input(bed_filename);
    if(!input.good()){
        cout << "[Error] Read bed file error" << endl;
        return;
    }

    string line;
    while( std::getline( input, line ).good() )
    {
    	if(line[0] == '#') continue;
    	vector<string> columns = split(line, '\t');
    	string chr_name = columns[0];
    	if(chrname_2_index.find(chr_name) == chrname_2_index.end()){
    		chrname_2_index[chr_name] = chr_num;
    		map<int, int> temp;
    		chr_end_start.push_back(temp);
    		chr_num++;
    	}
    	int chr_index = chrname_2_index[chr_name];
    	//cout << line << endl;
    	int startp = stoi(columns[1]);
    	int endp = stoi(columns[2]);
    	chr_end_start[chr_index][endp] = startp;
    }
    cout << "finish reading bed file" << endl;
    return;

}

void FilterVcfFile(string vcf_filename, 
	map<string, int> & chrname_2_index, 
	vector<map<int, int>> & chr_end_start){
	
	string filter_filename = vcf_filename + ".nhc.vcf";

	ifstream input(vcf_filename);
	if(!input.good()){
		cout << "[Error] Read vcf file " + vcf_filename + " error" << endl;
		return;
	}

	vector<string> output_lines;
	string line;
    while( std::getline( input, line ).good() )
    {
    	if(line[0] == '#'){
    		output_lines.push_back(line);
    		continue;
    	}
    	vector<string> columns = split(line, '\t');
    	string chr_name = columns[0];
    	if(chrname_2_index.find(chr_name) == chrname_2_index.end()){
    		output_lines.push_back(line);
    		continue;
    	}
    	int chr_index = chrname_2_index[chr_name];
    	int varp = stoi(columns[1])-1;
    	map<int, int> & end_start = chr_end_start[chr_index];
    	map<int, int>::iterator itlow, itup;

    	itlow = end_start.lower_bound(varp);
    	int startp = itlow->second;
    	int endp = itlow->first;
    	if(varp >= startp && varp< endp){
    		continue;
    	}
    	
    	output_lines.push_back(line);
    
    }

	ofstream filter_file;
    filter_file.open(filter_filename);
    for(auto line: output_lines){
    	filter_file << line << endl;
    }
    filter_file.close();

}

int main(int argc, char* argv[]){

    Args args;
    TclapParser(args, argc, argv);

    vector<map<int, int>> chr_end_start;
    map<string, int> chrname_2_index;

    ReadBedfile(args.bed_filename, chrname_2_index, chr_end_start);

    vector<string> vcf_filenames = args.vcf_filenames;

    for(auto vcf_filename: vcf_filenames){
    	FilterVcfFile(vcf_filename, chrname_2_index, chr_end_start);
    }

    return 0;
}
