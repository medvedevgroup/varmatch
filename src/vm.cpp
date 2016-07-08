// concurrent.cpp : Defines the entry point for the console application.
//
//#include "stdafx.h"
#include <iostream>
#include <thread>
#include <tclap/CmdLine.h>
#include "wholegenome.h"

using namespace std;

typedef struct Args {
	string ref_vcf_filename;
	string que_vcf_filename;
	string genome_seq_filename;
	string output_dir;
    string output_prefix;
	int thread_num;
	int score_unit;
	int match_mode;
	int score_scheme;
    bool detail_results;
    vector<string> query_file_list;
    bool pr_curves;
    bool direct_match;

//	bool direct_search;
//	string chr_name;
//	string stat_filename;
//	bool remove_duplicates;
//	string single_vcf_filename;
//	bool match_genotype;
//	bool normalization;
//	bool score_basepair;
//	bool overlap_match;
//	bool variant_check; // check if variant matches
//	bool whole_genome;
}Args;

bool TclapParser(Args & args, int argc, char** argv){
	string version = "0.9";

	try {
		std::string desc = "Please cite our paper if you are using this program in your research. \n";
		TCLAP::CmdLine cmd(desc, ' ', version);
		//TCLAP::ValueArg<std::string> arg_input_vcf_file("i", "i", "input VCF file", true, "", "file", cmd);
		TCLAP::ValueArg<std::string> arg_genome_seq_filename("g", "genome_sequence", "genome sequence FASTA file", true, "", "file");
		TCLAP::ValueArg<std::string> arg_baseline_vcf_filename("b", "baseline", "baseline variant VCF file", true, "", "file");
		TCLAP::MultiArg<std::string> arg_query_vcf_filename("q", "query", "query variant VCF file list", true, "file list");
		TCLAP::ValueArg<std::string> arg_output_dir("o", "output_dir", "output directory, default is current working directory", false, ".", "string");
		TCLAP::ValueArg<std::string> arg_output_prefix("p", "file_prefix", "output filename prefix, default is \"out\"", false, "out", "string");
        int thread_num = (int)thread::hardware_concurrency();


		int max_cores = (int)thread::hardware_concurrency();
		if(max_cores <= 0) max_cores = 1;

		string thread_string = "number of threads, default is the number of available cores (For this machine: " + to_string(max_cores) + ").\n"
                            "If larger than number of available cores or less than 1, automatically set to default value";
		TCLAP::ValueArg<int> arg_thread_num("t", "thread_num", thread_string, false, thread_num, "int");
		vector<int> allowed_two = {-1, 0,1};
		TCLAP::ValuesConstraint<int> allowedVals(allowed_two);

        string score_unit_string = "scoring function/score unit: (Default: -1)\n"
        "-1 : iterate both 0 and 1.\n"
        "0 : the score that a VCF entry contributes is 1.\n"
        "1 : the score that a VCF entry contributes is the edit distance between the new allele and the reference one.\n";

        TCLAP::ValueArg<int> arg_score_unit("u", "score_unit", score_unit_string, false, -1, &allowedVals);


        string match_mode_string = "matching mode: (Default: -1)\n"
        "-1 : iterate both 0 and 1.\n"
        "0 : a set of query entries match a set of baseline entries if, "
        "for each entry, we can select one of the alleles such that the inferred sequences are identical\n"
        "1 : a set of query entries match a set of baseline entries if there exist a phasing of each set such that "
        "the two inferred haplotypes from the query are equal to the two inferred haplotypes from the baseline.\n";

        TCLAP::ValueArg<int> arg_match_mode("m", "match_mode", match_mode_string, false, -1, &allowedVals);


        string score_scheme_string = "scoring scheme: (Default: -1)\n"
        "-1 : iterate 0, 1, and 2 (not including 3)\n"
        "0 : find two subsets of non-overlapping equivalent variants such that "
        "the score of the matched variants is maximized (Default)\n"
        "1 : find two subsets of non-overlapping equivalent variants such that"
        " the score of the chosen baseline variants is maximized\n"
        "2 : find a maximum scoring set of variants in the query such that"
        " each variant can be matched by a subset of the baseline variants\n"
        "3 : (1 to 1 direct match) find a maximum scoring set of entry pairs such that each entry pair contains"
        " one query and one baseline variant that result in the same sequence. In this scheme, different scoring functions and "
        "matching mode have no difference.\n";
        vector<int> allowed_four = {-1,0,1,2,3};
        TCLAP::ValuesConstraint<int> allowedFour(allowed_four);
        TCLAP::ValueArg<int> arg_score_scheme("s", "score_scheme", score_scheme_string, false, -1, &allowedFour);

        //string direct_match_string = "Direct Match. \n";
        //TCLAP::SwitchArg arg_direct_match("d", "direct_match", direct_match_string, cmd, false);

        string detail_results_string = "output detail matching results, by default do not output.\n"
        "filename in format PREFIX.PARAMETER.match\n"
        "The results present which variants in baseline match which variants in query.";

        TCLAP::SwitchArg arg_detail_results("e","detail_results", detail_results_string, cmd, false);

        string precision_recall_string = "Disable Precision-Recall curves. \n";
        TCLAP::SwitchArg arg_disable_curves("C", "disable_curves", precision_recall_string, cmd, false);

        cmd.add(arg_score_scheme);
        cmd.add(arg_match_mode);
        cmd.add(arg_score_unit);
        cmd.add(arg_thread_num);
        //cmd.add(arg_output_prefix);
        cmd.add(arg_output_dir);
        cmd.add(arg_query_vcf_filename);
        cmd.add(arg_baseline_vcf_filename);
        cmd.add(arg_genome_seq_filename);

		cmd.parse(argc, argv);

		args.genome_seq_filename = arg_genome_seq_filename.getValue();
		args.ref_vcf_filename = arg_baseline_vcf_filename.getValue();
		args.query_file_list = arg_query_vcf_filename.getValue();
        //args.que_vcf_filename = arg_query_vcf_filename.getValue();
        args.output_dir = arg_output_dir.getValue();
		args.output_prefix = arg_output_prefix.getValue();
		args.thread_num = arg_thread_num.getValue();
		if(args.thread_num <= 0 || args.thread_num > max_cores) args.thread_num = max_cores;
		args.score_unit = arg_score_unit.getValue();
		args.match_mode = arg_match_mode.getValue();
		args.score_scheme = arg_score_scheme.getValue();
        args.detail_results = arg_detail_results.getValue();
        args.pr_curves = ! arg_disable_curves.getValue();
        //args.direct_match = arg_direct_match.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		std::cerr << "error: " << e.error() << " for arg " << e.argId() << "\n";
		abort();
	}
	return true;
}

int usage(char* command) {
	cout << "\n";
	cout << "\tPlease cite our paper if you are using this program in your research." << endl;
    cout << endl;
	cout << "Usage: " << endl;
    cout << command << " -g genome file path(FASTA format)" << endl;
    cout << "\t-r reference VCF file path" << endl;
    cout << "\t-q query VCF file path" << endl;
    cout << "\t-o output file prefix" << endl;
    cout << "\t[-t thread number]" << endl;
    cout << "\t[-n normalize VCF entries before comparing]" << endl;
    cout << "\t[-m single VCF file to remove duplicates]" << endl;
    cout << "\t[-G do not match genotype when match vcf records]" << endl;
    cout << endl;

	return 0;
}

int main(int argc, char* argv[])
{
//	dout << "Debug Mode" << endl;
//    WholeGenome wg(1);
//    wg.test();
//    return 0;

    Args args;
    TclapParser(args, argc, argv);

    //return 0;
    WholeGenome wg(args.thread_num,
                   args.output_dir,
                   args.pr_curves);

    // if(args.direct_match){
    //     for(int i = 0; i < args.query_file_list.size(); i++){
    //         string query_filename = args.query_file_list[i];
    //         wg.DirectMatch(args.ref_vcf_filename, query_filename, args.match_mode);
    //     }
    //     return 0;
    // }

    wg.ReadRef(args.genome_seq_filename, 
        args.ref_vcf_filename);

    // use a loop 
    for(int i = 0; i < args.query_file_list.size(); i++){
        string query_filename = args.query_file_list[i];

        wg.Compare(query_filename,
            "query"+to_string(i+1),
            args.detail_results,
            args.score_unit,
            args.match_mode,
            args.score_scheme);
    }

    return 0;

//
//    if(args.remove_duplicates){
//        RemoveDuplicate rd(args.thread_num);
//        rd.Deduplicate(args.single_vcf_filename,
//            args.genome_seq_filename,
//            args.direct_search,
//            args.output_filename);
//        return 0;
//	}
//
//	DiploidVCF dv(args.thread_num);
//    dv.Compare(args.ref_vcf_filename,
//		args.que_vcf_filename,
//		args.genome_seq_filename,
//		args.direct_search,
//		args.output_filename,
//		args.match_genotype,
//		args.normalization,
//		args.score_basepair,
//		args.overlap_match,
//		args.variant_check);
//	return 0;
//
//	VCF vcf(args.thread_num);
//	vcf.Compare(args.ref_vcf_filename,
//			args.que_vcf_filename,
//			args.genome_seq_filename,
//			args.direct_search,
//            args.output_filename,
//            args.match_genotype,
//			args.normalization);
    return 0;
}
