// concurrent.cpp : Defines the entry point for the console application.
//
#define DEBUG 1
#include "stdafx.h"
#include <iostream>
#include <thread>
#include "vcf.h"

using namespace std;

typedef struct Args {
	string ref_vcf_filename;
	string que_vcf_filename;
	string genome_seq_filename;
	string fp_filename;
	string fn_filename;
	string output_filename;
	bool direct_search;
	string chr_name;
	string stat_filename;
}Args;

bool ParserArgs(Args & args, int argc, char* argv[]) {
	//ez command parsing
	args.direct_search = false;

	for (int i = 1; i < argc - 1; i++)
	{
		if (!strcmp(argv[i], "-r")) {
			args.ref_vcf_filename = string(argv[++i]);
			if (!FileExists(args.ref_vcf_filename)) {
				cout << "[Error] ParserArgs: input reference vcf file not exist" << endl;
				return false;
			}
		}
		else if (!strcmp(argv[i], "-q")) {
			args.que_vcf_filename = string(argv[++i]);
			if (!FileExists(args.ref_vcf_filename)) {
				cout << "[Error] ParserArgs: input query vcf file not exist" << endl;
				return false;
			}
		}
		else if (!strcmp(argv[i], "-g")) {
			args.genome_seq_filename = string(argv[++i]);
			if (!FileExists(args.genome_seq_filename)) {
				cout << "[Error] ParserArgs: input genome sequence fasta file not exist" << endl;
				return false;
			}
		}
		else if (!strcmp(argv[i], "-p")) {
			args.fp_filename = string(argv[++i]);
		}
		else if (!strcmp(argv[i], "-n")) {
			args.fn_filename = string(argv[++i]);
		}
		else if (!strcmp(argv[i], "-o")) {
			args.output_filename = string(argv[++i]);
		}
		else if (!strcmp(argv[i], "-d")) {
			args.direct_search = true;
		}
		else if (!strcmp(argv[i], "-s")) {
			args.stat_filename = string(argv[++i]);
		}
	}
	return true;
	// end parsing command parameters ...
}

int usage(char* command) {
	cout << "\n";
	cout << "\tPlease cite our paper if you are using this program in your research." << endl;
	cout << endl;

	return 0;
}

int main(int argc, char* argv[])
{
	dout << "Debug Mode" << endl;
	if (argc < 2) {
		usage(argv[0]);
		//return 0;
	}
	Args args;
	if (!ParserArgs(args, argc, argv)) {
		usage(argv[0]);
		return 0;
	}

	args.ref_vcf_filename = "E:\\data\\CHM1.bt2.fb.norm.chr1.vcf";

	VCF vcf;
	dsptime();
	dout << " Read reference vcf file... " << endl;
	vcf.ReadRefVCF(args.ref_vcf_filename);
	dsptime();
	dout << " Finish reading reference vcf file." << endl;
	

    return 0;
}

