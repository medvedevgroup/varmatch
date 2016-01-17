// concurrent.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <thread>
#include "vcf.h"

using namespace std;

struct Args {
	string refFilename;
	string donorFilename;
	bool inputFastaRef;
	string refBedFilename;
	string donorBedFilename;
	bool inputBedRef;
	string snpFilename;
	bool snvModification;
	string delVcfFilename;
	bool delModification;
	bool delFromVcf;
	string delBedFilename;
	bool delFromBed;
	string insVcfFilename;
	bool insFromVcf;
	string insFastaFilename;
	bool insModification;
	bool insFromFasta;
	string insBedFilename;
	bool insFromBed;
	string allBedFilename;
	bool allFromBed;

};

bool ParserArgs(Args & args, int argc, char* argv[]) {
	//ez command parsing
	args.snvModification = false;
	args.delModification = false;
	args.insModification = false;
	args.delFromVcf = false;
	args.delFromBed = false;
	args.insFromVcf = false;
	args.insFromFasta = false;
	args.insFromBed = false;
	args.allFromBed = false;
	args.inputFastaRef = false;
	args.inputBedRef = false;

	for (int i = 1; i < argc - 1; i++)
	{
		if (!strcmp(argv[i], "-i")) {
			args.refFilename = string(argv[++i]);
			if (!FileExists(args.refFilename)) {
				cout << "Error: input fasta file not exist" << endl;
				return false;
			}
			args.inputFastaRef = true;
		}
		else if (!strcmp(argv[i], "-o")) {
			args.donorFilename = string(argv[++i]);
		}
		else if (!strcmp(argv[i], "-I")) {
			args.refBedFilename = string(argv[++i]);
			if (!FileExists(args.refBedFilename)) {
				cout << "Error: input bed reference file not exist" << endl;
				return false;
			}
			args.inputBedRef = true;
		}
		else if (!strcmp(argv[i], "-O")) {
			args.donorBedFilename = string(argv[++i]);
		}
		else if (!strcmp(argv[i], "-s")) {
			args.snpFilename = string(argv[++i]);
			if (!FileExists(args.snpFilename)) {
				cout << "Error: input snp file not exist" << endl;
				return false;
			}
			args.snvModification = true;
		}
		else if (!strcmp(argv[i], "-vd")) {
			args.delVcfFilename = string(argv[++i]);
			args.delModification = true;
			args.delFromVcf = true;
		}
		else if (!strcmp(argv[i], "-vi")) {
			args.insVcfFilename = string(argv[++i]);
			args.insModification = true;
			args.insFromVcf = true; //does not implement
		}
		else if (!strcmp(argv[i], "-fi")) {
			args.insFastaFilename = string(argv[++i]);
			args.insModification = true;
			args.insFromFasta = true;
		}
		else if (!strcmp(argv[i], "-bd")) {
			args.delBedFilename = string(argv[++i]);
			args.delModification = true;
			args.delFromBed = true;
			if (!FileExists(args.delBedFilename)) {
				cout << "Error: input del bed file not exist" << endl;
				return false;
			}
		}
		else if (!strcmp(argv[i], "-bi")) {
			args.insBedFilename = string(argv[++i]);
			args.insModification = true;
			args.insFromBed = true;
			if (!FileExists(args.insBedFilename)) {
				cout << "Error: input ins bed file not exist" << endl;
				return false;
			}
		}
		else if (!strcmp(argv[i], "-b")) {
			args.allBedFilename = string(argv[++i]);
			args.insModification = true;
			args.delModification = true;
			args.allFromBed = true;
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
	/*auto a = "a\tb";
	auto b = split(a, '\t');
	for (auto i = 0; i < b.size(); i++) {
		cout << b[i] << endl;
	}*/
	cout << "hello world" << endl;

	if (argc < 2) {
		usage(argv[0]);
		return 0;
	}
	Args args;
	if (!ParserArgs(args, argc, argv)) {
		usage(argv[0]);
		return 0;
	}
    return 0;
}

