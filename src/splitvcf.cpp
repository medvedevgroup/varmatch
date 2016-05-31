//#include "stdafx.h"
#include "splitvcf.h"

SplitVcf::SplitVcf(int argc, char* argv) {
	version = "0.9";

	try {
		std::string desc = "split vcf file according to chromosome. \n";
		TCLAP::CmdLine cmd(desc, ' ', version);
		//TCLAP::ValueArg<std::string> arg_input_vcf_file("i", "i", "input VCF file", true, "", "file", cmd);
		TCLAP::UnlabeledValueArg<std::string> arg_input_vcf_file("<in.vcf>", "input VCF file", true, "", "file", cmd);
		TCLAP::ValueArg<std::string> arg_genome_list_file("g", "g", "genome list file", true, "", "file", cmd);
		cmd.parse(argc, argv);
	}
	catch (TCLAP::ArgException &e)
	{
		std::cerr << "error: " << e.error() << " for arg " << e.argId() << "\n";
		abort();
	}
}

SplitVcf::~SplitVcf() {

}

int main(int argc, char* argv[]) {
	SplitVcf sf = new SplitVcf(argc, argv);
	sf->Split();
}

