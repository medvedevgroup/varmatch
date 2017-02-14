#include "formatconverter.h"
typedef struct Args {
	string ref_vcf_filename;
	string que_vcf_filename;
	string match_filename;
	string ref_miss_filename;
    string que_miss_filename;
    string output_ref_filename;
    string output_que_filename;
}Args;

bool TclapParser(Args & args, int argc, char** argv){
	string version = "0.9";

	try {
		std::string desc = "Please cite the following paper if you are using this program in your research: \n";
		desc += "@article{sun2016varmatch,\n";
		desc += "  title={VarMatch: robust matching of small variant datasets using flexible scoring schemes},\n";
		desc += "  author={Sun, Chen and Medvedev, Paul},\n";
		desc += "  journal={Bioinformatics},\n";
  		desc += "  pages={btw797},\n";
  		desc += "  year={2016},\n";
  		desc += "  publisher={Oxford Univ Press}\n";
		desc += "}\n";
		TCLAP::CmdLine cmd(desc, ' ', version);
		//TCLAP::ValueArg<std::string> arg_input_vcf_file("i", "i", "input VCF file", true, "", "file", cmd);
		TCLAP::ValueArg<std::string> arg_match_filename("m", "match", "Match file: filename.match", true, "", "file");
		
		TCLAP::ValueArg<std::string> arg_ref_vcf_filename("b", "baseline", "baseline variant VCF file", true, "", "file");
		TCLAP::ValueArg<std::string> arg_que_vcf_filename("q", "query", "query variant VCF file list", true, "", "file");

		TCLAP::ValueArg<std::string> arg_ref_miss_filename("B", "miss_baseline", "VCF entries not accessed in baseline variant VCF file", false, "", "file");
		TCLAP::ValueArg<std::string> arg_que_miss_filename("Q", "miss_query", "VCF entries not accessed in query variant VCF file list", false, "", "file");
		
		TCLAP::ValueArg<std::string> arg_output_ref_filename("o", "out_baseline", "intermediate baseline VCF filename", false, "baseline.intermediate.vcf", "file");
		TCLAP::ValueArg<std::string> arg_output_que_filename("u", "out_query", "intermediate query VCF filename", false, "query.intermediate.vcf", "file");

        cmd.add(arg_output_que_filename);
        cmd.add(arg_output_ref_filename);
        cmd.add(arg_que_miss_filename);
        cmd.add(arg_ref_miss_filename);

        cmd.add(arg_que_vcf_filename);
        cmd.add(arg_ref_vcf_filename);
        cmd.add(arg_match_filename);

		cmd.parse(argc, argv);

		args.match_filename = arg_match_filename.getValue();
		args.ref_vcf_filename = arg_ref_vcf_filename.getValue();
		args.que_vcf_filename = arg_que_vcf_filename.getValue();
		args.ref_miss_filename = arg_ref_miss_filename.getValue();
		args.que_miss_filename = arg_que_vcf_filename.getValue();
		args.output_ref_filename = arg_output_ref_filename.getValue();
		args.output_que_filename = arg_output_que_filename.getValue();
	}
	catch (TCLAP::ArgException &e)
	{
		std::cerr << "error: " << e.error() << " for arg " << e.argId() << "\n";
		abort();
	}
	return true;
}

void noMoreMemory()
{
	cerr << "Unable to satisfy request for memory\n";
	abort();
}

int main(int argc, char* argv[])
{
	
    Args args;
    TclapParser(args, argc, argv);

    // handle not enough memory error
    set_new_handler(noMoreMemory);

    // use smart ptr
    unique_ptr<FormatConverter> fc(new FormatConverter(args.match_filename));

    fc->convertVCF(args.ref_vcf_filename, args.output_ref_filename, args.ref_miss_filename);
    
    fc->convertVCF(args.que_vcf_filename, args.output_que_filename, args.que_miss_filename, true);

    return 0;
}