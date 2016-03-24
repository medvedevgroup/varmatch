#include <iostream>
#include "util.h"
#include <tclap/CmdLine.h>

class SplitVcf
{
private:
	std::string genome_list_filename;
	std::string vcf_filename;

public:
	SplitVcf(int argc, char* argv);
	~SplitVcf();
	bool Split();
};
