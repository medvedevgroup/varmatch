#include "diploid.h"


bool operator <(const DiploidVariant& x, const DiploidVariant& y) {
	return x.pos < y.pos;
}

bool operator ==(const DiploidVariant& x, const DiploidVariant& y) {
	if (x.pos == y.pos && x.var_type == y.var_type && x.heterozygous == y.heterozygous && x.multi_alts == y.multi_alts) {
		if (x.multi_alts && x.heterozygous) {
			int match_times = 0;
			for (int i = 0; i < 2; i++) {
				for (int j = 0; j < 2; j++) {
					if (x.alts[i] == y.alts[j])
						match_times++;
				}
			}
			if (match_times >= 2)
				return true;
		}
		else if(x.alts[0] == y.alts[0]){
			return true;
		}
	}
	return false;
}

DiploidVCF::DiploidVCF(int thread_num_)
{
    debug_f = 0;
	genome_sequence = "";
	boundries_decided = false;
	clustering_search = false;
	match_genotype = true;
    if (thread_num_ <= 0) {
		thread_num = 1;
	}
	else {
		thread_num = min(thread_num_, (int)thread::hardware_concurrency());
	}
	dout << "Thread Number: " << thread_num << endl;
    chromosome_name = ".";
}

DiploidVCF::~DiploidVCF()
{
}

bool DiploidVCF::NormalizeDiploidVariant(DiploidVariant & var) {
	pos = var.pos;
	parsimonious_ref = var.ref;
	parsimonious_alt0 = var.alts[0];
	parsimonious_alt1 = var.alts[0];
	if (var.heterozygous && var.multi_alts)
		parsimonious_alt1 = var.alts[1]

	int left_index = pos;
	if (genome_sequence.size() == 0) return false;
	if (parsimonious_ref.size() == 1 && parsimonious_alt0.size() == 1 && parsimonious_alt1.size() == 1) return true;
	if (toupper(genome_sequence[left_index]) != toupper(parsimonious_ref[0])) {
		dout << "[Error] genome sequence, subsequence, offset does not match." << endl;
		return false;
	}
	bool change_in_allels = true;
	while (change_in_allels) {
		change_in_allels = false;
		if (toupper(parsimonious_ref.back()) == toupper(parsimonious_alt0.back()) && toupper(parsimonious_ref.back()) == toupper(parsimonious_alt1.back())) {
			if ((parsimonious_ref.size() > 1 && parsimonious_alt0.size() > 1 && parsimonious_alt1.size() > 1) || left_index > 0) {
				parsimonious_ref.pop_back();
				parsimonious_alt0.pop_back();
				parsimonious_alt1.pop_back();
				change_in_allels = true;
			}
			else {
				return false;
			}
		}
		if (parsimonious_ref.length() == 0 || parsimonious_alt0.length() == 0 || parsimonious_alt1.length() == 0) {
			left_index--;
			char left_char = genome_sequence[left_index];
			parsimonious_ref = left_char + parsimonious_ref;
			parsimonious_alt0 = left_char + parsimonious_alt0;
			parsimonious_alt1 = left_char + parsimonious_alt1;
		}
	}
	while (toupper(parsimonious_ref[0]) == toupper(parsimonious_alt0[0]) && toupper(parsimonious_ref[0]) == toupper(parsimonious_alt1[0]) && parsimonious_ref.size() > 1 && parsimonious_alt0.size() > 1 && parsimonious_alt1.size() > 1) {
		parsimonious_ref.erase(0, 1);
		parsimonious_alt0.erase(0, 1);
		parsimonious_alt1.erase(0, 1);
	}
	var.pos = left_index;
	var.ref = parsimonious_ref;
	var.alts[0] = parsimonious_alt0;
	if (var.heterozygous && var.multi_alts)
		var.alts[1] = parsimonious_alt1;
	return true;
}

bool DiploidVCF::ReadDiploidVCF() {

}


