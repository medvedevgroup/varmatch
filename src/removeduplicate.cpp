#include "removeduplicate.h"

RemoveDuplicate::RemoveDuplicate(int thread_num_):VCF(thread_num_){}
RemoveDuplicate::~RemoveDuplicate(){}

int RemoveDuplicate::GetThreadIndex(int pos){
    for(int i = 0; i < pos_boundries.size(); i++){
        if(pos < pos_boundries[i]){
            return i;
        }
    }
}

int RemoveDuplicate::ReadVCFWithoutDup(string filename){
    if(!boundries_decided){
        cout << "[Error: RemoveDuplicate] ReadVCFWithoutDup can not read vcf file before read genome file" << endl;
        return -1;
    }

    ifstream vcf_file;
    vcf_file.open(filename.c_str());
    if (!vcf_file.good()) {
        cout << "[Error] RemoveDuplicate::ReadVCFWithoutDup can not open vcf file" << endl;
        return -1;
    }
    int var_num = 0;
    int nodup_var_num = 0;
    while(!vcf_file.eof()){
        string line;
        getline(vcf_file, line, '\n');
        if (line.length() <= 1) continue;
        if (line[0] == '#') continue;
        auto columns = split(line, '\t');
        if(chromosome_name == ".") chromosome_name = columns[0];
        auto pos = atoi(columns[1].c_str()) - 1;
        string ref = columns[3];
        string alt = columns[4];
        string quality = columns[6];

        vector<string> alt_list;
        if(alt.find(",") != string::npos){
            continue;
            // deal with multi alt
            alt_list = split(alt, ',');
        }else{
            alt_list.push_back(alt);
        }

        int thread_index = GetThreadIndex(pos);

        char snp_type;
        for (auto it = alt_list.begin(); it != alt_list.end(); ++it){
            snp_type = 'S';
            string a = *it;
            if(ref.length() > alt.length()){
                snp_type = 'D';
            }else if(ref.length() < alt.length()){
                snp_type = 'I';
            }
            var_num ++;
            string varid = to_string(pos) + "_" + ref + "_" + a;
            transform(varid.begin(), varid.end(), varid.begin(), ::toupper);
            if(nondup_vcfentry_hash_list[thread_index].find(varid) != nondup_vcfentry_hash_list[thread_index].end()){
                nodup_var_num ++;
                nondup_vcfentry_hash_list[thread_index][varid] = line;
                nondup_pos_snp_map_list[thread_index][pos].push_back(SNP(pos, snp_type, ref, a));
            }
        }
    }

    vcf_file.close();
    return var_num;
}
