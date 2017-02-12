// data structure for direct search
class DiploidVariant {
public:
    DiploidVariant(int pos_ = -1,
        string ref_ = "",
        vector<string> alts_ = {"",""},
        bool heterozygous_ = false,
        bool multi_alts_ = false,
        int mdl_ = 0,
        int mil_ = 0,
        bool flag_ = false,
        double qual_ = 0.0,
        bool zero_one_var_ = false) :
        pos(pos_),
        ref(ref_),
        alts(alts_),
        heterozygous(heterozygous_),
        multi_alts(multi_alts_),
        mdl(mdl_),
        mil(mil_),
        flag(flag_),
        qual(qual_),
        zero_one_var(zero_one_var_)
        { }

    int pos;
    string ref;
    vector<string> alts;
    bool heterozygous;
    bool multi_alts;
    bool zero_one_var; // which means the phasing should be 0/1 or 1/0, no matter if it contains multi_alts
    // i.e. multi_alts does not mean that it is 1/2 or 2/1
    int mdl;
    int mil;
    bool flag; //in DiploidVariant, flag = false is reference, flag = true is query
    // keep flag as int? not necessary
    double qual;
    long long matching_condition;
    int unique;

    // for giab pipeline, initialize as "."
    string GQ_val;
    string DP_val;
    string AD_val;

    int accumulated_vcf_num; // this is for debug

//    int get_pos() const{return pos};
//    string get_ref() const{return ref};
//    vector<string> get_alts() const{return alts};
//    bool get_heterozygous() const{return heterozygous};
//    bool get_multi_alts() const{return multi_alts};

    bool operator <(const DiploidVariant& y) const {
        return pos < y.pos;
    }

    // this is based on the assumption that all sequence are in upper case
    bool operator ==(const DiploidVariant& y) {
        if (pos == y.pos && ref == y.ref) {
            if(heterozygous == y.heterozygous && multi_alts == y.multi_alts){
                if (multi_alts && heterozygous) {
                    int match_times = 0;
                    for (int i = 0; i < 2; i++) {
                        for (int j = 0; j < 2; j++) {
                            if (alts[i] == y.alts[j])
                                match_times++;
                        }
                    }
                    if (match_times >= 2)
                        return true;
                }
                else if(alts[0] == y.alts[0]){
                    return true;
                }
            }
            if(multi_alts && zero_one_var && y.multi_alts && y.zero_one_var){
                int match_times = 0;
                for (int i = 0; i < 2; i++) {
                    for (int j = 0; j < 2; j++) {
                        if (alts[i] == y.alts[j])
                            match_times++;
                    }
                }
                if(match_times > 1) return true;
            }
        }
        return false;
    }

    bool DirectCompare(const DiploidVariant& y){
        if (pos == y.pos && ref == y.ref) {
            if (multi_alts && heterozygous && y.multi_alts && y.heterozygous) {
                int match_times = 0;
                for (int i = 0; i < 2; i++) {
                    for (int j = 0; j < 2; j++) {
                        if (alts[i] == y.alts[j])
                            match_times++;
                    }
                }
                if (match_times > 0)
                    return true;
            }
            else if(alts[0] == y.alts[0]){
                return true;
            }
        }
        return false;
    }

    bool CompareNoGenotype(const DiploidVariant & y){
        if(pos == y.pos && ref == y.ref){
            if(alts[0] == y.alts[0]) return true;
            if(multi_alts){
                if(alts[1] == y.alts[0]) return true;
                if(y.multi_alts && alts[1] == y.alts[1]){
                    return true;
                }
            }
            if(y.multi_alts && alts[0] == y.alts[1]){
                return true;
            }
        }
        return false;
    }

};
