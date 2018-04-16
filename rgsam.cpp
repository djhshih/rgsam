#include <fstream>
#include <iostream>
#include <string>
#include <stdexcept>
#include <cstring>
#include <list>
#include <set>
#include <map>

using namespace std;


struct fastq_entry {
    /// read name
    string qname;
    /// read sequence
    string seq;
    /// read quality scores
    string qual;
};

struct sam_opt_field {
    char tag[2];
    char type;
    string value;

    sam_opt_field(const char _tag[], char _type, const string& _value)
    : type(_type), value(_value) {
        tag[0] = _tag[0];
        tag[1] = _tag[1];
    }

    sam_opt_field(const string& x) {
        if (x.length() > 5) {
            tag[0] = x[0];
            tag[1] = x[1];
            if (x[2] == ':') {
                type = x[3];
                if (x[4] == ':') {
                    value = x.substr(5);
                }
            }
        }
    }

    void write(ostream& f) const {
        f << tag[0] << tag[1] << ':' << type << ':' << value;
    }
};

sam_opt_field read_group_field(const string& rg) {
    return sam_opt_field("RG", 'Z', rg);
}


struct sam_entry {
    // read name
    string qname;
    // bitwise flag
    unsigned short int flag;
    // reference sequence name
    string rname;
    // 1-based leftmost mapping position
    unsigned long int pos;
    // mapping quality
    unsigned char mapq;
    // CIGAR tring
    string cigar;
    // reference name of mate/next read
    string rnext;
    // position of the mate/next read
    unsigned long int pnext;
    // observed template length
    signed long long int tlen; 
    // sequence
    string seq;
    // ascii of phred-scaled base quality + 33
    string qual;
    // optional fields
    list<sam_opt_field> opts;
};

/// SAM entry where core fields are not parsed
struct raw_sam_entry {
    // raw core fields
    string core;
    // optional fields
    list<sam_opt_field> opts;
};

const size_t n_sam_core_fields = 11;
const char sam_delim = '\t';
const char* platform = "illumina";

/**
 * Natural string comparison, account for numbers within a string.
 * This is different from lexigraphical comparison (e.g. std::string::compare).
 * Taken from https://github.com/samtools/samtools/blob/develop/bam_sort.c#L13
 */
int strnum_cmp(const char *_a, const char *_b) {
    const unsigned char *a = (const unsigned char*)_a, *b = (const unsigned char*)_b;
    const unsigned char *pa = a, *pb = b;
    while (*pa && *pb) {
        if (isdigit(*pa) && isdigit(*pb)) {
            while (*pa == '0') ++pa;
            while (*pb == '0') ++pb;
            while (isdigit(*pa) && isdigit(*pb) && *pa == *pb) ++pa, ++pb;
            if (isdigit(*pa) && isdigit(*pb)) {
                int i = 0;
                while (isdigit(pa[i]) && isdigit(pb[i])) ++i;
                return isdigit(pa[i])? 1 : isdigit(pb[i])? -1 : (int)*pa - (int)*pb;
            } else if (isdigit(*pa)) return 1;
            else if (isdigit(*pb)) return -1;
            else if (pa - a != pb - b) return pa - a < pb - b? 1 : -1;
        } else {
            if (*pa != *pb) return (int)*pa - (int)*pb;
            ++pa; ++pb;
        }
    }
    return *pa? 1 : *pb? -1 : 0;
}

/**
 * Compare two read names by lexical or natural order.
 */
int compare_qnames(const string& a, const string& b, bool natural) {
    size_t a_end = a.length();
    size_t b_end = b.length();
    if (a.substr(a_end - 2) == "/1") {
        a_end -= 2;
    }
    if (b.substr(b_end - 2) == "/2") {
        b_end -= 2;
    }
    
    if (natural) {
        return strnum_cmp(a.substr(0, a_end).c_str(), b.substr(0, b_end).c_str());
    }
    return a.substr(0, a_end).compare(b.substr(0, b_end));
}

size_t find_in_string(const string& x, char c, size_t pos, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        pos = x.find(c, pos);
        if (pos == string::npos) break;
    }
    return pos;
}

void parse_opts(const string& x, list<sam_opt_field>& opts) {
    size_t pos = 0;
    while (pos != string::npos) {
        size_t end = x.find(sam_delim);
        opts.push_back(sam_opt_field(x.substr(pos, end - pos)));
    }
}

/**
 * Read a SAM entry but do not parse the core fields.
 */
bool extract_raw_sam_entry(const string& line, raw_sam_entry& x) {
    if (line.empty()) return false;

    // split entry into core string and parsed opts
    size_t pos = find_in_string(line, sam_delim, 0, n_sam_core_fields);
    x.core = line.substr(0, pos);
    parse_opts(line.substr(pos + 1), x.opts);
    
    return true;
}

/**
 * Write a raw SAM entry.
 */
bool write_raw_sam_entry(ostream& f, const raw_sam_entry& x) {
    f << x.core;
   
    for (list<sam_opt_field>::const_iterator it = x.opts.begin(); it != x.opts.end(); ++it) {
        f << sam_delim;
        it->write(f);
    }
    f << endl;

    return true;
}

/**
 * Functor for matching a tag.
 */
struct tag_match
{
    tag_match(const char _tag[]) {
        tag[0] = _tag[0];
        tag[1] = _tag[1];
    }

    bool operator()(const sam_opt_field& x) const {
        return x.tag[0] == tag[0] && x.tag[1] == tag[1];
    }

    char tag[2];
};

void replace_opt_field(list<sam_opt_field>& opts, const sam_opt_field& x) {
    opts.remove_if(tag_match(x.tag));
    opts.push_back(x);
}

/**
 * Infer read-group based on flowcell id and lane id,
 * assuming Illumin 1.8 read name format.
 */
string infer_read_group_illumina18(const string& qname) {
    // extract flowcell
    size_t start = find_in_string(qname, ':', 0, 2);
    if (start == string::npos) return "";
    ++start;
    size_t end = find_in_string(qname, ':', start, 1);
    if (end == string::npos) return "";
    string flowcell = qname.substr(start, end - start);

    // extract lane
    start = end + 1;
    end = find_in_string(qname, ':', start, 1);
    if (end == string::npos) return flowcell;
    string lane = qname.substr(start, end - start);

    return flowcell + '_' + lane;
}

/** 
 * Read one fastq entry from file.
 *
 * Assume that each sequence and quality score entry `r1.fq` and `r2.fq` is single-line.
 */
bool read_fastq_entry(istream& f, fastq_entry& x) {
    getline(f, x.qname);
    if (x.qname.empty()) return false;
    
    getline(f, x.seq);

    string marker;
    getline(f, marker);
    if (marker != "+") {
        throw runtime_error("fastq entry is malformed");
    }

    getline(f, x.qual);
    
    return true;
}

string get_qname_from_sam_core(const string& core) {
    return core.substr(0, core.find(sam_delim));
}

/** 
 * Write one fastq entry to file.
 */
void write_fastq_entry(ostream& f, const fastq_entry& x) {
    f << x.qname << endl;
    f << x.seq << endl;
    f << '+' << endl;
    f << x.qual << endl;
}

bool read_read_groups(ifstream& f, map<string, string>& rgs) {
    while(true) {
        string line;
        getline(f, line);
        if (line.empty()) break;

        if (line.find("@RG") != 0) break;

        size_t start = line.find("ID:", 4);
        if (start == string::npos) break;
        size_t end = line.find(sam_delim, start);
        if (end == string::npos) break;

        start += 3;
        string id = line.substr(start, end - start);
        rgs[id] = line;
    }

    return !rgs.empty();
}

void write_read_groups(ostream& f, const set<string> rgs, char* sample, char* library) {
    for (set<string>::iterator it = rgs.begin(); it != rgs.end(); ++it) {
        f   << "@RG" << sam_delim
            << "ID:" << *it << sam_delim
            << "PU:" << *it << sam_delim
            << "SM:" << sample << sam_delim
            << "LB:" << library << sam_delim
            << "PL:" << platform << endl;
    }
}

void collect_rg_from_sam(char* in_fname, char* sample, char* library, char* out_rg_fname) {
    // collect read-groups
    set<string> rgs;
    ifstream sam_f(in_fname);
    while (true) {
        string line;
        getline(sam_f, line);

        if (line.empty()) break;

        // skip header lines
        if (line[0] == '@') continue;
    
        // infer read-group
        raw_sam_entry x; 
        if (!extract_raw_sam_entry(line, x)) break;
        string rg = infer_read_group_illumina18(get_qname_from_sam_core(x.core));
        rgs.insert(rg);
    }
    sam_f.close();

    // write all read groups
    ofstream rg_f(out_rg_fname);
    write_read_groups(rg_f, rgs, sample, library);
    rg_f.close();
}

void collect_rg_from_fq(char* in_fname, char* sample, char* library, char* out_rg_fname) {
    // collect read-groups
    set<string> rgs;
    ifstream fq_f(in_fname);
    while (true) {
        fastq_entry x;
        if (!read_fastq_entry(fq_f, x)) break;
    
        // infer read-group
        string rg = infer_read_group_illumina18(x.qname);
        rgs.insert(rg);
    }
    fq_f.close();

    // write all read groups
    ofstream rg_f(out_rg_fname);
    write_read_groups(rg_f, rgs, sample, library);
    rg_f.close();
}

void tag_sam_with_rg(char* in_fname, char* rg_fname, char* out_sam_fname) {
    ifstream rg_f(rg_fname);
    map<string, string> rgs;
    read_read_groups(rg_f, rgs);
    rg_f.close();

    ifstream in_f(in_fname);
    ofstream out_f(out_sam_fname);

    string line;

    while (true) {
        getline(in_f, line);

        if (line.empty()) break;

        if (line[0] == '@') {
            // skip RG header line but copy other header lines verbatim
            if (line.find("@RG") != 0) {
                out_f << line;
            }
            continue;
        } else {
            break;
        }
    }
    
    // process SAM entries
    while (true) {
        raw_sam_entry x; 
        if (!extract_raw_sam_entry(line, x)) break;

        string rg = infer_read_group_illumina18(get_qname_from_sam_core(x.core));

        if (rgs.find(rg) == rgs.end()) {
            cerr << "Warning: read group ID " << rg << " is not found in input read-groups" << endl;
        }
        
        // tag read with inferred read group
        replace_opt_field(x.opts, read_group_field(rg));

        // write modified SAM entry to out file
        write_raw_sam_entry(out_f, x);

        // get next line
        getline(in_f, line);
        if (line.empty()) break;
    }

    in_f.close();
    out_f.close();
}


/**
 * Utility programs
 *
 * collect    collect read-group information from SAM file
 * collectfq  collect read-group information from FASTQ file
 * tag        add read-group field to reads
 * split      split SAM file based on read-group
 * splitfq    split FASTQ file based on read-group
 *
 * Read-groups identifier (ID) and platform unit (PU) are inferred from read 
 * names according to Illumina's read name format.
 * Platform (PL) is assumed to be "illumina".
 * Sample (SM) and library identifier (LB) must be given.
 *
 * Files with reads from more than one sample or library are not supported.
 *
 * To split BAM or SAM files with proper read-group information, use instead:
 * samtools view -r <rgid> <in.bam>
 */
int main(int argc, char* argv[]) {

    int nargs = argc - 1;

    if (nargs < 1) {
        cerr << "usage: rgsam collect <in.sam> <sample> <library> <rg.txt>" << endl;
        cout << "       rgsam collectfq <in.fq> <sample> <library> <rg.txt>" << endl;
        cout << "       rgsam tag <in.sam> <rg.txt> <out.sam>" << endl;
        cout << "       rgsam split <in.sam> <sample> <library>" << endl;
        cout << "       rgsam splitfq <in.fq> <sample> <library>" << endl;
        return 1;
    }

    if (strcmp(argv[1], "collect") == 0) {

        char* in_fname = argv[2];
        char* sample = argv[3];
        char* library = argv[4];
        char* out_rg_fname = argv[5];
        collect_rg_from_sam(in_fname, sample, library, out_rg_fname);

    } else if (strcmp(argv[1], "collectfq") == 0) {

        char* in_fname = argv[2];
        char* sample = argv[3];
        char* library = argv[4];
        char* out_rg_fname = argv[5];
        collect_rg_from_fq(in_fname, sample, library, out_rg_fname);

    } else if (strcmp(argv[1], "tag") == 0) {

        char* in_fname = argv[2];
        char* rg_fname = argv[3];
        char* out_sam_fname = argv[4];
        tag_sam_with_rg(in_fname, rg_fname, out_sam_fname);

    } else if (strcmp(argv[1], "split") == 0) {

        char* in_fname = argv[2];
        char* sample = argv[3];
        char* library = argv[4];

        cerr << "Unimplemented command: " << argv[1] << endl;

    } else if (strcmp(argv[1], "splitfq") == 0) {

        char* in_fname = argv[2];
        char* sample = argv[3];
        char* library = argv[4];

        cerr << "Unimplemented command: " << argv[1] << endl;

    } else {

        cerr << "Unsupported command: " << argv[1] << endl;
        return 1;

    }

    return 0;
}
