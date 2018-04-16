#include <fstream>
#include <iostream>
#include <string>
#include <stdexcept>
#include <cstring>
#include <list>
#include <set>
#include <map>

#include "rgsam/fastq.hpp"

using namespace std;

const char* platform = "illumina";

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

/**
 * Create a read-group optional field entry.
 */
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

/**
 * Find the n-th occurence of a character after offset in a string.
 */
size_t find_in_string(const string& x, char c, size_t pos, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        pos = x.find(c, pos + 1);
        if (pos == string::npos) break;
    }
    return pos;
}

/**
 * Parse optional fields from a string.
 */
void parse_opts(const string& x, list<sam_opt_field>& opts) {
    size_t pos = 0;
    while (true) {
        size_t end = x.find(sam_delim, pos);
        opts.push_back(sam_opt_field(x.substr(pos, end - pos)));
        if (end == string::npos) break;
        pos = end + 1;
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

/**
 * Replace an optional field.
 */
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
 * Get the read name from the string of the SAM core fields.
 */
string get_qname_from_sam_core(const string& core) {
    return core.substr(0, core.find(sam_delim));
}

/**
 * Read read groups from a SAM header file.
 *
 * Each line begins with `@RG`.
 */
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

/**
 * Write read groups into a SAM header file.
 */
void write_read_groups(ostream& f, const map<string, string>& rgs) {
    for (map<string, string>::const_iterator it = rgs.begin(); it != rgs.end(); ++it) {
        f << it->second << endl;
    }
}

/**
 * Write read groups into a SAM header file.
 *
 * @param rgs  raw read-group value strings.
 */
void write_read_groups(ostream& f, const set<string> rgs, const char* sample, 
        const char* library, const char*platform) {
    for (set<string>::const_iterator it = rgs.begin(); it != rgs.end(); ++it) {
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
    if (!sam_f.is_open()) throw runtime_error("Error: input SAM file could not be opened.");
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
    if (!rg_f.is_open()) throw runtime_error("Error: output read group file could not be opened.");
    write_read_groups(rg_f, rgs, sample, library, platform);
    rg_f.close();
}

void collect_rg_from_fq(char* in_fname, char* sample, char* library, char* out_rg_fname) {
    // collect read-groups
    set<string> rgs;
    ifstream fq_f(in_fname);
    if (!fq_f.is_open()) throw runtime_error("Error: input FASTQ file could not be opened.");
    while (true) {
        fastq::entry x;
        if (!fastq::read_entry(fq_f, x)) break;
    
        // infer read-group
        string rg = infer_read_group_illumina18(x.qname);
        rgs.insert(rg);
    }
    fq_f.close();

    // write all read groups
    ofstream rg_f(out_rg_fname);
    if (!rg_f.is_open()) throw runtime_error("Error: output read group file could not be opened.");
    write_read_groups(rg_f, rgs, sample, library, platform);
    rg_f.close();
}

void tag_sam_with_rg(char* in_fname, char* rg_fname, char* out_sam_fname) {
    ifstream rg_f(rg_fname);
    if (!rg_f.is_open()) throw runtime_error("Error: input read group file could not be opened.");
    map<string, string> rgs;
    read_read_groups(rg_f, rgs);
    rg_f.close();

    ifstream in_f(in_fname);
    if (!in_f.is_open()) throw runtime_error("Error: input SAM file could not be opened.");
    ofstream out_f(out_sam_fname);
    if (!out_f.is_open()) throw runtime_error("Error: output SAM file could not be opened.");

    string line;

    while (true) {
        getline(in_f, line);

        if (line.empty()) break;

        if (line[0] == '@') {
            // skip RG header line but copy other header lines verbatim
            if (line.find("@RG") != 0) {
                out_f << line << endl;
            }
            continue;
        } else {
            break;
        }
    }

    // write read-group header
    write_read_groups(out_f, rgs);
    
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
 * Utility programs.
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
