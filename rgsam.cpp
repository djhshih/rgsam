#include <fstream>
#include <iostream>
#include <string>
#include <stdexcept>
#include <cstring>
#include <list>
#include <set>
#include <map>

#include "rgsam/fastq.hpp"
#include "rgsam/sam.hpp"
#include "rgsam/string.hpp"

using namespace std;

const char* platform = "illumina";

/**
 * Infer read-group based on flowcell id and lane id,
 * assuming Illumina v1.0 read name format.
 */
void infer_read_group_illumina10(const string& qname, string& rg) {
    // extract flowcell
    size_t start = 0;
    size_t end = find_in_string(qname, '-', start, 1);
    if (end == string::npos) return;
    string flowcell = qname.substr(start, end - start);

    // extract lane
    start = find_in_string(qname, ':', end + 1, 1);
    if (start == string::npos) return;
    ++start;
    end = find_in_string(qname, ':', start, 1);
    if (end == string::npos) return;
    string lane = qname.substr(start, end - start);

    rg = flowcell + '_' + lane;
}

/**
 * Infer read-group based on flowcell id and lane id,
 * assuming Illumina v1.8 read name format.
 */
void infer_read_group_illumina18(const string& qname, string& rg) {
    // extract flowcell
    size_t start = find_in_string(qname, ':', 0, 2);
    if (start == string::npos) return;
    ++start;
    size_t end = find_in_string(qname, ':', start, 1);
    if (end == string::npos) return;
    string flowcell = qname.substr(start, end - start);

    // extract lane
    start = end + 1;
    end = find_in_string(qname, ':', start, 1);
    if (end == string::npos) return;
    string lane = qname.substr(start, end - start);

    rg = flowcell + '_' + lane;
}

/**
 * Infer read-group based on flowcell id and lane id,
 * assuming Broad v1.0 read name format.
 */
void infer_read_group_broad10(const string& qname, string& rg) {
    // extract flowcell
    string flowcell = qname.substr(0, 5);

    // extract lane
    size_t start = find_in_string(qname, ':', 5, 1);
    if (start == string::npos) return;
    ++start;
    size_t end = find_in_string(qname, ':', start, 1);
    if (end == string::npos) return;
    string lane = qname.substr(start, end - start);

    rg = flowcell + '_' + lane;
}

/**
 * Infer read-group based on flowcell id and lane id.
 */
void infer_read_group(const char* format, const string& qname, string& rg) {
    if (strcmp(format, "illumina-1.0") == 0) {
        infer_read_group_illumina10(qname, rg);
    } else if (strcmp(format, "illumina-1.8") == 0) {
        infer_read_group_illumina18(qname, rg);
    } else if (strcmp(format, "broad-1.0") == 0) {
        infer_read_group_broad10(qname, rg);
    } else {
        throw runtime_error("Unsupported read format");
    }
}


void collect_rg_from_sam(const char* format, const char* in_fname, const char* sample, const char* library, const char* out_rg_fname) {
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
        sam::raw_entry x; 
        if (!sam::extract_raw_entry(line, x)) break;
        string rg;
        infer_read_group(format, sam::get_qname_from_core(x.core), rg);
        rgs.insert(rg);
    }
    sam_f.close();

    // write all read groups
    ofstream rg_f(out_rg_fname);
    if (!rg_f.is_open()) throw runtime_error("Error: output read group file could not be opened.");
    sam::write_read_groups(rg_f, rgs, sample, library, platform);
    rg_f.close();
}

void collect_rg_from_fq(const char* format, const char* in_fname, const char* sample, const char* library, const char* out_rg_fname) {
    // collect read-groups
    set<string> rgs;
    ifstream fq_f(in_fname);
    if (!fq_f.is_open()) throw runtime_error("Error: input FASTQ file could not be opened.");
    while (true) {
        fastq::entry x;
        if (!fastq::read_entry(fq_f, x)) break;
    
        // infer read-group
        string rg;
        infer_read_group(format, x.qname.substr(1), rg);
        rgs.insert(rg);
    }
    fq_f.close();

    // write all read groups
    ofstream rg_f(out_rg_fname);
    if (!rg_f.is_open()) throw runtime_error("Error: output read group file could not be opened.");
    sam::write_read_groups(rg_f, rgs, sample, library, platform);
    rg_f.close();
}

void tag_sam_with_rg(const char* format, const char* in_fname, const char* rg_fname, const char* out_sam_fname) {
    ifstream rg_f(rg_fname);
    if (!rg_f.is_open()) throw runtime_error("Error: input read group file could not be opened.");
    map<string, string> rgs;
    sam::read_read_groups(rg_f, rgs);
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
    sam::write_read_groups(out_f, rgs);
    
    // process SAM entries
    while (true) {
        sam::raw_entry x; 
        if (!sam::extract_raw_entry(line, x)) break;

        string rg;
        infer_read_group(format, sam::get_qname_from_core(x.core), rg);

        if (rgs.find(rg) == rgs.end()) {
            cerr << "Warning: read group ID " << rg << " is not found in input read-groups" << endl;
        }
        
        // tag read with inferred read group
        sam::replace_opt_field(x.opts, sam::read_group_field(rg));

        // write modified SAM entry to out file
        sam::write_raw_entry(out_f, x);

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
        cerr << "usage: rgsam collect <format> <in.sam> <sample> <library> <rg.txt>" << endl;
        cout << "       rgsam collectfq <format> <in.fq> <sample> <library> <rg.txt>" << endl;
        cout << "       rgsam tag <format> <in.sam> <rg.txt> <out.sam>" << endl;
        cout << "       rgsam split <format> <in.sam> <sample> <library>" << endl;
        cout << "       rgsam splitfq <format> <in.fq> <sample> <library>" << endl;
        cout << "       rgsam formats" << endl;
        return 1;
    }

    if (strcmp(argv[1], "collect") == 0) {

        char* format = argv[2];
        char* in_fname = argv[3];
        char* sample = argv[4];
        char* library = argv[5];
        char* out_rg_fname = argv[6];
        collect_rg_from_sam(format, in_fname, sample, library, out_rg_fname);

    } else if (strcmp(argv[1], "collectfq") == 0) {

        char* format = argv[2];
        char* in_fname = argv[3];
        char* sample = argv[4];
        char* library = argv[5];
        char* out_rg_fname = argv[6];
        collect_rg_from_fq(format, in_fname, sample, library, out_rg_fname);

    } else if (strcmp(argv[1], "tag") == 0) {

        char* format = argv[2];
        char* in_fname = argv[3];
        char* rg_fname = argv[4];
        char* out_sam_fname = argv[5];
        tag_sam_with_rg(format, in_fname, rg_fname, out_sam_fname);

    } else if (strcmp(argv[1], "split") == 0) {

        char* format = argv[2];
        char* in_fname = argv[3];
        char* sample = argv[4];
        char* library = argv[5];

        cerr << "Unimplemented command: " << argv[1] << endl;

    } else if (strcmp(argv[1], "splitfq") == 0) {

        char* format = argv[2];
        char* in_fname = argv[3];
        char* sample = argv[4];
        char* library = argv[5];

        cerr << "Unimplemented command: " << argv[1] << endl;

    } else if (strcmp(argv[1], "formats") == 0) {

        cout << "illumina-1.0: @HWUSI-EAS100R:6:73:941:1973#0/1" << endl;
        cout << "illumina-1.8: @EAS139:136:FC706VJ:2:2104:15343:197393" << endl;
        cout << "broad-1.0: @H0164ALXX140820:2:1101:10003:23460" << endl;

    } else {

        cerr << "Unsupported command: " << argv[1] << endl;
        return 1;

    }

    return 0;
}

