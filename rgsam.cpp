#include <fstream>
#include <iostream>
#include <string>
#include <stdexcept>
#include <cstring>
#include <list>
#include <set>
#include <map>

#include "rgsam/arg.hpp"
#include "rgsam/fastq.hpp"
#include "rgsam/sam.hpp"
#include "rgsam/string.hpp"
#include "rgsam/file.hpp"

using namespace std;

const char* rgsam_version = "0.1";


namespace file_format {
    enum Format {
        SAM,
        FASTQ
    };
}

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


void collect_rg_from_sam(const char* format, const char* in_fname, const char* sample, const char* library, const char* platform, const char* out_rg_fname) {
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
        sam::raw_entry x; 
        if (!sam::extract_raw_entry(line, x)) break;
        string rg;
        infer_read_group(format, sam::get_qname_from_core(x.core), rg);
        rgs.insert(rg);
    }
    sam_f.close();

    // write all read groups
    ofstream rg_f(out_rg_fname);
    sam::write_read_groups(rg_f, rgs, sample, library, platform);
    rg_f << "@CO\t" << "QF:" << format << endl;
    rg_f.close();
}

void collect_rg_from_fq(const char* format, const char* in_fname, const char* sample, const char* library, const char* platform, const char* out_rg_fname) {
    // collect read-groups
    set<string> rgs;
    ifstream fq_f(in_fname);
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
    sam::write_read_groups(rg_f, rgs, sample, library, platform);
    rg_f << "@CO\t" << "QF:" << format << endl;
    rg_f.close();
}

void tag_sam_with_rg(const char* format, const char* in_fname, const char* rg_fname, const char* out_sam_fname) {
    ifstream rg_f(rg_fname);
    map<string, string> rgs;
    sam::read_read_groups(rg_f, rgs);
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
                out_f << line << endl;
            }
            continue;
        } else {
            break;
        }
    }

    // write read-group header
    sam::write_read_groups(out_f, rgs);
    out_f << "@CO\t" << "QF:" << format << endl;
    
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
 * Read-group identifier (ID) and platform unit (PU) are inferred from read 
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

    argc -= (argc > 0); argv += (argc > 0);  // skip program name if present

    if (argc < 1) {
        cerr << "usage: rgsam [command]" << endl << endl;
        cerr << "commands:" << endl;
        cerr << "  collect    collect read-group information from SAM or FASTQ file" << endl;
        cerr << "  split      split SAM or FASTQ file based on read-group" << endl;
        cerr << "  tag        tag reads in SAM file with read-group field" << endl;
        cout << "  qnames     list supported read name formats" << endl;
        cout << "  version    print version" << endl;
        return 1;
    }

    if (strcmp(argv[0], "collect") == 0) {

        --argc; ++argv;  // skip command

        enum optionIndex { UNKNOWN, HELP, INPUT, OUTPUT, FORMAT, QNFORMAT, SAMPLE, LIBRARY, PLATFORM };
        const option::Descriptor usage[] =
        {
          { UNKNOWN, 0, "", "", Arg::None, "usage: rgsam collect [options]\n\noptions:" },
          { INPUT, 0, "i", "input", Arg::InFile,     "  --input     SAM file" },
          { OUTPUT, 0, "o", "output", Arg::OutFile,  "  --output    read-group header file" },
          { FORMAT, 0, "f", "format", Arg::Some,     "  --format    input file format [sam, fastq]" },
          { QNFORMAT, 0, "q", "qnformat", Arg::Some, "  --qnformat  read name format" },
          { SAMPLE, 0, "s", "sample", Arg::Some,     "  --sample    sample name" },
          { LIBRARY, 0, "l", "library", Arg::Some,   "  --library   library name" },
          { PLATFORM, 0, "p", "plaform", Arg::Some,  "  --platform  sequencing platform [default: illumina]" },
          { HELP, 0, "h", "help", Arg::None,         "  --help      print usage and exit" },
          { 0, 0, 0, 0, 0, 0 }
        };

        option::Stats stats(usage, argc, argv);
        option::Option options[stats.options_max], buffer[stats.buffer_max];
        option::Parser parse(usage, argc, argv, options, buffer);

        if (parse.error()) return 1;

        if (argc == 0) {
            option::printUsage(cerr, usage);
            return 1;
        }

        if (options[HELP]) {
            option::printUsage(cerr, usage);
            return 0;
        }

        const char* qnformat;
        if (options[QNFORMAT].arg == NULL) {
            cerr << "Warning: read name format is not specified; assume `illumina-1.8`" << endl;
            qnformat = "illumina-1.8"; } else {
            qnformat = options[QNFORMAT].arg;
        }

        const char* input;
        if (options[INPUT].arg == NULL || strcmp(options[INPUT].arg, "-") == 0) {
            cerr << "Info: reading from stdin" << endl;
            input = "/dev/stdin";
        } else {
            input = options[INPUT].arg;
        }

        const char* output;
        if (options[OUTPUT].arg == NULL || strcmp(options[OUTPUT].arg, "-") == 0) {
            cerr << "Info: writing to stdout" << endl;
            output = "/dev/stdout";
        } else {
            output = options[OUTPUT].arg;
        }

        const char* sample = options[SAMPLE].arg;
        if (sample == NULL) {
            if (options[INPUT].arg == NULL) {
                cerr << "Error: sample name must be specified if input name is not specified" << endl;
                return 1;
            }
            string stem;
            get_file_stem(options[INPUT].arg, stem);
            sample = stem.c_str();
        }
        
        const char* library = options[LIBRARY].arg;
        if (library == NULL) {
            library = sample;
        }

        const char* platform;
        if (options[PLATFORM].arg == NULL) {
            platform = "illumina";
        } else {
            platform = options[PLATFORM].arg;
        }

        enum file_format::Format format = file_format::SAM;
        if (options[FORMAT].arg == NULL) {
            // attempt to infer format from input file name
            if (options[INPUT].arg == NULL) {
                cerr << "Warning: input file name and format are not specified; assume `SAM`" << endl;
            } else {
                string ext;
                get_file_ext(options[INPUT].arg, ext);
                to_lower(ext);
                if (ext.empty()) {
                    cerr << "Warning: input file format could not be inferred; assume `SAM`" << endl;
                } else if (ext == "fastq" || ext == "fq") {
                    format = file_format::FASTQ;
                } else if (ext == "sam") {
                    format = file_format::SAM;
                } else {
                    cerr << "Error: unsupported input file format `" << ext << '`' << endl;
                }
            }
        } else {
            string format_str = options[FORMAT].arg;
            if (format_str == "fastq" || format_str == "fq") {
                format = file_format::FASTQ;
            } else if (format_str == "sam") {
                format = file_format::SAM;
            } else {
                cerr << "Error: unsupported input file format `" << format_str << '`' << endl;
            }
        }

        switch (format) {
            case file_format::SAM:
                collect_rg_from_sam(qnformat, input, sample, library, platform, output);
                break;
            case file_format::FASTQ:
                collect_rg_from_fq(qnformat, input, sample, library, platform, output);
                break;
        }

    } else if (strcmp(argv[0], "split") == 0) {

        --argc; ++argv;  // skip command

        // TODO this command would be more practical if compression is
        // possible on output file

        cerr << "Unimplemented command: " << argv[0] << endl;

    } else if (strcmp(argv[0], "tag") == 0) {

        --argc; ++argv;  // skip command

        enum optionIndex { UNKNOWN, HELP, INPUT, INPUT_RG, OUTPUT, QNFORMAT };
        const option::Descriptor usage[] =
        {
            { UNKNOWN, 0, "", "", Arg::None, "usage: rgsam collect [options]\n\noptions:" },
            { INPUT, 0, "i", "input", Arg::InFile,     "  --input     input SAM file" },
            { INPUT_RG, 0, "r", "rg", Arg::InFile,     "  --rg        input read-group header file" },
            { OUTPUT, 0, "o", "output", Arg::OutFile,  "  --output    output SAM file" },
            { QNFORMAT, 0, "q", "qnformat", Arg::Some, "  --qnformat  read name format" },
            { HELP, 0, "h", "help", Arg::None,         "  --help      print usage and exit" },
            { 0, 0, 0, 0, 0, 0 }
        };

        option::Stats stats(usage, argc, argv);
        option::Option options[stats.options_max], buffer[stats.buffer_max];
        option::Parser parse(usage, argc, argv, options, buffer);

        if (parse.error()) return 1;

        if (argc == 0) {
            option::printUsage(cerr, usage);
            return 1;
        }

        if (options[HELP]) {
            option::printUsage(cerr, usage);
            return 0;
        }

        const char* qnformat;
        if (options[QNFORMAT].arg == NULL) {
            cerr << "Warning: read name format is not specified; assume `illumina-1.8`" << endl;
            qnformat = "illumina-1.8"; } else {
            qnformat = options[QNFORMAT].arg;
        }

        const char* input;
        if (options[INPUT].arg == NULL || strcmp(options[INPUT].arg, "-") == 0) {
            cerr << "Info: reading from stdin" << endl;
            input = "/dev/stdin";
        } else {
            input = options[INPUT].arg;
        }

        const char* input_rg = options[INPUT_RG].arg;
        if (input_rg == NULL) {
            cerr << "Error: read-group header file is required and can be"
                    "acquired by running `rgsam collect`" << endl;
            return 1;
        }

        const char* output;
        if (options[OUTPUT].arg == NULL || strcmp(options[OUTPUT].arg, "-") == 0) {
            cerr << "Info: writing to stdout" << endl;
            output = "/dev/stdout";
        } else {
            output = options[OUTPUT].arg;
        }

        tag_sam_with_rg(qnformat, input, input_rg, output);

    } else if (strcmp(argv[0], "qnames") == 0) {

        --argc; ++argv;  // skip command

        cout << "{" << endl
             << "  \"illumina-1.0\": {" << endl
             << "    \"format\": \"@{flowcell}-{instrument}:{lane}:{tile}:{x}:{y}#{sample}/{pair}\"," << endl
             << "    \"example\": \"@HWUSI-EAS100R:6:73:941:1973#0/1\"" << endl
             << "  }," << endl
             << "  \"illumina-1.8\": {" << endl
             << "    \"format\": \"@{flowcell}:{run}:{flowcell}:{lane}:{tile}:{x}:{y}\"," << endl
             << "    \"exaample\": \"@EAS139:136:FC706VJ:2:2104:15343:197393\"" << endl 
             << "  }," << endl
             << "  \"broad-1.0\": {" << endl
             << "    \"format\": \"@{flowcell,5}:{barcode}:{lane}:{tile}:{x}:{y}\"," << endl
             << "    \"example\": \"@H0164ALXX140820:2:1101:10003:23460\"" << endl
             << "  }" << endl
             << "}" << endl;

    } else if (strcmp(argv[0], "version") == 0) {

        cout << rgsam_version << endl;

    } else {

        cerr << "Unsupported command: " << argv[0] << endl;
        return 1;

    }

    return 0;
}

