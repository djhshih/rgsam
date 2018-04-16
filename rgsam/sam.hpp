#ifndef _RGSAM_SAM_HPP_
#define _RGSAM_SAM_HPP_

#include <string>
#include <set>
#include <map>
#include <list>
#include <fstream>

#include "string.hpp"

namespace sam {

using namespace std;

struct opt_field {
    char tag[2];
    char type;
    string value;

    opt_field(const char _tag[], char _type, const string& _value)
    : type(_type), value(_value) {
        tag[0] = _tag[0];
        tag[1] = _tag[1];
    }

    opt_field(const string& x) {
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
opt_field read_group_field(const string& rg) {
    return opt_field("RG", 'Z', rg);
}


/**
 * SAM entry.
 */
struct entry {
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
    list<opt_field> opts;
};

/// SAM entry where core fields are not parsed
struct raw_entry {
    // raw core fields
    string core;
    // optional fields
    list<opt_field> opts;
};

const size_t n_core_fields = 11;
const char delim = '\t';

/**
 * Parse optional fields from a string.
 */
void parse_opts(const string& x, list<opt_field>& opts) {
    size_t pos = 0;
    while (true) {
        size_t end = x.find(delim, pos);
        opts.push_back(opt_field(x.substr(pos, end - pos)));
        if (end == string::npos) break;
        pos = end + 1;
    }
}

/**
 * Read a SAM entry but do not parse the core fields.
 */
bool extract_raw_entry(const string& line, raw_entry& x) {
    if (line.empty()) return false;

    // split entry into core string and parsed opts
    size_t pos = find_in_string(line, delim, 0, n_core_fields);
    x.core = line.substr(0, pos);
    parse_opts(line.substr(pos + 1), x.opts);
    
    return true;
}

/**
 * Write a raw SAM entry.
 */
bool write_raw_entry(ostream& f, const raw_entry& x) {
    f << x.core;
   
    for (list<opt_field>::const_iterator it = x.opts.begin(); it != x.opts.end(); ++it) {
        f << delim;
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

    bool operator()(const opt_field& x) const {
        return x.tag[0] == tag[0] && x.tag[1] == tag[1];
    }

    char tag[2];
};

/**
 * Replace an optional field.
 */
void replace_opt_field(list<opt_field>& opts, const opt_field& x) {
    opts.remove_if(tag_match(x.tag));
    opts.push_back(x);
}

/**
 * Get the read name from the string of the SAM core fields.
 */
string get_qname_from_core(const string& core) {
    return core.substr(0, core.find(delim));
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
        size_t end = line.find(delim, start);
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
        f   << "@RG" << delim
            << "ID:" << *it << delim
            << "PU:" << *it << delim
            << "SM:" << sample << delim
            << "LB:" << library << delim
            << "PL:" << platform << endl;
    }
}

}  // namespace sam

#endif  // _RGSAM_SAM_HPP_

