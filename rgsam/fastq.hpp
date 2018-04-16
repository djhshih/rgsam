#ifndef _RGSAM_FASTQ_HPP_
#define _RGSAM_FASTQ_HPP_

#include <string>
#include <stdexcept>
#include <fstream>

namespace fastq {

using namespace std;

/**
 * Fastq entry.
 */
struct entry {
    /// read name
    string qname;
    /// read sequence
    string seq;
    /// read quality scores
    string qual;
};

/** 
 * Read one fastq entry from file.
 *
 * Assume that each sequence and each quality score entry is single-line.
 */
bool read_entry(istream& f, entry& x) {
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

/** 
 * Write one fastq entry to file.
 */
void write_entry(ostream& f, const entry& x) {
    f << x.qname << endl;
    f << x.seq << endl;
    f << '+' << endl;
    f << x.qual << endl;
}

}  // namespace fastq

#endif  // _RGSAM_FASTQ_HPP_

