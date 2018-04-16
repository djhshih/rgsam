#ifndef _RGSAM_FILE_HPP_
#define _RGSAM_FILE_HPP_

#include <fstream>
#include <string>

inline bool file_exists (const char* fname) {
    std::ifstream f(fname);
    return f.good();
}

inline bool file_writable (const char* fname) {
    std::ofstream f(fname);
    return f.good();
}

/**
 * Get file stem.
 *
 * Assumes POSIX path.
 */
void get_file_stem(const std::string& fname, std::string& stem) {
    size_t start = fname.rfind("/");
    if (start == std::string::npos) {
        start = 0;
    } else {
        ++start;
    }

    size_t end = fname.rfind(".");
    stem = fname.substr(start, end - start);
}

/**
 * Get file extension.
 */
void get_file_ext(const std::string& fname, std::string& ext) {
    size_t start = fname.rfind(".");
    if (start == std::string::npos || start == fname.length() - 1) {
        ext = "";
        return;
    }
    ext = fname.substr(start + 1);
}

#endif  // _RGSAM_FILE_HPP_
