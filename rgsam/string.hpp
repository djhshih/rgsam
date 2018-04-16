#ifndef _RGSAM_STRING_HPP_
#define _RGSAM_STRING_HPP_

#include <string>
#include <algorithm>

/**
 * Find the n-th occurence of a character after offset in a string.
 */
size_t find_in_string(const std::string& x, char c, size_t pos, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        pos = x.find(c, pos + 1);
        if (pos == std::string::npos) break;
    }
    return pos;
}

void to_lower(std::string& x) {
    std::transform(x.begin(), x.end(), x.begin(), ::tolower);
}

void to_upper(std::string& x) {
    std::transform(x.begin(), x.end(), x.begin(), ::toupper);
}

#endif  // _RGSAM_STRING_HPP_
