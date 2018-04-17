#ifndef _RGSAM_ARG_HPP_
#define _RGSAM_ARG_HPP_

#include <iostream>
#include <cstring>

#include "optionparser.hpp"
#include "file.hpp"

struct Arg: public option::Arg
{
    static option::ArgStatus Some(const option::Option& option, bool msg) {
        if (option.arg != NULL) return option::ARG_OK;
        if (msg) std::cerr << "Option `" << option.name << "` requires an argument" << std::endl;
        return option::ARG_ILLEGAL;
    }

    static option::ArgStatus InFile(const option::Option& option, bool msg) {
        if (option.arg != NULL) {
            if (std::strcmp(option.arg, "-") == 0 || file_exists(option.arg)) {
                return option::ARG_OK;
            }
        }
        if (msg) std::cerr << "Error: option `" << option.name << "` requires a valid input file" << std::endl;
        return option::ARG_ILLEGAL;
    }

    static option::ArgStatus OutFile(const option::Option& option, bool msg) {
        if (option.arg != NULL) {
            if (std::strcmp(option.arg, "-") == 0 || file_writable(option.arg)) {
                return option::ARG_OK;
            }
        }
        if (msg) std::cerr << "Error: option `" << option.name << "` requires a writable output path" << std::endl;
        return option::ARG_ILLEGAL;
    }
};

#endif  // _RGSAM_ARG_HPP_

