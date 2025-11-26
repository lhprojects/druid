#include "Options.h"
#include <iostream>
#include <cstring>
#include <vector>

bool ends_with(const std::string &value, const std::string &ending) {
    if (ending.size() > value.size()) return false;
    return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

bool ends_with_any(const std::string &value, const std::vector<std::string> &endings) {
    for (const auto &ending : endings) {
        if (ends_with(value, ending)) {
            return true;
        }
    }
    return false;
}


void printVersion() {
    std::cout << "Druid version 1.0.0\n";
}

void printHelp() {
    std::cout << "Usage: druid [options]\n"
                << "Options:\n"
                << "  -h   Show this help message\n"
                << "  -v    Show version information\n"
                << "  -MCPtCut    MC Particle Pt Cut (GeV) (default=0.1)\n";
}


bool read_double(int argc, char** argv, int &index, char const *arg, double &value) {
    if (strcmp(argv[index], arg) != 0) return false;

    // Check if the next argument is available and not another option
    if (index + 1 < argc) {
        try {
            value = std::stod(argv[index+1]);
            index += 2;
            return true;
        } catch (const std::invalid_argument&) {
            std::cerr << "Error: Invalid value for " << arg << ": " << argv[index] << "\n";
            exit(1);
        }
    } else {
        std::cerr << "Error: " << arg << " option requires a value.\n";
        exit(1);
    }
}
bool read_int(int argc, char** argv, int &index, char const *arg, int &value) {
    if (strcmp(argv[index], arg) != 0) return false;

    // Check if the next argument is available and not another option
    if (index + 1 < argc) {
        try {
            value = std::stoi(argv[index+1]);
            index += 2;
            return true;
        } catch (const std::invalid_argument&) {
            std::cerr << "Error: Invalid value for " << arg << ": " << argv[index+1] << "\n";
            exit(1);
        }
    } else {
        std::cerr << "Error: " << arg << " option requires a value.\n";
        exit(1);
    }
}

bool read_bool(int argc, char** argv, int &index, char const *arg, bool &value) {
    if (strcmp(argv[index], arg) != 0) return false;

    // If the argument is present, set the value to true
    value = true;
    index++;
    return true;
}

bool add_string(int argc, char** argv, int &index, char const *arg, std::vector<std::string> &value) {
    if (strcmp(argv[index], arg) != 0) return false;

    // Check if the next argument is available and not another option
    if (index + 1 < argc) {
        try {
            value.push_back(argv[index+1]);
            index += 2;
            return true;
        } catch (const std::invalid_argument&) {
            std::cerr << "Error: Invalid value for " << arg << ": " << argv[index+1] << "\n";
            exit(1);
        }
    } else {
        std::cerr << "Error: " << arg << " option requires a value.\n";
        exit(1);
    }
}

void Options::parse(int &argc, char** &argv) {
    static std::vector<char *> argv_;

    MCPtCut = 0.1;
    printHelp=false;
    printVersion=false;

    argv_.push_back(argv[0]);

    for (int i = 1; i < argc;) {
        if(read_double(argc, argv, i, "-MCPtCut", MCPtCut)) {
        } else if (read_bool(argc, argv, i, "-h", printHelp)) {
        } else if (read_bool(argc, argv, i, "-v", printVersion)) {
        } else if(add_string(argc, argv, i, "-coll.caloHit.filterOutSuffix", coll_caloHit_filterOutSuffixes)) {
        } else {
            if(argv[i][0] == '-') {
                std::cerr << "Error: Unknown option " << argv[i] << "\n";
                exit(1);
            }
            argv_.push_back(argv[i]);
            i++;
        }
    }
    if(argc == 1) {
        printHelp = true;
    }

    if(printHelp) {
        ::printHelp();
    }
    if(printVersion) {
        ::printVersion();
    }

    argc = argv_.size();
    argv = argv_.data();


}


Options gOptions;
