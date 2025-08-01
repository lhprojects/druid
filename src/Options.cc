#include "Options.h"
#include <iostream>
#include <cstring>

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

void Options::parse(int &argc, char** &argv) {
    static char * argv_[100];
    int argc_ = 0;

    MCPtCut = 0.1;
    printHelp=false;
    printVersion=false;

    argv_[argc_++] = argv[0];

    for (int i = 1; i < argc;) {
        if(read_double(argc, argv, i, "-MCPtCut", MCPtCut)) {
        } else if (read_bool(argc, argv, i, "-h", printHelp)) {
        } else if (read_bool(argc, argv, i, "-v", printVersion)) {
        } else {
            if(argv[i][0] == '-') {
                std::cerr << "Error: Unknown option " << argv[i] << "\n";
                exit(1);
            }
            argv_[argc_++] = argv[i];
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
    argc = argc_;
    argv = argv_;


}


Options gOptions;
