#include "Options.h"
#include "Rtypes.h"
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

bool equals_any(const std::string &value, const std::vector<std::string> &candidates) {
    for (const auto &candidate : candidates) {
        if (value == candidate) {
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
                << "  -MCPtCut    MC Particle Pt Cut (GeV) (default=0.1)\n"
                << "  -BField     Magnetic field strength in Tesla (default=3.5)\n";
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

bool read_string(int argc, char** argv, int &index, char const *arg, std::string &value) {
    if (strcmp(argv[index], arg) != 0) return false;

    // Check if the next argument is available and not another option
    if (index + 1 < argc) {
        value = argv[index+1];
        index += 2;
        return true;
    } else {
        std::cerr << "Error: " << arg << " option requires a value.\n";
        exit(1);
    }
}

Options::Options() 
{
}

void Options::parse(int &argc, char** &argv) {
    static std::vector<char *> argv_;

    MCPtCut = 0.1;
    BField = 3.5;
    logLevel = 0;  // 0 means no debug logging
    tpc_innerRadius = -1.0;  // negative means not set, use default from geometry
    tpc_outerRadius = -1.0;
    tpc_halfZ = -1.0;
    draw_tpc_cylinder = false;
    tpc_innerBarrelColor = kCyan;
    tpc_outerBarrelColor = kMagenta;
    printHelp=false;
    printVersion=false;
    exe_path = "";  // Will be set from argv[0] if not provided via option
    coll_recoRelation_collections = "";  // Empty by default, must be explicitly set

    argv_.push_back(argv[0]);

    for (int i = 1; i < argc;) {
        if(read_double(argc, argv, i, "-MCPtCut", MCPtCut)) {
        } else if(read_double(argc, argv, i, "-BField", BField)) {
        } else if(read_int(argc, argv, i, "-logLevel", logLevel)) {
        } else if(read_double(argc, argv, i, "-tpc.innerRadius", tpc_innerRadius)) {
        } else if(read_double(argc, argv, i, "-tpc.outerRadius", tpc_outerRadius)) {
        } else if(read_double(argc, argv, i, "-tpc.halfZ", tpc_halfZ)) {
        } else if (read_bool(argc, argv, i, "-draw.tpc.cylinder", draw_tpc_cylinder)) {
        } else if(read_int(argc, argv, i, "-draw.tpc.innerBarrelColor", tpc_innerBarrelColor)) {
        } else if(read_int(argc, argv, i, "-draw.tpc.outerBarrelColor", tpc_outerBarrelColor)) {
        } else if (read_bool(argc, argv, i, "-h", printHelp)) {
        } else if (read_bool(argc, argv, i, "-v", printVersion)) {
        } else if(add_string(argc, argv, i, "-coll.caloHit.filterOutSuffix", coll_caloHit_filterOutSuffixes))  {
        } else if(add_string(argc, argv, i, "-coll.simCaloHit.filterOutSuffix", coll_simCaloHit_filterOutSuffixes))  {
        } else if(add_string(argc, argv, i, "-coll.MCP.add", coll_MCP_collections)) 
        {
        } else if(add_string(argc, argv, i, "-coll.track.add", coll_track_collections)) 
        {
        } else if(add_string(argc, argv, i, "-coll.cluster.add", coll_cluster_collections)) 
        {
        } else if(read_string(argc, argv, i, "-coll.recoRelation", coll_recoRelation_collections)) 
        {
        } else if(read_string(argc, argv, i, "--exe-path", exe_path)) 
        {
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

    if(coll_MCP_collections.size() == 0) {
        coll_MCP_collections.push_back("MCParticle");
        coll_MCP_collections.push_back("MCParticles");
    }

    // If exe_path was not set via --exe-path option, use argv[0]
    if(exe_path.empty()) {
        exe_path = argv_.at(0);
    }

    argc = argv_.size();
    argv = argv_.data();


}


Options gOptions;
