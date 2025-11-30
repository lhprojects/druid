#ifndef DRUID_ARGS_H
#define DRUID_ARGS_H

#include <string>
#include <vector>
struct Options
{

    double MCPtCut;
    double BField;
    int logLevel;
    double tpc_innerRadius;
    double tpc_outerRadius;
    double tpc_halfZ;
    bool draw_tpc_cylinder;
    int tpc_innerBarrelColor;
    int tpc_outerBarrelColor;
    bool printVersion;
    bool printHelp;
    std::vector<std::string> coll_caloHit_filterOutSuffixes;
    std::vector<std::string> coll_simCaloHit_filterOutSuffixes;
    std::vector<std::string> coll_MCP_collections;
    std::vector<std::string> coll_track_collections;
    
    Options();
    void parse(int &argc, char** &argv);

};

bool ends_with(const std::string &value, const std::string &ending);
bool ends_with_any(const std::string &value, const std::vector<std::string> &endings);
bool equals_any(const std::string &value, const std::vector<std::string> &candidates);

extern Options gOptions;


#endif // DRUID_ARGS_H

