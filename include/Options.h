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
    bool draw_axis;
    double draw_camera_center_x;
    double draw_camera_center_y;
    double draw_camera_center_z;
    double draw_axis_length;
    bool printVersion;
    bool printHelp;
    std::vector<std::string> coll_caloHit_filterOutSuffixes;
    std::vector<std::string> coll_simCaloHit_filterOutSuffixes;
    std::vector<std::string> coll_MCP_collections;
    std::vector<std::string> coll_track_collections;
    std::vector<std::string> coll_cluster_collections;
    std::vector<std::string> coll_pfo_collections;  // Reconstructed particle collections
    std::string coll_recoRelation_collections;  // Reco relation collection name for cluster connections
    bool show_reco_diff;  // Show only connections where reco differs from truth
    std::string exe_path;  // Path to executable (from argv[0] or --exe-path option)
    
    Options();
    void parse(int &argc, char** &argv);

};

bool ends_with(const std::string &value, const std::string &ending);
bool ends_with_any(const std::string &value, const std::vector<std::string> &endings);
bool equals_any(const std::string &value, const std::vector<std::string> &candidates);

extern Options gOptions;


#endif // DRUID_ARGS_H

