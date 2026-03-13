#ifndef CONNECTION_TITLE_UTILS_H
#define CONNECTION_TITLE_UTILS_H

#include <iomanip>
#include <sstream>
#include <string>
#include <utility>

namespace ConnectionTitleUtils {

inline const char *getConnectionObjectType(int motherType)
{
    if (motherType == 1) return "Track";
    if (motherType == 2) return "Cluster";
    return "Unknown";
}

inline std::pair<std::string, float> decodeRecoRelationWeight(float weight)
{
    if (weight >= 0.f && weight <= 1.f) {
        return {"normal", weight};
    }
    if (weight >= 2.f && weight < 4.f) {
        return {"neutral fragment", weight - 2.f};
    }
    if (weight >= 4.f && weight < 6.f) {
        return {"charged fragment", weight - 4.f};
    }
    return {"excluded", weight - 6.f};
}

inline std::string formatWeight(float weight)
{
    std::ostringstream os;
    os << std::fixed << std::setprecision(3) << weight;
    return os.str();
}

inline std::string buildTruthConnectionTitle(const std::string &motherID,
                                             int motherType,
                                             const std::string &daughterID)
{
    std::ostringstream os;
    os << "Connection\n"
       << "Mode: truth\n"
       << "From: " << motherID << " (" << getConnectionObjectType(motherType) << ")\n"
       << "To: " << daughterID;
    return os.str();
}

inline std::string buildRecoConnectionTitle(const std::string &motherID,
                                            int motherType,
                                            const std::string &daughterID,
                                            float decodedWeight)
{
    std::ostringstream os;
    os << "Connection\n"
       << "Mode: reco\n"
       << "From: " << motherID << " (" << getConnectionObjectType(motherType) << ")\n"
       << "To: " << daughterID << "\n"
       << "Weight: " << formatWeight(decodedWeight);
    return os.str();
}

inline std::string buildNoRecoConnectionTitle(const std::string &label,
                                              const std::string &daughterID,
                                              float decodedWeight,
                                              const std::string &weightType)
{
    std::ostringstream os;
    os << "Connection\n"
       << "Mode: reco\n"
       << "From: None (Unknown)\n"
       << "Mother type: " << weightType << "\n"
       << "To: " << daughterID << "\n"
       << "Weight: " << formatWeight(decodedWeight) << "\n";
    return os.str();
}

inline std::string buildNoTruthConnectionTitle(const std::string &label,
                                               const std::string &daughterID)
{
    std::ostringstream os;
    os << "Connection\n"
       << "Mode: truth\n"
       << "From: None\n"
       << "To: " << daughterID << "\n";
    return os.str();
}
}

#endif