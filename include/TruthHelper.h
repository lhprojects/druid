#ifndef TRUTHHELPER_H_
#define TRUTHHELPER_H_

#include "EVENT/MCParticle.h"
#include "EVENT/Cluster.h"
#include "EVENT/CalorimeterHit.h"
#include "EVENT/SimCalorimeterHit.h"
#include "EVENT/SimTrackerHit.h"
#include "EVENT/Track.h"
#include "EVENT/LCEvent.h"
#include <vector>
#include <map>

struct TrackerHitInfo
{
    double position[3];  // x, y, z position
    float time;          // time of the hit
    
    TrackerHitInfo(const double* pos, float t) : time(t)
    {
        position[0] = pos[0];
        position[1] = pos[1];
        position[2] = pos[2];
    }
};

struct LCObjectConnection
{

    EVENT::LCObject *m_mother = nullptr;
    int m_motherType = 0; // 1 = Trajck, 2 = Cluster
    double m_motherX;
    double m_motherY;
    double m_motherZ;
    double m_daughterX;
    double m_daughterY;
    double m_daughterZ;

};

struct TruthHelper
{
    std::vector<std::string> m_MCPCollNames;
    void ResetMCTruth(EVENT::LCEvent *evt);

    EVENT::MCParticle *GetMainMCP(EVENT::Cluster *cluster);
    EVENT::MCParticle *GetMainMCP(EVENT::CalorimeterHit *caloHit);
    EVENT::MCParticle *GetMainMCP(EVENT::SimCalorimeterHit *caloHit);
    EVENT::MCParticle *GetMainMCP(EVENT::Track *track);

    std::string GetStringID(EVENT::MCParticle *mcp);
    std::string GetStringID(EVENT::Track *track);

    LCObjectConnection GetTracsterConnection(EVENT::LCObject *tracster);

    // Get tracker hits for a given MCParticle, sorted by time
    const std::vector<TrackerHitInfo>& GetTrackerHits(EVENT::MCParticle *mcp);

private:
    std::map<EVENT::MCParticle*, std::vector<TrackerHitInfo>> m_mcpTrackerHits;
    std::vector<TrackerHitInfo> m_emptyTrackerHits;  // For returning empty reference

};  // struct TruthHelper

extern TruthHelper gTruthHelper;

#endif
