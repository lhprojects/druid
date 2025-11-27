#ifndef TRUTHHELPER_H_
#define TRUTHHELPER_H_

#include "EVENT/MCParticle.h"
#include "EVENT/Cluster.h"
#include "EVENT/CalorimeterHit.h"
#include "EVENT/SimCalorimeterHit.h"
#include "EVENT/LCEvent.h"

struct TruthHelper
{
    std::vector<std::string> m_MCPCollNames;
    void ResetMCTruth(EVENT::LCEvent *evt);
    EVENT::MCParticle *GetMainMCP(EVENT::Cluster *cluster);
    EVENT::MCParticle *GetMainMCP(EVENT::CalorimeterHit *caloHit);
    EVENT::MCParticle *GetMainMCP(EVENT::SimCalorimeterHit *caloHit);

    std::string GetStringID(EVENT::MCParticle *mcp);

};  // struct TruthHelper

extern TruthHelper gTruthHelper;

#endif
