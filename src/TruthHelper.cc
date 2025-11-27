#include "TruthHelper.h"
#include "EVENT/CalorimeterHit.h"
#include "EVENT/SimCalorimeterHit.h"
#include "EVENT/Cluster.h"
#include <iostream>
#include <MLPFA/LCIO_MLPFA.h>
#include "Options.h"
using namespace MLPFA;

MLDataReader_LCIO mlpfa_reader;
MLInputData gMLPFAInputData;
MLPFA::MLMetaData gMLPFAMetaData;
PFAA::LCIOData gLCIOData;
std::map<EVENT::Cluster*, EVENT::MCParticle*> mainMCParticles;

// in PFAA, I didnot record information of this level yet
std::map<EVENT::SimCalorimeterHit *, EVENT::MCParticle* > simCaloHitMainMCP;
TruthHelper gTruthHelper;

void TruthHelper::ResetMCTruth(EVENT::LCEvent *evt)
{
    simCaloHitMainMCP.clear();
    mainMCParticles.clear();

    mlpfa_reader.fillInputData(evt, gLCIOData, gMLPFAInputData, gMLPFAMetaData);
    MLPFA::MLGeom::instance().setBField(gOptions.BField);

    for(int iobj = 0; iobj < gMLPFAInputData.m_objects.size(); ++iobj)
    {
        MLMothered* mlMothered = gMLPFAInputData.m_objects[iobj];
        if(mlMothered->isCluster())
        {
            MLCluster *mlCluster = static_cast<MLCluster*>(mlMothered);
            MLMCPart* mlPart = mlCluster->getTrueMainMCPart();
            int index = mlPart->getIndex();
            void *ptr = ML_at(gMLPFAMetaData.m_MCParts, index);
            EVENT::MCParticle* mcp = static_cast<EVENT::MCParticle*>(ptr);
            void *originCluster = ML_at(gMLPFAMetaData.m_objects, iobj);
            EVENT::Cluster *cluster = static_cast<EVENT::Cluster *>(originCluster);
            mainMCParticles[cluster] = mcp;
        }
    }
    for(int ilink = 0; ilink < gLCIOData.m_caloHitLinks.size(); ++ilink)
    {
        EVENT::LCRelation *rel = gLCIOData.m_caloHitLinks[ilink];
    }
}
EVENT::MCParticle *TruthHelper::GetMainMCP(EVENT::CalorimeterHit *caloHit)
{
    return nullptr;
}

EVENT::MCParticle *TruthHelper::GetMainMCP(EVENT::SimCalorimeterHit *caloHit)
{
    if(simCaloHitMainMCP.find(caloHit) != simCaloHitMainMCP.end())
    {
        return simCaloHitMainMCP[caloHit];
    }

    std::map<EVENT::MCParticle *, float> energyMap;
    for(int ic = 0; ic < caloHit->getNMCContributions(); ++ic)
    {
        EVENT::MCParticle *part = caloHit->getParticleCont(ic);
        float energy = caloHit->getEnergyCont(ic);
        energyMap[part] += energy;
    }

    EVENT::MCParticle *mainPart = nullptr;
    float maxEnergy = -1.0f;
    for(auto &it : energyMap)
    {
        if(it.second > maxEnergy)
        {
            maxEnergy = it.second;
            mainPart = it.first;
        }
    }
    simCaloHitMainMCP[caloHit] = mainPart;

    return mainPart;
}

EVENT::MCParticle *TruthHelper::GetMainMCP(EVENT::Cluster *cluster)
{
    return mainMCParticles[cluster];    
}

#include "../thirdparty/MLPFA/core/src/MLPFA.cc"
#include "../thirdparty/MLPFA/core/src/LCIO_MLPFA.cc"
#include "../thirdparty/MLPFA/core/src/MarlinUtil/ClusterShapes.cc"
#include "../thirdparty/MLPFA/core/src/Helix.cc"
#include "../thirdparty/MLPFA/core/src/MLInputData.cc"
#include "../thirdparty/MLPFA/core/src/MLResolution.cc"

#include "../thirdparty/MLPFA/thirdparty/PFAA/src/Implementations.cc"
#include "../thirdparty/MLPFA/thirdparty/PFAA/src/Traits_LCIO.cc"
#include "../thirdparty/MLPFA/thirdparty/PFAA/src/MCTruth_LCIO.cc"
#include "../thirdparty/MLPFA/thirdparty/PFAA/src/MCTruth.cc"


