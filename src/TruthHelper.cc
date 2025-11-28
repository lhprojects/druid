#include "TruthHelper.h"
#include "EVENT/CalorimeterHit.h"
#include "EVENT/SimCalorimeterHit.h"
#include "EVENT/LCRelation.h"
#include "EVENT/Cluster.h"
#include "lcio.h"
#include "IMPL/LCCollectionVec.h"
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
std::map<EVENT::CalorimeterHit *, EVENT::MCParticle* > caloHitMainMCP;
std::map<EVENT::MCParticle*, MLMCPart*> mcpartIDs;
std::map<EVENT::CalorimeterHit*, EVENT::SimCalorimeterHit *> caloHitSimCaloHitMap;


TruthHelper gTruthHelper;

void TruthHelper::ResetMCTruth(EVENT::LCEvent *evt)
{
    gMLPFAInputData.reset();
    gMLPFAMetaData.reset();
    gLCIOData.clear();

    simCaloHitMainMCP.clear();
    caloHitMainMCP.clear();
    mainMCParticles.clear();
    mcpartIDs.clear();
    caloHitSimCaloHitMap.clear();


    MLPFA::MLGeom::instance().setBField(gOptions.BField);
    
    // Set TPC geometry if provided via command line options
    if (gOptions.tpc_innerRadius > 0) {
        MLPFA::MLGeom::instance().setTPCInnerRadius(gOptions.tpc_innerRadius);
    }
    if (gOptions.tpc_outerRadius > 0) {
        MLPFA::MLGeom::instance().setTPCOuterRadius(gOptions.tpc_outerRadius);
    }
    
    mlpfa_reader.m_MCParticleCollectionNames = m_MCPCollNames;
    
    // Set up relation collection names for CaloHit-SimCaloHit links
    mlpfa_reader.m_RelCaloHitCollectionNames.clear();
    
    // Try to find calorimeter relation collections in the event
    const std::vector<std::string> *collNames = evt->getCollectionNames();
    for(std::string const &name : *collNames)
    {
        try
        {
            EVENT::LCCollection *col = evt->getCollection(name);
            if(col->getTypeName() == "LCRelation" && col->getNumberOfElements() > 0)
            {
                // Test if this is a CaloHit-SimCaloHit relation by checking first element
                EVENT::LCRelation *rel = dynamic_cast<EVENT::LCRelation*>(col->getElementAt(0));
                if(rel)
                {
                    EVENT::CalorimeterHit *caloHit = dynamic_cast<EVENT::CalorimeterHit*>(rel->getFrom());
                    EVENT::SimCalorimeterHit *simHit = dynamic_cast<EVENT::SimCalorimeterHit*>(rel->getTo());
                    
                    // If both casts succeed, this is a CaloHit-SimCaloHit relation
                    if(caloHit && simHit)
                    {
                        mlpfa_reader.m_RelCaloHitCollectionNames.push_back(name);
                        std::cout << "Found CaloHit-SimCaloHit relation collection: " << name << std::endl;
                    }
                }
            }
        }
        catch(...) {}
    }

    mlpfa_reader.fillInputData(evt, gLCIOData, gMLPFAInputData, gMLPFAMetaData);

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
        EVENT::CalorimeterHit *recoHit = dynamic_cast<EVENT::CalorimeterHit*>(rel->getFrom());
        EVENT::SimCalorimeterHit *simHit = dynamic_cast<EVENT::SimCalorimeterHit*>(rel->getTo());        
        caloHitSimCaloHitMap[recoHit] = simHit;
    }
    
    std::cout << "Loaded " << gLCIOData.m_caloHitLinks.size() << " CaloHit-SimCaloHit relations" << std::endl;
    std::cout << "caloHitSimCaloHitMap size: " << caloHitSimCaloHitMap.size() << std::endl;

    std::cout << evt->getCollection("MCParticle")->getNumberOfElements() << " MCParticles loaded." << std::endl;
    std::cout << " gLCIOData.m_mcParts.size()" << gLCIOData.m_mcParts.size() << std::endl;
    for(int ipart = 0; ipart < (int)gLCIOData.m_mcParts.size(); ++ipart)
    {
        EVENT::MCParticle *mcPart = ML_at(gLCIOData.m_mcParts, ipart);
        MLMCPart *mlPart = ML_at(gMLPFAInputData.m_MCParts, ipart);
        mcpartIDs[mcPart] = mlPart;
    }

}

std::string TruthHelper::GetStringID(EVENT::MCParticle *mcp)
{
    if(mcp == nullptr) {
        return "MCP@Null";
    }
    auto iter = mcpartIDs.find(mcp);
    if(iter ==  mcpartIDs.end()) {
        std::cout << "MCParticle not found in mcpartIDs map: " << mcp << std::endl;
        std::cout << "MCParticle not found in mcpartIDs map: " << mcp << std::endl;
        std::cout << "MCParticle not found in mcpartIDs map: " << mcp << std::endl;
        std::cout << "MCParticle not found in mcpartIDs map: " << mcp << std::endl;
        std::cout << "MCParticle not found in mcpartIDs map: " << mcp << std::endl;
        std::cout << "MCParticle not found in mcpartIDs map: " << mcp << std::endl;
        std::cout << "MCParticle not found in mcpartIDs map: " << mcp << std::endl;
        std::cout << "MCParticle not found in mcpartIDs map: " << mcp << std::endl;

        return "MCP@Unknown";
    }
    return iter->second->getStringID();
}

EVENT::MCParticle *TruthHelper::GetMainMCP(EVENT::CalorimeterHit *caloHit)
{
    // Check cache first
    if(caloHitMainMCP.find(caloHit) != caloHitMainMCP.end())
    {
        return caloHitMainMCP[caloHit];
    }

    // CalorimeterHit is a reconstructed hit, need to find corresponding SimCalorimeterHit via relations
    EVENT::SimCalorimeterHit *simHit = nullptr;
    
    auto iter = caloHitSimCaloHitMap.find(caloHit);
    if(iter != caloHitSimCaloHitMap.end())
    {
        simHit = iter->second;
    }
    else
    {
        static int warnCount = 0;
        if(warnCount < 3)
        {
            std::cout << "Warning: CalorimeterHit " << caloHit << " not found in caloHitSimCaloHitMap (map size: " 
                      << caloHitSimCaloHitMap.size() << ")" << std::endl;
            warnCount++;
        }
    }
    
    EVENT::MCParticle *mainPart = nullptr;
    if(simHit != nullptr)
    {
        // Use the SimCalorimeterHit method to get the main MCParticle
        mainPart = GetMainMCP(simHit);
    }
    
    // Cache the result (even if nullptr)
    caloHitMainMCP[caloHit] = mainPart;
    
    return mainPart;
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


