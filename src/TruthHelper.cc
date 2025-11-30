#include "TruthHelper.h"
#include "GlobalDefs.hh"
#include "EVENT/CalorimeterHit.h"
#include "EVENT/SimCalorimeterHit.h"
#include "EVENT/SimTrackerHit.h"
#include "EVENT/LCRelation.h"
#include "EVENT/Cluster.h"
#include "lcio.h"
#include "IMPL/LCCollectionVec.h"
#include <iostream>
#include <algorithm>
#include <MLPFA/LCIO_MLPFA.h>
#include "Options.h"
#include "PFAA/JERDInternal.h"

using namespace MLPFA;

MLDataReader_LCIO mlpfa_reader;

MLInputData gMLPFAInputData;
MLPFA::MLMetaData gMLPFAMetaData;
PFAA::LCIOData gLCIOData;
PFAA::MCTruth_Ans mcTruthAns;

std::map<EVENT::Cluster*, EVENT::MCParticle*> mainMCParticles;
// in PFAA, I didnot record information of this level yet
std::map<EVENT::SimCalorimeterHit *, EVENT::MCParticle* > simCaloHitMainMCP;
std::map<EVENT::CalorimeterHit *, EVENT::MCParticle* > caloHitMainMCP;
std::map<EVENT::Track *, EVENT::MCParticle* > trackMainMCP;
std::map<EVENT::MCParticle*, MLMCPart*> mcpartIDs;
std::map<EVENT::CalorimeterHit*, EVENT::SimCalorimeterHit *> caloHitSimCaloHitMap;
std::map<EVENT::Track*, MLPFA::MLTrack*> m_track2MLTrack;
std::map<EVENT::LCObject*, LCObjectConnection> m_tracsterConnections;

TruthHelper gTruthHelper;

void init_tracsterConnections()
{
    for(int iobj = 0; iobj < gMLPFAInputData.m_objects.size(); ++iobj)
    {
        MLMothered* mlMothered = gMLPFAInputData.m_objects[iobj];
        void *originObj = ML_at(gMLPFAMetaData.m_objects, iobj);
        EVENT::LCObject *lcObj = static_cast<EVENT::LCObject *>(originObj);
        MLMothered* mlMother = mlMothered->getMother();
        LCObjectConnection conn;
        if(mlMother)
        {
            int motherIndex = mlMother->getFlatIndex();
            EVENT::LCObject *lcMom = (EVENT::LCObject *)ML_at(gMLPFAMetaData.m_objects, motherIndex);
            conn.m_mother = lcMom;
            conn.m_motherType = mlMother->isTrack() ? 1 : (mlMother->isCluster() ? 2 : 0);

            MLPFA::Vect3f motherPos = mlMothered->getMotherConnectionPoint();
            conn.m_motherX = motherPos.x;
            conn.m_motherY = motherPos.y;
            conn.m_motherZ = motherPos.z;
            MLPFA::Vect3f thisPos = mlMothered->getThisConnectionPoint();
            conn.m_daughterX = thisPos.x;
            conn.m_daughterY = thisPos.y;
            conn.m_daughterZ = thisPos.z;
            
        } else {
            conn.m_mother = nullptr;
            conn.m_motherType = 0;
            // For primary objects without mother, use their own start point as daughter position
            if(mlMothered->isCluster()) {
                MLPFA::MLCluster* cluster = static_cast<MLPFA::MLCluster*>(mlMothered);
                MLPFA::Vect3f startPos = cluster->getTrueStartPoint();
                conn.m_daughterX = startPos.x;
                conn.m_daughterY = startPos.y;
                conn.m_daughterZ = startPos.z;
            } else if(mlMothered->isTrack()) {
                MLPFA::MLTrack* track = static_cast<MLPFA::MLTrack*>(mlMothered);
                MLPFA::Vect3f startPos = track->m_startPoint;
                conn.m_daughterX = startPos.x;
                conn.m_daughterY = startPos.y;
                conn.m_daughterZ = startPos.z;
            }
        }
        m_tracsterConnections[lcObj] = conn;
    }
}


void TruthHelper::ResetMCTruth(EVENT::LCEvent *evt)
{
    gMLPFAInputData.reset();
    gMLPFAMetaData.reset();
    gLCIOData.clear();

    simCaloHitMainMCP.clear();
    caloHitMainMCP.clear();
    trackMainMCP.clear();
    mainMCParticles.clear();
    mcpartIDs.clear();
    caloHitSimCaloHitMap.clear();
    m_mcpTrackerHits.clear();
    m_track2MLTrack.clear();


    MLPFA::MLGeom::instance().setBField(gOptions.BField);

    if (gOptions.logLevel > 0) {
        PFAA::setLogLevel(gOptions.logLevel);
    }
    
    // Set TPC geometry if provided via command line options
    if (gOptions.tpc_innerRadius > 0) {
        MLPFA::MLGeom::instance().setTPCInnerRadius(gOptions.tpc_innerRadius);
    }
    if (gOptions.tpc_outerRadius > 0) {
        MLPFA::MLGeom::instance().setTPCOuterRadius(gOptions.tpc_outerRadius);
    }
    
    mlpfa_reader.m_MCParticleCollectionNames = m_MCPCollNames;
    
    // Get collection names from the event
    const std::vector<std::string> *collNames = evt->getCollectionNames();
    
    // Set up track collection names using the helper function
    mlpfa_reader.m_TrackCollectionNames = get_track_collections_to_use(evt);
    for(std::string const &name : mlpfa_reader.m_TrackCollectionNames)
    {
        try
        {
            EVENT::LCCollection *col = evt->getCollection(name);
            std::cout << "Found Track collection: " << name << " with " 
                      << col->getNumberOfElements() << " tracks" << std::endl;
        }
        catch(...) {}
    }
    
    // Set up cluster collection names using the helper function
    mlpfa_reader.m_ClusterCollectionName = get_cluster_collections_to_use(evt);
    for(std::string const &name : mlpfa_reader.m_ClusterCollectionName)
    {
        try
        {
            EVENT::LCCollection *col = evt->getCollection(name);
            std::cout << "Found Cluster collection: " << name << " with " 
                      << col->getNumberOfElements() << " clusters" << std::endl;
        }
        catch(...) {}
    }
    
    // Set up relation collection names for CaloHit-SimCaloHit links
    mlpfa_reader.m_RelCaloHitCollectionNames.clear();    
    // Set up relation collection names for TrackerHit-SimTrackerHit links
    mlpfa_reader.m_RelTrackCollectionNames.clear();
    
    // Try to find calorimeter relation collections in the event
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
                    else
                    {
                        // Try TrackerHit-SimTrackerHit relation
                        EVENT::TrackerHit *trackerHit = dynamic_cast<EVENT::TrackerHit*>(rel->getFrom());
                        EVENT::SimTrackerHit *simTrackerHit = dynamic_cast<EVENT::SimTrackerHit*>(rel->getTo());
                        
                        if(trackerHit && simTrackerHit)
                        {
                            mlpfa_reader.m_RelTrackCollectionNames.push_back(name);
                            std::cout << "Found TrackerHit-SimTrackerHit relation collection: " << name << std::endl;
                        }
                    }
                }
            }
        }
        catch(...) {}
    }

    mlpfa_reader.fillInputData(evt, gLCIOData, mcTruthAns, gMLPFAInputData, gMLPFAMetaData);

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
    for(int itrack = 0; itrack < (int)gLCIOData.m_tracks.size(); ++itrack)
    {
        EVENT::Track *track = ML_at(gLCIOData.m_tracks, itrack);
        MLTrack *mltrack = ML_at(gMLPFAInputData.m_tracks, itrack);
        m_track2MLTrack[track] = mltrack;
    }


    for (auto &pair : mcTruthAns.m_track2MainPart)
    {
        PFAA::Track *track = pair.first;
        PFAA::MCPart *mcp = pair.second;
        int index = mcp->getIndex();
        int trackIndex = track->getIndex();
        EVENT::MCParticle *mcPart = ML_at(gLCIOData.m_mcParts, index);
        EVENT::Track *lcTrack =  ML_at(gLCIOData.m_tracks, trackIndex);
        trackMainMCP[lcTrack] = mcPart;
    }

    // Read all SimTrackerHit collections and organize by MCParticle
    for(std::string const &name : *collNames)
    {
        try
        {
            EVENT::LCCollection *col = evt->getCollection(name);
            if(col->getTypeName() == lcio::LCIO::SIMTRACKERHIT)
            {
                int nHits = col->getNumberOfElements();
                std::cout << "Reading " << nHits << " SimTrackerHits from collection: " << name << std::endl;
                
                for(int i = 0; i < nHits; ++i)
                {
                    EVENT::SimTrackerHit *hit = dynamic_cast<EVENT::SimTrackerHit*>(col->getElementAt(i));
                    if(hit && hit->getMCParticle())
                    {
                        EVENT::MCParticle *mcp = hit->getMCParticle();
                        const double *pos = hit->getPosition();
                        float time = hit->getTime();
                        if(time < 1E-4) {
                            continue; // Skip hits with zero time
                        }
                        
                        // Check if this hit is too far from the last hit for this MCParticle
                        std::vector<TrackerHitInfo> &mcpHits = m_mcpTrackerHits[mcp];                        
                        if(0) {
                            // Print SimTrackerHit information
                            std::cout << "  SimTrackerHit " << i << " at (" 
                                    << pos[0] << ", " << pos[1] << ", " << pos[2] << ") mm, "
                                    << "Time=" << time << " ns, "
                                    << "MCP=" << mcp << std::endl;
                            
                        }
                        m_mcpTrackerHits[mcp].emplace_back(pos, time);
                    }
                }
            }
        }
        catch(...) {}
    }

    init_tracsterConnections();
    
    // Sort tracker hits by time for each MCParticle
    int totalHits = 0;
    for(auto &pair : m_mcpTrackerHits)
    {
        std::vector<TrackerHitInfo> &hits = pair.second;
        bool debug = false;

        std::sort(hits.begin(), hits.end(), 
                  [](const TrackerHitInfo &a, const TrackerHitInfo &b) { return a.time < b.time; });

        // Split hits into segments based on distance threshold (1000 mm = 100 cm)
        std::vector<TrackerHitInfo> filteredHits;
        const double maxDistance = 30.0; // mm
        const double minInterSegmentDistance = 500.0; // mm
        const int minSegmentSize = 5; // Minimum hits per segment (changed from 54)
        
        std::vector<TrackerHitInfo> currentSegment;
        
        for(size_t i = 0; i < hits.size(); ++i)
        {
            if(currentSegment.empty())
            {
                // Starting new segment - check distance from previous segment's last hit
                if(!filteredHits.empty())
                {
                    const TrackerHitInfo &lastFilteredHit = filteredHits.back();
                    double dx = hits[i].position[0] - lastFilteredHit.position[0];
                    double dy = hits[i].position[1] - lastFilteredHit.position[1];
                    double dz = hits[i].position[2] - lastFilteredHit.position[2];
                    double distance = std::sqrt(dx*dx + dy*dy + dz*dz);
                    
                    // Skip this hit if too far from previous segment
                    if(distance > minInterSegmentDistance)
                    {
                        if(debug) {
                            std::cout << "  Skipping hit - too far from previous segment: " << distance << " mm" << std::endl;
                        }
                        continue;
                    }
                }
                currentSegment.push_back(hits[i]);
            }
            else
            {
                // Calculate distance from last hit in current segment
                const TrackerHitInfo &lastHit = currentSegment.back();
                double dx = hits[i].position[0] - lastHit.position[0];
                double dy = hits[i].position[1] - lastHit.position[1];
                double dz = hits[i].position[2] - lastHit.position[2];
                double distance = std::sqrt(dx*dx + dy*dy + dz*dz);
                
                if(distance > maxDistance)
                {
                    if (debug)
                    {
                        std::cout << "  Distance " << distance << " mm too large between hits at times "
                                  << lastHit.time << " ns and " << hits[i].time << " ns. "
                                  << "Ending current segment of size " << currentSegment.size() << "." << std::endl;
                    }
                    // Distance too large - end current segment and start new one
                    if(currentSegment.size() >= minSegmentSize)
                    {
                        // Add current segment to filtered hits
                        filteredHits.insert(filteredHits.end(), currentSegment.begin(), currentSegment.end());
                    }
                    else
                    {
                        if(debug) {
                            std::cout << "  Ignoring segment with only " << currentSegment.size() 
                                    << " hits (< " << minSegmentSize << ")" << std::endl;
                        }
                    }
                    // Start new segment - will be empty and checked on next iteration
                    currentSegment.clear();
                }
                else
                {
                    // Distance OK - add to current segment
                    currentSegment.push_back(hits[i]);
                }
            }
        }
        
        // Don't forget the last segment
        if(!currentSegment.empty())
        {
            if(currentSegment.size() >= minSegmentSize)
            {
                filteredHits.insert(filteredHits.end(), currentSegment.begin(), currentSegment.end());
            }
            else
            {
                std::cout << "  Ignoring final segment with only " << currentSegment.size() 
                          << " hits (< " << minSegmentSize << ")" << std::endl;
            }
        }
        
        // Replace hits with filtered hits
        hits = filteredHits;
        
        totalHits += hits.size();
    }
    
    std::cout << "Organized " << totalHits << " SimTrackerHits for " 
              << m_mcpTrackerHits.size() << " MCParticles" << std::endl;

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

std::string TruthHelper::GetStringID(EVENT::Track *track)
{
    if(track == nullptr) {
        return "Track@Null";
    }
    auto iter = m_track2MLTrack.find(track);
    if(iter ==  m_track2MLTrack.end()) {
        std::cout << "Track not found in m_track2MLTrack map: " << track << std::endl;
        std::cout << "Track not found in m_track2MLTrack map: " << track << std::endl;
        std::cout << "Track not found in m_track2MLTrack map: " << track << std::endl;
        std::cout << "Track not found in m_track2MLTrack map: " << track << std::endl;
        std::cout << "Track not found in m_track2MLTrack map: " << track << std::endl;
        std::cout << "Track not found in m_track2MLTrack map: " << track << std::endl;
        std::cout << "Track not found in m_track2MLTrack map: " << track << std::endl;
        std::cout << "Track not found in m_track2MLTrack map: " << track << std::endl;

        return "Track@Unknown";
    }
    return iter->second->getStringID();
}

LCObjectConnection TruthHelper::GetTracsterConnection(EVENT::LCObject *tracster)
{
    auto iter = m_tracsterConnections.find(tracster);
    if(iter !=  m_tracsterConnections.end()) {
        return iter->second;
    }
    return LCObjectConnection();
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

EVENT::MCParticle *TruthHelper::GetMainMCP(EVENT::Track *track)
{
    if(!track) return nullptr;
    
    // Check cache first
    if(trackMainMCP.find(track) != trackMainMCP.end())
    {
        return trackMainMCP[track];
    }    
    return nullptr;
}

const std::vector<TrackerHitInfo>& TruthHelper::GetTrackerHits(EVENT::MCParticle *mcp)
{
    auto iter = m_mcpTrackerHits.find(mcp);
    if(iter != m_mcpTrackerHits.end())
    {
        return iter->second;  // Already sorted by time
    }
    return m_emptyTrackerHits;  // Empty vector reference if not found
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


