///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//    Build TEve Tracks from EVENT::Track collection                        //
//    Visualizes track trajectories from first to last hit                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TEveManager.h"
#include "TEveElement.h"
#include "TEveTrack.h"
#include "TEveTrackPropagator.h"
#include "TEveVSDStructs.h"
#include "TEvePathMark.h"
#include "TVector3.h"
#include "lcio.h"
#include "EVENT/LCCollection.h"
#include "EVENT/Track.h"
#include "EVENT/TrackerHit.h"
#include "GlobalDefs.hh"
#include "Options.h"
#include "MLPFA/MLPFA.h"
#include "TruthHelper.h"
#include <iostream>
#include <map>

using namespace lcio;
using namespace EVENT;
using namespace std;

extern int GlobalRandomColorIndex;

TEveElementList* BuildTracks(LCCollection* col, string name)
{
    // Filter track collections based on command line options
    if(gOptions.coll_track_collections.size() > 0) {
        if(!equals_any(name, gOptions.coll_track_collections)) {
            cout << "  Skipping track collection: " << name << " (not in -coll.track.add list)" << endl;
            return nullptr;
        }
    }
    
    cout << "  Track collection: " << name.c_str() 
         << ". Number of Tracks: " << col->getNumberOfElements() << endl;

    TEveTrackList* trackList = new TEveTrackList();
    trackList->SetName(name.c_str());
    trackList->SetMainColor(kBlue);

    // Create track propagator for charged tracks
    TEveTrackPropagator* propsetCharged = new TEveTrackPropagator();
    propsetCharged->SetMagFieldObj(new TEveMagFieldDuo(350, -gOptions.BField, 2.0));
    propsetCharged->SetMaxR(2500);
    propsetCharged->SetMaxZ(5000);
    propsetCharged->SetMaxOrbs(1.2);

    int nTracks = col->getNumberOfElements();
    
    for(int itrk = 0; itrk < nTracks; itrk++)
    {
        Track* track = dynamic_cast<Track*>(col->getElementAt(itrk));
        if(!track) continue;

        TrackerHitVec const &hits = track->getTrackerHits();
        if(hits.size() < 2) continue; // Need at least 2 hits to draw a track

        // Get first and last hit positions
        const double* firstPos = hits.front()->getPosition();
        const double* lastPos = hits.back()->getPosition();


        // Get track parameters
        double omega = track->getOmega();
        double tanLambda = track->getTanLambda();
        double phi0 = track->getPhi();
        double d0 = track->getD0();
        double z0 = track->getZ0();
        
        // Calculate momentum using MLPFA Helix
        MLPFA::Helix helix(phi0, d0, z0, omega, tanLambda, gOptions.BField);
        MLPFA::Vect3f mom = helix.m_momentum;
        double momentum = sqrt(mom.GetX() * mom.GetX() + 
                              mom.GetY() * mom.GetY() + 
                              mom.GetZ() * mom.GetZ());
        
        int charge = (omega > 0) ? 1 : -1;

        // print out all tracker hit
        if (itrk == 21 || itrk == 22)
        {
            for (int iHit = 0; iHit < static_cast<int>(hits.size()); iHit++)
            {
                const double *pos = hits[iHit]->getPosition();
                const double radius = sqrt(pos[0]*pos[0] + pos[1]*pos[1]);
                std::cout << "   " << iHit
                          << ": (x,y,z) = (" << pos[0] << ", " << pos[1] << ", " << pos[2] << ") " << " r = " << radius << " mm"
                          << " time: " << hits[iHit]->getTime()
                          << std::endl;
            }
        }
        
        const TrackState *trackState = nullptr;
        if ((lastPos[2] - firstPos[2]) * mom.GetZ() < 0)
        {
            std::swap(firstPos, lastPos);
            trackState = track->getTrackState(TrackState::AtLastHit);
        } else {
            trackState = track->getTrackState(TrackState::AtFirstHit);
        }
        // we need to use momentum at firstPos
        MLPFA::Helix helix2(trackState->getPhi(),
                    trackState->getD0(),
                    trackState->getZ0(),
                    trackState->getOmega(),
                    trackState->getTanLambda(),
                    gOptions.BField);
        mom = helix2.m_momentum;




        // Convert to Eve units (cm) and create vectors
        TEveVector vtx(firstPos[0] / 10.0, firstPos[1] / 10.0, firstPos[2] / 10.0);
        TEveVector end(lastPos[0] / 10.0, lastPos[1] / 10.0, lastPos[2] / 10.0);
        // Create TEveRecTrack
        TEveRecTrack* recTrack = new TEveRecTrack();
        recTrack->fV.Set(vtx);
        recTrack->fP.Set(mom.GetX(), mom.GetY(), mom.GetZ());
        recTrack->fSign = charge;


        // Create TEveTrack
        TEveTrack* eveTrack = new TEveTrack(recTrack, propsetCharged);
        
        // Add end point
        TEvePathMark* pm = new TEvePathMark(TEvePathMark::kDecay);
        pm->fV.Set(end);
        eveTrack->AddPathMark(*pm);
        
        // Set track properties
        eveTrack->SetLineWidth(2);
        eveTrack->SetSmooth(kTRUE);
        
        // Set default color (will be updated by UpdateTrackColorsByMCParticle)
        eveTrack->SetLineColor(kGray);
        
        // Store track in map for later color updates
        gGUIManager._RecoTracks[track] = eveTrack;

        // Get track string ID for display
        std::string trackID = gTruthHelper.GetStringID(track);
        EVENT::MCParticle *mainMCP = gTruthHelper.GetMainMCP(track);
        std::string mcpID = gTruthHelper.GetStringID(mainMCP);
        double mcPx = mainMCP ? mainMCP->getMomentum()[0] : 0.0;
        double mcPy = mainMCP ? mainMCP->getMomentum()[1] : 0.0;
        double mcPz = mainMCP ? mainMCP->getMomentum()[2] : 0.0;
        double mcMomentum = sqrt(mcPx * mcPx + mcPy * mcPy + mcPz * mcPz);


        // Set track name and info
        eveTrack->SetName(Form("%s: p=%.2f GeV/c", trackID.c_str(), momentum));
        eveTrack->SetTitle(Form("%s\n"
                                "Charge: %d\n"
                                "Momentum: %.3f GeV/c, (%.2f, %.2f, %.2f)\n"
                                "MCP %s, Momentum: %.3f GeV/c, (%.2f, %.2f, %.2f)\n"
                                "First Hit: (%.2f, %.2f, %.2f) cm\n"
                                "Last Hit: (%.2f, %.2f, %.2f) cm\n"
                                "Tracker hits: %lu",
                                trackID.c_str(),
                                charge, 
                                momentum, mom.GetX(), mom.GetY(), mom.GetZ(),
                                mcpID.c_str(), mcMomentum, mcPx, mcPy, mcPz,
                                vtx.fX, vtx.fY, vtx.fZ,
                                end.fX, end.fY, end.fZ,
                                hits.size()));

        trackList->AddElement(eveTrack);
    }

    trackList->MakeTracks();  
    trackList->SetRnrSelfChildren(kTRUE, kTRUE);  
    return trackList;
}

// Function to update reconstructed Track colors to match MCParticle colors
void UpdateTrackColorsByMCParticle() {
    std::cout << "Updating reconstructed Track colors to match MCParticle colors..." << std::endl;
    
    int updatedCount = 0;
    int grayCount = 0;
    
    // Iterate through all stored reconstructed tracks
    for (auto& trackPair : gGUIManager._RecoTracks) {
        EVENT::Track* track = trackPair.first;
        TEveTrack* eveTrack = trackPair.second;
        
        if (!track || !eveTrack) continue;
        
        // Get the main MCParticle for this track
        EVENT::MCParticle* mainMCP = gTruthHelper.GetMainMCP(track);
        int newColor = kGray;
        bool mcpIsVisible = false;
        std::cout << "Track " << gTruthHelper.GetStringID(track) << " associated MCParticle: " 
                  << (mainMCP ? gTruthHelper.GetStringID(mainMCP) : "None") << std::endl;
        
        if (mainMCP && gGUIManager._MCPTracks.find(mainMCP) != gGUIManager._MCPTracks.end()) {
            TEveTrack* mcpTrack = gGUIManager._MCPTracks[mainMCP];
            
            // Check if the MCParticle track is visible
            bool trackVisible = mcpTrack && mcpTrack->GetRnrSelf();
            
            // Check if the particle's group is showing
            bool groupVisible = true;
            auto groupIt = gGUIManager._MCPGroups.find(mainMCP);
            if (groupIt != gGUIManager._MCPGroups.end()) {
                TEveElement* group = groupIt->second;
                if (group) {
                    groupVisible = group->GetRnrSelf() && group->GetRnrChildren();
                }
            }
            
            if (trackVisible && groupVisible) {
                mcpIsVisible = true;
                // Use the same color as the MCParticle track
                newColor = mcpTrack->GetLineColor();
                updatedCount++;
            }
        }
        
        if (!mcpIsVisible) {
            newColor = kGray;
            grayCount++;
        }
        
        // Update the track color
        eveTrack->SetLineColor(newColor);
    }
    
    std::cout << "Updated " << updatedCount << " tracks with MCP colors, " 
              << grayCount << " tracks set to gray" << std::endl;
}
