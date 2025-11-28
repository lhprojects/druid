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
#include <iostream>

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
    
    for(int i = 0; i < nTracks; i++)
    {
        Track* track = dynamic_cast<Track*>(col->getElementAt(i));
        if(!track) continue;

        TrackerHitVec const &hits = track->getTrackerHits();
        if(hits.size() < 2) continue; // Need at least 2 hits to draw a track

        // Get first and last hit positions
        const double* firstPos = hits.front()->getPosition();
        const double* lastPos = hits.back()->getPosition();

        // Convert to Eve units (cm) and create vectors
        TEveVector vtx(firstPos[0] / 10.0, firstPos[1] / 10.0, firstPos[2] / 10.0);
        TEveVector end(lastPos[0] / 10.0, lastPos[1] / 10.0, lastPos[2] / 10.0);

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
        
        // Color based on charge
        if(charge > 0) {
            eveTrack->SetLineColor(kRed);
        } else {
            eveTrack->SetLineColor(kBlue);
        }

        // Set track name and info
        eveTrack->SetName(Form("Track %d: p=%.2f GeV/c", i, momentum));
        eveTrack->SetTitle(Form("Track Collection: %s\n"
                                "Track #%d\n"
                                "Charge: %d\n"
                                "Momentum: %.3f GeV/c\n"
                                "First Hit: (%.2f, %.2f, %.2f) cm\n"
                                "Last Hit: (%.2f, %.2f, %.2f) cm\n"
                                "Omega: %.6f mm^-1\n"
                                "tan(lambda): %.4f\n"
                                "Tracker hits: %lu",
                                name.c_str(), i, charge, momentum,
                                vtx.fX, vtx.fY, vtx.fZ,
                                end.fX, end.fY, end.fZ,
                                omega, tanLambda, hits.size()));

        trackList->AddElement(eveTrack);
    }

    trackList->MakeTracks();  
    trackList->SetRnrSelfChildren(kTRUE, kTRUE);  

    return trackList;
}
