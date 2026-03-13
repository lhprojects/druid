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
#include "TEveArrow.h"
#include "TEveGeoShape.h"
#include "TGeoSphere.h"
#include "TVector3.h"
#include "lcio.h"
#include "EVENT/LCCollection.h"
#include "EVENT/LCEvent.h"
#include "EVENT/Track.h"
#include "EVENT/TrackerHit.h"
#include "GlobalDefs.hh"
#include "Options.h"
#include "MLPFA/MLPFA.h"
#include "TruthHelper.h"
#include "ConnectionTitleUtils.h"
#include <iostream>
#include <map>
#include <algorithm>

using namespace lcio;
using namespace EVENT;
using namespace std;

extern int GlobalRandomColorIndex;

std::tuple<TEveElementList*, TEveElementList*, TEveElementList*> BuildTracks(LCEvent* evt, LCCollection* col, string name)
{
    // Filter track collections based on command line options
    if(gOptions.coll_track_collections.size() > 0) {
        if(!equals_any(name, gOptions.coll_track_collections)) {
            cout << "  Skipping track collection: " << name << " (not in -coll.track.add list)" << endl;
            return std::make_tuple(nullptr, nullptr, nullptr);
        }
    }

    cout << "  Track collection: " << name.c_str()
         << ". Number of Tracks: " << col->getNumberOfElements() << endl;

    TEveTrackList* trackList = new TEveTrackList();
    trackList->SetName((name).c_str());
    trackList->SetMainColor(kBlue);

    TEveElementList* connectionList = new TEveElementList();
    connectionList->SetName((name + "_true_conn").c_str());
    connectionList->SetMainColor(kOrange);

    TEveElementList* recoConnectionList = new TEveElementList();
    recoConnectionList->SetName((name + "_reco_conn").c_str());
    recoConnectionList->SetMainColor(gRecoArrowColor);

    // Build a map of reco relations for fast lookup
    std::map<LCObject*, RecoRelationData> recoRelationMap = buildRecoRelationMap(evt, gOptions.coll_recoRelation_collections);

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

        // print out all tracker hit
        bool print_track = false;
        if (0) // example: print hits for track 16
        {
            print_track = false;
        }

        if (print_track)
        {
            std::string trackID = gTruthHelper.GetStringID(track);
            std::cout << "Track " << trackID << " has " << hits.size() << " hits:" << std::endl;
            for (int iHit = 0; iHit < static_cast<int>(hits.size()); iHit++)
            {
                const double *pos = hits[iHit]->getPosition();
                const double radius = sqrt(pos[0] * pos[0] + pos[1] * pos[1]);
                std::cout << "   " << iHit
                            << ": (x,y,z) = (" << pos[0] << ", " << pos[1] << ", " << pos[2] << ") " << " r = " << radius << " mm"
                            << " time: " << hits[iHit]->getTime()
                            << std::endl;
            }

            std::cout << "Global:" << "omega: " << track->getOmega()
                      << " tanLambda: " << track->getTanLambda()
                      << " phi: " << track->getPhi() << " d0: " << track->getD0()
                      << " z0: " << track->getZ0() << std::endl;

            auto print_trackState = [&](int loc, const char* locName)
             {
                 const TrackState *trackState = track->getTrackState(loc);
                 if (trackState) {
                     std::cout << locName << ":"
                               << " omega: " << trackState->getOmega()
                               << " tanLambda: " << trackState->getTanLambda()
                               << " phi: " << trackState->getPhi()
                               << " d0: " << trackState->getD0()
                               << " z0: " << trackState->getZ0() << std::endl;
                 } else {
                     std::cout << locName << ": No track state available" << std::endl;
                 }
             };

            print_trackState(TrackState::AtIP, "AtIP");
            print_trackState(TrackState::AtFirstHit, "AtFirstHit");
            print_trackState(TrackState::AtLastHit, "AtLastHit");
        }

        const TrackState *trackState = nullptr;
        double tanLambda = track->getTanLambda();
        if ((lastPos[2] - firstPos[2]) * tanLambda < 0)
        {
            std::swap(firstPos, lastPos);
            trackState = track->getTrackState(TrackState::AtLastHit);
            std::cout << "Track " << gTruthHelper.GetStringID(track) << ": Swapped first and last hit based on z and tanLambda sign." << std::endl;
        } else {
            trackState = track->getTrackState(TrackState::AtFirstHit);
        }
        // we need to use momentum at firstPos
        float const *referencePoint = trackState->getReferencePoint();
        MLPFA::Helix helix2(trackState->getPhi(),
                    trackState->getD0(),
                    trackState->getZ0(),
                    trackState->getOmega(),
                    trackState->getTanLambda(),
                    gOptions.BField);
        MLPFA::Vect3f mom = helix2.m_momentum;
        double const momentum = sqrt(mom.GetX() * mom.GetX() +
                              mom.GetY() * mom.GetY() +
                              mom.GetZ() * mom.GetZ());
        // Get track parameters
        double omega = trackState->getOmega();
        int const charge = (omega > 0) ? 1 : -1;

        if(print_track) {
            std::cout << "Track " << gTruthHelper.GetStringID(track) << ": Charge = " << charge
                      << ", Momentum = " << momentum << " GeV/c"
                      << ",\n First Hit = (" << firstPos[0] << ", " << firstPos[1] << ", " << firstPos[2] << ") mm"
                      << ",\n Last Hit = (" << lastPos[0] << ", " << lastPos[1] << ", " << lastPos[2] << ") mm"
                      << ",\n Reference Point = (" << referencePoint[0] << ", " << referencePoint[1] << ", " << referencePoint[2] << ") mm"
                      << std::endl;
        }

        // Convert to Eve units (cm) and create vectors
        TEveVector vtx(referencePoint[0] / 10.0, referencePoint[1] / 10.0, referencePoint[2] / 10.0);
        TEveVector end(lastPos[0] / 10.0, lastPos[1] / 10.0, lastPos[2] / 10.0);

        MLPFA::Vect3f endPointHelix;
        MLPFA::Vect3f endPoint0 = MLPFA::Vect3f(lastPos[0] - referencePoint[0],
                                  lastPos[1] - referencePoint[1],
                                  lastPos[2] - referencePoint[2]);
        MLPFA::Vect3f distance;
        helix2.getPointOnHelix(endPoint0, endPointHelix, distance);
        end = TEveVector((endPointHelix.GetX() + referencePoint[0]) / 10.0,
                         (endPointHelix.GetY() + referencePoint[1]) / 10.0,
                         (endPointHelix.GetZ() + referencePoint[2]) / 10.0);

        if (print_track)
        {
            std::cout << "End point on helix: (" << endPointHelix.GetX() + referencePoint[0]
                      << ", " << endPointHelix.GetY() + referencePoint[1]
                      << ", " << endPointHelix.GetZ() + referencePoint[2] << ") mm" << std::endl;
            std::cout << "distance" << ": (xy, z, 3D) = (" << distance.GetX() << ", " << distance.GetY() << ", " << distance.GetZ() << ") mm"
                      << std::endl;
        }

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

        // Draw mother-daughter connection arrow for truth
        LCObjectConnection conn = gTruthHelper.GetTracsterConnection(track);
        TEveArrow* connArrow = createConnectionArrow(conn);
        if (connArrow) {
            std::string motherID = "None";
            if (conn.m_mother) {
                if (conn.m_motherType == 1) {
                    motherID = gTruthHelper.GetStringID(dynamic_cast<Track*>(conn.m_mother));
                } else {
                    motherID = gTruthHelper.GetStringID(dynamic_cast<Cluster*>(conn.m_mother));
                }
            }
            connArrow->SetName(Form("Conn %s->%s p=%.3f", motherID.c_str(), trackID.c_str(), momentum));
            connArrow->SetPickable(kTRUE);
            connArrow->SetTitle(ConnectionTitleUtils::buildTruthConnectionTitle(motherID,
                                                                                 conn.m_motherType,
                                                                                 trackID).c_str());
            connectionList->AddElement(connArrow);
        }

        // Draw mother-daughter connection for reco using the pre-built map
        if(!recoRelationMap.empty())
        {
            auto it = recoRelationMap.find(track);
            if(it != recoRelationMap.end())
            {
                // Found the relation for this track
                LCObject *fromObj = it->second.fromObj;
                float weight = it->second.weight;

                LCObjectConnection recoConn = createConnectionFromRelation(fromObj, track);

                // Check if we should show this arrow
                bool showArrow = true;
                if (gOptions.show_reco_diff)
                {
                    // Only show if reco differs from truth
                    LCObjectConnection truthConn = gTruthHelper.GetTracsterConnection(track);
                    // Compare mother objects - show if they differ
                    showArrow = (recoConn.m_mother != truthConn.m_mother) || weight > 1;
                }

                if (showArrow)
                {
                    TEveArrow *recoArrow = createConnectionArrow(recoConn);
                    if (recoArrow)
                    {
                        auto [weightType, fixedWeight] = ConnectionTitleUtils::decodeRecoRelationWeight(weight);

                        // Use red for diff mode, yellow for normal reco
                        recoArrow->SetMainColor(gOptions.show_reco_diff ? gRecoDiffArrowColor : gRecoArrowColor);
                        recoArrow->SetMainAlpha(gRecoArrowAlpha);
                        recoArrow->SetPickable(kTRUE);
                        std::string motherID = "None";
                        if (recoConn.m_mother) {
                            if (recoConn.m_motherType == 1) {
                                motherID = gTruthHelper.GetStringID(dynamic_cast<Track*>(recoConn.m_mother));
                            } else {
                                motherID = gTruthHelper.GetStringID(dynamic_cast<Cluster*>(recoConn.m_mother));
                            }
                        }
                        recoArrow->SetName(Form("Conn %s->%s p=%.3f", motherID.c_str(), trackID.c_str(), momentum));
                        recoArrow->SetTitle(ConnectionTitleUtils::buildRecoConnectionTitle(motherID,
                                                                                             recoConn.m_motherType,
                                                                                             trackID,
                                                                                             fixedWeight).c_str());
                        recoConnectionList->AddElement(recoArrow);
                    } else {
                        // draw ball when arrow can't be created (no mother connection point)
                        // Calculate fixed weight for display (same logic as above)
                        auto [weightType, fixedWeight] = ConnectionTitleUtils::decodeRecoRelationWeight(weight);

                        TGeoSphere *sphere = new TGeoSphere(0, 2.0); // 2 cm radius ball
                        TEveGeoShape *marker = new TEveGeoShape(Form("NoConn %s", trackID.c_str()));
                        marker->SetShape(sphere);
                        marker->SetMainColor(kRed);
                        marker->SetMainAlpha(0.9);
                        marker->RefMainTrans().SetPos(0.1*recoConn.m_daughterX,
                                                       0.1*recoConn.m_daughterY,
                                                       0.1*recoConn.m_daughterZ);
                        marker->SetPickable(kTRUE);
                        marker->SetTitle(ConnectionTitleUtils::buildNoRecoConnectionTitle("Track",
                                                                                           trackID,
                                                                                           fixedWeight,
                                                                                           weightType).c_str());
                        recoConnectionList->AddElement(marker);
                    }
                }
            }
        }
    }

    trackList->MakeTracks();
    trackList->SetRnrSelfChildren(kTRUE, kTRUE);

    connectionList->SetRnrSelfChildren(kTRUE, kTRUE);
    recoConnectionList->SetRnrSelfChildren(kTRUE, kTRUE);

    return std::make_tuple(trackList, connectionList, recoConnectionList);
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

void loadTracks(LCEvent *evt, string coltype)
{
    gGUIManager._RecoTracks.clear();

    std::map<string, int> DisplayFlag;
    bool Draw = true;
    bool ChildDraw = true;

    extern std::map<string, TEveElementList*> collectionClasses;

    if (collectionClasses.find(coltype) != collectionClasses.end() &&
        collectionClasses[coltype])
    {
        Draw = collectionClasses[coltype]->GetRnrSelf();
        ChildDraw = collectionClasses[coltype]->GetRnrChildren();
        for (TEveElement::List_i itt = collectionClasses[coltype]->BeginChildren();
             itt != collectionClasses[coltype]->EndChildren(); itt++)
        {
            std::string colname = (*itt)->GetElementName();
            if (DisplayFlag.find(colname) == DisplayFlag.end())
            {
                DisplayFlag[colname] = (*itt)->GetRnrSelf();
            }
        }
        collectionClasses[coltype]->DestroyElements();
        collectionClasses[coltype]->Destroy();
    }

    collectionClasses[coltype] = new TEveElementList();
    collectionClasses[coltype]->SetRnrSelfChildren(Draw, ChildDraw);
    collectionClasses[coltype]->SetRnrSelf(Draw);
    collectionClasses[coltype]->SetName(coltype.c_str());

    // Get track collections to use from helper function
    std::vector<std::string> trackCollNames = get_track_collections_to_use(evt);

    for (std::string const &collName : trackCollNames)
    {
        LCCollection *col = nullptr;
        try
        {
            col = evt->getCollection(collName);
        }
        catch(...) {}
        auto result = BuildTracks(evt, col, collName);
        if (std::get<0>(result))
        {
            bool FlagDraw = DisplayFlag.find(collName) != DisplayFlag.end() ? DisplayFlag[collName] : true;
            std::get<0>(result)->SetRnrSelfChildren(FlagDraw, FlagDraw);
            collectionClasses[coltype]->AddElement(std::get<0>(result));

            if (std::get<1>(result)) {
                std::get<1>(result)->SetRnrSelfChildren(FlagDraw, FlagDraw);
                collectionClasses[coltype]->AddElement(std::get<1>(result));
            }

            if (std::get<2>(result)) {
                std::get<2>(result)->SetRnrSelfChildren(FlagDraw, FlagDraw);
                collectionClasses[coltype]->AddElement(std::get<2>(result));
            }
        }
    }
    UpdateTrackColorsByMCParticle();
}

// Helper function to create arrow for mother-daughter connections
TEveArrow* createConnectionArrow(LCObjectConnection const &conn)
{
	if (conn.m_mother == nullptr) return nullptr;

	// Arrow from mother connection point to daughter connection point
	float dx = conn.m_daughterX - conn.m_motherX;
	float dy = conn.m_daughterY - conn.m_motherY;
	float dz = conn.m_daughterZ - conn.m_motherZ;
	float arrowLength = std::sqrt(dx*dx + dy*dy + dz*dz);

	// Make radius inversely proportional to length (as fraction of length)
	// Target: ~2mm absolute radius, so fraction = 2mm / length
	float targetRadius = 5.0;  // Target absolute radius in mm
	float tubeR = targetRadius / arrowLength;  // Fraction of arrow length
	float coneR = 2.0 * tubeR;  // Cone 2x larger than tube
	float coneL = 0.2;  // Cone length as fraction of arrow length

	TEveArrow* connArrow = new TEveArrow(0.1*dx, 0.1*dy, 0.1*dz,
	                                      0.1*conn.m_motherX, 0.1*conn.m_motherY, 0.1*conn.m_motherZ);
	connArrow->SetTubeR(tubeR);
	connArrow->SetConeR(coneR);
	connArrow->SetConeL(coneL);
	connArrow->SetMainColor(gTruthArrowColor);
	connArrow->SetMainAlpha(0.8);
	connArrow->SetPickable(kTRUE);
	connArrow->SetTitle(Form("Mother-Daughter Connection\nMother type: %s",
	                         conn.m_motherType == 1 ? "Track" : "Cluster"));
	return connArrow;
}
