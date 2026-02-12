
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//    Build Clustered CaloHits into TEveBox,								 //
//    Color Specified to Index.												 //
//    Option to be add as trace back to Reco PID if needed		             //
//                                                                           //
//    Date:   03 Dec 2010                                                    //
//    Author: Manqi Ruan (LLR)                                               //
//                                                                           //
//    Last Modified, 17, Nov 2011, Cleanning                                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TRint.h"
#include "TEveBox.h"
#include "lcio.h"
#include "IMPL/LCTOOLS.h"
#include "EVENT/LCCollection.h"
#include "EVENT/Cluster.h"
#include "EVENT/LCRelation.h"
#include "EVENT/Track.h"
#include "TEveElement.h"
#include "TSystem.h"
#include "TEvePointSet.h"
#include "TVector3.h"
#include "TEveArrow.h"
#include "TEveGeoShape.h"
#include "TGeoSphere.h"
#include "TRandom3.h"

#include <map>
#include <limits>
#include "trajectory.h"
#include "fitting_root.h"
#include "geometry.h"
#include "Options.h"
#include "GlobalDefs.hh"
#include "point.h"
#include "segment3.h"
#include "GlobalDefs.hh"
#include "TruthHelper.h"

using namespace lcio;
using namespace std;
using namespace EVENT;

float ClusterHitSize = 1.0;
float gRecoArrowAlpha = 0.8;  // Alpha transparency for reco arrows
int gRecoArrowColor = kYellow;  // Color for reco arrows (yellow by default)
int gTruthArrowColor = kGreen;  // Color for truth arrows
int gRecoDiffArrowColor = kRed;  // Color for reco diff arrows (showing errors)
extern int flagdetectortype;
const float HCALBarrelLengthILD = 309.3;
const float HCALBarrelLengthSID = 320.0;
extern TEveBox *BoxPhi(TVector3 &HitPos, TVector3 &Scale, int Type, int SegOrStaveNumber, float HitEnergy);
extern int GlobalRandomColorIndex;
extern bool Flag_AttachTextToHit;
TRandom *r0 = new TRandom();

// Build a map of LCRelations for fast lookup
// Returns a map from daughter (to) object to mother (from) object
std::map<LCObject *, LCObject *> buildRecoRelationMap(LCEvent *evt, const std::string &relationName)
{
    std::map<LCObject *, LCObject *> recoRelationMap;

    if (relationName.empty())
    {
        return recoRelationMap;
    }

    try
    {
        LCCollection *rel = evt->getCollection(relationName);
        for (int i = 0; i < int(rel->getNumberOfElements()); ++i)
        {
            LCRelation *relation = dynamic_cast<LCRelation *>(rel->getElementAt(i));
            if (!relation)
                continue;

            LCObject *fromObj = relation->getFrom();
            LCObject *toObj = relation->getTo();
            recoRelationMap[fromObj] = toObj;
            //std::cout << "  Reco relation: fromObj=" << fromObj << " toObj=" << toObj << std::endl;
        }
        std::cout << "  Built reco relation map with " << recoRelationMap.size() << " entries from " << relationName
                  << " with " << rel->getNumberOfElements() << std::endl;
    }
    catch (...)
    {
        // Collection not found or wrong type - silently skip
        std::cout << "  Reco relation collection " << relationName << " not found or invalid" << std::endl;
    }

    return recoRelationMap;
}

// Helper function to get the position of the hit nearest to the center (origin)
static void getNearestHitToCenter(Cluster *cluster, float &x, float &y, float &z)
{
    CalorimeterHitVec hits = cluster->getCalorimeterHits();
    if (hits.empty())
    {
        // Fallback to cluster position if no hits
        const float *pos = cluster->getPosition();
        x = pos[0];
        y = pos[1];
        z = pos[2];
        return;
    }

    // Find hit with minimum distance to origin
    float minDist2 = std::numeric_limits<float>::max();
    int minIndex = 0;
    for (int i = 0; i < hits.size(); ++i)
    {
        const float *hitPos = hits[i]->getPosition();
        float dist2 = hitPos[0] * hitPos[0] + hitPos[1] * hitPos[1] + hitPos[2] * hitPos[2];
        if (dist2 < minDist2)
        {
            minDist2 = dist2;
            minIndex = i;
        }
    }

    const float *nearestHitPos = hits[minIndex]->getPosition();
    x = nearestHitPos[0];
    y = nearestHitPos[1];
    z = nearestHitPos[2];
}

// Create LCObjectConnection from a relation mother object and daughter object
LCObjectConnection createConnectionFromRelation(LCObject *motherObj, LCObject *daughterObj)
{
    LCObjectConnection conn;
    conn.m_mother = motherObj;

    // Get mother position (could be Track or Cluster)
    Track *motherTrack = dynamic_cast<Track *>(motherObj);
    Cluster *motherCluster = dynamic_cast<Cluster *>(motherObj);

    if (motherTrack)
    {
        conn.m_motherType = 1; // Track
        // Track end is at last hit
        const TrackState *trackState = motherTrack->getTrackState(TrackState::AtLastHit);
        if (trackState)
        {
            const float *refPoint = trackState->getReferencePoint();
            conn.m_motherX = refPoint[0];
            conn.m_motherY = refPoint[1];
            conn.m_motherZ = refPoint[2];
        }
    }
    else if (motherCluster)
    {
        conn.m_motherType = 2; // Cluster
        float x, y, z;
        getNearestHitToCenter(motherCluster, x, y, z);
        conn.m_motherX = x;
        conn.m_motherY = y;
        conn.m_motherZ = z;
    }

    // Get daughter position (could be Cluster or Track)
    Cluster *daughterCluster = dynamic_cast<Cluster *>(daughterObj);
    Track *daughterTrack = dynamic_cast<Track *>(daughterObj);

    if (daughterCluster)
    {
        float x, y, z;
        getNearestHitToCenter(daughterCluster, x, y, z);
        conn.m_daughterX = x;
        conn.m_daughterY = y;
        conn.m_daughterZ = z;
    }
    else if (daughterTrack)
    {
        // Track begin is at first hit
        const TrackState *trackState = daughterTrack->getTrackState(TrackState::AtFirstHit);
        if (trackState)
        {
            const float *refPoint = trackState->getReferencePoint();
            conn.m_daughterX = refPoint[0];
            conn.m_daughterY = refPoint[1];
            conn.m_daughterZ = refPoint[2];
        }
    }

    return conn;
}

TEveElementList *ClusterHits(LCEvent *evt, LCCollection *col, string name)
{
    // Filter cluster collections based on command line options
    if (gOptions.coll_cluster_collections.size() > 0)
    {
        if (!equals_any(name, gOptions.coll_cluster_collections))
        {
            cout << "  Skipping cluster collection: " << name << " (not in -coll.cluster.add list)" << endl;
            return nullptr;
        }
    }

    std::cout << "  Cluster collection: " << name.c_str() << ". Number of Cluster: " << col->getNumberOfElements() << endl;
    TEveElementList *CaloCluster = new TEveElementList;
    CaloCluster->SetName(name.c_str());
    CaloCluster->SetMainColor(5);

    TEveElementList *truthConnectionList = new TEveElementList;
    truthConnectionList->SetName((name + "_true_conn").c_str());
    truthConnectionList->SetMainColor(kOrange);

    TEveElementList *recoConnectionList = new TEveElementList;
    recoConnectionList->SetName((name + "_reco_conn").c_str());
    recoConnectionList->SetMainColor(gRecoArrowColor);

    float HitEn, ClusterEnergy;
    HitEn = 0;
    ClusterEnergy = 0;
    int nCluster = col->getNumberOfElements();
    int ClusterPID = 0;
    float HBLength = 0;
    int LocalRandomIndex = 0;
    int CluSize = 0;
    TVector3 cluPos, cluBeginPos, cluEndPos, cluDir, cluDir_B, cluDir_E;

    if (flagdetectortype == 3)
    {
        HBLength = HCALBarrelLengthSID;
    }
    else if (flagdetectortype < 2)
    {
        HBLength = HCALBarrelLengthILD;
    }

    std::cout << "HCAL Barrel Length " << HBLength << std::endl;

    TVector3 HitScale(0.5 * ClusterHitSize, 0.5 * ClusterHitSize, 0.5 * ClusterHitSize);

    // Build a map of relations for fast lookup
    std::map<LCObject *, LCObject *> recoRelationMap = buildRecoRelationMap(evt, gOptions.coll_recoRelation_collections);

    for (int iC = 0; iC < nCluster; iC++)
    {
        TEveElementList *Recocluster = new TEveElementList;
        Recocluster->SetMainColor(iC);

        Cluster *acluster = dynamic_cast<Cluster *>(col->getElementAt(iC));
        cluPos = acluster->getPosition();
        float cluPosMag = cluPos.Mag();
        if (cluPosMag > 0)
        {
            cluDir = 10.0 / cluPosMag * cluPos;
        }
        else
        {
            cluDir.SetXYZ(0, 0, 0);
        }
        CalorimeterHitVec Hits = acluster->getCalorimeterHits();
        CluSize = Hits.size();
        ClusterEnergy = acluster->getEnergy();
        ClusterPID = acluster->getType();
        std::string cluStrID = gTruthHelper.GetStringID(acluster);
        Recocluster->SetName(Form("%s, En=%.3f", cluStrID.c_str(), ClusterEnergy));
        // std::cout<<"CLUSTERPIDTYPE "<<ClusterPID<<std::endl;

        LocalRandomIndex = int(100 * r0->Rndm(iC));

        if (true)
        {
            for (int j = 0; j < CluSize; j++)
            {
                TVector3 HitPosition = Hits[j]->getPosition();
                HitEn = Hits[j]->getEnergy();
                HitPosition *= 0.1;

                // Create a cube with size based on ClusterHitSize (in cm)
                float cubeHalfSize = 0.5 * ClusterHitSize;
                TEveBox *q = new TEveBox();
                q->SetVertex(0, HitPosition.X() - cubeHalfSize, HitPosition.Y() - cubeHalfSize, HitPosition.Z() - cubeHalfSize);
                q->SetVertex(1, HitPosition.X() + cubeHalfSize, HitPosition.Y() - cubeHalfSize, HitPosition.Z() - cubeHalfSize);
                q->SetVertex(2, HitPosition.X() + cubeHalfSize, HitPosition.Y() + cubeHalfSize, HitPosition.Z() - cubeHalfSize);
                q->SetVertex(3, HitPosition.X() - cubeHalfSize, HitPosition.Y() + cubeHalfSize, HitPosition.Z() - cubeHalfSize);
                q->SetVertex(4, HitPosition.X() - cubeHalfSize, HitPosition.Y() - cubeHalfSize, HitPosition.Z() + cubeHalfSize);
                q->SetVertex(5, HitPosition.X() + cubeHalfSize, HitPosition.Y() - cubeHalfSize, HitPosition.Z() + cubeHalfSize);
                q->SetVertex(6, HitPosition.X() + cubeHalfSize, HitPosition.Y() + cubeHalfSize, HitPosition.Z() + cubeHalfSize);
                q->SetVertex(7, HitPosition.X() - cubeHalfSize, HitPosition.Y() + cubeHalfSize, HitPosition.Z() + cubeHalfSize);

                if (gGUIManager.ClusterHitColourType == 0)
                {
                    q->SetMainColor(5);
                    q->SetLineColor(5);
                }
                else if (gGUIManager.ClusterHitColourType == 1)
                {
                    int colorindex = (iC * 11 + GlobalRandomColorIndex * 13 + LocalRandomIndex) % 50 + 51;
                    q->SetMainColor(colorindex);
                    q->SetLineColor(colorindex);
                }
                else if (gGUIManager.ClusterHitColourType == 2)
                { // Uniform Color: Red for PFO used
                    q->SetMainColor(2);
                    q->SetLineColor(2);
                }
                if (Flag_AttachTextToHit)
                {
                    q->SetTitle(Form("CluserHit%d, En = %.3f keV\n"
                                     "PosX = %.3f mm, PosY = %.3f mm, PosZ = %.3f mm\n"
                                     "Cluster %s\n"
                                     "Cluster Energy = %f GeV\n",
                                     j,
                                     HitEn * 1E6, 10 * HitPosition[0], 10 * HitPosition[1], 10 * HitPosition[2],
                                     cluStrID.c_str(), ClusterEnergy));
                }
                q->SetPickable(kTRUE);
                if (1)
                {
                    Recocluster->AddElement(q);
                }
            }
        }

        CaloCluster->AddElement(Recocluster);

        // Draw mother-daughter connection for truth
        if (true)
        {
            LCObjectConnection conn = gTruthHelper.GetTracsterConnection(acluster);
            //std::cout << "a cluster" << cluStrID << " from " << conn.m_mother << std::endl;
            //std::cout << conn.m_daughterX << " " << conn.m_daughterY << " " << conn.m_daughterZ << std::endl;
            //std::cout << conn.m_motherX << " " << conn.m_motherY << " " << conn.m_motherZ << std::endl;

            TEveArrow *connArrow = createConnectionArrow(conn);
            if (connArrow)
            {
                if (1)
                {
                    std::string strid = Form("Conn %s E=%.3f", cluStrID.c_str(), ClusterEnergy);
                    connArrow->SetName(strid.c_str());
                    truthConnectionList->AddElement(connArrow);
                }
            }
            else
            {
                if (1)
                {
                    // No mother connection - draw a sphere at cluster start position
                    TGeoSphere *sphere = new TGeoSphere(0, 2.0); // 2 cm radius sphere
                    TEveGeoShape *marker = new TEveGeoShape("Primary");
                    marker->SetShape(sphere);
                    marker->SetMainColor(kYellow);
                    marker->SetMainTransparency(20);
                    marker->RefMainTrans().SetPos(0.1 * conn.m_daughterX, 0.1 * conn.m_daughterY, 0.1 * conn.m_daughterZ);
                    std::string strid = Form("Primary cluster %s E=%.3f", cluStrID.c_str(), ClusterEnergy);
                    marker->SetName(strid.c_str());
                    truthConnectionList->AddElement(marker);
                }
            }
        }

        // Draw mother-daughter connection for reco using the pre-built map
        if (!recoRelationMap.empty())
        {
            auto it = recoRelationMap.find(acluster);

            if (it != recoRelationMap.end())
            {
                // Found the relation for this cluster
                LCObject *fromObj = it->second;
                LCObjectConnection recoConn = createConnectionFromRelation(fromObj, acluster);

                // Check if we should show this arrow
                bool showArrow = true;
                if (gOptions.show_reco_diff)
                {
                    // Only show if reco differs from truth
                    LCObjectConnection truthConn = gTruthHelper.GetTracsterConnection(acluster);
                    // Compare mother objects - show if they differ
                    showArrow = (recoConn.m_mother != truthConn.m_mother);
                }

                if (showArrow)
                {
                    TEveArrow *connArrow = createConnectionArrow(recoConn);
                    if (connArrow)
                    {
                        // Use red for diff mode, yellow for normal reco
                        connArrow->SetMainColor(gOptions.show_reco_diff ? gRecoDiffArrowColor : gRecoArrowColor);
                        connArrow->SetMainAlpha(gRecoArrowAlpha);
                        std::string strid = Form("Conn %s E=%.3f", cluStrID.c_str(), ClusterEnergy);
                        connArrow->SetName(strid.c_str());
                        recoConnectionList->AddElement(connArrow);
                    }
                    else
                    {
                        TGeoSphere *sphere = new TGeoSphere(0, 2.0); // 2 cm radius sphere
                        TEveGeoShape *marker = new TEveGeoShape("Primary_Reco");
                        marker->SetShape(sphere);
                        // Use red for diff mode, yellow for normal reco
                        marker->SetMainColor(gOptions.show_reco_diff ? gRecoDiffArrowColor : gRecoArrowColor);
                        marker->SetMainTransparency(20);
                        marker->RefMainTrans().SetPos(0.1 * recoConn.m_daughterX, 0.1 * recoConn.m_daughterY, 0.1 * recoConn.m_daughterZ);
                        std::string strid = Form("Primary cluster %s E=%.3f", cluStrID.c_str(), ClusterEnergy);
                        marker->SetName(strid.c_str());
                        recoConnectionList->AddElement(marker);
                    }
                }
            }
        }

    } // GlobalRandomColorIndex++;

    CaloCluster->AddElement(truthConnectionList);
    CaloCluster->AddElement(recoConnectionList);

    return CaloCluster;
}
