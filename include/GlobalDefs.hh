#ifndef GLOBALDEFS_H_
#define GLOBALDEFS_H_

#include <string>
#include "TEveVector.h"
#include "TEvePathMark.h"

//class TEveVector;
//class TEvePathMark;  
TEvePathMark* PathMarkEndTrack2DClu(TEveVector &Vtx, TEveVector &End);
TEvePathMark* PathMarkEndTrackDecay(TEveVector &Vtx, TEveVector &End);
bool IsNeutrino(int PID);

#include "Rtypes.h"

void BuildGeoGDMLRoot(std::string gdmlroot);
void make_gui();
void load_event(int EventNum);

#include "TEveCompound.h"
#include "TEveStraightLineSet.h"
#include "TEveBox.h"

#ifndef __CINT__
#include "lcio.h"
using namespace lcio;
#include "EVENT/LCEvent.h"
#include "EVENT/LCCollection.h"
using namespace EVENT;
#else
class LCEvent;
class LCCollection;
#endif

void load_collections(LCEvent* evt, std::string coltype);

TEveElementList* CaloHits( LCCollection* col, std::string hh);
TEveElementList* TrackerHits( LCCollection* col, std::string hh);
TEveElementList* TrackAssignedHits( LCCollection* col, std::string hh );
TEveElementList* ClusterHits( LCCollection* col, std::string hh);
TEveElementList* Vertex(LCCollection* col, std::string hh);
TEveElementList* BuildMCParticles( LCEvent* evt );
TEveElementList* BuildPFOs( LCCollection* col, std::string hh);
TEveCompound* ConnectTrees( LCCollection* col, std::string hh );
TEveCompound* RecoJets( LCCollection* col, std::string name);
TEveBox* BoxPhi( TVector3 &HitPos, TVector3 &Scale, int Type, int SegOrStaveNumber, float HitEnergy );

//------------------ For Prototype ------------------
TVector3 GetScEcalHitPos(int LayerIDs, int ChipIDs, int ChannelIDs);
TVector3 GetAhcalHitPos(int LayerIDs, int ChipIDs, int ChannelIDs);

#endif //GLOBALDEFS_H_
