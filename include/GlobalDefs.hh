#ifndef GLOBALDEFS_H_
#define GLOBALDEFS_H_

#include <string>
#include <map>
#include <vector>

#include "TEveVector.h"
#include "TEvePathMark.h"
#include "TGNumberEntry.h"
#include "TEveCompound.h"
#include "TEveElement.h"
#include "TEveTrack.h"
#include "IMPL/LCTOOLS.h"
#include "EVENT/LCCollection.h"
#include "EVENT/SimTrackerHit.h"
#include "EVENT/SimCalorimeterHit.h"
#include "EVENT/MCParticle.h"
#include "EVENT/LCEvent.h"
#include <string>
#include <vector>

struct EventState
{
    int eventNumber = -1;
    int runNumber = -1;

};
extern EventState gEventState;

struct GUIManager
{
    TGNumberEntry *_EventNumberEntry = nullptr;
    TEveElementList* _MCPList = nullptr;
    std::map<EVENT::MCParticle*, TEveTrack*> _MCPTracks;
    std::map<EVENT::MCParticle*, bool> _MCPShowing;
    std::map<std::string, int> _MCParticleDisplayFlag;


    std::map<EVENT::MCParticle*, TEveTrack*> getMCPTracks();
    bool isMCPShowing( EVENT::MCParticle* part );
    bool setIsMCPShowing( EVENT::MCParticle* part, bool isShowing );


};

extern GUIManager gGUIManager;


struct DisplayState
{
    int eventNumber = -1;
    bool eventNumberSet = false;
    int runNumber = -1;
    bool runNumberSet = false;

    void setRunNumber( int rn ) {
        runNumber = rn;
        runNumberSet = true;
    }
    bool isRunNumberSet() const {
        return runNumberSet;
    }
    int getRunNumber() const {
        return runNumber;
    }

    void setEventNumber( int en ) {
        eventNumber = en;
        eventNumberSet = true;
    }
    bool isEventNumberSet() const {
        return eventNumberSet;
    }
    int getEventNumber() const {
        return eventNumber;
    }
};


extern DisplayState gDisplayState;


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
TEveElementList* BuildMCParticles( LCEvent* evt, std::vector<std::string> const &mcpartcollnames);
TEveElementList* BuildPFOs( LCCollection* col, std::string hh);
TEveCompound* ConnectTrees( LCCollection* col, std::string hh );
TEveCompound* RecoJets( LCCollection* col, std::string name);
TEveBox* BoxPhi( TVector3 &HitPos, TVector3 &Scale, int Type, int SegOrStaveNumber, float HitEnergy );

//------------------ For Prototype ------------------
TVector3 GetScEcalHitPos(int LayerIDs, int ChipIDs, int ChannelIDs);
TVector3 GetAhcalHitPos(int LayerIDs, int ChipIDs, int ChannelIDs);


#endif //GLOBALDEFS_H_
