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
#include "TEveBox.h"
#include "IMPL/LCTOOLS.h"
#include "EVENT/LCCollection.h"
#include "EVENT/SimTrackerHit.h"
#include "EVENT/SimCalorimeterHit.h"
#include "EVENT/CalorimeterHit.h"
#include "EVENT/MCParticle.h"
#include "EVENT/LCEvent.h"
#include <string>
#include <vector>
#include <cctype>

// String utility functions
// Helper function to check if string starts with a prefix
inline bool starts_with(const std::string& str, const std::string& prefix) {
    return str.size() >= prefix.size() && str.compare(0, prefix.size(), prefix) == 0;
}

// Helper function to check if string contains a substring
inline bool contains(const std::string& str, const std::string& substr) {
    return str.find(substr) != std::string::npos;
}

// Helper function to convert string to lowercase
inline std::string to_lower(const std::string& str) {
    std::string result = str;
    for (size_t i = 0; i < result.length(); ++i) {
        result[i] = std::tolower(static_cast<unsigned char>(result[i]));
    }
    return result;
}

// Forward declarations
class EventNavigator;

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
    std::map<EVENT::MCParticle*, TEveElement*> _MCPGroups;
    std::map<EVENT::MCParticle*, bool> _MCPShowing;
    std::map<std::string, int> _MCParticleDisplayFlag;


    std::map<TEveBox*, EVENT::SimCalorimeterHit*> _SimCaloHitBoxes;
    std::map<TEveBox*, EVENT::CalorimeterHit*> _CaloHitBoxes;
    std::vector<std::string> m_SimCaloHitCollNames;
    std::vector<std::string> m_CaloHitCollNames;
    std::map<TEveBox*, int> _SimCaloHitType;
    std::map<TEveBox*, int> _CaloHitType;

    int HitColourType = 0;
    int PFOHitColourType = 0;
    int ClusterHitColourType = 1;

    std::map<EVENT::MCParticle*, TEveTrack*> getMCPTracks();


};

extern GUIManager gGUIManager;

// Global EventNavigator instance for signal connections
extern EventNavigator gEventNavigator;

// Function to update SimCalorimeterHit colors to match MCParticle visibility and colors
void UpdateHitColorsToMatchMCParticles();


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
void DrawTPCCylinder(double innerRadius, double outerRadius, double halfZ);
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


typedef EVENT::MCParticle MCParticle;
TEveElementList* createSimCaloHits(LCCollection* col, std::string name);
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
TEveBox* createBox( TVector3 &HitPos, TVector3 &Scale, int Type, int SegOrStaveNumber, float HitEnergy );


//------------------ For Prototype ------------------
TVector3 GetScEcalHitPos(int LayerIDs, int ChipIDs, int ChannelIDs);
TVector3 GetAhcalHitPos(int LayerIDs, int ChipIDs, int ChannelIDs);


#endif //GLOBALDEFS_H_
