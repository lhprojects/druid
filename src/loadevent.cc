
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//    Select the collections to be displayed.                                //
//    Stauts inherited from last event.
//    //
//																			 //
//    Date:   03 Dec 2010                                                    //
//    Author: Manqi Ruan, Daniel Jeans (LLR)                                 //
//                                                                           //
//    Last Modified: 17, Nov 2011, cleaning                                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TGLAnnotation.h>
#include <TGLViewer.h>

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "EVENT/LCCollection.h"
#include "EVENT/LCEvent.h"
#include "EVENT/MCParticle.h"
#include "GlobalDefs.hh"
#include "EventNavigator.hh"
#include "TruthHelper.h"
#include "IO/LCReader.h"
#include "MultiView.hh"
#include "TEveCompound.h"
#include "TEveElement.h"
#include "TEveManager.h"
#include "TEveBoxSet.h"
#include "TEveRGBAPalette.h"
#include "TEveRGBAPaletteOverlay.h"
#include "TEveText.h"
#include "TEveTextGL.h"
#include "TEveTrans.h"
#include "lcio.h"
#include "Options.h"

extern TEveCompound* allmarkers(int markercolorindex);
TEveCompound* marks = 0;

using namespace lcio;
using namespace std;

int last_event_id = -100;
int HenriCount;
LCEvent* evt;
LCEvent* lastevt;
// TEveText *t;

extern TEveManager* gEve;
bool extern FlagMultiView;
extern MultiView* gMultiView;
// MultiView * gMultiView;

extern LCReader* lcReader;
extern int GlobalRandomColorIndex;
extern int flagdetectortype;
extern bool flagslcio;
extern bool DefaultCollectionFlag;


std::map<float, int>
    randomColor;  // used to give an random color to each MCParticle
std::map<MCParticle*, int> OriginColor;  // used to give ... to Origin Particle

std::map<string, TEveElementList*> collectionClasses;

TEveRGBAPalette* p = new TEveRGBAPalette(0, 10);
TEveRGBAPaletteOverlay* po = new TEveRGBAPaletteOverlay(p, 0.15, 0.1, 0.8, 0.05);
TEveBoxSet* cal_shell = new TEveBoxSet("Shell");
bool cal_shell_shown = false;

// adding for Annotation

class TimeAnnotation : public TGLAnnotation {
 public:
  TimeAnnotation(TGLViewer* parent) : TGLAnnotation(parent, "ff", 0.1, 0.9) {}
  ~TimeAnnotation() {}

  void UpdateMyText(const char* x) {
    fText = x;
    fParent->RequestDraw();
  }
};

TimeAnnotation* ann;

int AnnoTime = 0;
int MarkerTime = 0;

/*
void OriginParticle()
{
        LCCollection *MCPColl = evt->getCollection("MCParticle");
        int NMCP = MCP->getNumberOfElements();
        std::vector<MCParticle* > candiMCP;
        std::vector<MCParticle* > candiOriP;
        for(int i = 0; i < NMCP; i++)
        {
                MCParticle * a_MCP =
dynamic_cast<MCParticle*>(MCPColl->getElementAt(i)); if(a_MCP->getPDG() == 92 &&
a_MCP->getParents().size() > 4 && a_MCP->getDaughters().size() > 33)
        }

}
*/

void load_event(int EventNum) {
  OriginColor.clear();

  HenriCount = 0;

  if (!flagslcio) {
    cout << "SLCIO file not available, Please check!" << endl;
    return;
  }

  GlobalRandomColorIndex = 0;

  cout << endl
       << "********************************************************************"
          "***********"
       << endl
       << endl;
  ;
  cout << "Start to display event " << EventNum << endl;
  cout << "LCIO Flag " << flagslcio << endl;

  std::cout << "Run number : " << gDisplayState.getRunNumber() << 
                ", Event number : " << EventNum << std::endl;

  try {

      if (!gDisplayState.isRunNumberSet())
      {
          EVENT::LCEvent *firstEvent = lcReader->readNextEvent();
          firstEvent = lcReader->readNextEvent(); // I dont know why need read twice

          int runNumber = firstEvent->getRunNumber();
          gDisplayState.setRunNumber(runNumber);
          std::cout << "    Reset Run Number: " << gDisplayState.getRunNumber() << std::endl;
          if (!gDisplayState.isEventNumberSet())
          {
              gDisplayState.setEventNumber(firstEvent->getEventNumber());
              EventNum = firstEvent->getEventNumber();
              std::cout << "    Reset Event Number: " << gDisplayState.getEventNumber() << std::endl;
          }
      }

      evt = nullptr;
      evt = lcReader->readEvent(gDisplayState.getRunNumber(), EventNum);

      if (!evt)
      {
          cout << "    Event Not Found!" << endl;
      }
      else
      {
          gDisplayState.setEventNumber(evt->getEventNumber());
          std::cout << "    Actual Event Number: " << gDisplayState.getEventNumber() << std::endl;

          string k =
              evt->getDetectorName();             // Initial Detector Type from slcio header;
          if (k == "ILD_00" || k == "ILD_o1_v05") // ILD00_AHCAL
          {
              flagdetectortype = 0;
          }
          else if (k == "ILD_00Dhcal" || k == "ILD_01pre00" ||
                   k == "ILD_o2_v05" || k == "ILD_o2_v06" || k == "CEPC_v0" ||
                   k == "CEPC_v1") // ILD00_DHCAL
          {
              flagdetectortype = 1;
          }
          else if (k.find("TB") != string::npos)
          {
              flagdetectortype = 2; // Simulated test beam events; currently only
                                    // considering about AHCAL Calibrations
          }
          else if (k == "sidloi3") // SID_LOI
          {
              flagdetectortype = 3;
          }
          else if (k == "clic_sid_cdr_b") // SID_CLIC
          {
              flagdetectortype = 5;
          }
          else
          {
              flagdetectortype =
                  10; // Risky... currently used as for the real data of CALICE TB
          }
          // As initial Value = -1, thus = 4 means has data file.

          cout << "    Event Located at : " << evt << endl;
          cout << endl
               << "****************************************************************"
                  "***************"
               << endl
               << endl;
          cout << "    Event Statistics: " << endl
               << endl;

          // Reset truth helper when loading new event

          load_collections(evt, "");

          last_event_id = EventNum;
          last_event_id = gDisplayState.getEventNumber();

          cout << "****************************************************************"
                  "*********************"
               << endl;
          cout << "   Display of event " << gDisplayState.getEventNumber() << " finished... " << endl
               << endl;
      }

    // 3D Annotations...

    TGLViewer* v = gEve->GetDefaultGLViewer();

    //	if(ann) v->DeleteOverlayAnnotations();
    if (AnnoTime == 0) {
      ann = new TimeAnnotation(v);
      AnnoTime = 1;
    }
    if (ann) {
      ann->SetAllowClose(0);
      ann->SetTextSize(0.03);
      ann->SetText(
          Form("DRUID, RunNum = %d, EventNum = %d", gDisplayState.getRunNumber(), gDisplayState.getEventNumber()));
    }

    if (MarkerTime == 0) {
      MarkerTime = 1;
      //		if(marks) marks->Destroy();
      marks = allmarkers(3);
      marks->SetMainColor(3);
      marks->SetRnrSelfChildren(0, 0);  // Hidden by default

      gEve->AddGlobalElement(marks);
      gMultiView->ImportGeomRPhi(marks);
      gMultiView->ImportGeomRhoZ(marks);
    }

  } catch (lcio::DataNotAvailableException zero) {
    cout << "    error... Event not found!" << endl;
  }
}

void load_collections(LCEvent* evt, string coltype) {
  if (flagslcio) {
    // collection types we do not draw explicitly for the moment...

    std::vector<std::string> ignoredTypes;
    // ignoredTypes.push_back(LCIO::LCRELATION);
    ignoredTypes.push_back(LCIO::LCFLOATVEC);
    ignoredTypes.push_back(LCIO::MCPARTICLE);
    ignoredTypes.push_back(LCIO::LCGENERICOBJECT);
    // ignoredTypes.push_back(LCIO::RECONSTRUCTEDPARTICLE);

    if (DefaultCollectionFlag) {  // Only Display Simulated Hits and
                                  // Reconstructed PFO (if exist)
      // ignoredTypes.push_back(LCIO::CLUSTER);
      ignoredTypes.push_back(LCIO::TRACK);
      ignoredTypes.push_back(LCIO::VERTEX);
      ignoredTypes.push_back(LCIO::TRACKERHIT);
      if (flagdetectortype != 10) {
        ignoredTypes.push_back(LCIO::CALORIMETERHIT);
      }
    }

    bool Draw = true;
    bool ChildDraw = true;
    bool MCPDraw = true;
    bool MCPChildDraw = true;

    std::map<string, int> DisplayFlag;
    std::map<string, int> CollectionScan;

    const std::vector<std::string> *strVec_ = evt->getCollectionNames();
    std::vector<std::string> collNames = *strVec_;
    std::sort(collNames.begin(), collNames.end());


    if (coltype == "") {
      for (std::map<string, TEveElementList*>::iterator AA =
               collectionClasses.begin();
           AA != collectionClasses.end(); AA++) {
        CollectionScan[AA->first] = 0;
      }
    }

    for (std::string const &collName : collNames) {
      LCCollection* col = evt->getCollection(collName);
      std::string ct = col->getTypeName();
      if (std::find(ignoredTypes.begin(), ignoredTypes.end(), ct) !=
          ignoredTypes.end())
        continue;

      if (coltype == "" || coltype == ct) {
        if (collectionClasses.find(ct) != collectionClasses.end() &&
            collectionClasses[ct]) {
          Draw = collectionClasses[ct]->GetRnrSelf();
          ChildDraw = collectionClasses[ct]->GetRnrChildren();

          for (TEveElement::List_i itt = collectionClasses[ct]->BeginChildren();
               itt != collectionClasses[ct]->EndChildren(); itt++) {
            std::string colname = (*itt)->GetElementName();
            if (DisplayFlag.find(colname) == DisplayFlag.end()) {
              DisplayFlag[colname] = (*itt)->GetRnrSelf();
            }
          }
          collectionClasses[ct]->DestroyElements();
          collectionClasses[ct]->Destroy();
        }
        collectionClasses[ct] = new TEveElementList();
        collectionClasses[ct]->SetRnrSelfChildren(Draw, ChildDraw);
        collectionClasses[ct]->SetRnrSelf(Draw);
        collectionClasses[ct]->SetName(ct.c_str());
        CollectionScan[ct] = 1;
      }
    }

    if (coltype == "") {
      for (std::map<string, TEveElementList*>::iterator AA =
               collectionClasses.begin();
           AA != collectionClasses.end(); AA++) {
        if (CollectionScan[AA->first] == 0 &&
            collectionClasses[AA->first]) {  //&& AA->second){		//Exist in last
                                             //evt but not this evt
          std::cout << "    Remove the collection exist in previos display "
                       "setting or last event but not this evt"
                    << std::endl;
          collectionClasses[AA->first]->DestroyElements();
          collectionClasses[AA->first]->Destroy();
          collectionClasses[AA->first] = 0;
        }
      }
    }

    bool FlagDraw = true;

    /*
       gMultiView->DestroyEventRPhi();
       gMultiView->DestroyEventRhoZ();
       */

    // Build MC Particles only when loading all collections or specifically MCParticle type
    if (coltype == "" || coltype == LCIO::MCPARTICLE)
    {
        std::vector<std::string> MCPartCollNames;
        for (std::string const &collName : collNames)
        {
            if (equals_any(collName, gOptions.coll_MCP_collections))
            {
                MCPartCollNames.push_back(collName);
            }
        }

        // Always rebuild MCParticles (unconditionally)
        gGUIManager._MCPList = BuildMCParticles(evt, MCPartCollNames);
        gEve->AddElement(gGUIManager._MCPList);

        if (FlagMultiView)
        {
            gMultiView->ImportEventRPhi(gGUIManager._MCPList);
            gMultiView->ImportEventRhoZ(gGUIManager._MCPList);
        }
        gTruthHelper.m_MCPCollNames = MCPartCollNames;
        
        for(std::string const &name : MCPartCollNames)
        {
            std::cout << " MCParticle collection for TruthHelper: " << name << std::endl;
        }
        gTruthHelper.ResetMCTruth(evt);
    }

    if (coltype == "" || coltype == LCIO::SIMCALORIMETERHIT)
    {
        gGUIManager._SimCaloHitBoxes.clear();
        gGUIManager._SimCaloHitType.clear();
        gGUIManager.m_CollNames.clear();
    }

    for (std::string const &collName : collNames) {
      LCCollection* col = evt->getCollection(collName);

      std::string CollHead(collName, 0, 5);
      string ct = col->getTypeName();

      if (col->getNumberOfElements() == 0) continue;

      if (DisplayFlag.find(collName) == DisplayFlag.end()) {
        // if ( ct==LCIO::SIMCALORIMETERHIT ||
        if (flagdetectortype == 10) {
          FlagDraw = true;
        } else {
          FlagDraw = false;
        }
      } else
        FlagDraw = DisplayFlag[collName];

      if (coltype != "" && coltype != ct) continue;

      if (equals_any(collName, gOptions.coll_MCP_collections))
      {
        // Already displayed in BuildMCParticles
        cout << "  Collection <" << collName
             << "  already displayed in MCParticles form" << endl
             << endl;
      } else if (ct == LCIO::MCPARTICLE) {
        cout << "  Collection <" << collName
             << "  skiped" << endl
             << endl;
      } else if (find(ignoredTypes.begin(), ignoredTypes.end(), ct) !=
                 ignoredTypes.end()) {
        cout << "  Collection <" << collName << "> will not be displayed. " << endl
             << endl;
      } else if (ct == LCIO::VERTEX) {
        TEveElementList* temp = Vertex(col, collName);
        temp->SetRnrSelfChildren(FlagDraw, FlagDraw);
        collectionClasses[ct]->AddElement(temp);
      } else if (ct == LCIO::RAWCALORIMETERHIT) {
        TEveElementList* temp = CaloHits(col, collName);
        temp->SetRnrSelfChildren(FlagDraw, FlagDraw);
        collectionClasses[ct]->AddElement(temp);
      } else if (ct == LCIO::CALORIMETERHIT) {
          if (ends_with_any(collName, gOptions.coll_caloHit_filterOutSuffixes))
          {
              std::cout << "  Collection <" << collName
                   << "> filtered out" << endl
                   << endl;
              continue;
          }

        TEveElementList* temp = CaloHits(col, collName);
        temp->SetRnrSelfChildren(FlagDraw, FlagDraw);
        collectionClasses[ct]->AddElement(temp);
      } else if (ct == LCIO::SIMCALORIMETERHIT) {
          if (ends_with_any(collName, gOptions.coll_simCaloHit_filterOutSuffixes))
          {
              std::cout << "  SimCalorimeterHit Collection <" << collName
                   << "> filtered out" << endl
                   << endl;
              continue;
          }
        TEveElementList* temp = createSimCaloHits(col, collName);
        temp->SetRnrSelfChildren(FlagDraw, FlagDraw);
        collectionClasses[ct]->AddElement(temp);
      } else if (ct == LCIO::TRACKERHIT) {
        TEveElementList* temp = TrackerHits(col, collName);
        temp->SetRnrSelfChildren(FlagDraw, FlagDraw);
        collectionClasses[ct]->AddElement(temp);
      } else if (ct == LCIO::SIMTRACKERHIT) {
        TEveElementList* temp = TrackerHits(col, collName);
        temp->SetRnrSelfChildren(FlagDraw, FlagDraw);
        collectionClasses[ct]->AddElement(temp);
      } else if (ct == LCIO::TRACK) {
        TEveElementList* temp = TrackAssignedHits(col, collName);
        temp->SetRnrSelfChildren(FlagDraw, FlagDraw);
        collectionClasses[ct]->AddElement(temp);
      } else if (ct == LCIO::RECONSTRUCTEDPARTICLE) {
        if (collName == "PandoraPFOs" || collName == "PandoraPFANewPFOs" ||
            CollHead == "Arbor") {  // Supposed to be modified if user needed...
          TEveElementList* temp = BuildPFOs(col, collName);
          temp->SetRnrSelfChildren(FlagDraw, FlagDraw);
          collectionClasses[ct]->AddElement(temp);
        } else if (collName == "Durham_4Jets" || collName == "Durham_6Jets") {
          TEveElementList* temp = RecoJets(col, collName);
          temp->SetRnrSelfChildren(FlagDraw, FlagDraw);
          collectionClasses[ct]->AddElement(temp);
        } else {
          cout << "   " << collName
               << " is currently considering as Reconstructed Particles... "
                  "will skip "
               << endl;
        }
      } else if (ct == LCIO::LCRELATION &&
                 (CollHead == "Henri" || CollHead == "InitH" ||
                  CollHead == "InitE" || CollHead == "Links")) {
        HenriCount++;
        TEveCompound* temp = ConnectTrees(col, collName);
        temp->SetRnrSelfChildren(FlagDraw, FlagDraw);
        collectionClasses[ct]->AddElement(temp);
      } else if (ct == LCIO::CLUSTER) {
        TEveElementList* temp = ClusterHits(col, collName);
        temp->SetRnrSelfChildren(FlagDraw, FlagDraw);
        collectionClasses[ct]->AddElement(temp);
      } else {
        cout << "  Unknown collection type " << col->getTypeName()
             << " for collection " << collName << endl
             << endl;
      }
    }

    if (coltype == "") {
      for (std::map<string, TEveElementList*>::iterator ff =
               collectionClasses.begin();
           ff != collectionClasses.end(); ff++) {
        // std::cout<<ff->first<<std::endl;
        if (ff->second) {
          gEve->AddElement(ff->second);
          if (FlagMultiView) {
            gMultiView->ImportEventRPhi(ff->second);
            gMultiView->ImportEventRhoZ(ff->second);
          }
        }
      }
    } else {
      if (collectionClasses.find(coltype) != collectionClasses.end() &&
          collectionClasses[coltype]) {
        gEve->AddElement(collectionClasses[coltype]);

        // Yuzhi Che: date 2023-06-03, try to add an overlay color palatte.
        /*
        cal_shell->Reset(TEveBoxSet::kBT_FreeBox, kFALSE, 64);
        TGLViewer* v = gEve->GetDefaultGLViewer();
        v->RemoveOverlayElement(&po);
        p.SetupColorArray();
        v->AddOverlayElement(&po);
        v->UpdateScene();
        */
      }

      // std::cout<<collectionClasses[coltype]<<std::endl;

      if (FlagMultiView && collectionClasses[coltype]) {
        gMultiView->ImportEventRPhi(collectionClasses[coltype]);
        gMultiView->ImportEventRhoZ(collectionClasses[coltype]);
      }
    }
    
    // Update SimCalorimeterHit colors after ALL collections are loaded
    // Only for color scheme 7 (MCParticle colors)
    if ((coltype == LCIO::SIMCALORIMETERHIT || coltype == "") && 
        gGUIManager.HitColourType == 7) {
      UpdateHitColorsToMatchMCParticles();
    }
  } else {
    std::cout
        << " Slcio data file is not available! skip loading event information. "
        << std::endl;
  }

  return;
}
