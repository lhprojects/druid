
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

extern TEveCompound* allmarkers(int markercolorindex);
TEveCompound* marks = 0;

using namespace lcio;
using namespace std;

int last_event_id = -100;
int HenriCount;
LCEvent* evt;
LCEvent* lastevt;
// TEveText *t;

extern int event_id;
extern TEveManager* gEve;
bool extern FlagMultiView;
extern MultiView* gMultiView;
// MultiView * gMultiView;

extern int runNumber;
extern LCReader* lcReader;
extern int GlobalRandomColorIndex;
extern int flagdetectortype;
extern bool flagslcio;
extern bool DefaultCollectionFlag;

TEveElementList* pMCParticle = 0;

std::map<float, int>
    randomColor;  // used to give an random color to each MCParticle
std::map<MCParticle*, int> OriginColor;  // used to give ... to Origin Particle
std::map<string, int> MCParticleDisplayFlag;

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
  cout << "    Start to display event " << EventNum << endl;
  cout << "    Last event " << last_event_id << endl;
  cout << "	 Event Flag " << flagslcio << endl;

  try {
    if (runNumber == -1) {
      lcReader->skipNEvents(EventNum);
      evt = lcReader->readNextEvent();
      runNumber = evt->getRunNumber();  // This trick only works for the LCIO
                                        // file with one single run number
      if (evt) event_id = evt->getEventNumber();
    } else {
      evt = lcReader->readEvent(runNumber, EventNum);
    }

    if (!evt) {
      cout << "    Event Not Found!" << endl;
    } else {
      event_id = evt->getEventNumber();

      string k =
          evt->getDetectorName();  // Initial Detector Type from slcio header;
      if (k == "ILD_00" || k == "ILD_o1_v05")  // ILD00_AHCAL
      {
        flagdetectortype = 0;
      } else if (k == "ILD_00Dhcal" || k == "ILD_01pre00" ||
                 k == "ILD_o2_v05" || k == "ILD_o2_v06" || k == "CEPC_v0" ||
                 k == "CEPC_v1")  // ILD00_DHCAL
      {
        flagdetectortype = 1;
      } else if (k.find("TB") != string::npos) {
        flagdetectortype = 2;  // Simulated test beam events; currently only
                               // considering about AHCAL Calibrations
      } else if (k == "sidloi3")  // SID_LOI
      {
        flagdetectortype = 3;
      } else if (k == "clic_sid_cdr_b")  // SID_CLIC
      {
        flagdetectortype = 5;
      } else {
        flagdetectortype =
            10;  // Risky... currently used as for the real data of CALICE TB
      }
      // As initial Value = -1, thus = 4 means has data file.

      cout << "    Event Located at : " << evt << endl;
      cout << endl
           << "****************************************************************"
              "***************"
           << endl
           << endl;
      cout << "    Event Statistics: " << endl << endl;

      load_collections(evt, "");

      last_event_id = EventNum;
      last_event_id = event_id;

      cout << "****************************************************************"
              "*********************"
           << endl;
      cout << "   Display of event " << event_id << " finished... " << endl
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
          Form("DRUID, RunNum = %d, EventNum = %d", runNumber, event_id));
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

    const std::vector<std::string>* strVec = evt->getCollectionNames();

    std::vector<std::string>::const_iterator name;

    if (coltype == "") {
      for (std::map<string, TEveElementList*>::iterator AA =
               collectionClasses.begin();
           AA != collectionClasses.end(); AA++) {
        CollectionScan[AA->first] = 0;
      }
    }

    for (name = strVec->begin(); name != strVec->end(); name++) {
      LCCollection* col = evt->getCollection(*name);
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

    for (name = strVec->begin(); name != strVec->end(); name++) {
      LCCollection* col = evt->getCollection(*name);

      std::string CollHead(*name, 0, 5);
      string ct = col->getTypeName();

      if (col->getNumberOfElements() == 0) continue;

      if (DisplayFlag.find(*name) == DisplayFlag.end()) {
        // if ( ct==LCIO::SIMCALORIMETERHIT ||
        if (flagdetectortype == 10) {
          FlagDraw = true;
        } else {
          FlagDraw = false;
        }
      } else
        FlagDraw = DisplayFlag[*name];

      if (coltype != "" && coltype != ct) continue;

      if (*name == "MCParticle") {
        if (pMCParticle) {
          MCPDraw = pMCParticle->GetRnrSelf();
          MCPChildDraw = pMCParticle->GetRnrChildren();

          for (TEveElement::List_i itt = pMCParticle->BeginChildren();
               itt != pMCParticle->EndChildren(); itt++) {
            std::string colname = (*itt)->GetElementName();
            if (MCParticleDisplayFlag.find(colname) ==
                MCParticleDisplayFlag.end())
              MCParticleDisplayFlag[colname] = (*itt)->GetRnrSelf();
          }

          pMCParticle->DestroyElements();
          pMCParticle->Destroy();
        }
        pMCParticle = BuildMCParticles(evt);
        pMCParticle->SetRnrSelfChildren(MCPDraw, MCPChildDraw);
        pMCParticle->SetName("MCPARTICLE");
        gEve->AddElement(pMCParticle);

        if (FlagMultiView) {
          gMultiView->ImportEventRPhi(pMCParticle);
          gMultiView->ImportEventRhoZ(pMCParticle);
        }

      } else if (ct == LCIO::MCPARTICLE) {
        cout << "  Collection <" << *name
             << "> already displayed in MCParticles form" << endl
             << endl;
      } else if (find(ignoredTypes.begin(), ignoredTypes.end(), ct) !=
                 ignoredTypes.end()) {
        cout << "  Collection <" << *name << "> will not be displayed. " << endl
             << endl;
      } else if (ct == LCIO::VERTEX) {
        TEveElementList* temp = Vertex(col, *name);
        temp->SetRnrSelfChildren(FlagDraw, FlagDraw);
        collectionClasses[ct]->AddElement(temp);
      } else if (ct == LCIO::RAWCALORIMETERHIT) {
        TEveElementList* temp = CaloHits(col, *name);
        temp->SetRnrSelfChildren(FlagDraw, FlagDraw);
        collectionClasses[ct]->AddElement(temp);
      } else if (ct == LCIO::CALORIMETERHIT) {
        TEveElementList* temp = CaloHits(col, *name);
        temp->SetRnrSelfChildren(FlagDraw, FlagDraw);
        collectionClasses[ct]->AddElement(temp);
      } else if (ct == LCIO::SIMCALORIMETERHIT) {
        TEveElementList* temp = CaloHits(col, *name);
        temp->SetRnrSelfChildren(FlagDraw, FlagDraw);
        collectionClasses[ct]->AddElement(temp);
      } else if (ct == LCIO::TRACKERHIT) {
        TEveElementList* temp = TrackerHits(col, *name);
        temp->SetRnrSelfChildren(FlagDraw, FlagDraw);
        collectionClasses[ct]->AddElement(temp);
      } else if (ct == LCIO::SIMTRACKERHIT) {
        TEveElementList* temp = TrackerHits(col, *name);
        temp->SetRnrSelfChildren(FlagDraw, FlagDraw);
        collectionClasses[ct]->AddElement(temp);
      } else if (ct == LCIO::TRACK) {
        TEveElementList* temp = TrackAssignedHits(col, *name);
        temp->SetRnrSelfChildren(FlagDraw, FlagDraw);
        collectionClasses[ct]->AddElement(temp);
      } else if (ct == LCIO::RECONSTRUCTEDPARTICLE) {
        if (*name == "PandoraPFOs" || *name == "PandoraPFANewPFOs" ||
            CollHead == "Arbor") {  // Supposed to be modified if user needed...
          TEveElementList* temp = BuildPFOs(col, *name);
          temp->SetRnrSelfChildren(FlagDraw, FlagDraw);
          collectionClasses[ct]->AddElement(temp);
        } else if (*name == "Durham_4Jets" || *name == "Durham_6Jets") {
          TEveElementList* temp = RecoJets(col, *name);
          temp->SetRnrSelfChildren(FlagDraw, FlagDraw);
          collectionClasses[ct]->AddElement(temp);
        } else {
          cout << "   " << *name
               << " is currently considering as Reconstructed Particles... "
                  "will skip "
               << endl;
        }
      } else if (ct == LCIO::LCRELATION &&
                 (CollHead == "Henri" || CollHead == "InitH" ||
                  CollHead == "InitE" || CollHead == "Links")) {
        HenriCount++;
        TEveCompound* temp = ConnectTrees(col, *name);
        temp->SetRnrSelfChildren(FlagDraw, FlagDraw);
        collectionClasses[ct]->AddElement(temp);
      } else if (ct == LCIO::CLUSTER) {
        TEveElementList* temp = ClusterHits(col, *name);
        temp->SetRnrSelfChildren(FlagDraw, FlagDraw);
        collectionClasses[ct]->AddElement(temp);
      } else {
        cout << "  Unknown collection type " << col->getTypeName()
             << " for collection " << *name << endl
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
  } else {
    std::cout
        << " Slcio data file is not available! skip loading event information. "
        << std::endl;
  }

  return;
}
