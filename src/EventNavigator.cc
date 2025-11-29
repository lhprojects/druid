
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//    GUI function defined.
//    //
//                                                                           //
//    Date:   03 Dec 2010                                                    //
//    Author: Manqi Ruan (LLR)                                               //
//                                                                           //
//    Last Modified: 17, Nov 2011, Cleaning //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "EventNavigator.hh"

#include <cstdlib>
#include <iostream>
#include <string>

#include "EVENT/LCEvent.h"
#include "GlobalDefs.hh"
#include "MultiView.hh"
#include "Rtypes.h"
#include "TGNumberEntry.h"
#include "TEveManager.h"
#include "TGLViewer.h"
#include "TEveBoxSet.h"
#include "TEveRGBAPalette.h"
#include "TEveRGBAPaletteOverlay.h"
#include "Options.h"

extern LCReader* lcReader;
extern LCEvent* evt;

extern bool HitSizeLog;
extern bool Flag_AttachTextToHit;
extern float PFOHitSize;
extern float ClusterHitSize;
extern float cellfactor;
extern double cellEnergyThresh;
extern float HitColourFactor;
extern float HitColourUnitFactor;
extern int GlobalRandomColorIndex;
extern float HitEnergyCut;

// Global EventNavigator instance
EventNavigator gEventNavigator;

// for TEveRGBAPalette powered coloring style
extern TEveRGBAPalette* p;
extern TEveRGBAPaletteOverlay* po;
extern TEveBoxSet* cal_shell;
extern bool cal_shell_shown;
extern float HitColourOverflowLimit;
extern float HitColourUnderflowLimit;

extern GUIManager gGUIManager;

extern TGNumberEntry* _cellSizeEntry;
extern TGNumberEntry* _SimucellSizeEntry;
extern TGNumberEntry* _PFOcellSizeEntry;
extern TGNumberEntry* _ClustercellSizeEntry;
extern TGNumberEntry* _cellEnergyThrEntry;
extern TGNumberEntry* _cellEnergyScaleEntry;
extern TGNumberEntry* _PTCutEntry;
extern TGNumberEntry* _CellECutEntry;
// Yuzhi CHE: RGBAPalette
extern TGNumberEntry* _cellColorOverflowLimit;
extern TGNumberEntry* _cellColorUnderflowLimit;

extern bool DefaultCollectionFlag;  // kTRUE, default display collections;
                                    // kFALSE; all collections
extern bool FlagMultiView;
extern bool HiddenLowESimCell;
extern bool HiddenPFOCluster;

extern TEveManager * gEve;

int EventNumberToDisplay;

void EventNavigator::Fwd()
{
    int eventNumber = gDisplayState.getEventNumber();
    if (eventNumber < 1000000)
    {
        ++eventNumber;
        while (lcReader->readEvent(gDisplayState.getRunNumber(), eventNumber) == nullptr)
        {
            ++eventNumber;
        }
        load_event(eventNumber);
    }
    else
    {
        printf("Already at last event.\n");
    }
}

void EventNavigator::Bck() {
    int eventNumber = gDisplayState.getEventNumber();
  if (eventNumber > 0) {
    --eventNumber;
    while (lcReader->readEvent(gDisplayState.getRunNumber(), eventNumber) == nullptr &&
           eventNumber > 0)
    {
        --eventNumber;
    }
    load_event(eventNumber);
  } else {
    printf("Already at first event.\n");
  }
}

/*
   void EventNavigator::setCollection() {
   DefaultCollectionFlag = !DefaultCollectionFlag;
   if(DefaultCollectionFlag) printf("Display default collections. \n");
   else printf("Display all collections. \n");
   load_collections(evt,"");
   return;
   }
   */

void EventNavigator::setCollection()
{
    printf(" Update colors \n");

    if (gGUIManager.HitColourType == 7)
    {
        UpdateHitColorsToMatchMCParticles();
        UpdateTrackColorsByMCParticle();
    }
    else
    {
        load_collections(evt, LCIO::TRACK);
        load_collections(evt, LCIO::SIMCALORIMETERHIT);
        load_collections(evt, LCIO::CALORIMETERHIT);
    }
}

void EventNavigator::MultiViewSwitch() {
  std::cout << "about to be added!" << std::endl;
  FlagMultiView = !FlagMultiView;
  load_event(gDisplayState.getEventNumber());
}

void EventNavigator::setCellColour(int ds) {
  std::cout << "SetCellColor: " << ds << std::endl;
  gGUIManager.HitColourType = ds;
  if (gGUIManager.HitColourType < 0) gGUIManager.HitColourType = 0;
  if (gGUIManager.HitColourType > 7) gGUIManager.HitColourType = 0;
  
  if(gGUIManager.HitColourType == 7) {
    // For color scheme 7, just update the colors to match current MCParticle visibility
    UpdateHitColorsToMatchMCParticles();
  } else {
    // For other color schemes (0-6), need to reload to apply new coloring
    load_collections(evt, LCIO::SIMCALORIMETERHIT);
  }
  
  return;
}

void EventNavigator::HideLowSimE() {
  HiddenLowESimCell = !HiddenLowESimCell;
  load_collections(evt, LCIO::SIMCALORIMETERHIT);
  return;
}

void EventNavigator::HidePFOClu() {
  HiddenPFOCluster = !HiddenPFOCluster;
  load_collections(evt, LCIO::RECONSTRUCTEDPARTICLE);
  return;
}

void EventNavigator::setPFOCellColour(int ds) {
  gGUIManager.PFOHitColourType = ds;

  if (gGUIManager.PFOHitColourType < 0 || gGUIManager.PFOHitColourType > 2) {
    gGUIManager.PFOHitColourType = 1;
  }

  load_collections(evt, LCIO::RECONSTRUCTEDPARTICLE);
  return;
}

void EventNavigator::colorReroll() {
  if (gGUIManager.HitColourType != 2)
    gGUIManager.HitColourType = 3;  // if Energy not of Origin type, automatically swithched
                        // to according to Random index
  gGUIManager.PFOHitColourType = 1;
  gGUIManager.ClusterHitColourType = 1;
  GlobalRandomColorIndex++;
  std::cout << "    GlobalRandomColorIndex: " << GlobalRandomColorIndex
            << std::endl;

  load_collections(evt, LCIO::SIMCALORIMETERHIT);
  std::cout << "    Reload SimCaloHit Finish" << std::endl;
  load_collections(evt, LCIO::RECONSTRUCTEDPARTICLE);
  std::cout << "    Reload RECOParticle Finish" << std::endl;
  load_collections(evt, LCIO::CLUSTER);
  std::cout << "    Reload Cluster Finsh" << std::endl;
  load_collections(evt, LCIO::TRACK);
  std::cout << "    Reload TRACK Finsh" << std::endl;
  load_collections(evt, LCIO::VERTEX);
  std::cout << "    Reload VERTEX Finsh" << std::endl;
  return;
}

void EventNavigator::ScaleModify() {
  HitSizeLog = !HitSizeLog;
  std::cout << "Scale the Size of Digitized CaloHits according to logrithm of "
               "its Energy! "
            << std::endl;
  load_collections(evt, LCIO::CALORIMETERHIT);
  return;
}

void EventNavigator::HitTextAttach() {
  Flag_AttachTextToHit = !Flag_AttachTextToHit;
  std::cout << "HitAttachFlagValue = " << Flag_AttachTextToHit << std::endl;
  load_collections(evt, LCIO::CALORIMETERHIT);
  load_collections(evt, LCIO::RAWCALORIMETERHIT);
  load_collections(evt, LCIO::SIMCALORIMETERHIT);
  load_collections(evt, LCIO::TRACKERHIT);
  load_collections(evt, LCIO::SIMTRACKERHIT);
  load_collections(evt, LCIO::CLUSTER);
  return;
}

void EventNavigator::PTCutModify() {
  gOptions.MCPtCut = _PTCutEntry->GetNumber();
  std::cout << "Change the PTCut for the MCParticle Collection: " << gOptions.MCPtCut
            << std::endl;
  load_collections(evt, LCIO::MCPARTICLE);
  return;
}

void EventNavigator::CellECutModify() {
  HitEnergyCut = _CellECutEntry->GetNumber();
  std::cout << "Change the Energy Cut to Calo Hits: " << HitEnergyCut
            << " Mips." << std::endl;
  load_collections(evt, LCIO::CALORIMETERHIT);
  load_collections(evt, LCIO::SIMCALORIMETERHIT);
  return;
}

void EventNavigator::SizeModify() {
  cellfactor = _cellSizeEntry->GetNumber();
  std::cout << "Size of the Digitized Hit in following collection will be "
               "scaled by a factor of: "
            << cellfactor << std::endl;
  load_collections(evt, LCIO::CALORIMETERHIT);
  return;
}

void EventNavigator::SimuHitSizeModify() {
  cellfactor = _SimucellSizeEntry->GetNumber();
  std::cout << "Cellfactor " << cellfactor << std::endl;
  load_collections(evt, LCIO::SIMCALORIMETERHIT);
  return;
}

void EventNavigator::PFOHitSizeModify() {
  PFOHitSize = _PFOcellSizeEntry->GetNumber();
  std::cout << "  Modify PFO Associated Hit Size: " << std::endl;
  std::cout << "  Cellfactor " << PFOHitSize << std::endl;
  load_collections(evt, LCIO::RECONSTRUCTEDPARTICLE);
  return;
}

void EventNavigator::ClusterHitSizeModify() {
  ClusterHitSize = _ClustercellSizeEntry->GetNumber();
  std::cout << "  Modify Cluster Hit Size: " << std::endl;
  std::cout << "  Cellfactor " << PFOHitSize << std::endl;
  load_collections(evt, LCIO::CLUSTER);
  load_collections(evt, LCIO::RECONSTRUCTEDPARTICLE);
  return;
}

void EventNavigator::GotoEvent() {
  EventNumberToDisplay = int(gGUIManager._EventNumberEntry->GetNumber());
  std::cout << "  Let's go to display " << EventNumberToDisplay << "th Event!"
            << std::endl;
  gDisplayState.setEventNumber(EventNumberToDisplay);
  load_event(EventNumberToDisplay);
  return;
}

void EventNavigator::setEnergyScale() {
  gGUIManager.HitColourType = 1;
  cellEnergyThresh = _cellEnergyThrEntry->GetNumber();
  //  cellEnergyThresh *= 0.00001;
  std::cout << std::endl
            << "DHCAL Threshold Changing! Current Thresholds: " << std::endl
            << " Blue: from " << cellEnergyThresh * 1000 << " MeV to "
            << cellEnergyThresh * 10000 << " MeV" << std::endl
            << " Green: from " << cellEnergyThresh * 10000 << " MeV to "
            << cellEnergyThresh * 100000 << " MeV" << std::endl
            << " Red: > " << cellEnergyThresh * 100000 << " MeV" << std::endl
            << std::endl;
  load_collections(evt, LCIO::CALORIMETERHIT);
  load_collections(evt, LCIO::SIMCALORIMETERHIT);
  return;
}

void EventNavigator::setGlobalEnergyScale() {
  HitColourFactor = 10;
  HitColourFactor = _cellEnergyScaleEntry->GetNumber();

  std::cout << std::endl
            << "Global Energy Scale Changing! " << HitColourFactor << std::endl;
  std::cout << "Energy Scale = n Means 10/n Mip has Color Blue, and > 100/n "
               "Mips has clor Red!"
            << std::endl
            << std::endl;
  ;

  load_collections(evt, LCIO::SIMCALORIMETERHIT);
  load_collections(evt, LCIO::CALORIMETERHIT);

  return;
}

void EventNavigator::setColorOverflowLimit() {
  HitColourOverflowLimit = 10;
  HitColourOverflowLimit = _cellColorOverflowLimit->GetNumber();

  std::cout << std::endl
            << "Color Overflow Limit Changed! " << HitColourFactor << std::endl;
  std::cout << "automatically turn on the alternative coloring style"
            << std::endl;

  cal_shell->Reset(TEveBoxSet::kBT_FreeBox, kFALSE, 64);
  load_collections(evt, LCIO::SIMCALORIMETERHIT);
  load_collections(evt, LCIO::CALORIMETERHIT);

  // Turn on the alternative coloring style powered by TEveRGBAPalette
  p->SetFixColorRange(true);
  p->SetLimits(HitColourUnderflowLimit * HitColourUnitFactor,
                       HitColourOverflowLimit * HitColourUnitFactor);
  if (cal_shell_shown) {
    gEve->EditElement(cal_shell);
  } else {
    gEve->AddElement(cal_shell);
    cal_shell_shown = true;
  }
  TGLViewer* v = gEve->GetDefaultGLViewer();
  p->SetupColorArray();
  v->RemoveOverlayElement(po);
  v->AddOverlayElement(po);
  v->UpdateScene();

  return;
}

void EventNavigator::setColorUnderflowLimit() {
  HitColourUnderflowLimit = 0;
  HitColourUnderflowLimit = _cellColorUnderflowLimit->GetNumber();

  std::cout << std::endl
            << "Color Underflow Limit Changed! " << HitColourFactor << std::endl;
  std::cout << "automatically turn on the alternative coloring style"
            << std::endl;

  cal_shell->Reset(TEveBoxSet::kBT_FreeBox, kFALSE, 64);
  load_collections(evt, LCIO::SIMCALORIMETERHIT);
  load_collections(evt, LCIO::CALORIMETERHIT);

  // Turn on the alternative coloring style powered by TEveRGBAPalette
  p->SetFixColorRange(true);
  p->SetLimits(HitColourUnderflowLimit * HitColourUnitFactor,
                       HitColourOverflowLimit * HitColourUnitFactor);
  if (cal_shell_shown) {
    gEve->EditElement(cal_shell);
  } else {
    gEve->AddElement(cal_shell);
    cal_shell_shown = true;
  }
  TGLViewer* v = gEve->GetDefaultGLViewer();
  p->SetupColorArray();
  v->AddOverlayElement(po);
  v->UpdateScene();

  return;
}

void EventNavigator::OnTrackVisibilityChanged()
{
  std::cout << "=== Track Visibility Changed ===" << std::endl;

  // Update SimCalorimeterHit colors only if using color scheme 7 (MCParticle colors)
  if (gGUIManager.HitColourType == 7) {
    UpdateHitColorsToMatchMCParticles();
  }

	std::cout << "=== Update Complete ===" << std::endl;
}

void EventNavigator::OnElementVisibilityChanged()
{
	OnTrackVisibilityChanged();
}