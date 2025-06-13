
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//    Druid option Pannel													 //
//                                                                           //
//    Date:   03 Dec 2010                                                    //
//    Author: Manqi Ruan (LLR)                                               //
//                                                                           //
//    Last Modified: 17, Nov 2011, Cleaning									 //
//                                                                           // 
///////////////////////////////////////////////////////////////////////////////


// 定义一个函数，用于创建GUI

// 创建GUI

#include "GlobalDefs.hh"
#include "EventNavigator.hh"
#include "TEveBrowser.h"
#include "TGNumberEntry.h"
#include "TGButton.h"
#include "TGLabel.h"
#include "TStyle.h"
#include "TGComboBox.h"
#include "TEveManager.h"
#include "TSystem.h"
#include "TGLViewer.h"
#include "TEveViewer.h"
#include "TGLPerspectiveCamera.h"
#include "TGLCamera.h"

extern TEveManager * gEve;
extern TSystem * gSystem;

TGNumberEntry *_cellEnergyThrEntry;
TGNumberEntry *_EventNumberEntry;
TGNumberEntry *_cellSizeEntry;
TGNumberEntry *_SimucellSizeEntry;
TGNumberEntry *_PTCutEntry;
TGNumberEntry *_CellECutEntry;
TGNumberEntry *_PFOcellSizeEntry;
TGNumberEntry *_ClustercellSizeEntry;
TGNumberEntry *_cellEnergyScaleEntry;	//Tune EF factor for different color of Energy
TGNumberEntry *_cellColorOverflowLimit;	
TGNumberEntry *_cellColorUnderflowLimit;

void make_gui()
{
  TEveBrowser* browser = gEve->GetBrowser();
	browser->StartEmbedding(TRootBrowser::kLeft);
	TGMainFrame* frmMain = new TGMainFrame(gClient->GetRoot(), 1000, 600);
	frmMain->SetWindowName("XX GUI");
	frmMain->SetCleanup(kDeepCleanup);

	TGLViewer* glviewer = gEve->GetDefaultViewer()->GetGLViewer();
	TGLViewer* v = gEve->GetDefaultGLViewer();
	EventNavigator* fh = new EventNavigator;

	//event navigation   
	{
		TGGroupFrame *frmEvent = 
			new TGGroupFrame(frmMain, "Event Navigation", kHorizontalFrame);   
		TGHorizontalFrame* hf = new TGHorizontalFrame(frmMain);

		frmEvent->AddFrame(hf);

		TString icondir( Form("%s/icons/", gSystem->Getenv("DRUIDDIR")) );

		TGPictureButton* b = 0;

		b = new TGPictureButton(hf, gClient->GetPicture(icondir + "GoBack.gif"));     
		hf->AddFrame(b);
		b->Connect("Clicked()", "EventNavigator", fh, "Bck()");

		b = new TGPictureButton(hf, gClient->GetPicture(icondir + "GoForward.gif"));     
		hf->AddFrame(b);
		b->Connect("Clicked()", "EventNavigator", fh, "Fwd()");

		TGHorizontalFrame *EventNrFrame = new TGHorizontalFrame(frmEvent);

		TGLabel *EventNrLabel = new TGLabel(EventNrFrame, " Go to \n  Evt" );

		_EventNumberEntry = new TGNumberEntry(EventNrFrame, 0, 5, -1,
				TGNumberFormat::kNESInteger,
				TGNumberFormat::kNEAAnyNumber,
				TGNumberFormat::kNELLimitMin,
				0, 100);

		EventNrFrame->AddFrame(EventNrLabel);
		EventNrFrame->AddFrame(_EventNumberEntry);

		_EventNumberEntry->Connect("ValueSet(Long_t)", "EventNavigator", fh, "GotoEvent()");

		frmEvent->AddFrame(EventNrFrame);
		frmMain->AddFrame(frmEvent);
	}	


	{
		TGGroupFrame *frmRotateColor =
			new TGGroupFrame(frmMain, "Rotation Center, Hits Color", kHorizontalFrame);

		TString icondir( Form("%s/icons/", gSystem->Getenv("DRUIDDIR")) );

		TGPictureButton* d = 0;
		d = new TGPictureButton(frmRotateColor, gClient->GetPicture(icondir + "RotationCenter.png"));
		d->Connect("Clicked()", "TGLViewer", glviewer, "PickCameraCenter()");
		frmRotateColor->AddFrame(d);

		d = new TGPictureButton(frmRotateColor, gClient->GetPicture(icondir + "ColorSelection.png"));
		d->Connect("Clicked()", "EventNavigator", fh, "colorReroll()");
		frmRotateColor->AddFrame(d);

		d = new TGPictureButton(frmRotateColor, gClient->GetPicture(icondir + "ReloadPage.gif"));
		d->Connect("Clicked()", "EventNavigator", fh, "setCollection()");
		frmRotateColor->AddFrame(d);

		d = new TGPictureButton(frmRotateColor, gClient->GetPicture(icondir + "HitText.gif"));
		d->Connect("Clicked()", "EventNavigator", fh, "HitTextAttach()");
		frmRotateColor->AddFrame(d);

		frmMain->AddFrame(frmRotateColor);
	}

	//MCParticle, PTCut
	{
		TGGroupFrame *frmPTCut = new TGGroupFrame(frmMain, "En Cut (MCParticle); E Cut (Cell) ");
		frmMain->AddFrame(frmPTCut, new TGLayoutHints(kLHintsNormal, 2, 2, 0, 0));

		TGHorizontalFrame *cellSize = new TGHorizontalFrame(frmPTCut);
		TGLabel *cellLabel = new TGLabel(cellSize, "GeV");

		cellSize->AddFrame(cellLabel);

		_PTCutEntry = new TGNumberEntry(frmPTCut, 1.5, 5, -1,
				TGNumberFormat::kNESRealOne,
				TGNumberFormat::kNEAPositive,
				TGNumberFormat::kNELNoLimits,
				1, 100);

		frmPTCut->AddFrame(_PTCutEntry);
		frmPTCut->AddFrame(cellSize);
		_PTCutEntry->Connect("ValueSet(Long_t)", "EventNavigator", fh, "PTCutModify()");

		TGHorizontalFrame *cellECut = new TGHorizontalFrame(frmPTCut);
		TGLabel *ECutLabel = new TGLabel(cellECut, "Mip");

		cellECut->AddFrame(ECutLabel);
		_CellECutEntry = new TGNumberEntry(frmPTCut, 0.2, 5, -1,
				TGNumberFormat::kNESRealOne,
				TGNumberFormat::kNEAPositive,
				TGNumberFormat::kNELNoLimits,
				1, 100);

		frmPTCut->AddFrame(_CellECutEntry);
		frmPTCut->AddFrame(cellECut);
		_CellECutEntry->Connect("ValueSet(Long_t)", "EventNavigator", fh, "CellECutModify()");
	}

	//Simulated CaloHits Style;
	{
		TGGroupFrame *frmHitColour = new TGGroupFrame(frmMain, "SimCaloHit Colur: ");   
		frmMain->AddFrame(frmHitColour, new TGLayoutHints(kLHintsNormal, 2, 2, 0, 0));

		TGComboBox* menuHitColour = new TGComboBox(frmHitColour);
		frmHitColour->AddFrame(menuHitColour, new TGLayoutHints(kLHintsNormal, 1, 1, 1, 1));
		//        menuHitColour->AddEntry("Hit Energy", 0);
		menuHitColour->AddEntry("PDG of Track", 0);
		menuHitColour->AddEntry("Hit Energy", 1);
		menuHitColour->AddEntry("PDG of Mother", 2);
		menuHitColour->AddEntry("Naive index", 3);
		menuHitColour->AddEntry("Uniform Yellow", 4);
		menuHitColour->AddEntry("EM, Had & Neutron", 5);
		menuHitColour->AddEntry("Timing/ns", 6);
		menuHitColour->Select(0, kFALSE);		//Default choice: Color with PID
		//menuHitColour->Connect("Selected(Int_t)", "EventNavigator", fh, "setCellColour(Int_t)");
		menuHitColour->Connect("Selected(int)", "EventNavigator", fh, "setCellColour(int)");
		menuHitColour->Resize(150, 20);

		TGCheckButton * CellHideButton = new TGCheckButton(frmHitColour, "Show Cell E < Cut");
		frmHitColour->AddFrame(CellHideButton);
		CellHideButton->Connect("Clicked()", "EventNavigator", fh, "HideLowSimE()");

		TGHorizontalFrame *cellSize = new TGHorizontalFrame(frmHitColour);

		TGLabel *cellLabel = new TGLabel(cellSize, " Size (cm. *3 for AHCAL) ");

		cellSize->AddFrame(cellLabel);

		_SimucellSizeEntry = new TGNumberEntry(frmHitColour, 1, 5, -1,
				TGNumberFormat::kNESRealOne,
				TGNumberFormat::kNEAPositive,
				TGNumberFormat::kNELNoLimits,
				1, 100);

		frmHitColour->AddFrame(_SimucellSizeEntry);
		frmHitColour->AddFrame(cellSize);
		_SimucellSizeEntry->Connect("ValueSet(Long_t)", "EventNavigator", fh, "SimuHitSizeModify()");

	}

	//PFO CaloHit Style
	{
		TGGroupFrame *frmPFOHitColour = new TGGroupFrame(frmMain, "PFOCaloHit Colur: ");
		frmMain->AddFrame(frmPFOHitColour, new TGLayoutHints(kLHintsNormal, 2, 2, 0, 0));

		TGComboBox* menuPFOHitColour = new TGComboBox(frmPFOHitColour);
		frmPFOHitColour->AddFrame(menuPFOHitColour, new TGLayoutHints(kLHintsNormal, 1, 1, 1, 1));
		menuPFOHitColour->AddEntry("PDG of Track", 0);
		menuPFOHitColour->AddEntry("Naive index", 1);
		menuPFOHitColour->AddEntry("Uniform Green", 2);
		menuPFOHitColour->Select(0, kFALSE);		//Default Choice: color with PID
		menuPFOHitColour->Connect("Selected(Int_t)", "EventNavigator", fh, "setPFOCellColour(Int_t)");
		menuPFOHitColour->Resize(150, 20);

		TGCheckButton * PFOCluHideButton = new TGCheckButton(frmPFOHitColour, "Hide PFO Cluster");
                frmPFOHitColour->AddFrame(PFOCluHideButton);
                PFOCluHideButton->Connect("Clicked()", "EventNavigator", fh, "HidePFOClu()");

		TGHorizontalFrame *PFOcellSize = new TGHorizontalFrame(frmPFOHitColour);

		TGLabel *PFOcellLabel = new TGLabel(PFOcellSize, "Size of PFO hits (cm) ");

		PFOcellSize->AddFrame(PFOcellLabel);

		_PFOcellSizeEntry = new TGNumberEntry(frmPFOHitColour, 1, 5, -1,
				TGNumberFormat::kNESRealOne,
				TGNumberFormat::kNEAPositive,
				TGNumberFormat::kNELNoLimits,
				1, 100);

		/*
		   frmPFOHitColour->AddFrame(_PFOcellSizeEntry, new TGLayoutHints(kLHintsTop |
		   kLHintsLeft |
		   kLHintsCenterY,
		   16, 0, 1, 1));
		   */

		frmPFOHitColour->AddFrame(_PFOcellSizeEntry);
		frmPFOHitColour->AddFrame(PFOcellSize);
		_PFOcellSizeEntry->Connect("ValueSet(Long_t)", "EventNavigator", fh, "PFOHitSizeModify()");

		/*
		   frmPFOHitColour->AddFrame(cellSize, new TGLayoutHints(kLHintsTop |
		   kLHintsLeft |
		   kLHintsCenterY, 0, 0, 0, 0));
		   */
	}

	{
		TGGroupFrame *frmClusterHit = new TGGroupFrame(frmMain, "Cluster Hit Size: ");
		frmMain->AddFrame(frmClusterHit, new TGLayoutHints(kLHintsNormal, 2, 2, 0, 0));

		//	TGHorizontalFrame *ClustercellSize = new TGHorizontalFrame(frmClusterHit);
		TGVerticalFrame *ClustercellSize = new TGVerticalFrame(frmClusterHit);
		TGLabel *ClustercellLabel = new TGLabel(ClustercellSize, "Size of \ncluster hits: ");

		ClustercellSize->AddFrame(ClustercellLabel, new TGLayoutHints(kLHintsTop |
					kLHintsLeft |
					kLHintsCenterY, 0, 0, 0, 0));

		_ClustercellSizeEntry = new TGNumberEntry(frmClusterHit, 1, 5, -1,
				TGNumberFormat::kNESRealOne,
				TGNumberFormat::kNEAPositive,
				TGNumberFormat::kNELNoLimits,
				1, 100);

		//        frmClusterHit->AddFrame(ClustercellSize);

		_ClustercellSizeEntry->Connect("ValueSet(Long_t)", "EventNavigator", fh, "ClusterHitSizeModify()");

		frmClusterHit->AddFrame(_ClustercellSizeEntry);

	}

	//For cell size
	{
		TGGroupFrame *HitSizeFrame =
			new TGGroupFrame(frmMain, "Calorimeter Hit Size", kVerticalFrame);
		frmMain->AddFrame(HitSizeFrame, new TGLayoutHints(kLHintsNormal, 2, 2, 0, 0));

		TGHorizontalFrame *cellSize = new TGHorizontalFrame(HitSizeFrame);

		TGLabel *cellLabel = new TGLabel(cellSize, "Size (cm, *3 for AHCAL)""\n");

		cellSize->AddFrame(cellLabel);

		_cellSizeEntry = new TGNumberEntry(HitSizeFrame, 1, 5, -1,
				TGNumberFormat::kNESRealOne,
				TGNumberFormat::kNEAPositive,
				TGNumberFormat::kNELNoLimits,
				1, 100);

		HitSizeFrame->AddFrame(_cellSizeEntry);
		HitSizeFrame->AddFrame(cellSize);
		_cellSizeEntry->Connect("ValueSet(Long_t)", "EventNavigator", fh, "SizeModify()");

		TGVerticalFrame *manualFrame = new TGVerticalFrame(HitSizeFrame);
		HitSizeFrame->AddFrame(manualFrame);

		TGHorizontalFrame *energyFrame = new TGHorizontalFrame(manualFrame);
		manualFrame->AddFrame(energyFrame);

		TGLabel *energyLabel = new TGLabel(energyFrame, "DHCAL Thres(mips): \n0.2, 5, 10/SF:");		//Mips
		energyFrame->AddFrame(energyLabel,
				new TGLayoutHints(kLHintsLeft |
					kLHintsCenterY)
				);

		_cellEnergyThrEntry = new TGNumberEntry(energyFrame, 0, 5, -1,
				TGNumberFormat::kNESInteger,
				TGNumberFormat::kNEAAnyNumber,
				TGNumberFormat::kNELLimitMinMax,
				0, 100);

		energyFrame->AddFrame(_cellEnergyThrEntry);

		_cellEnergyThrEntry->Connect("ValueSet(Long_t)", "EventNavigator", fh, "setEnergyScale()");
	}

	{	
		TGGroupFrame *GlobalenergyScale =
			new TGGroupFrame(frmMain, "Global SF for Calo Hit", kVerticalFrame);
		frmMain->AddFrame(GlobalenergyScale, new TGLayoutHints(kLHintsNormal, 2, 2, 0, 0));

		TGVerticalFrame *SSmanualFrame = new TGVerticalFrame(GlobalenergyScale);
		GlobalenergyScale->AddFrame(SSmanualFrame);

		TGHorizontalFrame *GlobalenergyFrame = new TGHorizontalFrame(SSmanualFrame);

		SSmanualFrame->AddFrame(GlobalenergyFrame,
				new TGLayoutHints(kLHintsTop |
					kLHintsLeft |
					kLHintsCenterY, 0, 0, 0, 0)
				);

		TGLabel *energyLabel = new TGLabel(GlobalenergyFrame,
				"Color Scale:\n"
				);

		GlobalenergyFrame->AddFrame(energyLabel,
				new TGLayoutHints(kLHintsLeft |
					kLHintsCenterY)
				);

		_cellEnergyScaleEntry = new TGNumberEntry(GlobalenergyFrame, 10, 5, -1,
				TGNumberFormat::kNESInteger,
				TGNumberFormat::kNEAAnyNumber,
				TGNumberFormat::kNELLimitMinMax,
				0, 2000);

		GlobalenergyFrame->AddFrame(_cellEnergyScaleEntry,
				new TGLayoutHints(kLHintsTop |
					kLHintsLeft |
					kLHintsCenterY,
					16, 0, 1, 1));

		_cellEnergyScaleEntry->Connect("ValueSet(Long_t)", "EventNavigator", fh, "setGlobalEnergyScale()");

    // Overflow limit
		TGVerticalFrame *overflow_vert_frame = new TGVerticalFrame(GlobalenergyScale);
		GlobalenergyScale->AddFrame(overflow_vert_frame);
		TGHorizontalFrame *overflow_frame = new TGHorizontalFrame(overflow_vert_frame);
    overflow_vert_frame->AddFrame(overflow_frame);
		TGLabel *overflow_label = new TGLabel(overflow_frame,
				"Color Overflow Limit:\n"
				);

		overflow_frame->AddFrame(overflow_label,
				new TGLayoutHints(kLHintsLeft |
					kLHintsCenterY)
				);

		_cellColorOverflowLimit = new TGNumberEntry(overflow_frame, 10, 5, -1,
				TGNumberFormat::kNESRealOne,
				TGNumberFormat::kNEAAnyNumber,
				TGNumberFormat::kNELLimitMinMax,
				0, 2000);

		overflow_frame->AddFrame(_cellColorOverflowLimit,
				new TGLayoutHints(kLHintsTop |
					kLHintsLeft |
					kLHintsCenterY,
					16, 0, 1, 1));

		_cellColorOverflowLimit->Connect("ValueSet(Long_t)", "EventNavigator", fh, "setColorOverflowLimit()");

    // Overflow limit
		TGVerticalFrame *underflow_vert_frame = new TGVerticalFrame(GlobalenergyScale);
		GlobalenergyScale->AddFrame(underflow_vert_frame);
		TGHorizontalFrame *underflow_frame = new TGHorizontalFrame(underflow_vert_frame);
    underflow_vert_frame->AddFrame(underflow_frame);
		TGLabel *underflow_label = new TGLabel(underflow_frame,
				"Color underflow Limit:\n"
				);

		underflow_frame->AddFrame(underflow_label,
				new TGLayoutHints(kLHintsLeft |
					kLHintsCenterY)
				);

		_cellColorUnderflowLimit = new TGNumberEntry(underflow_frame, 0, 5, -1,
				TGNumberFormat::kNESRealOne,
				TGNumberFormat::kNEAAnyNumber,
				TGNumberFormat::kNELLimitMinMax,
				0, 2000);

		underflow_frame->AddFrame(_cellColorUnderflowLimit,
				new TGLayoutHints(kLHintsTop |
					kLHintsLeft |
					kLHintsCenterY,
					16, 0, 1, 1));

		_cellColorUnderflowLimit->Connect("ValueSet(Long_t)", "EventNavigator", fh, "setColorUnderflowLimit()");

    // Underflow limit
    /*
		TGLabel *energyLabel = new TGLabel(GlobalenergyFrame,
				"Color Underflow limit:\n"
				);

		GlobalenergyFrame->AddFrame(energyLabel,
				new TGLayoutHints(kLHintsLeft |
					kLHintsCenterY)
				);

		_cellEnergyScaleEntry = new TGNumberEntry(GlobalenergyFrame, 10, 5, -1,
				TGNumberFormat::kNESInteger,
				TGNumberFormat::kNEAAnyNumber,
				TGNumberFormat::kNELLimitMinMax,
				0, 2000);

		GlobalenergyFrame->AddFrame(_cellEnergyScaleEntry,
				new TGLayoutHints(kLHintsTop |
					kLHintsLeft |
					kLHintsCenterY,
					16, 0, 1, 1));

		_cellEnergyScaleEntry->Connect("ValueSet(Long_t)", "EventNavigator", fh, "setGlobalEnergyScale()");
    */

	}


	frmMain->MapSubwindows();
	frmMain->Resize();
	frmMain->MapWindow();
	browser->StopEmbedding();
	browser->SetTabTitle("Options", 0);
}


