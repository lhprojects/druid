
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//    Build TEveTracks according to MCParticle Collection					 //
//    Colors and Style are determined according to PID						 //
//    Default Cut at minimal PT(or Energy) = 1.5GeV --						 //
//							can be adjusted from CUI   						 //
//                                                                           //
//    Date:   03 Dec 2010                                                    //
//    Author: Manqi Ruan (LLR)                                               //
//                                                                           //
//    Last Modified: 17, Dec 2011, cleaning					                 //
//                                                                           // 
///////////////////////////////////////////////////////////////////////////////

#include "TEveManager.h"
#include "TEveElement.h"
#include "TEveTrack.h"
#include "TEveArrow.h"
#include "TEveTrackPropagator.h"
#include "TDatabasePDG.h"
#include "TEveVSDStructs.h"
#include "lcio.h"
#include "IMPL/LCTOOLS.h"
#include "EVENT/LCCollection.h"
#include "EVENT/SimTrackerHit.h"
#include "EVENT/SimCalorimeterHit.h"
#include "EVENT/MCParticle.h"
#include "EVENT/LCEvent.h"
#include "Options.h"
#include "GlobalDefs.hh"
#include "EventNavigator.hh"
#include "TruthHelper.h"

using namespace lcio;
using namespace EVENT;
using namespace std;


//float PTCut = 0.1; //GeV; Tracks with PT less than this threshold will not be displayed;

TEvePathMark * PathMarkEndTrackDecay(TEveVector &/*Vtx*/, TEveVector &End){
	TEveVector Mark = End;
	TEvePathMark* pm = new TEvePathMark(TEvePathMark::kDecay);
	pm->fV.Set(Mark);
	return pm;
}

bool IsNeutrino(int PID){
	int pid = abs(PID);
	if( pid==12 || pid==14 || pid==16 ) return true;
	return false;
}

// Get color variation for a particle based on its index
// pid is the PDG ID (for future use)
// index is the original particle index in the collection (before any filtering)
// baseColor is the ROOT color constant (kBlue, kRed, etc.)
// Returns a color from the color wheel (e.g., kBlue-7, kRed+2, etc.)
int GetParticleColor(int pid, int index, int baseColor) {
	// Each index increases color by 2
	// ROOT colors are valid in ranges like kBlue-9 to kBlue+10
	int colorOffset = (index * 2) % 20 - 9;
	int color = baseColor + colorOffset;
	// Clamp color to safe range: baseColor-9 to baseColor+9
	if (color > baseColor + 9) color = baseColor + 9;
	if (color < baseColor - 9) color = baseColor - 9;
	// Ensure color is positive
	if (color < 1) color = baseColor;
	return color;
}


TEveElementList* BuildMCParticles( LCEvent *evt, std::vector<std::string> const &mcpartcollnames)
{
    gGUIManager._MCPTracks.clear();
    gGUIManager._MCPGroups.clear();
    gGUIManager._MCPShowing.clear();

	// Preserve visibility states from previous MCPList
	bool MCPDraw = true;
	bool MCPChildDraw = true;
	if (gGUIManager._MCPList)
	{
		MCPDraw = gGUIManager._MCPList->GetRnrSelf();
		MCPChildDraw = gGUIManager._MCPList->GetRnrChildren();

		std::cout << " Saving MC Particle visibility states: " << std::endl;
		std::cout << "  MCPList RnrSelf: " << MCPDraw << ", RnrChildren: " << MCPChildDraw << std::endl;
		std::cout << "  Particle Type Visibility States: " << std::endl;

		for (TEveElement::List_i itt = gGUIManager._MCPList->BeginChildren();
		     itt != gGUIManager._MCPList->EndChildren(); ++itt)
		{
			std::string colname = (*itt)->GetElementName();
			// Always save current visibility state (overwrite if exists)
			gGUIManager._MCParticleDisplayFlag[colname] = (*itt)->GetRnrSelf();

			std::cout << "   " << colname << ": " << gGUIManager._MCParticleDisplayFlag[colname] << std::endl;
		}
		
		// Destroy old list
		gGUIManager._MCPList->DestroyElements();
		gGUIManager._MCPList->Destroy();
		gGUIManager._MCPList = nullptr;
	}

	std::cout<<"  Start to build MC Tracks collection: "<<endl;

	TEveElementList  *MCTracks = new TEveElementList();

	MCTracks->SetMainColor(kRed);
	MCTracks->SetName("MCParticles");

	TEveTrackPropagator* propsetNeutral = new TEveTrackPropagator();
	TEveTrackPropagator* propsetCharged = new TEveTrackPropagator();
	//  TEveTrackPropagator* propsetLowE = new TEveTrackPropagator();

	propsetCharged->SetMagFieldObj(new TEveMagFieldDuo(350, -gOptions.BField, 2.0));
	propsetCharged->SetName("Track propagator for charged particles");
	propsetCharged->SetMaxR(1000);
	propsetCharged->SetMaxZ(1000);
	propsetCharged->SetMaxOrbs(10.0);
	propsetCharged->SetDelta(0.01);
	//	propsetCharged->SetStepper(TEveTrackPropagator::kRungeKutta);

	propsetNeutral->SetMagFieldObj(new TEveMagFieldConst(0., 0., 0.));
	propsetNeutral->SetName("Track propagator for neutral particles");
	propsetNeutral->SetMaxR(1000);
	propsetNeutral->SetMaxZ(1000);
	propsetNeutral->SetMaxOrbs(1.0);

	//Hand Put ini
	float MCPartUnit = 0.1;
	double MCTracksMinLength = 0.5; //cm
	double MCTracksLowEThresh = 0.01;

	enum ETrType { kAucune=0, kChargedLepton, kChargedHadron, kPhoton, kNeutralHad, kNeutrino, kLowE, kLast};

	struct MCTrackDisplay {
		const char * Name;
		int Color;
		int Width;
		int Style;
	};

	MCTrackDisplay MCTParams[kLast] = {
		{"None          ",   0, 0, 0},
		{"ChargedLepton ",   kBlue, 2, 1},
		{"ChargedHadron ",   kRed, 2, 1},
		{"Photon        ",   kYellow, 1, 1},
		{"NeutralHadron ",   kGreen, 1, 1},
		{"Neutrino      ",   kGray, 1, 1},
		{"LowE          ",   15, 1, 1}
	};


	int PID, ParentNum, DaughterNum, EventNr, MotherPID, OriginPID;
	int displayedMCParticle = 0;
	int skippedMCParticle = 0;
	float energy, px, py, pz, mass, MotherEnergy, OriginEnergy, PT;
	float Vx, Vy, Vz;       //vertex position
	float Ex, Ey, Ez;       //endpoint position
	float charge;
	float KineticE, GenRadius, EndRadius;

	//Arrow to show the initial Mother Particle
	TEveTrackList* cpdLowE = new TEveTrackList("LowE");
	cpdLowE->SetMainColor(15);

	TEveTrackList* cpdChargedLepton = new TEveTrackList("Charged Lepton");
	cpdChargedLepton->SetMainColor(MCTParams[kChargedLepton].Color);

	TEveTrackList* cpdChargedHadron = new TEveTrackList("Charged Hadron");
	cpdChargedHadron->SetMainColor(MCTParams[kChargedHadron].Color);

	TEveTrackList* cpdPhoton = new TEveTrackList("Photon");
	cpdPhoton->SetMainColor(MCTParams[kPhoton].Color);

	TEveTrackList* cpdNeutralHad = new TEveTrackList("Neutral Hadron");
	cpdNeutralHad->SetMainColor(MCTParams[kNeutralHad].Color);

	TEveTrackList* cpdNeutrinos = new TEveTrackList(MCTParams[kNeutrino].Name);
	cpdNeutrinos->SetMainColor(MCTParams[kNeutrino].Color);

	//Fix missing PIDs.
	TDatabasePDG *pdgDB = TDatabasePDG::Instance();
	Int_t ionCode = 1000010020;
	if(!pdgDB->GetParticle(ionCode)){
		pdgDB->AddParticle("Deuteron","Deuteron",2+8.071e-3,kTRUE,0,1,"Ion",ionCode);
	}


	for(std::string const &collName : mcpartcollnames){
		LCCollection* col = evt->getCollection( collName ) ;
		EventNr = evt->getEventNumber();
		if(true)
		{
			int nMCParticle =  col->getNumberOfElements();
			cout<<"  Number of MCParticle: "<<nMCParticle<<endl;
			TEveTrackList * currCompound = 0;
			float EMCMax = -0.1;
			int countMother = 0; // used to identify Whizard event & Pythia event...
			for(int i=0; i<nMCParticle; i++)
			{
				MCParticle* part =  dynamic_cast<MCParticle*>( col->getElementAt( i ) ) ;
				if(part->getParents().size()==0 && part->getEnergy()>EMCMax)
				{
					if(part->getDaughters().size()==0)
					{countMother++;}
					EMCMax = part->getEnergy();
				}
			}

			for(int i=0; i<nMCParticle; i++)
			{
				MCParticle* part =  dynamic_cast<MCParticle*>( col->getElementAt( i ) ) ;

				PID=0; ParentNum=0; DaughterNum=0;
				charge=0; mass=0; energy=0;
				px=0; py=0; pz=0; Vx=0; Vy=0; Vz=0; Ex=0; Ey=0; Ez=0;
				KineticE = 0; GenRadius = 0; EndRadius = 0;
				MotherPID = 0; OriginPID = 0;
				MotherEnergy = 0; OriginEnergy = 0;
				PT = 0;

				px=part->getMomentum()[0];
				py=part->getMomentum()[1];
				pz=part->getMomentum()[2];
				PID=part->getPDG();
				Vx=part->getVertex()[0];
				Vy=part->getVertex()[1];
				Vz=part->getVertex()[2];
				Ex=part->getEndpoint()[0];
				Ey=part->getEndpoint()[1];
				Ez=part->getEndpoint()[2];
				charge=part->getCharge();
				mass=part->getMass();
				energy=part->getEnergy();	
				ParentNum=part->getParents().size();
				DaughterNum=part->getDaughters().size();
				EndRadius = sqrt(Ex*Ex+Ey*Ey);
				GenRadius = sqrt(Vx*Vx+Vy*Vy);
				TEveVector Vtx(Vx, Vy, Vz);
				TEveVector End(Ex, Ey, Ez);
				//	PT = sqrt(px*px+py*py);
				PT = energy;					//tmplate usage...

				if(PID == 22)	//Only show the information for particle/gamma with energy > 0.5GeV
				{	
					MCParticleVec mother = part->getParents();
					if (mother.size() > 0)
					{
						MotherPID = mother[0]->getPDG();
						MotherEnergy = mother[0]->getEnergy();
					}
				}

				Vtx *= MCPartUnit;
				End *= MCPartUnit;
				KineticE = sqrt(px*px+py*py+pz*pz);

				TEveTrack* track = NULL;
				ETrType TrType = kAucune;
				TEveArrow* a1 = NULL;
				int trackColor = kGray;

				float Length = Vtx.Distance(End);

				if( PT < gOptions.MCPtCut || Length<MCTracksMinLength) skippedMCParticle++;
				else displayedMCParticle++;


				if( Length<MCTracksMinLength ) continue; // Skip small tracks
				if( Length<=0) continue; // Protectin against bad parameters = 0      ??
				//  if( PID>=1000010020 ) continue;  //Mute the heavy hygen nuclea and so on.
				if( PT < gOptions.MCPtCut) continue;

				printf("MCP: charge %d, pdg %d, (%f,%f,%f)\n", int(charge), PID, px, py, pz);
				
				// First check for low energy particles (both charged and neutral)
				if(KineticE < MCTracksLowEThresh) {
					TrType = kLowE;
					currCompound = cpdLowE;
					trackColor = MCTParams[kLowE].Color;
				}
				// Charged particles
				else if(charge != 0) {
					// Classify into simplified categories
					int absPID = abs(PID);
					
					// Charged Leptons: electrons, muons, taus
					if(absPID == 11 || absPID == 13 || absPID == 15) {
						TrType = kChargedLepton;
						currCompound = cpdChargedLepton;
						trackColor = GetParticleColor(PID, i, kBlue);
					}
					// Charged Hadrons: pions, kaons, protons, and other charged particles
					else {
						TrType = kChargedHadron;
						currCompound = cpdChargedHadron;
						trackColor = GetParticleColor(PID, i, kRed);
					}
				}
			// Neutral particles
			else {
				// Neutrinos
				if(IsNeutrino(PID)) {
					TrType = kNeutrino;
					currCompound = cpdNeutrinos;
					trackColor = GetParticleColor(PID, i, kGray);
				}
				// Photons
				else if(abs(PID) == 22) {
					TrType = kPhoton;
					currCompound = cpdPhoton;
					trackColor = GetParticleColor(PID, i, kYellow);
				}
				// Neutral Hadrons: neutrons, K_L, and all other neutral hadrons
				else {
					TrType = kNeutralHad;
					currCompound = cpdNeutralHad;
					trackColor = GetParticleColor(PID, i, kGreen);
				}
			}				// Create track for charged particles
				if(TrType != kAucune && charge != 0)
				{
					propsetCharged->RefPMAtt().SetMarkerColor(kYellow);
					propsetCharged->RefPMAtt().SetMarkerStyle(kCircle);
					propsetCharged->RefPMAtt().SetMarkerSize(1.0);

					TEveRecTrack* ChargedTrack = new TEveRecTrack();
					ChargedTrack->fV.Set(Vtx);
					ChargedTrack->fP.Set(px, py, pz);
					ChargedTrack->fSign = int(charge);

					track = new TEveTrack(ChargedTrack, propsetCharged);

					// Get tracker hits for this MCParticle from TruthHelper
					const std::vector<TrackerHitInfo>& trackerHits = gTruthHelper.GetTrackerHits(part);
					
					// Add tracker hits as path marks (sorted by time), but only if they are far enough apart
					TEveVector lastPos = Vtx;  // Start from vertex
					for(const auto& hitInfo : trackerHits)
					{
						// Convert mm to cm (LCIO uses mm, ROOT uses cm)
						TEveVector hitPos(hitInfo.position[0]/10.0, hitInfo.position[1]/10.0, hitInfo.position[2]/10.0);
						
						// Only add hit if it's at least 10 cm away from the last position
						float distance = lastPos.Distance(hitPos);
						if(distance >= 1.0)
						{
							TEvePathMark* pm = new TEvePathMark(TEvePathMark::kDaughter);
							pm->fV.Set(hitPos);
							track->AddPathMark(*pm);
							lastPos = hitPos;  // Update last position
						}
					}

					// Add endpoint as decay mark
					TEvePathMark* pmEnd = new TEvePathMark(TEvePathMark::kDecay);
					pmEnd->fV.Set(End);
					track->AddPathMark(*pmEnd);
				} // charged
				// Create track for neutral particles
				else if(TrType != kAucune)
				{ // neutral
					TEveRecTrack* NeutralTrack = new TEveRecTrack();
					NeutralTrack->fV.Set(Vtx);
					NeutralTrack->fP.Set(px, py, pz);
					NeutralTrack->fSign = int(charge);

					track = new TEveTrack(NeutralTrack, propsetNeutral);

					TEvePathMark *pm = PathMarkEndTrackDecay(Vtx, End);
					track->AddPathMark(*pm);
				}

				if(track){
					track->SetName(Form("%s", gTruthHelper.GetStringID(part).c_str()));
					track->SetLineWidth(MCTParams[TrType].Width);
					track->SetLineColor(trackColor);  // Use specific particle color
					track->SetLineStyle(MCTParams[TrType].Style);
					track->SetSmooth(kTRUE);
					if(PID == 22){
						track->SetTitle(Form("%s\n"
									"Charge=%.3f, PID=%d, Energy=%.3f\n"
									"(Vx, Vy, Vz) = (%.3f, %.3f, %.3f)\n"
									"(Ex, Ey, Ez) = (%.3f, %.3f, %.3f)\n"
									"(Px, Py, Pz) = (%.3f, %.3f, %.3f)\n"
									"MotherPID = %d, MotherEnergy = %.3f",
									gTruthHelper.GetStringID(part).c_str(),
									charge, PID, energy,
									Vx, Vy, Vz, Ex, Ey, Ez, px, py, pz, MotherPID, MotherEnergy));
					}else{
						track->SetTitle(Form("%s\n"
									"Charge=%.3f, PID=%d, Energy=%.3f\n"
									"(Vx, Vy, Vz) = (%.3f, %.3f, %.3f)\n"
									"(Ex, Ey, Ez) = (%.3f, %.3f, %.3f)\n"
									"(Px, Py, Pz) = (%.3f, %.3f, %.3f)\n",
									gTruthHelper.GetStringID(part).c_str(),
									charge, PID, energy,
									Vx, Vy, Vz, Ex, Ey, Ez, px, py, pz));
				}
				if ( currCompound ) {
					currCompound->AddElement(track);
					gGUIManager._MCPTracks[part] = track;
					gGUIManager._MCPGroups[part] = currCompound;
				}
			}

			currCompound->MakeTracks();			}

			//        currCompound->CloseCompound();

			//        MCTracks->AddElement(cpdMother);
			MCTracks->AddElement(cpdLowE);
			MCTracks->AddElement(cpdNeutrinos);
			MCTracks->AddElement(cpdPhoton);
			MCTracks->AddElement(cpdNeutralHad);
			MCTracks->AddElement(cpdChargedLepton);
			MCTracks->AddElement(cpdChargedHadron);


			cpdLowE->SetRnrSelfChildren(kFALSE, kFALSE);
			cpdNeutralHad->SetRnrSelfChildren(kFALSE, kFALSE);

		} // if MCTrack collection

	} // loop over collections

	bool FlagDraw;
	for (TEveElement::List_i itt=MCTracks->BeginChildren(); itt!=MCTracks->EndChildren(); itt++){
		std::string colname = (*itt)->GetElementName();

		// Restore saved visibility state from previous event
		if (gGUIManager._MCParticleDisplayFlag.find(colname) != gGUIManager._MCParticleDisplayFlag.end())
		{
			FlagDraw = gGUIManager._MCParticleDisplayFlag[colname];
		}
		else
		{
			// Default visibility: false for LowE/NeutralHad/Neutrons, true for others
			if (colname == "LowE" || colname == "NeutralHad" || colname == "Neutrons")
			{
				FlagDraw = false;
			}
			else
			{
				FlagDraw = true;
			}
		}
		(*itt)->SetRnrSelfChildren(FlagDraw, FlagDraw);
	}
	// Apply saved top-level visibility states from previous event
	MCTracks->SetRnrSelfChildren(MCPDraw, MCPChildDraw);

	std::cout<<"  With current gOptions.MCPtCut "<<gOptions.MCPtCut<<" GeV, "<<displayedMCParticle<<" MCparticle has been displayed, and "<<skippedMCParticle<<" particles has been skipped"<<std::endl<<std::endl<<std::endl;
	return MCTracks;

}
