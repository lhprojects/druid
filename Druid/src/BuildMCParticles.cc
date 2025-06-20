
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

using namespace lcio;
using namespace EVENT;
using namespace std;

extern std::map <string, bool> MCParticleDisplayFlag;

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


TEveElementList* BuildMCParticles( LCEvent *evt )
{
	std::cout<<"  "<<endl;
	std::cout<<"  Start to build MC Tracks collection: "<<endl;

	TEveElementList  *MCTracks = new TEveElementList();

	MCTracks->SetMainColor(kRed);
	MCTracks->SetName("MC Particles");

	TEveTrackPropagator* propsetNeutral = new TEveTrackPropagator();
	TEveTrackPropagator* propsetCharged = new TEveTrackPropagator();
	//  TEveTrackPropagator* propsetLowE = new TEveTrackPropagator();

	propsetCharged->SetMagFieldObj(new TEveMagFieldDuo(350, -3.5, 2.0));
	propsetCharged->SetName("Track propagator for charged particles");
	propsetCharged->SetMaxR(1000);
	propsetCharged->SetMaxZ(1000);
	propsetCharged->SetMaxOrbs(10.0);
	propsetCharged->SetDelta(0.01);
	//	propsetCharged->SetStepper(TEveTrackPropagator::kRungeKutta);

	propsetNeutral->SetMagFieldObj(new TEveMagFieldConst(0., 0., -3.5));
	propsetNeutral->SetName("Track propagator for neutral particles");
	propsetNeutral->SetMaxR(1000);
	propsetNeutral->SetMaxZ(1000);
	propsetNeutral->SetMaxOrbs(1.0);

	//Hand Put ini
	float MCPartUnit = 0.1;
	double MCTracksMinLength = 0.5; //cm
	double MCTracksLowEThresh = 0.01;

	enum ETrType { kAucune=0, kElectron, kPositron, kMuonP, kMuonN, kPionP, kPionN, kKaonP, kKaonN, kProton, kNeutron, kKlong, kGamma, kIon, kNeutralHad, kNeutrino, kLowE, kLast};

	struct MCTrackDisplay {
		const char * Name;
		int Color;
		int Width;
		int Style;
	};

	MCTrackDisplay MCTParams[kLast] = {
		{"None       ",   0, 0, 0},
		{"Electron   ",   98, 2, 1},
		{"Positron   ",   53, 2, 1},
		{"Muon+      ",   2, 2, 1},
		{"Muon-      ",   4, 2, 1},
		{"Pion+      ",   96, 2, 1},
		{"Pion-      ",   66, 2, 1},
		{"Kaon+      ",   6, 2, 1},
		{"Kaon-      ",   7, 2, 1},
		{"Proton     ",   6, 2, 1},
		{"Neutron    ",   7, 1, 1},
		{"Klong     ",   3, 1, 1},
		{"Gamma     ",   5, 1, 1},
		{"Ion       ",   15, 1, 1},
		{"NeutralHad",   5, 1, 1},
		{"Neutrino  ",   7, 1, 1},
		{"LowE      ",   15, 1, 1}
	};


	std::string MCTrackName;
	MCTrackName="MCParticle";
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

	TEveTrackList* cpdMuons = new TEveTrackList("Muons");
	cpdMuons->SetMainColor(MCTParams[kMuonP].Color);

	TEveTrackList* cpdPions = new TEveTrackList("Pions");
	cpdPions->SetMainColor(MCTParams[kKaonP].Color);

	TEveTrackList* cpdElectrons = new TEveTrackList("Electrons");
	cpdElectrons->SetMainColor(MCTParams[kElectron].Color);

	TEveTrackList* cpdChargedKaons = new TEveTrackList("Charged Kaons");
	cpdChargedKaons->SetMainColor(MCTParams[kKaonP].Color);

	TEveTrackList* cpdProtons = new TEveTrackList("Protons");
	cpdProtons->SetMainColor(MCTParams[kProton].Color);

	TEveTrackList* cpdNeutrons = new TEveTrackList("Neutrons");
	cpdNeutrons->SetMainColor(MCTParams[kNeutron].Color);

	TEveTrackList* cpdKlongs = new TEveTrackList("Klong");
	cpdKlongs->SetMainColor(MCTParams[kKlong].Color);

	TEveTrackList* cpdGamma = new TEveTrackList("Gamma");
	cpdGamma->SetMainColor(MCTParams[kGamma].Color);

	TEveTrackList* cpdIon = new TEveTrackList("Ion");       //default: charged tracks besides the defined ones
	cpdIon->SetMainColor(MCTParams[kIon].Color);

	TEveTrackList* cpdNeutralHad = new TEveTrackList("NeutralHad");
	cpdNeutralHad->SetMainColor(MCTParams[kNeutralHad].Color);

	TEveTrackList* cpdNeutrinos = new TEveTrackList(MCTParams[kNeutrino].Name);
	cpdNeutrinos->SetMainColor(MCTParams[kNeutrino].Color);

	//Fix missing PIDs.
	TDatabasePDG *pdgDB = TDatabasePDG::Instance();
	Int_t ionCode = 1000010020;
	if(!pdgDB->GetParticle(ionCode)){
		pdgDB->AddParticle("Deuteron","Deuteron",2+8.071e-3,kTRUE,0,1,"Ion",ionCode);
	}

	std::vector<std::string>::const_iterator name;

	const std::vector< std::string >* strVec = evt->getCollectionNames() ;

	for( name = strVec->begin() ; name != strVec->end() ; name++){
		LCCollection* col = evt->getCollection( *name ) ;
		EventNr = evt->getEventNumber();
		if(*name == MCTrackName)
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

				if(PID == 22 && energy > 0.5)	//Only show the information for particle/gamma with energy > 0.5GeV
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

				float Length = Vtx.Distance(End);

				if( PT < gOptions.MCPtCut || Length<MCTracksMinLength) skippedMCParticle++;
				else displayedMCParticle++;


				if( Length<MCTracksMinLength ) continue; // Skip small tracks
				if( Length<=0) continue; // Protectin against bad parameters = 0      ??
				//  if( PID>=1000010020 ) continue;  //Mute the heavy hygen nuclea and so on.
				if( PT < gOptions.MCPtCut) continue;

				if(charge!=0 && KineticE >= MCTracksLowEThresh){

					switch(PID){
						case 11:
							TrType = kElectron;
							currCompound = cpdElectrons;
							break;

						case -11:
							TrType = kPositron;
							currCompound = cpdElectrons;
							break;

						case 13:
							TrType = kMuonN;
							currCompound = cpdMuons;
							break;

						case -13:
							TrType = kMuonP;
							currCompound = cpdMuons;
							break;

						case 211:
							TrType = kPionP;
							currCompound = cpdPions;
							break;

						case -211:
							TrType = kPionN;
							currCompound = cpdPions;
							break;

						case 321:
							TrType = kKaonP;
							currCompound = cpdChargedKaons;
							break;

						case -321:
							TrType = kKaonN;
							currCompound = cpdChargedKaons;
							break;

						case 2212:
							TrType = kProton;
							currCompound = cpdProtons;
							break;

						default:
							TrType = kIon; 
							currCompound = cpdIon;
							break;
					}

					propsetCharged->RefPMAtt().SetMarkerColor(kYellow);
					propsetCharged->RefPMAtt().SetMarkerStyle(kCircle);
					propsetCharged->RefPMAtt().SetMarkerSize(1.0);

					TEveRecTrack* ChargedTrack = new TEveRecTrack();
					ChargedTrack->fV.Set(Vtx);
					ChargedTrack->fP.Set(px, py, pz);
					ChargedTrack->fSign = int(charge);

					track = new TEveTrack(ChargedTrack, propsetCharged);

					TEvePathMark* pm1 = new TEvePathMark(TEvePathMark::kDaughter);
					TEvePathMark* pm2 = new TEvePathMark(TEvePathMark::kDaughter);

					TEvePathMark* pm3 = new TEvePathMark(TEvePathMark::kDecay);


					if( (Vz<2350 && Vz>-2350 && GenRadius<1810) && (Ez>2350 || Ez<-2350 || EndRadius>1810) )   // if cross the board of TPC
					{
						std::string SETHitCollection = "SETCollection";
						try{
							LCCollection* col = evt->getCollection( SETHitCollection ) ;
							int nHits = col->getNumberOfElements();
							int count = 0;
							for(int j=0; j<nHits; j++)
							{
								SimTrackerHit* hit = dynamic_cast<SimTrackerHit*>( col->getElementAt(j) );
								MCParticle* hitMCPart = dynamic_cast<MCParticle*>( hit->getMCParticle());
								if(hitMCPart==part && count==0)
								{
									TEveVector SetHit(hit->getPosition()[0]/10.0, hit->getPosition()[1]/10.0, hit->getPosition()[2]/10.0);
									pm1->fV.Set(SetHit);
									track->AddPathMark(*pm1);
									count=1;
								}
							}
						}catch (lcio::DataNotAvailableException zero) { }

					}


					if ( (Vz<3381.6 && Vz>-3381.6 && GenRadius<3973.6) && (Ez>3381.6 || Ez<-3381.6 || EndRadius > 3973.6) )   // if end outside the Calo
					{
						float MuonCaloHitDis = 0;
						float EndPointDisMax = 0;
						float Xmax = 0;
						float Ymax = 0;
						float Zmax = 0;

						const std::vector< std::string >* strVec = evt->getCollectionNames() ;

						for( std::vector<std::string>::const_iterator  name2 = strVec->begin() ; name2 != strVec->end() ; name2++){
							try{
								LCCollection* col = evt->getCollection( *name2 ) ;
								string SubD (*name2, 0, 4);
								if ( SubD=="Muon" )
								{
									int nMuonHits = col->getNumberOfElements();
									for(int i = 0; i<nMuonHits; i++)
									{
										SimCalorimeterHit* hit11 = dynamic_cast<SimCalorimeterHit*>( col->getElementAt(i) );
										MCParticle* hitMCPart11 = dynamic_cast<MCParticle*>( hit11->getParticleCont(0));
										if( hitMCPart11==part )
										{ 
											MuonCaloHitDis = sqrt(hit11->getPosition()[0]*hit11->getPosition()[0]+hit11->getPosition()[1]*hit11->getPosition()[1]+hit11->getPosition()[2]*hit11->getPosition()[2]); 
											if(MuonCaloHitDis>EndPointDisMax)
											{
												EndPointDisMax = MuonCaloHitDis;
												Xmax = hit11->getPosition()[0]/10.0;
												Ymax = hit11->getPosition()[1]/10.0;
												Zmax = hit11->getPosition()[2]/10.0;
											}
										}
									}
								}
							}catch (lcio::DataNotAvailableException zero) { }

						}

						TEveVector MuonFarHit(Xmax, Ymax, Zmax);
						if(EndPointDisMax>10)
						{
							pm2->fV.Set(MuonFarHit);
							track->AddPathMark(*pm2);
						}
					}

					pm3->fV.Set(End);
					track->AddPathMark(*pm3);


				} // charged
				else 
				{ // neutral
					/*
					   propsetNeutral->SetMagFieldObj(new TEveMagFieldConst(0., 0., -3.5));
					   propsetNeutral->SetName("Track propagator for neutral particles");
					   propsetNeutral->SetMaxR(1000);
					   propsetNeutral->SetMaxZ(1000);
					   propsetNeutral->SetMaxOrbs(1.0);
					   */
					if( KineticE < MCTracksLowEThresh && ! IsNeutrino(PID) )
					{
						TrType = kLowE;
						currCompound = cpdLowE;
					}else{

						switch( abs(PID) ){
							case  12:   // pass through
							case  14:   // pass through
							case  16:   // pass through  //Neutrinos
									 TrType = kNeutrino;
									 currCompound = cpdNeutrinos;
									 break;

							case 22:    //Gammas
									 TrType = kGamma;
									 currCompound = cpdGamma;
									 std::cout<<"  Displaying Gamma with energy: "<<energy<<std::endl;
									 break;

							case 2112:
									 TrType = kNeutron;
									 currCompound = cpdNeutrons;
									 break;

							case 130:
									 TrType = kKlong;
									 currCompound = cpdKlongs;
									 break;

							default:    //All neutral hadrons
									 TrType = kNeutralHad;
									 currCompound = cpdNeutralHad;
									 break;
						}
					}
					if( TrType != kAucune )
					{

						TEveRecTrack* NeutralTrack = new TEveRecTrack();
						NeutralTrack->fV.Set(Vtx);
						NeutralTrack->fP.Set(px, py, pz);
						NeutralTrack->fSign = int(charge);

						track = new TEveTrack(NeutralTrack, propsetNeutral);

						TEvePathMark *pm = PathMarkEndTrackDecay(Vtx, End);
						track->AddPathMark(*pm);
					}

				}

				if(track){
					track->SetName(Form("Track %d", i));    // i = tracknum
					track->SetLineWidth(MCTParams[TrType].Width);
					track->SetLineColor(MCTParams[TrType].Color);
					track->SetLineStyle(MCTParams[TrType].Style);
					track->SetSmooth(kTRUE);
					if(PID == 22){
						track->SetTitle(Form("MCParticles: \n"
									"EventNr=%d, Track No.=%d\n""Charge=%.3f, PID=%d, Energy=%.3f\n"
									"(Vx, Vy, Vz) = (%.3f, %.3f, %.3f)\n"
									"(Ex, Ey, Ez) = (%.3f, %.3f, %.3f)\n"
									"(Px, Py, Pz) = (%.3f, %.3f, %.3f)\n"
									"MotherPID = %d, MotherEnergy = %.3f",
									EventNr, i, charge, PID, energy,
									Vx, Vy, Vz, Ex, Ey, Ez, px, py, pz, MotherPID, MotherEnergy));
					}else{
						track->SetTitle(Form("MCParticles: \n"
									"EventNr=%d, Track No.=%d\n""Charge=%.3f, PID=%d, Energy=%.3f\n"
									"(Vx, Vy, Vz) = (%.3f, %.3f, %.3f)\n"
									"(Ex, Ey, Ez) = (%.3f, %.3f, %.3f)\n"
									"(Px, Py, Pz) = (%.3f, %.3f, %.3f)\n",
									EventNr, i, charge, PID, energy,
									Vx, Vy, Vz, Ex, Ey, Ez, px, py, pz));
					}
					if ( currCompound ) {
						currCompound->AddElement(track);
					}
				}

				currCompound->MakeTracks();

			}

			//        currCompound->CloseCompound();

			//        MCTracks->AddElement(cpdMother);
			MCTracks->AddElement(cpdLowE);
			MCTracks->AddElement(cpdNeutrinos);
			MCTracks->AddElement(cpdGamma);
			MCTracks->AddElement(cpdNeutralHad);
			MCTracks->AddElement(cpdMuons);
			MCTracks->AddElement(cpdPions);
			MCTracks->AddElement(cpdElectrons);
			MCTracks->AddElement(cpdChargedKaons);
			MCTracks->AddElement(cpdProtons);
			MCTracks->AddElement(cpdNeutrons);
			MCTracks->AddElement(cpdKlongs);
			MCTracks->AddElement(cpdIon);

			cpdLowE->SetRnrSelfChildren(kFALSE, kFALSE);
			cpdNeutralHad->SetRnrSelfChildren(kFALSE, kFALSE);
			cpdNeutrons->SetRnrSelfChildren(kFALSE, kFALSE);

		} // if MCTrack collection

	} // loop over collections

	bool FlagDraw;
	for (TEveElement::List_i itt=MCTracks->BeginChildren(); itt!=MCTracks->EndChildren(); itt++){
		std::string colname = (*itt)->GetElementName();
		if(colname!="LowE" && colname!="NeutralHad" && colname!="Neutrons"){
			if(MCParticleDisplayFlag.find(colname)!=MCParticleDisplayFlag.end()) FlagDraw=MCParticleDisplayFlag[colname];
			else FlagDraw = true;
			(*itt)->SetRnrSelfChildren(FlagDraw, FlagDraw);
		}
	}

	std::cout<<"  With current PTCut "<<gOptions.MCPtCut<<" GeV, "<<displayedMCParticle<<" MCparticle has been displayed, and "<<skippedMCParticle<<" particles has been skipped"<<std::endl<<std::endl<<std::endl;
	return MCTracks;

}
