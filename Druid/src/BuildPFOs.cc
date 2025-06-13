
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//    Build TEve Tracks on Reconstructed Particle Collection                 //
//    Color & Style according to PID                                         //
//    Unlike MCParticle, every reconstructed particle will be displayed      //
//                                                                           //
//    Date:   03 Dec 2010                                                    //
//    Author: Manqi Ruan (LLR)                                               //
//                                                                           //
//    Last Modified: 17, Dec 2011, cleaning                                  //
//                                                                           // 
///////////////////////////////////////////////////////////////////////////////


#include "TEveManager.h"
#include "TEveElement.h"
#include "TEveTrack.h"
#include "TEveArrow.h"
#include "TEveTrackPropagator.h"
#include "TEveBrowser.h"
#include "TDatabasePDG.h"
#include "TEveBox.h"
#include "TEveVSDStructs.h"

#include "lcio.h"
#include "IMPL/LCTOOLS.h"
#include "EVENT/LCCollection.h"
#include "EVENT/SimCalorimeterHit.h"
#include "EVENT/ReconstructedParticle.h"
#include "EVENT/CalorimeterHit.h"
#include "EVENT/LCEvent.h"

using namespace lcio;
using namespace EVENT;
using namespace std;

extern TEvePathMark * PathMarkEndTrackDecay(TEveVector &/*Vtx*/, TEveVector &End);
extern TEveBox* BoxPhi( TVector3 &HitPos, TVector3 &Scale, int Type, int SegOrStaveNumber, float HitEnergy );

extern int GlobalRandomColorIndex;
extern float HCALBarrelLength;
extern int flagdetectortype;
extern float ClusterHitSize;
int PFOHitColourType = 0;
float PFOHitSize = 1.0;
bool HiddenPFOCluster = 1.0;
extern int event_id;


TEveElementList* BuildPFOs( LCCollection* col, string name )
{

	std::cout<<"  Reconstructed particle collection: "<<name.c_str()<<std::endl;
	std::cout<<"  Number of PFO: "<<col->getNumberOfElements()<<std::endl;
	std::cout<<std::endl;
	//	cout<<"HCALBarrelLength: "<<HCALBarrelLength<<endl;

	TEveElementList  *RecoTracks = new TEveElementList();
	RecoTracks->SetMainColor(kRed);
	RecoTracks->SetName(name.c_str());
	TEveElementList  *RecoClus = new TEveElementList();
	RecoClus->SetMainColor(kRed);
	RecoClus->SetName("ClusterAttached");

	//    	RecoTracks->OpenCompound();

	float HBLength = 0;
	if(flagdetectortype == 3)
        {
                HBLength = 309.3;
        }
        else if(flagdetectortype < 2)
        {
                HBLength = 320.0;
        }


	TEveTrackPropagator* propsetNeutral = new TEveTrackPropagator();
	TEveTrackPropagator* propsetCharged = new TEveTrackPropagator();
	//	TEveTrackPropagator* propsetLowE = new TEveTrackPropagator();

	propsetCharged->SetMagFieldObj(new TEveMagFieldDuo(350, -3.5, 2.0));
	propsetCharged->SetName("Track propagator for charged particles");
	propsetCharged->SetMaxR(1000);
	propsetCharged->SetMaxZ(1000);
	propsetCharged->SetMaxOrbs(1.0);
	propsetCharged->SetDelta(0.01); // Step

	propsetCharged->RefPMAtt().SetMarkerColor(kYellow);
	propsetCharged->RefPMAtt().SetMarkerStyle(kCircle);
	propsetCharged->RefPMAtt().SetMarkerSize(1.0);

	propsetNeutral->SetMagFieldObj(new TEveMagFieldConst(0., 0., -3.5));
	propsetNeutral->SetName("Track propagator for neutral particles");
	propsetNeutral->SetMaxR(1000);
	propsetNeutral->SetMaxZ(1000);
	propsetNeutral->SetMaxOrbs(1.0);

	float MCPartUnit = 0.1;
	double MCTracksMinLength = 0.5; //cm
	double MCTracksLowEThresh = 0.01;

	enum ERecType { kRecAucune=0, kElectron, kPositron, kMuonP, kMuonN, kPionP, kPionN, kKaonP, kKaonN, kProton, kNeutron, kKlong, kGamma, kIonP, kIonN, kNeutralHad, kLowE, kLowETrack, kRecLast};

	struct PFODisplay {
		const char * Name;
		int Color;
		int Width;
		int Style;
		float CaloHitColor;
	};

	PFODisplay PFOParams[kRecLast] = {
		{"None       ", 0, 0, 0, 0},
		{"Electron   ", 53, 2, 9, 10},
		{"Positron   ", 98, 2, 9, 90},
		{"Muon+      ", 2, 2, 9, 95},
		{"Muon-      ", 4, 2, 9, 5},
		{"Pion+      ", 96, 2, 9, 85},
		{"Pion-      ", 66, 2, 9, 15},
		{"Kaon+      ", 6, 2, 9, 75},
		{"Kaon-      ", 7, 2, 9, 25},
		{"Proton     ", 6, 2, 9, 80},
		{"Neutron    ", 7, 1, 2, 30},
		{"Klong     ", 3, 1, 2, 40},
		{"Gamma     ", 5, 1, 2, 70},
		{"Ion+      ", 99, 2, 9, 75},
		{"Ion-	    ", 58, 4, 9, 75},
		{"NeutralHad", 5, 2, 2, 35},
		{"LowE      ", 15, 2, 2, 33},
		{"LowETrack ", 6,  2, 2, 20}
	};

	int PID, ParentNum, DaughterNum;
	float energy, px, py, pz, mass;
	float charge;
	float KineticE, GenRadius, EndRadius;
	// float clusterE;
	int Ncluster = 0; 
	int NTrack = 0; 
	int NTrackHits = 0; 

	//	TEveCompound* cpdLowE = new TEveCompound(PFOParams[kLowE].Name, "All low E tracks");
	TEveTrackList* cpdLowE = new TEveTrackList(PFOParams[kLowE].Name);
	cpdLowE->SetMainColor(PFOParams[kLowE].Color);

	TEveTrackList* cpdMuons = new TEveTrackList("Muons");
	cpdMuons->SetMainColor(PFOParams[kMuonP].Color);

	TEveTrackList* cpdPions = new TEveTrackList("Pions");
	cpdPions->SetMainColor(PFOParams[kKaonP].Color);

	TEveTrackList* cpdElectrons = new TEveTrackList("Electrons");
	cpdElectrons->SetMainColor(PFOParams[kElectron].Color);

	TEveTrackList* cpdChargedKaons = new TEveTrackList("Charged Kaons");
	cpdChargedKaons->SetMainColor(PFOParams[kKaonP].Color);

	TEveTrackList* cpdProtons = new TEveTrackList("Protons");
	cpdProtons->SetMainColor(PFOParams[kProton].Color);

	TEveTrackList* cpdNeutrons = new TEveTrackList("Neutrons");
	cpdNeutrons->SetMainColor(PFOParams[kNeutron].Color);

	TEveTrackList* cpdKlongs = new TEveTrackList("Klong");
	cpdKlongs->SetMainColor(PFOParams[kKlong].Color);

	TEveTrackList* cpdRecGamma = new TEveTrackList("Gamma");
	cpdRecGamma->SetMainColor(PFOParams[kGamma].Color);

	TEveTrackList* cpdIonP = new TEveTrackList("Undef_chargedP");	//default: charged tracks besides the defined ones
	cpdIonP->SetMainColor(PFOParams[kIonP].Color);

	TEveTrackList* cpdIonN = new TEveTrackList("Undef_chargedN");    
	cpdIonN->SetMainColor(PFOParams[kIonN].Color);

	TEveTrackList* cpdNeutralHad = new TEveTrackList("NeutralHad");
	cpdNeutralHad->SetMainColor(PFOParams[kNeutralHad].Color);

	TEveTrackList* LowETrack = new TEveTrackList("LowETrack");	//Reco charged PFO without cluster
	LowETrack->SetMainColor(PFOParams[kLowETrack].Color);

	TEveVector End(0.0, 0.0, 0.0);
	TEveVector Vtx(0.0, 0.0, 0.0);

	try{

		int nMCParticle =  col->getNumberOfElements();

		TEveTrackList * currCompound = 0;

		for(int i=0; i<nMCParticle; i++)
		{
			ReconstructedParticle* part =  dynamic_cast<ReconstructedParticle*>( col->getElementAt( i ) ) ;

			PID=0; ParentNum=0; DaughterNum=0;
			charge=0; mass=0; energy=0;
			px=0; py=0; pz=0; 
			KineticE = 0; GenRadius = 0; EndRadius = 0;

			Ncluster = part->getClusters().size();
			NTrack = part->getTracks().size();
			// clusterE =  ;
			px=part->getMomentum()[0];
			py=part->getMomentum()[1];
			pz=part->getMomentum()[2];
			KineticE = sqrt(px*px+py*py+pz*pz);

			if(part->getParticleIDs().size()>0)
			{
				PID = part->getParticleIDs()[0]->getPDG(); //First one of the PID lists
			}
			else
			{
				PID = -99; 
				PID = part->getType();
				//cout<<"PID not identified. Please check your data file."<<endl;	
			}

			charge=part->getCharge();
			mass=part->getMass();
			energy=part->getEnergy();
			Vtx = part->getReferencePoint();

			/*
			   if(Ncluster>0)
			   {
			   End = part->getClusters()[0]->getPosition();	//Test 
			   }
			   else 
			 */
			if(NTrack>0)
			{
				Track * a_trk = part->getTracks()[0];
				NTrackHits = a_trk->getTrackerHits().size();
				if(NTrackHits>0){ 
					TrackerHit * last_hit = a_trk->getTrackerHits()[NTrackHits - 1];
					End = last_hit->getPosition();
				}else{
					End = part->getMomentum() ;
					End *= 3000.0/KineticE ;
				}
			}
			else		// reserved for where cluster is dropped
			{
				End = part->getMomentum() ;
				End *= 3000.0/KineticE ;
			}

			// add the part to display calohits;

			if(Ncluster > 0 && HiddenPFOCluster)
			{
	
				TVector3 HitScale(0.5*PFOHitSize, 0.5*PFOHitSize, 0.5*PFOHitSize);

				int colorindex = 0;

				if( PFOHitColourType == 0)
				{
					if(PID == 22)
						colorindex = 5;
					else if(PID == 11)	
						colorindex = 53;
					else if(PID == -11)
						colorindex = 98;
					else if(PID == 13)
						colorindex = 2;
					else if(PID == -13)
						colorindex = 4;
					else if(PID == 211)
						colorindex = 96;
					else if(PID == -211)
						colorindex = 66; 
					else if(PID == 310 || PID == 130 || PID == 2112)
						colorindex = 3;
					else
						colorindex = 10; 
				}
				else if( PFOHitColourType == 2 )
				{
					colorindex = 3;
				}
				else if ( PFOHitColourType == 1 )
				{
					colorindex = (i*11 + GlobalRandomColorIndex*13)%50 + 51;
				}

				for(int i = 0; i < Ncluster; i++)
				{
					Cluster * a_clu = part->getClusters()[i];
					int CluSize = a_clu->getCalorimeterHits().size();

					TEveElementList* a_pfocluster = new TEveElementList;
					a_pfocluster->SetMainColor(i);
					a_pfocluster->SetName(Form("PFOClu, En = %.6f", a_clu->getEnergy()));

					for(int j = 0; j<CluSize; j++)
					{
						TVector3 HitPosition = a_clu->getCalorimeterHits()[j]->getPosition();
						float HitEn = a_clu->getCalorimeterHits()[j]->getEnergy();
						HitPosition *= 0.1;

						TEveBox* q = new TEveBox();

						if(HitPosition[2]>HBLength || HitPosition[2]<-1*HBLength)
						{
							q = BoxPhi(HitPosition, HitScale, -1, 0, HitEn );
						}
						else
						{
							q = BoxPhi(HitPosition, HitScale, 1, 8, HitEn );
						}

						q->SetTitle(Form( "Type = %d, Cluster En = %.6f GeV, PFO En = %.6f GeV", PID, a_clu->getEnergy(), energy ));

                                        	q->SetPickable(kTRUE);

						// int colorindex = (j*11 + GlobalRandomColorIndex*13)%50 + 51;
						q->SetMainColor(colorindex);
						q->SetLineColor(colorindex);
						a_pfocluster->AddElement(q);

					}

					RecoClus->AddElement(a_pfocluster);
				}
			}

			int PFOshowhit = 1;

			Vtx *= MCPartUnit;
			End *= MCPartUnit;

			TEveTrack* track = 0;

			ERecType TrType = kRecAucune;

			//	//used for which has not pass PID

			if(PID == 0) 
			{
				if(charge > 0) PID = 211; 
				else if(charge < 0) PID = -211; 
				else PID = 22; 
			}

			if(charge!=0 && KineticE >= MCTracksLowEThresh && Ncluster){

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
						TrType = kIonP;
						if(charge > 0 ) currCompound = cpdIonP;
						else if(charge < 0) currCompound = cpdIonN;
						break;

				}

				TEveRecTrack* ChargedTrack = new TEveRecTrack();
				ChargedTrack->fV.Set(Vtx);
				ChargedTrack->fP.Set(px, py, pz);
				ChargedTrack->fSign = int(charge);

				track = new TEveTrack(ChargedTrack, propsetCharged);

				if(currCompound != cpdMuons)	//Any track besides Muon track will be end at the first calorimeter hit it corresponding to 
				{
					TEvePathMark* pm = new TEvePathMark(TEvePathMark::kDecay);
					pm->fV.Set(End);
					track->AddPathMark( *pm );
				}

			} 
			else 
			{

				if( KineticE < MCTracksLowEThresh )
				{
					TrType = kLowE;
					currCompound = cpdLowE;
				}
				else if( Ncluster == 0)
				{
					TrType = kLowETrack;
					currCompound = LowETrack; 
				}
				else
				{

					switch( abs(PID) )
					{
						case  12:; case  14:; case  16:;    //Neutrinos	actually never be reconstructed in ILD I guess
							 TrType = kRecAucune;
							 currCompound = cpdNeutralHad;
							 break;

						case 22:    
							 TrType = kGamma;
							 currCompound = cpdRecGamma;
							 break;

						case 2112:
							 TrType = kNeutron;
							 currCompound = cpdNeutrons;		
							 break;

						case 130:
							 TrType = kKlong;
							 currCompound = cpdKlongs;
							 break;

						default:    
							 TrType = kNeutralHad;	
							 currCompound = cpdNeutralHad;
							 break;
					}
				}
				if( TrType != kRecAucune )
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

			if(track && currCompound)
			{
				track->SetName(Form("Track %d, Energy=%f\n", i, energy));    // i = tracknum
				track->SetLineWidth(PFOParams[TrType].Width);
				track->SetLineColor(PFOParams[TrType].Color);
				track->SetLineStyle(PFOParams[TrType].Style);
				track->SetSmooth(kTRUE);
				track->SetTitle(Form("Reconstructed PFOs: %s\n"
							"EventNr=%d, Track No.=%d\n""Charge=%f, PID=%d\n"
							"Energy=%f\n"
							"Vtx position= (%.3f, %.3f, %.3f)\n"
							"Cluster pos = (%.3f, %.3f, %.3f)\n"
							"3-momentum = (%.3f, %.3f, %.3f)",
							name.c_str(), event_id, i, charge, PID, energy,
							10*Vtx[0], 10*Vtx[1], 10*Vtx[2], 10*End[0], 10*End[1], 10*End[2], px, py, pz));

				currCompound->AddElement(track);
			}

			currCompound->IncDenyDestroy();
			currCompound->MakeTracks();

		}

		RecoTracks->AddElement(LowETrack);
		RecoTracks->AddElement(cpdLowE);
		RecoTracks->AddElement(cpdRecGamma);
		RecoTracks->AddElement(cpdNeutralHad);
		RecoTracks->AddElement(cpdMuons);
		RecoTracks->AddElement(cpdPions);
		RecoTracks->AddElement(cpdElectrons);
		RecoTracks->AddElement(cpdChargedKaons);
		RecoTracks->AddElement(cpdProtons);
		RecoTracks->AddElement(cpdKlongs);
		RecoTracks->AddElement(cpdNeutrons);
		RecoTracks->AddElement(cpdIonP);
		RecoTracks->AddElement(cpdIonN);
		RecoTracks->AddElement(RecoClus);

		return RecoTracks;

	}
	catch(lcio::DataNotAvailableException zero) { }

}

