
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
#include "GlobalDefs.hh"

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


TVector3 computeVtx(const TrackState* state) {
	double refX = state->getReferencePoint()[0];
	double refY = state->getReferencePoint()[1];
	double refZ = state->getReferencePoint()[2];
	double phi = state->getPhi();
	double d0 = state->getD0();
	float x0 = refX - d0 * sin(phi);
	float y0 = refY + d0 * cos(phi);
	float z0 = refZ + z0;
	return TVector3(x0, y0, z0);
}


TVector3 computeMomentum(const TrackState* state, double Bz = 3.5 /* T */)
{
	const double kAlpha = 2.99792458e-4;     
    if (!state) return TVector3();              // protect against null

    const double omega      = state->getOmega();          // 1/mm
    if (omega == 0.0) return TVector3();        // pT → ∞ : undefined

    const double pT         = kAlpha * fabs(Bz) / fabs(omega);
    const double phi        = state->getPhi();
    const double tanLambda  = state->getTanLambda();

    const double px = pT * std::cos(phi);
    const double py = pT * std::sin(phi);
    const double pz = pT * tanLambda;

    return TVector3(px, py, pz);
}

inline int chargeSign(const TrackState* st, double Bz = 3.5 /* T */)
{
    if (!st)                    return 0;        // null-check
    const double omega = st->getOmega();
    const int sOmega = (omega > 0) ? +1 : -1;
    return sOmega;                         // +1 or –1
}



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

	enum ERecType { kRecAucune=0, kElectron, kPositron, kMuonP, kMuonN, kPionP, kPionN, kKaonP, kKaonN, 
		kProton, kNeutron, kKlong, kGamma, kIonP, kIonN, kNeutralHad, kLowE, kLowETrack, kSubTrack, kRecLast};

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
		{"LowETrack ", 6,  2, 2, 20},
		{"SubTrack", 53, 2, 9, 10}
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

	TEveTrackList* SubTrack = new TEveTrackList("SubTrack");	//Reco charged PFO without cluster
	SubTrack->SetMainColor(PFOParams[kSubTrack].Color);

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
			px=0; py=0; pz=0;                       // particle momentum
			double trkPxAtVtx = 0, trkPyAtVtx = 0, trkPzAtVtx = 0; // track momentum
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
			if(NTrack==1 || (charge && NTrack>1))
			{
				// FIXME: for charge and NTrack>1, first is the parent, is this what we need?
				Track * a_trk = part->getTracks()[0];
				NTrackHits = a_trk->getTrackerHits().size();
				if(NTrackHits>0){ 
					TrackerHit * last_hit = a_trk->getTrackerHits()[NTrackHits - 1];
					End = last_hit->getPosition();
				}else{
					End = part->getMomentum() ;
					End *= 3000.0/KineticE ;
					End = Vtx + End;
				}
			}
			else		// reserved for where cluster is dropped
			{
				End = part->getMomentum() ;
				End *= 3000.0/KineticE ;
				End = Vtx + End;
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

			printf("PFO %d, PID=%d, Charge=%f, Energy=%f, mom=(%f,%f,%f),"
					"KineticE=%f, Vtx=(%.3f, %.3f, %.3f), End=(%.3f, %.3f, %.3f)\n",
			       i, PID, charge, energy, 
					px, py, pz,
				   KineticE,
			       Vtx[0], Vtx[1], Vtx[2],
			       End[0], End[1], End[2]);

			for(Track *trk : part->getTracks())
			{
				if(trk->getTrackStates().size() > 0)
				{					
					TrackState const *s = trk->getTrackState(TrackState::AtFirstHit);
					double ch = chargeSign(s, -3.5);
					TVector3 trkMom = computeMomentum(s, -3.5);
					printf("Track %d: charge %f, momentum (%f,%f,%f)\n", i, ch, trkMom.X(), trkMom.Y(), trkMom.Z());
				}

				TrackStateVec const &trkStateVec = trk->getTrackStates();
				bool VtxSet = false;
				bool EndSet = false;
				double trkPx = px;
				double trkPy = py;
				double trkPz = pz;
				TEveVector trkVtx = Vtx;
				TEveVector trkEnd = End;
				double trkCharge = 0;

				printf("Global Vtx %f %f %f\n", Vtx.fX, Vtx.fY, Vtx.fZ);
				printf("Global End %f %f %f\n", End.fX, End.fY, End.fZ);
				printf("Global Momentum %f(%f,%f,%f)\n", sqrt(px*px + py*py + pz*pz), 
						px, py, pz);

				for (TrackState *trkState : trkStateVec)
				{
					TVector3 trkMom = computeMomentum(trkState, -3.5);
					int trkstateCharge = chargeSign(trkState, -3.5);
					double mom = trkMom.Mag();
					if (trkState->getLocation() == TrackState::AtFirstHit)
					{
						TVector3 v = computeVtx(trkState);
						trkVtx = TEveVector(v.X(), v.Y(), v.Z());
						trkVtx = trkVtx * MCPartUnit;
						trkPx = trkMom.X();
						trkPy = trkMom.Y();
						trkPz = trkMom.Z();
						trkCharge = trkstateCharge;
						printf("TrackState at first hit: %f %f %f %f\n", trkCharge, trkVtx.fX, trkVtx.fY, trkVtx.fZ);
						printf("momentum: %f(%f,%f,%f)\n", mom, trkMom.X(), trkMom.Y(), trkMom.Z());
						VtxSet = true;
					}
					else if (trkState->getLocation() == TrackState::AtLastHit)
					{
						trkEnd = trkState->getReferencePoint();
						trkEnd = trkEnd * MCPartUnit;
						printf("TrackState at last hit: %f %f %f %f\n", trkCharge, trkEnd.fX, trkEnd.fY, trkEnd.fZ);
						printf("momentum: %f(%f,%f,%f)\n", mom, trkMom.X(), trkMom.Y(), trkMom.Z());
						EndSet = true;
					}
					else if (!VtxSet && trkState->getLocation() == TrackState::AtIP)
					{
						TVector3 v = computeVtx(trkState);
						trkVtx = TEveVector(v.X(), v.Y(), v.Z());
						trkVtx = trkVtx * MCPartUnit;
						trkPx = trkMom.X();
						trkPy = trkMom.Y();
						trkPz = trkMom.Z();
						trkCharge = trkstateCharge;
						printf("TrackState at IP: %f %f %f %f\n", trkCharge, trkVtx.fX, trkVtx.fY, trkVtx.fZ);
						printf("momentum: %f(%f,%f,%f)\n", mom, trkMom.X(), trkMom.Y(), trkMom.Z());
						VtxSet = true;
					}
					else if (!EndSet && trkState->getLocation() == TrackState::AtCalorimeter)
					{
						trkEnd = trkState->getReferencePoint();
						trkEnd = trkEnd * MCPartUnit;
						printf("TrackState at Calorimeter: %f %f %f %f\n", trkCharge, trkEnd.fX, trkEnd.fY, trkEnd.fZ);
						//printf("D0 Z0: %f %f\n", trkState->getD0(), trkState->getZ0());
						printf("momentum: %f(%f,%f,%f)\n", mom, trkMom.X(), trkMom.Y(), trkMom.Z());
						EndSet = true;
					}
				}
			}

			if (NTrack >= 1)
			{
				TrackVec const &trks = part->getTracks();
				int trkIdx = 0;
				for (auto *trk : trks)
				{

					TrackStateVec const &trkStateVec = trk->getTrackStates();
					bool VtxSet = false;
					bool EndSet = false;
					double trkPx = px;
					double trkPy = py;
					double trkPz = pz;
					TEveVector trkVtx = Vtx;
					TEveVector trkEnd = End;
					double trkCharge = 0;

					printf("Global Vtx %f %f %f\n", Vtx.fX, Vtx.fY, Vtx.fZ);
					printf("Global End %f %f %f\n", End.fX, End.fY, End.fZ);
					printf("Global Momentum %f(%f,%f,%f)\n", sqrt(px*px + py*py + pz*pz), 
							px, py, pz);

					for (TrackState *trkState : trkStateVec)
					{
						TVector3 trkMom = computeMomentum(trkState, -3.5);
						int trkstateCharge = chargeSign(trkState, -3.5);
						TVector3 location = computeVtx(trkState);						
						double mom = trkMom.Mag();
						if (trkState->getLocation() == TrackState::AtFirstHit)
						{
							trkVtx = TEveVector(location.X(), location.Y(), location.Z());
							trkVtx = trkVtx * MCPartUnit;
							trkPx = trkMom.X();
							trkPy = trkMom.Y();
							trkPz = trkMom.Z();
							trkCharge = trkstateCharge;
							printf("TrackState at first hit: %f %f %f %f\n", trkCharge, trkVtx.fX, trkVtx.fY, trkVtx.fZ);
							//printf("D0 Z0: %f %f\n", trkState->getD0(), trkState->getZ0());
							//printf("momentum: %f(%f,%f,%f)\n", mom, trkPx, trkPy, trkPz);
							VtxSet = true;
						}
						else if (trkState->getLocation() == TrackState::AtLastHit)
						{
							trkEnd = TEveVector(location.X(), location.Y(), location.Z());
							trkEnd = trkEnd * MCPartUnit;
							printf("TrackState at last hit: %f %f %f %f\n", trkCharge, trkEnd.fX, trkEnd.fY, trkEnd.fZ);
							//printf("D0 Z0: %f %f\n", trkState->getD0(), trkState->getZ0());
							//printf("momentum: %f(%f,%f,%f)\n", mom, trkPx, trkPy, trkPz);
							EndSet = true;
						}
						else if (!VtxSet && trkState->getLocation() == TrackState::AtIP)
						{
							trkVtx = TEveVector(location.X(), location.Y(), location.Z());
							trkVtx = trkVtx * MCPartUnit;
							trkPx = trkMom.X();
							trkPy = trkMom.Y();
							trkPz = trkMom.Z();
							trkCharge = trkstateCharge;
							printf("TrackState at IP: %f %f %f %f\n", trkCharge, trkVtx.fX, trkVtx.fY, trkVtx.fZ);
							//printf("D0 Z0: %f %f\n", trkState->getD0(), trkState->getZ0());
							//printf("momentum: %f(%f,%f,%f)\n", mom, trkPx, trkPy, trkPz);
							VtxSet = true;
						}
						else if (!EndSet && trkState->getLocation() == TrackState::AtCalorimeter)
						{
							trkEnd = TEveVector(location.X(), location.Y(), location.Z());
							trkEnd = trkEnd * MCPartUnit;

							printf("TrackState at Calorimeter: %f %f %f %f\n", trkCharge, trkEnd.fX, trkEnd.fY, trkEnd.fZ);
							//printf("D0 Z0: %f %f\n", trkState->getD0(), trkState->getZ0());
							//printf("momentum: %f(%f,%f,%f)\n", mom, trkPx, trkPy, trkPz);
							EndSet = true;
						}
					}

					printf("---------------------\n");
					printf("Track %d: trkCharge %f, Vtx %f %f %f, End %f %f %f, Momentum %f(%f,%f,%f)\n",
							trkIdx, trkCharge, trkVtx.fX, trkVtx.fY, trkVtx.fZ,
							trkEnd.fX, trkEnd.fY, trkEnd.fZ,
							sqrt(trkPx*trkPx + trkPy*trkPy + trkPz*trkPz),
							trkPx, trkPy, trkPz);

					if (charge && trkIdx == 0)
					{  // skip the first track if it is charged, because it is the PFO itself
						Vtx = trkVtx;
						End = trkEnd;
						trkPxAtVtx = trkPx;
						trkPyAtVtx = trkPy;
						trkPzAtVtx = trkPz;
						++trkIdx;
						//continue;
					}

					TEveRecTrack* ChargedTrack = new TEveRecTrack();
					ChargedTrack->fV.Set(trkVtx);
					ChargedTrack->fP.Set(trkPx, trkPy, trkPz);
					ChargedTrack->fSign = int(trkCharge);

					track = new TEveTrack(ChargedTrack, propsetCharged);

					TEvePathMark* pm = new TEvePathMark(TEvePathMark::kDecay);
					pm->fV.Set(trkEnd);
					//track->AddPathMark( *pm );

					TrType = kSubTrack;
					track->SetName(Form("Sub Track %d:%d, Energy=%f\n", i, trkIdx, energy));    // i = tracknum
					track->SetLineWidth(PFOParams[TrType].Width);
					track->SetLineColor(PFOParams[TrType].Color);
					track->SetLineStyle(PFOParams[TrType].Style);
					track->SetSmooth(kTRUE);
					track->SetTitle(Form("Sub Track: %s, Energy=%f\n\n"
								"EventNr=%d, Track No.=%d\n"
								"Charge=%f, PID=%d\n"
								"Energy=%f\n"
								"Vtx = (%.3f, %.3f, %.3f) cm\n"
								"End = (%.3f, %.3f, %.3f) cm\n"
								"Momentum = (%.3f, %.3f, %.3f)",
								name.c_str(), energy, gDisplayState.getEventNumber(), i, (double)trkCharge, 
								PID, energy,
								10*trkVtx[0], 10*trkVtx[1], 10*trkVtx[2],
								10*trkEnd[0], 10*trkEnd[1], 10*trkEnd[2],
								trkPx, trkPy, trkPz));

					SubTrack->AddElement(track);
					SubTrack->IncDenyDestroy();
					SubTrack->MakeTracks();
					++trkIdx;

				}

			}


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
						TrType = kIonP;
						if(charge > 0 ) currCompound = cpdIonP;
						else if(charge < 0) currCompound = cpdIonN;
						break;

				}


				TEveRecTrack* ChargedTrack = new TEveRecTrack();
				ChargedTrack->fV.Set(Vtx);
				ChargedTrack->fP.Set(trkPxAtVtx, trkPyAtVtx, trkPzAtVtx);
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
				else if( Ncluster == 0 && false)
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
							name.c_str(), gDisplayState.getEventNumber(), i, charge, PID, energy,
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
		RecoTracks->AddElement(SubTrack);

		return RecoTracks;

	}
	catch(lcio::DataNotAvailableException zero) { }

}

