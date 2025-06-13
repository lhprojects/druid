#include "TRint.h"
#include "lcio.h"
#include "IO/LCReader.h"
#include "IMPL/LCTOOLS.h"
#include "EVENT/LCCollection.h"
#include "EVENT/SimCalorimeterHit.h"
#include "EVENT/CalorimeterHit.h"
#include "EVENT/LCRelation.h"
#include "EVENT/ReconstructedParticle.h"
#include "TEveElement.h"
#include "TEveStraightLineSet.h"
#include "TEveCompound.h"
#include "TEvePointSet.h"
#include "TEveBox.h"
#include "TSystem.h"
#include "TEveManager.h"
#include "TEveArrow.h"

extern int GlobalRandomColorIndex;
extern float ClusterHitSize;
extern TEveBox* BoxPhi( TVector3 &HitPos, TVector3 &Scale, int Type, int SegOrStaveNumber, float HitEnergy );

using namespace lcio;
using namespace std;

TEveCompound* RecoJets( LCCollection* col, string name)
{
	TEveCompound *Jets = new TEveCompound();
	Jets -> SetName(name.c_str());

	int NJets = col->getNumberOfElements();
	int NPart = 0;
	int NTrk = 0; 
	int NClu = 0; 
	float HitX = 0;
	float HitY = 0;
	float HitZ = 0;
	int colorindex = 0;
	TVector3 HitScale(0.5*ClusterHitSize, 0.5*ClusterHitSize, 0.5*ClusterHitSize);

	for(int i = 0; i < NJets; i++)
	{
		ReconstructedParticle * currJet = dynamic_cast<ReconstructedParticle*>( col->getElementAt(i) );
		NPart = currJet->getParticles().size();
//		colorindex = ((i%2)*50+5*i+GlobalRandomColorIndex*31)%105;
//		colorindex = 51 + (i*9 + GlobalRandomColorIndex*15)%50 ;
//		colorindex = i*9 + 51;
		colorindex = i + 2; 
		if(colorindex == 5) colorindex = 94;
		TEveCompound *a_Jet = new TEveCompound();
		a_Jet->SetName("aJet");
		a_Jet->SetMainColor(colorindex);

		for(int j = 0; j < NPart; j++)
		{
			ReconstructedParticle * currRecoP = currJet->getParticles()[j];
			NTrk = currRecoP -> getTracks().size();
			NClu = currRecoP -> getClusters().size();
			TEveElementList* TrackAsshits = new TEveElementList;
			TrackAsshits->SetName("PFOTrk");

			TEveElementList* RecoPCluster = new TEveElementList;
			RecoPCluster->SetName("PFOClu");

			for(int k = 0; k < NTrk; k++)
			{
				Track* atrack = currRecoP -> getTracks()[k];
				int AssoHitNum = atrack->getTrackerHits().size();
				float Omega = atrack->getOmega();
				float TanL = atrack->getTanLambda();
				cout<<"Omega "<<Omega<<endl;
				if(fabs(Omega) < 1e-3 || fabs(TanL) > 1)
				{
					for(int l(0); l<AssoHitNum; l++)
					{
						TrackerHit* hit = dynamic_cast<TrackerHit*>( atrack->getTrackerHits()[l] );
						HitX=hit->getPosition()[0];
						HitY=hit->getPosition()[1];
						HitZ=hit->getPosition()[2];

						TEvePointSet* q1 = new TEvePointSet(1);
						q1->SetName("Calorimeter Hit ");
						q1->SetMarkerStyle(3);
						q1->SetPoint(0, 0.1*HitX, 0.1*HitY, 0.1*HitZ);
						q1->SetMarkerColor(colorindex);

						q1->SetMarkerSize(0.1);
						TrackAsshits->AddElement(q1);
					}
				}
			}

			for(int k2 = 0; k2 < NClu; k2++)
			{
				Cluster* aclu = currRecoP -> getClusters()[k2];
				CalorimeterHitVec Hits = aclu->getCalorimeterHits();
				for(int j2 = 0; j2<Hits.size(); j2++)
				{
					TVector3 HitPosition = Hits[j2]->getPosition();
					float HitEn = Hits[j2]->getEnergy();
					HitPosition *= 0.1;
					TEveBox* q = new TEveBox();
					q = BoxPhi(HitPosition, HitScale, -1, 0, HitEn );
					q->SetMainColor(colorindex);
					q->SetLineColor(colorindex);
					RecoPCluster->AddElement(q);
				}
			}

			a_Jet->AddElement(TrackAsshits);
			a_Jet->AddElement(RecoPCluster);
		}

		Jets->AddElement(a_Jet);
	}

	return Jets; 

}
