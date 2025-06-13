
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//    Track: Displayed as set of Trackerhits			             //
//                                                                           //
//    Date:   03 Dec 2010                                                    //
//    Author: Manqi Ruan (LLR)                                               //
//                                                                           //
//    Last Modified: 1, Dec 2013, group into tracks			     //
//                                                                           // 
///////////////////////////////////////////////////////////////////////////////

#include "TRint.h"
#include "lcio.h"
#include "EVENT/LCCollection.h"
#include "EVENT/SimTrackerHit.h"
#include "EVENT/TrackerHit.h"
#include "EVENT/Track.h"
#include "TEveElement.h"
#include "TSystem.h"
#include "TEveManager.h"
#include "TEvePointSet.h"
#include "GlobalDefs.hh"

using namespace lcio;
using namespace std;

extern int GlobalRandomColorIndex;
extern bool Flag_AttachTextToHit;

TEveElementList* TrackAssignedHits( LCCollection* col, string name)
{

	string TT (name, 0, 15);

	cout<<"  Track collection: "<<name.c_str()<<". Number of Tracks: "<<col->getNumberOfElements()<<endl;
	cout<<endl;
	TEveElementList* TrackAsshitsColl = new TEveElementList;
	TrackAsshitsColl->SetName(name.c_str());

	TrackAsshitsColl->SetMainColor(5);

	float HitX, HitY, HitZ, dEdx, time;
	HitX = 0;
	HitY = 0;
	HitZ = 0;
	dEdx = 0;
	time = 0;
	int MCPID = 0;

	int nTrack = col->getNumberOfElements();
	int AssoHitNum = 0;	

	float Omega = 0; 
	float TanL = 0; 
	float TrkP = 0; 
	int colorindex; 

	for(int i(0); i<nTrack; i++)
	{
		Track* atrack = dynamic_cast<Track*>( col->getElementAt(i) );
		AssoHitNum = atrack->getTrackerHits().size();
		Omega = atrack->getOmega();
		TanL = atrack->getTanLambda();
		TrkP = 0.00105*sqrt(TanL*TanL + 1)/Omega;

		if(Omega > 0) colorindex = 2; 
		else colorindex = 3; 

		//if(fabs(Omega) < 1e-2 || fabs(TanL) > 1)	//Eq to High Q Tracks...
		// {

		TEveElementList * a_trk = new TEveElementList; 
		a_trk->SetName(Form( "TrkP = %.6f", TrkP ));

		for(int j(0); j<AssoHitNum; j++)
		{
			TrackerHit* hit = dynamic_cast<TrackerHit*>( atrack->getTrackerHits()[j] );
			HitX=hit->getPosition()[0];
			HitY=hit->getPosition()[1];
			HitZ=hit->getPosition()[2];
			dEdx=hit->getEDep();
			time=hit->getTime();

			TEvePointSet* q1 = new TEvePointSet(1);
			q1->SetName("Calorimeter Hit ");
			q1->SetMarkerStyle(3);
			q1->SetPoint(0, 0.1*HitX, 0.1*HitY, 0.1*HitZ);
			if(Flag_AttachTextToHit)
			{
				q1->SetTitle(Form("Associated Tracker Hit\n"
							"TrackerCollection = %s, TrackID = %d\n"
							"(HitX, HitY, HitZ) = (%.3f, %.3f, %.3f)\n"
							"Time = %f, dEdx = %E",
							TT.c_str(), i, HitX, HitY, HitZ,  time, dEdx));
			}
			// q1->SetMarkerColor(((i%2)*50+5*i+GlobalRandomColorIndex*31)%105);
			q1->SetMarkerColor(colorindex);

			q1->SetMarkerSize(0.1);

			a_trk->AddElement(q1);

			//TrackAsshitsColl->AddElement(q1);
		}

		 TrackAsshitsColl->AddElement(a_trk);

		// }

	}

	TrackAsshitsColl->SetRnrSelfChildren(kFALSE, kFALSE);

	return TrackAsshitsColl;
}



