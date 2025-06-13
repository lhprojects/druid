
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//    (Sim) Tracker Hits: Displayed as points.                               //
//                                                                           //
//    Date:   03 Dec 2010                                                    //
//    Author: Manqi Ruan (LLR)                                               //
//                                                                           //
//    Last Modified: 17, Nov 2011, cleaning                                  //
//                                                                           // 
///////////////////////////////////////////////////////////////////////////////

#include "lcio.h"
#include "IMPL/LCTOOLS.h"
#include "EVENT/LCCollection.h"
#include "EVENT/SimTrackerHit.h"
#include "EVENT/TrackerHit.h"
#include "TEveElement.h"
#include "TEvePointSet.h"

using namespace lcio;
using namespace std;

extern bool Flag_AttachTextToHit;

TEveElementList* TrackerHits( LCCollection* col, string name)
{
	bool isSimHit = col->getTypeName() == LCIO::SIMTRACKERHIT;	

	cout<<"  Tracker Hits collection: "<<name.c_str()<<". Number of Hits: "<<col->getNumberOfElements()<<endl;
	cout<<endl;
	TEveElementList* TrackerhitsColl = new TEveElementList;
	TrackerhitsColl->SetName(name.c_str());

	if(isSimHit)
	{
		TrackerhitsColl->SetMainColor(4);
	} 
	else
	{
		TrackerhitsColl->SetMainColor(2);
	}

	float HitX, HitY, HitZ, dEdx, time;
	HitX = 0;
	HitY = 0;
	HitZ = 0;
	dEdx = 0;
	time = 0;
	int MCPID = 0;
	float MCPartEnergy = 0;

	string TT (name, 0, 15);

	if (isSimHit){			//SimTrackerHit
		if(col->getTypeName() == LCIO::SIMTRACKERHIT)
		{

			int nHits = col->getNumberOfElements();
			for(int i=0; i<nHits; i++)
			{
				SimTrackerHit* hit = dynamic_cast<SimTrackerHit*>( col->getElementAt(i) );
				HitX=hit->getPosition()[0];
				HitY=hit->getPosition()[1];
				HitZ=hit->getPosition()[2];
				dEdx=hit->getEDep();
				time=hit->getTime();

				MCParticle* hitMCPart = dynamic_cast<MCParticle*>( hit->getMCParticle());
				if(hitMCPart){
					MCPID=hitMCPart->getPDG();
					MCPartEnergy = hitMCPart->getEnergy();
				}else {
					MCPID=-99;
					MCPartEnergy = -99999;
				}

				TEvePointSet* q1 = new TEvePointSet(1);
				q1->SetName("Calorimeter Hit ");
				q1->SetMarkerStyle(3);
				q1->SetPoint(0, 0.1*HitX, 0.1*HitY, 0.1*HitZ);
				if(Flag_AttachTextToHit)
				{
					q1->SetTitle(Form("Simulated Tracker Hit\n"
								"SubDetector=%s\n"
								"(HitX, HitY, HitZ) = (%.3f, %.3f, %.3f)\n"
								"MCPID=%d, MCPEn=%E, Time = %f, dEdx = %E",
								TT.c_str(), HitX, HitY, HitZ, MCPID, MCPartEnergy, time, dEdx));
				}
				if(hitMCPart){
					if(hitMCPart->getCharge()==1)
					{q1->SetMarkerColor(3);}
					else if(hitMCPart->getCharge()==-1)
					{q1->SetMarkerColor(2);}   //Show different color to the MCParticle Track
					else
					{q1->SetMarkerColor(3);}   //Just in case: Unexpected Hits
				}else{
					q1->SetMarkerColor(4);
				}


				q1->SetMarkerSize(0.1);
				TrackerhitsColl->AddElement(q1);
			}
		}
	}
	else if(!isSimHit){
		if(col->getTypeName() == LCIO::TRACKERHIT)
		{

			int nHits = col->getNumberOfElements();
			for(int i=0; i<nHits; i++)
			{
				TrackerHit* hit = dynamic_cast<TrackerHit*>( col->getElementAt(i) );
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
					q1->SetTitle(Form("Reconstructed Tracker Hit\n"
								"SubDetector=%s\n"
								"(HitX, HitY, HitZ) = (%.3f, %.3f, %.3f)\n"
								"Time = %f, dEdx = %E",
								TT.c_str(), HitX, HitY, HitZ,  time, dEdx));
				}
				q1->SetMarkerColor(3);   

				q1->SetMarkerSize(0.1);
				TrackerhitsColl->AddElement(q1);
			}
		}
	}

	TrackerhitsColl->SetRnrSelfChildren(kFALSE, kFALSE);

	return TrackerhitsColl;
}



