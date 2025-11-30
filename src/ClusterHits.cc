
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//    Build Clustered CaloHits into TEveBox,								 //
//    Color Specified to Index.												 //
//    Option to be add as trace back to Reco PID if needed		             //
//                                                                           //
//    Date:   03 Dec 2010                                                    //
//    Author: Manqi Ruan (LLR)                                               //
//                                                                           //
//    Last Modified, 17, Nov 2011, Cleanning                                 //
//                                                                           // 
///////////////////////////////////////////////////////////////////////////////

#include "TRint.h"
#include "TEveBox.h"
#include "lcio.h"
#include "IMPL/LCTOOLS.h"
#include "EVENT/LCCollection.h"
#include "EVENT/Cluster.h"
#include "TEveElement.h"
#include "TSystem.h"
#include "TEvePointSet.h"
#include "TVector3.h"
#include "TEveArrow.h"
#include "TRandom3.h"

#include "trajectory.h"
#include "fitting_root.h"
#include "geometry.h"
#include "Options.h"
#include "GlobalDefs.hh"
#include "point.h"
#include "segment3.h"
#include "GlobalDefs.hh"
#include "TruthHelper.h"

using namespace lcio;
using namespace std;
using namespace EVENT;

float ClusterHitSize = 1.0;
extern int flagdetectortype;
const float HCALBarrelLengthILD = 309.3;
const float HCALBarrelLengthSID = 320.0; 
extern TEveBox* BoxPhi( TVector3 &HitPos, TVector3 &Scale, int Type, int SegOrStaveNumber, float HitEnergy );
extern int GlobalRandomColorIndex;
extern bool Flag_AttachTextToHit;
TRandom *r0 = new TRandom();

TEveElementList* ClusterHits( LCCollection* col, string name)
{
	// Filter cluster collections based on command line options
	if(gOptions.coll_cluster_collections.size() > 0) {
		if(!equals_any(name, gOptions.coll_cluster_collections)) {
			cout << "  Skipping cluster collection: " << name << " (not in -coll.cluster.add list)" << endl;
			return nullptr;
		}
	}

	cout<<"  Cluster collection: "<<name.c_str()<<". Number of Cluster: "<<col->getNumberOfElements()<<endl;
	cout<<endl;
	TEveElementList* CaloCluster = new TEveElementList;
	CaloCluster->SetName(name.c_str());
	CaloCluster->SetMainColor(5);

	float HitEn, ClusterEnergy ;
	HitEn = 0;
	ClusterEnergy = 0;
	int nCluster = col->getNumberOfElements();
	int ClusterPID = 0;
	float HBLength = 0;
	int LocalRandomIndex = 0;
	int CluSize = 0;
	TVector3 cluPos, cluBeginPos, cluEndPos, cluDir, cluDir_B, cluDir_E; 

	if(flagdetectortype == 3) 
	{
		HBLength = HCALBarrelLengthSID;
	}
	else if(flagdetectortype < 2) 
	{
		HBLength = HCALBarrelLengthILD;  
	}

	std::cout<<"HCAL Barrel Length "<<HBLength<<std::endl;

	TVector3 HitScale(0.5*ClusterHitSize, 0.5*ClusterHitSize, 0.5*ClusterHitSize);

	for(int iC = 0; iC<nCluster; iC++)
	{
		TEveElementList* Recocluster = new TEveElementList;
		Recocluster->SetMainColor(iC);

		Cluster* acluster  = dynamic_cast<Cluster*>( col->getElementAt(iC) );
		cluPos = acluster->getPosition();
		cluDir = 10.0/cluPos.Mag()*cluPos; 
		CalorimeterHitVec Hits = acluster->getCalorimeterHits();
		CluSize = Hits.size();
		ClusterEnergy = acluster->getEnergy();
		ClusterPID = acluster->getType();
		Recocluster->SetName(Form("Clu En=%.3f/", ClusterEnergy));
		// std::cout<<"CLUSTERPIDTYPE "<<ClusterPID<<std::endl;

		LocalRandomIndex = int(100*r0->Rndm(iC));

		if(true)
		{
			for(int j = 0; j<CluSize; j++)
			{
			TVector3 HitPosition = Hits[j]->getPosition();
			HitEn = Hits[j]->getEnergy();
			HitPosition *= 0.1;

			// Create a cube with size based on ClusterHitSize (in cm)
			float cubeHalfSize = 0.5 * ClusterHitSize;
			TEveBox* q = new TEveBox();
			q->SetVertex(0, HitPosition.X() - cubeHalfSize, HitPosition.Y() - cubeHalfSize, HitPosition.Z() - cubeHalfSize);
			q->SetVertex(1, HitPosition.X() + cubeHalfSize, HitPosition.Y() - cubeHalfSize, HitPosition.Z() - cubeHalfSize);
			q->SetVertex(2, HitPosition.X() + cubeHalfSize, HitPosition.Y() + cubeHalfSize, HitPosition.Z() - cubeHalfSize);
			q->SetVertex(3, HitPosition.X() - cubeHalfSize, HitPosition.Y() + cubeHalfSize, HitPosition.Z() - cubeHalfSize);
			q->SetVertex(4, HitPosition.X() - cubeHalfSize, HitPosition.Y() - cubeHalfSize, HitPosition.Z() + cubeHalfSize);
			q->SetVertex(5, HitPosition.X() + cubeHalfSize, HitPosition.Y() - cubeHalfSize, HitPosition.Z() + cubeHalfSize);
			q->SetVertex(6, HitPosition.X() + cubeHalfSize, HitPosition.Y() + cubeHalfSize, HitPosition.Z() + cubeHalfSize);
			q->SetVertex(7, HitPosition.X() - cubeHalfSize, HitPosition.Y() + cubeHalfSize, HitPosition.Z() + cubeHalfSize);

			if(gGUIManager.ClusterHitColourType == 0)
				{
					q->SetMainColor(5);
					q->SetLineColor(5);
				}
				else if(gGUIManager.ClusterHitColourType == 1)
				{
					int colorindex = (iC*11 + GlobalRandomColorIndex*13 + LocalRandomIndex)%50 + 51;
					q->SetMainColor(colorindex);
					q->SetLineColor(colorindex);
				}
				else if(gGUIManager.ClusterHitColourType == 2)
				{        //Uniform Color: Red for PFO used
					q->SetMainColor(2);
					q->SetLineColor(2);
				}
				if(Flag_AttachTextToHit)
				{
					q->SetTitle(Form( "CLUSTER (%s) Calo Hit, EventNr = %d\n"
								"Hit Energy=%.3f keV\n"
								"PosX = %.3f mm, PosY = %.3f mm, PosZ = %.3f mm\n"
								"PFOPDG = %d, PFOCharge = %f, PFOEnergy = %f\n"
								"Cluster Energy = %f GeV\n"
								"Cluster PosX = %.3f mm, PosY = %.3f mm, PosZ = %.3f mm\n", 
								name.c_str(), gDisplayState.getEventNumber(), HitEn*1000000, 10*HitPosition[0], 10*HitPosition[1], 10*HitPosition[2],
								0, 0., 0., ClusterEnergy, cluPos[0], cluPos[1], cluPos[2]));
				}
				q->SetPickable(kTRUE);
				Recocluster->AddElement(q);
			}
		}

		CaloCluster->AddElement(Recocluster);
		
		// Draw mother-daughter connection arrow
		LCObjectConnection conn = gTruthHelper.GetTracsterConnection(acluster);
		TEveArrow* connArrow = createConnectionArrow(conn);
		if (connArrow) {
			Recocluster->AddElement(connArrow);
		}
	}	

	// GlobalRandomColorIndex++;

	return CaloCluster;
}



