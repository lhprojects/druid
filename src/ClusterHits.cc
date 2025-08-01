
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
#include "point.h"
#include "segment3.h"

using namespace lcio;
using namespace std;
using namespace EVENT;

float ClusterHitSize = 1.0;
int ClusterHitColourType = 1; 
extern int flagdetectortype;
const float HCALBarrelLengthILD = 309.3;
const float HCALBarrelLengthSID = 320.0; 
extern int event_id;
extern TEveBox* BoxPhi( TVector3 &HitPos, TVector3 &Scale, int Type, int SegOrStaveNumber, float HitEnergy );
extern int GlobalRandomColorIndex;
extern bool Flag_AttachTextToHit;
TRandom *r0 = new TRandom();

TEveElementList* ClusterHits( LCCollection* col, string name)
{

	string TT (name, 0, 15);

	string NameHeader(name, 0, 5);	//Used to identify if Pandora Cluster Or ... Arbor Cluster...

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

		if((name == "HcalSortBranches" || name == "EcalSortBranches") && CluSize > 5 )
		{	
			double cone_crit = 30.; // Criteria for cone segmentation (in degrees)
			double chi2_crit = 50.;

			fit_type ft = HEL; // HEL (helix) sometimes is unstable for few points. Better use CIRC or SVD
			//point axis = point(0,0,1);
			point axis = point( cluPos.X(), cluPos.Y(), cluPos.Z() ); 

			trajectory t(ft, axis, cone_crit, chi2_crit);
			vector <point> points;
			for(unsigned j0 = 0; j0 < CluSize; j0++) // Now we add hits to t one by one
			{
				CalorimeterHit *a_hit = acluster->getCalorimeterHits()[j0];
				points.push_back( point( (double) a_hit->getPosition()[0],
							(double) a_hit->getPosition()[1],
							(double) a_hit->getPosition()[2]) );

				if(j0 == 0)
					cluEndPos = a_hit->getPosition();
				if(j0 == CluSize - 1)
					cluBeginPos = a_hit->getPosition();
			}

			unsigned inv;

			for(unsigned j1 = 0; j1 < CluSize; j1++)
			{
				inv = points.size() - j1 - 1;
				t.add_point( points[ inv ] );
			}

			point RefD = t.dir_begin();
			point EndD = t.dir_end(); 

			cluDir_B.SetXYZ(RefD.getX(), RefD.getY(), RefD.getZ());
			cluDir_E.SetXYZ(EndD.getX(), EndD.getY(), EndD.getZ());

			cluDir_B = 1.0/cluDir_B.Mag() * cluDir_B; 
			cluDir_E = 1.0/cluDir_E.Mag() * cluDir_E; 

			TEveArrow* a0 = new TEveArrow(10*cluDir_B.X(), 10*cluDir_B.Y(), 10*cluDir_B.Z(), 0.1*cluBeginPos.X(), 0.1*cluBeginPos.Y(), 0.1*cluBeginPos.Z());
			a0->SetTubeR(0.01);
			a0->SetConeR(0.03);
			a0->SetConeL(0.2);
			a0->SetMainColor(5);
			CaloCluster->AddElement(a0);

//			TEveArrow* a2 = new TEveArrow(10*EndD.getX(), 10*EndD.getY(), 10*EndD.getZ(), 0.1*cluEndPos.X(), 0.1*cluEndPos.Y(), 0.1*cluEndPos.Z());

			TEveArrow* a2 = new TEveArrow(10*cluDir_E.X(), 10*cluDir_E.Y(), 10*cluDir_E.Z(), 0.1*cluEndPos.X(), 0.1*cluEndPos.Y(), 0.1*cluEndPos.Z());
			a2->SetTubeR(0.01);
			a2->SetConeR(0.03);
			a2->SetConeL(0.2);
			a2->SetMainColor(5);
			CaloCluster->AddElement(a2);


		}

		if(CluSize > 20)
		{
			TEveArrow* a1 = new TEveArrow(cluDir.X(), cluDir.Y(), cluDir.Z(), 0.1*cluPos.X(), 0.1*cluPos.Y(), 0.1*cluPos.Z());
			a1->SetTubeR(0.01);
			a1->SetConeR(0.03);
			a1->SetConeL(0.2);
			a1->SetMainColor(3);
			CaloCluster->AddElement(a1);
		}

		LocalRandomIndex = int(100*r0->Rndm(iC));

		if(( (name == "HcalSortBranches" || name == "EcalSortBranches") && CluSize > 5 ) || (name != "HcalSortBranches" && name != "EcalSortBranches"))
		{
			for(int j = 0; j<CluSize; j++)
			{
				TVector3 HitPosition = Hits[j]->getPosition();
				HitEn = Hits[j]->getEnergy();
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

				if(ClusterHitColourType == 0)
				{
					q->SetMainColor(5);
					q->SetLineColor(5);
				}
				else if(ClusterHitColourType == 1)
				{
					int colorindex = (iC*11 + GlobalRandomColorIndex*13 + LocalRandomIndex)%50 + 51;
					q->SetMainColor(colorindex);
					q->SetLineColor(colorindex);
				}
				else if(ClusterHitColourType == 2)
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
								name.c_str(), event_id, HitEn*1000000, 10*HitPosition[0], 10*HitPosition[1], 10*HitPosition[2],
								0, 0., 0., ClusterEnergy, cluPos[0], cluPos[1], cluPos[2]));
				}
				q->SetPickable(kTRUE);
				Recocluster->AddElement(q);
			}
		}

		CaloCluster->AddElement(Recocluster);
	}	

	// GlobalRandomColorIndex++;

	return CaloCluster;
}



