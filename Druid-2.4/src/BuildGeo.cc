///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//    Build Geometry Volumes according to GDML file, with Projection         //
//                                                                           //
//    Date:   03 Dec 2010                                                    //
//    Author: Manqi Ruan (LLR)                                               //
//                                                                           //
//    Last Modified: 17, Nov 2011, Cleaning					                 //
//                                                                           // 
///////////////////////////////////////////////////////////////////////////////


#include "GlobalDefs.hh"
#include "TEveManager.h"
#include "TEveCompound.h"
#include "TEveGeoNode.h"
#include "TEvePointSet.h"
#include "TFile.h"
#include "TGeoManager.h"
#include "TGeoTube.h"
#include "TEveTrans.h"
#include "TGeoManager.h"
#include "TGeoNode.h"
#include "TGeoVolume.h"
#include "TGeoShape.h"
#include "TGeoMatrix.h"
#include "TEveGeoShape.h"
#include "MultiView.hh"

#include <string>
#include <cstdlib>
#include <fstream>
#include <sstream>

extern bool FlagMultiView;
extern MultiView * gMultiView;
extern TGeoManager * gGeoManager;
extern bool flagslcio; 
extern int flagdetectortype;	//after loadevt...

using namespace std;

TEveCompound* allmarkers(int markercolorindex)
{
	TEveCompound* mark = new TEveCompound;
	TEvePointSet* marker = new TEvePointSet(8);
	mark->SetName("Origin marker");
	mark->SetMainColor(2);
	marker->SetMarkerColor(markercolorindex);
	marker->SetMarkerStyle(3);
	marker->SetMarkerSize(2);

	Float_t a = 338.0;    //Unit in CM
	Float_t b = 385.0;
	marker->SetPoint(0, a,  +a, +b);
	marker->SetPoint(1, a,  -a, +b);
	marker->SetPoint(2, -a, -a, +b);
	marker->SetPoint(3, -a, +a, +b);
	marker->SetPoint(4, +a, +a, -b);
	marker->SetPoint(5, +a, -a, -b);
	marker->SetPoint(6, -a, +a, -b);
	marker->SetPoint(7, -a, -a, -b);
	mark->AddElement(marker);

	return mark;
}

void BuildGeoGDMLRoot(std::string gdmlroot){
	TFile* geom = TFile::Open(gdmlroot.c_str());
	if(!geom) return;

	gEve->RegisterGeometryAlias("whateverdetector",gdmlroot.c_str());
	gGeoManager = gEve->GetGeometryByAlias("whateverdetector");

	TGeoNode* node1 = gGeoManager->GetTopNode();

	//    TEveGeoTopNode* its = new TEveGeoTopNode(gGeoManager, node1);
	//    its->SetVisLevel(1);
	//    its->SetVisOption(1);

	TEveGeoTopNode* its = new TEveGeoTopNode(gGeoManager, node1, 1, 1, 10000);

	int k_nodes = node1->GetNdaughters();
	std::cout<<endl<<"   Number of Daughters (Geometry elements): "<<k_nodes<<std::endl<<std::endl;;    
	if(k_nodes>100) std::cout<<"   Too much first generation daughters! Will only print partial of their names."<<std::endl<<std::endl;
	for(int i=0; i<k_nodes; i++)
	{
		if( k_nodes<100 || i%20 == 0 )
			std::cout<<i<<"\t"<<node1->GetDaughter(i)->GetName()<<std::endl;
		TGeoNode * CurrentNode = node1->GetDaughter(i);

		CurrentNode->GetMedium()->GetMaterial()->SetTransparency(40);	//Default is 70
		//      CurrentNode->GetVolume()->SetLineColor(i%10);
	}

	//TEveCompound * markers = allmarkers( 2 );
	//gEve->AddGlobalElement(markers);
	gEve->AddGlobalElement(its);

	bool detectormatching = kFALSE; //Check, if geometry root file and detector type determined from lcio data file matchs
	if( (flagdetectortype == 0 && gdmlroot.find("ild_AHCAL") != string::npos) || ( flagdetectortype == 1 && gdmlroot.find("ild_DHCAL") !=string::npos ) ||((flagdetectortype == 2 || flagdetectortype == 4) && gdmlroot.find( "TB" ) != string::npos) || (flagdetectortype == 3 && gdmlroot.find( "loi3" ) != string::npos) || (flagdetectortype == 5 && gdmlroot.find("clic_sid_cdr_b") != string::npos))
		detectormatching = kTRUE;

	if( !detectormatching && flagslcio ) std::cout<<std::endl<<"   Detector type from gdml root file and slcio data file got conflict! Pls check! Will not make a geometry projection "<<std::endl<<std::endl;
	if( !flagslcio )
	{
		if(gdmlroot.find("ild_AHCAL") != string::npos) flagdetectortype = 0;
		else if(gdmlroot.find("ild_DHCAL") != string::npos) flagdetectortype = 1;
		else if(gdmlroot.find("TB") != string::npos) flagdetectortype = 2;				//Only Simu...?
		else if(gdmlroot.find("clic") != string::npos) flagdetectortype = 5;			//both, clic and sid with sid...
		else if(gdmlroot.find("sid") != string::npos) flagdetectortype = 3;			
	}


	if(FlagMultiView){
		//gMultiView->ImportGeomRPhi(markers);
		//gMultiView->ImportGeomRhoZ(markers);

		std::vector<std::string> VolProjection;
		std::vector<std::string> VolName;

		// if(ild)	
		if(detectormatching || !flagslcio)
		{
			if(flagdetectortype == 0)		//AHCAL: ILD00
			{
				VolProjection.push_back(std::string("/WorldLogical_1/TPCLog_11841"));
				VolProjection.push_back(std::string("/WorldLogical_1/CoilLogical_11894"));
				VolProjection.push_back(std::string("/WorldLogical_1/YokeBarrelLog_11895"));
				VolProjection.push_back(std::string("/WorldLogical_1/YokeEndcapLog_11896"));
				VolProjection.push_back(std::string("/WorldLogical_1/YokeEndcapLog_11897"));

				//HCAL Barrel Sides (RhoZ)
				VolProjection.push_back(std::string("/WorldLogical_1/barrelHcalModule_11763"));
				VolProjection.push_back(std::string("/WorldLogical_1/barrelHcalModule_11764"));
				//HCAL Barrel LOOP (RPhi)


				VolProjection.push_back(std::string("/WorldLogical_1/barrelHcalModule_11751"));
				VolProjection.push_back(std::string("/WorldLogical_1/barrelHcalModule_11753"));
				VolProjection.push_back(std::string("/WorldLogical_1/barrelHcalModule_11755"));
				VolProjection.push_back(std::string("/WorldLogical_1/barrelHcalModule_11757"));
				VolProjection.push_back(std::string("/WorldLogical_1/barrelHcalModule_11759"));
				VolProjection.push_back(std::string("/WorldLogical_1/barrelHcalModule_11761"));
				VolProjection.push_back(std::string("/WorldLogical_1/barrelHcalModule_11765"));


				//ECAL Barrel Sides (RhoZ)
				VolProjection.push_back(std::string("/WorldLogical_1/EnvLog_11717"));
				VolProjection.push_back(std::string("/WorldLogical_1/EnvLog_11718"));
				VolProjection.push_back(std::string("/WorldLogical_1/EnvLog_11719"));
				VolProjection.push_back(std::string("/WorldLogical_1/EnvLog_11720"));
				VolProjection.push_back(std::string("/WorldLogical_1/EnvLog_11721"));
				//ECAL Barrel LOOP (RPhi)
				VolProjection.push_back(std::string("/WorldLogical_1/EnvLog_11707"));
				VolProjection.push_back(std::string("/WorldLogical_1/EnvLog_11712"));
				VolProjection.push_back(std::string("/WorldLogical_1/EnvLog_11722"));
				VolProjection.push_back(std::string("/WorldLogical_1/EnvLog_11727"));
				VolProjection.push_back(std::string("/WorldLogical_1/EnvLog_11732"));
				VolProjection.push_back(std::string("/WorldLogical_1/EnvLog_11737"));
				VolProjection.push_back(std::string("/WorldLogical_1/EnvLog_11742"));

				//ECAL EndCap
				VolProjection.push_back(std::string("/WorldLogical_1/EndCapLog_11747"));
				VolProjection.push_back(std::string("/WorldLogical_1/EndCapLog_11748"));

				VolName.push_back(std::string("TPC"));
				VolName.push_back(std::string("Coil"));
				VolName.push_back(std::string("YokeBarrel"));
				VolName.push_back(std::string("YokeEndCap1"));
				VolName.push_back(std::string("YokeEndCap2"));

				VolName.push_back(std::string("HcalBarrel1"));
				VolName.push_back(std::string("HcalBarrel2"));


				VolName.push_back(std::string("HcalBarrel3"));
				VolName.push_back(std::string("HcalBarrel4"));
				VolName.push_back(std::string("HcalBarrel5"));
				VolName.push_back(std::string("HcalBarrel6"));
				VolName.push_back(std::string("HcalBarrel7"));
				VolName.push_back(std::string("HcalBarrel8"));
				VolName.push_back(std::string("HcalBarrel9"));


				VolName.push_back(std::string("EcalBarrel1"));
				VolName.push_back(std::string("EcalBarrel2"));
				VolName.push_back(std::string("EcalBarrel3"));
				VolName.push_back(std::string("EcalBarrel4"));
				VolName.push_back(std::string("EcalBarrel5"));

				VolName.push_back(std::string("EcalBarrel6"));
				VolName.push_back(std::string("EcalBarrel7"));
				VolName.push_back(std::string("EcalBarrel8"));
				VolName.push_back(std::string("EcalBarrel9"));
				VolName.push_back(std::string("EcalBarrel10"));
				VolName.push_back(std::string("EcalBarrel11"));
				VolName.push_back(std::string("EcalBarrel12"));

				VolName.push_back(std::string("EcalEndCap1"));
				VolName.push_back(std::string("EcalEndCap2"));
			}
			else if(flagdetectortype == 1)	//DHCAL: ILD00_dhcal
			{
				VolProjection.push_back(std::string("/WorldLogical_1/TPCLog_14062"));
				VolProjection.push_back(std::string("/WorldLogical_1/CoilLogical_14115"));
				VolProjection.push_back(std::string("/WorldLogical_1/YokeBarrelLog_14116"));
				VolProjection.push_back(std::string("/WorldLogical_1/YokeEndcapLog_14117"));
				VolProjection.push_back(std::string("/WorldLogical_1/YokeEndcapLog_14118"));
				//			VolProjection.push_back(std::string("/WorldLogical_1/HcalBarrelEnvLog_13993"));	//Memory malloc problem from projection of composed TEve objects
				//			VolProjection.push_back(std::string("/WorldLogical_1/EndCapLogical_13994"));
				//			VolProjection.push_back(std::string("/WorldLogical_1/EndCapLogical_13995"));

				VolName.push_back(std::string("TPC"));
				VolName.push_back(std::string("Coil"));
				VolName.push_back(std::string("YokeBarrel"));
				VolName.push_back(std::string("YokeEndCap1"));
				VolName.push_back(std::string("YokeEndCap2"));
				//			VolName.push_back(std::string("HCALBarrel"));
				//			VolName.push_back(std::string("HCALEndCap1"));
				//			VolName.push_back(std::string("HCALEndCap2"));
			}
			else if(flagdetectortype == 3){	//SiD
				// Tubes
				VolProjection.push_back(std::string("/world_volume_1/tracking_volume_12890"));
				VolProjection.push_back(std::string("/world_volume_1/SolenoidCoilBarrel_layer0_volume_12885"));
				VolProjection.push_back(std::string("/world_volume_1/SolenoidCoilBarrel_layer1_volume_12884"));
				VolProjection.push_back(std::string("/world_volume_1/SolenoidCoilBarrel_layer2_volume_12883"));
				//Polygons
				VolProjection.push_back(std::string("/world_volume_1/EcalBarrel_envelope_12866"));
				VolProjection.push_back(std::string("/world_volume_1/EcalEndcap_volume_12888"));
				VolProjection.push_back(std::string("/world_volume_1/EcalEndcap_volume_12889"));

				VolProjection.push_back(std::string("/world_volume_1/HcalBarrel_envelope_12887"));
				VolProjection.push_back(std::string("/world_volume_1/HcalEndcap_volume_12879"));
				VolProjection.push_back(std::string("/world_volume_1/HcalEndcap_volume_12880"));

				VolProjection.push_back(std::string("/world_volume_1/MuonBarrel_envelope_12886"));	//Muon Detector is not a tube...
				VolProjection.push_back(std::string("/world_volume_1/MuonEndcap_volume_12869"));
				VolProjection.push_back(std::string("/world_volume_1/MuonEndcap_volume_12870"));


				VolName.push_back(std::string("Tracker"));
				VolName.push_back(std::string("Coil0"));
				VolName.push_back(std::string("Coil1"));
				VolName.push_back(std::string("Coil2"));
				VolName.push_back(std::string("ECALBarrel"));
				VolName.push_back(std::string("ECALEndCap1"));
				VolName.push_back(std::string("ECALEndCap2"));
				VolName.push_back(std::string("HCALBarrel"));
				VolName.push_back(std::string("HCALEndCap1"));
				VolName.push_back(std::string("HCALEndCap2"));
				VolName.push_back(std::string("MuonBarrel"));
				VolName.push_back(std::string("MuonEndCap1"));
				VolName.push_back(std::string("MuonEndCap2"));
			}else if(flagdetectortype == 5){
				VolProjection.push_back(std::string("/world_volume_1/tracking_volume_13152"));
				VolProjection.push_back(std::string("/world_volume_1/SolenoidCoilBarrel_layer0_volume_13144"));
				VolProjection.push_back(std::string("/world_volume_1/SolenoidCoilBarrel_layer1_volume_13143"));
				VolProjection.push_back(std::string("/world_volume_1/SolenoidCoilBarrel_layer2_volume_13142"));
				VolProjection.push_back(std::string("/world_volume_1/EcalBarrel_envelope_13116"));
				VolProjection.push_back(std::string("/world_volume_1/EcalEndcap_volume_13150"));
				VolProjection.push_back(std::string("/world_volume_1/EcalEndcap_volume_13151"));

				VolProjection.push_back(std::string("/world_volume_1/HcalBarrel_envelope_13149"));
				VolProjection.push_back(std::string("/world_volume_1/HcalEndcap_volume_13138"));
				VolProjection.push_back(std::string("/world_volume_1/HcalEndcap_volume_13139"));

				VolProjection.push_back(std::string("/world_volume_1/MuonBarrel_envelope_13145"));
				VolProjection.push_back(std::string("/world_volume_1/MuonEndcap_volume_13121"));
				VolProjection.push_back(std::string("/world_volume_1/MuonEndcap_volume_13122"));

				VolName.push_back(std::string("Tracker"));
				VolName.push_back(std::string("Coil0"));
				VolName.push_back(std::string("Coil1"));
				VolName.push_back(std::string("Coil2"));
				VolName.push_back(std::string("ECALBarrel"));
				VolName.push_back(std::string("ECALEndCap1"));
				VolName.push_back(std::string("ECALEndCap2"));

				VolName.push_back(std::string("HCALBarrel"));
				VolName.push_back(std::string("HCALEndCap1"));
				VolName.push_back(std::string("HCALEndCap2"));
				VolName.push_back(std::string("MuonBarrel"));
				VolName.push_back(std::string("MuonEndCap1"));
				VolName.push_back(std::string("MuonEndCap2"));

			}else{
				std::cout<<std::endl<<"   Projection of Current detector is not supported. Contact the author at ruan@llr.in2p3.fr if you really need... "<<std::endl<<std::endl;
			}

			int numVolPro = VolProjection.size();

			if(VolProjection.size() != VolName.size() ) 
			{
				std::cout<<std::endl<<"   Vol of Projection and Name doesn't match! "<<std::endl<<std::endl;
			}else if(VolName.size()>0){

				std::cout<<std::endl<<"   Start Geometry Volume Projection! "<<std::endl<<std::endl;

				for(int i=0; i<numVolPro; ++i)
				{

					std::cout<<"   Volume Number "<< i<< ", "<<VolName[i].c_str()<<std::endl;

					gGeoManager -> cd(VolProjection[i].c_str());
					string ShapeName = VolName[i].c_str();
					TEveGeoShape * pmd = new TEveGeoShape(VolName[i].c_str(), VolName[i].c_str());

					pmd->SetNSegments(100);

					pmd->RefMainTrans().SetFrom(*gGeoManager->GetCurrentMatrix());
					pmd->SetShape((TGeoShape*) gGeoManager->GetCurrentVolume()->GetShape()->Clone());

					pmd->SetMainColor(i%40 + 11);
					pmd->SetLineColor(1);
					pmd->SetLineWidth(1);
					pmd->IncDenyDestroy();

					if(gGeoManager -> cd(VolProjection[i].c_str())){
						if( ShapeName.find("EndCap") == string::npos ) gMultiView->ImportGeomRPhi(pmd);
						gMultiView->ImportGeomRhoZ(pmd);
					}
				}

				std::cout<<std::endl<<"   Projection Finished. "<<std::endl;
			}

		}	
		//	std::cout<<" GEOMETRY MULTIVIEW LOADED FINISHED"<<std::endl;
	}
}

