
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//    Vertex: Displayed as Points			                                 //
//                                                                           //
//    Date:   03 Dec 2010                                                    //
//    Author: Manqi Ruan (LLR)                                               //
//                                                                           //
//    Last Modified: 17, Nov 2011, cleaning                                  //
//                                                                           // 
///////////////////////////////////////////////////////////////////////////////

#include "TRint.h"
#include "lcio.h"
#include "IO/LCReader.h"
#include "IMPL/LCTOOLS.h"
#include "EVENT/LCCollection.h"
#include "EVENT/Vertex.h"
#include "TEveElement.h"
#include "TSystem.h"
#include "TEveManager.h"
#include "TEvePointSet.h"

using namespace lcio;
using namespace std;

extern int GlobalRandomColorIndex;
extern int event_id;

TEveElementList* Vertex( LCCollection* col, string name)
{

	string TT (name, 0, 15);

	cout<<"  Vertex collection: "<<name.c_str()<<". Number of Vertex: "<<col->getNumberOfElements()<<endl;
	cout<<endl;
	TEveElementList* VertexPoint = new TEveElementList;
	VertexPoint->SetName(name.c_str());

	VertexPoint->SetMainColor(7);

	float HitX, HitY, HitZ, chi2;
	HitX = 0;
	HitY = 0;
	HitZ = 0;
	chi2 = 0;

	int nVtx = col->getNumberOfElements();

	for(int i(0); i<nVtx; i++)
	{
		EVENT::Vertex* avtx = dynamic_cast<EVENT::Vertex*>( col->getElementAt(i) );

		HitX = avtx->getPosition()[0]*0.1;
		HitY = avtx->getPosition()[1]*0.1;
		HitZ = avtx->getPosition()[2]*0.1;
		chi2 = avtx->getChi2();

		TEvePointSet* qv = new TEvePointSet(1);
		qv->SetName("Vertex ");
		qv->SetMarkerStyle(3);
		qv->SetPoint(0, 0.1*HitX, 0.1*HitY, 0.1*HitZ);
		qv->SetTitle(Form("Vertex of Event %d\n"
					"Vertex Collection = %s, Vertex ID = %d\n"
					"(X, Y, Z) = (%.3f, %.3f, %.3f)\n"
					"Chi2 = %f",
					event_id, TT.c_str(), i, HitX, HitY, HitZ, chi2));

		qv->SetMarkerColor(((i%2)*50+15*i+GlobalRandomColorIndex*31)%105);

		qv->SetMarkerSize(1.0);
		VertexPoint->AddElement(qv);
	}


	return VertexPoint;
}



