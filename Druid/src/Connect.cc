#include "TRint.h"
#include "lcio.h"
#include "IO/LCReader.h"
#include "IMPL/LCTOOLS.h"
#include "EVENT/LCCollection.h"
#include "EVENT/SimCalorimeterHit.h"
#include "EVENT/CalorimeterHit.h"
#include "EVENT/LCRelation.h"
#include "TEveElement.h"
#include "TEveStraightLineSet.h"
#include "TEveCompound.h"
#include "TSystem.h"
#include "TEveManager.h"
#include "TEveArrow.h"

using namespace lcio;
using namespace std;

extern int HenriCount;

//TEveStraightLineSet* ConnectTrees( LCCollection* col, string name)
TEveCompound* ConnectTrees( LCCollection* col, string name)
{

//	TEveStraightLineSet* ls = new TEveStraightLineSet();
	TEveCompound *ls = new TEveCompound();
	ls -> SetName(name.c_str());

	float X1, Y1, Z1;
	float X2, Y2, Z2;	//Coordinates...

		cout<<"NumofLinks "<<col->getNumberOfElements()<<endl;

		for(int k = 0; k<col->getNumberOfElements(); ++k)
		{
			LCRelation *rel = dynamic_cast<LCRelation*>( col->getElementAt(k) );
		
			CalorimeterHit *hit1 = dynamic_cast<CalorimeterHit*>(rel->getTo());
			CalorimeterHit *hit2 = dynamic_cast<CalorimeterHit*>(rel->getFrom());
		
		 /*
			SimCalorimeterHit *hit1 = dynamic_cast<SimCalorimeterHit*>(rel->getTo());
			SimCalorimeterHit *hit2 = dynamic_cast<SimCalorimeterHit*>(rel->getFrom());
		 */
			X1 = hit1->getPosition()[0];
			Y1 = hit1->getPosition()[1];
			Z1 = hit1->getPosition()[2];
			X2 = hit2->getPosition()[0];
			Y2 = hit2->getPosition()[1];
			Z2 = hit2->getPosition()[2];

			TEveArrow* a1 = new TEveArrow(0.1*(X1-X2), 0.1*(Y1-Y2), 0.1*(Z1-Z2), 0.1*X2, 0.1*Y2, 0.1*Z2);
			a1->SetTubeR(0.01);
			a1->SetConeR(0.03);
			a1->SetConeL(0.2);
			a1->SetMainColor(HenriCount+1);
			ls->AddElement(a1);
			//			ls->SetLineColor(3);
			//			ls->AddLine(0.1*X1, 0.1*Y1, 0.1*Z1, 0.1*X2, 0.1*Y2, 0.1*Z2);
		}

		return ls; 

}


