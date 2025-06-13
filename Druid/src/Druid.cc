
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//    Welcome to Druid.														 //
//    Main Program.															 //
//    Reading Parameters from steering file or Arguments				     //
//                                                                           //
//    Date:   17 Nov 2011                                                    //
//    Author: Manqi Ruan (LLR)                                               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "GlobalDefs.hh"
#include "TSystem.h"
#include "TEveManager.h"
#include "TTree.h"
#include "TFile.h"
#include "TGLViewer.h"
#include "TObjString.h"
#include "TROOT.h"
#include "TRootBrowser.h"
#include "TStyle.h"
#include "TRint.h"
#include <cstdlib>
#include <string>
#include <fstream>
#include <sstream>

using namespace lcio ;
using namespace std ;

#include "TGeoManager.h"
extern TGeoManager * gGeoManager;

int runNumber = 0;
int event_id = 0;

#include "lcio.h"
#include "IO/LCReader.h"
using namespace lcio ;

LCReader* lcReader;
string lcioFileName;
string gearFileName;
string detectortype;
string gdmlrootfileName;
int flagdetectortype = -1;
bool flagslcio = kFALSE;
bool DefaultCollectionFlag = kTRUE;
//bool FlagMultiView = kTRUE;
bool FlagMultiView = kFALSE;

#include "MultiView.hh"
MultiView* gMultiView = 0;

void SteerSetup(char *steeringFileName)
{

	ifstream steeringFile(steeringFileName);
	string parName, parValue;

	if(!steeringFile)
	{
		cout<<"Steering File Not Found!"<<endl;
		exit(2);
	}

	std::string line;

	while(getline(steeringFile, line)) {
		istringstream in(line);
		in >> parName >> parValue;
		if(parName.find('#') != string::npos || parValue=="") 
			continue;
		if(parName == "GEARFILE")
			gearFileName = parValue;
		if(parName == "LCIOFILE")
			lcioFileName = parValue;
		if(parName == "RUNNUMBER")
			runNumber = atoi (parValue.c_str());
		if(parName == "EVTNUMBER")
			event_id = atoi (parValue.c_str());
		if(parName == "DETECTOR")
			detectortype = parValue;
		if(parName == "GDMLROOTFILE")
			gdmlrootfileName = parValue;
	}
}

TString ExtendName(TString inputString){
	TObjArray *bb = inputString.Tokenize(".");
	int NameSize = bb->GetEntries();
	TString ExtName = ((TObjString*) bb->At(NameSize-1))->GetString();
	return ExtName;
}

int main(int argc, char *argv[])
{

	if(argc == 1) 
	{
		std::cout<<std::endl<<"Specify your steering file or input parameters with either of following format: "<<std::endl<<std::endl;
		std::cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl<<std::endl;
		std::cout<<"(1), bin/Druid *.steer    Parameters could be specified in steering file "<<std::endl<<std::endl;
		std::cout<<"(2), bin/Druid *.slcio    Browser the slcio file, with display the first event as default"<<std::endl<<std::endl;
		std::cout<<"(3), bin/Druid *.slcio n   Number n used to specify the first event to be displayed"<<std::endl<<std::endl;
		std::cout<<"(4), bin/Druid *gdml.root  or  Druid *geom.xml    Browser detector geometry according to xml file or root file converted from gdml file"<<std::endl<<std::endl;
		std::cout<<"(5), bin/Druid *.slcio *gdml.root/*geom.xml    Display the first event in lcio file and the detector geometry "<<std::endl<<std::endl;
		std::cout<<"(6), bin/Druid *.slcio *gdml.root/*geom.xml n     Display the event n in lcio file and the detector geometry "<<std::endl<<std::endl;
		std::cout<<"(7), bin/Druid *.slcio *gdml.root/*geom.xml n1 n2   Display event with RunNumber = n1 and EventNum = n2 in lcio file and the detector geometry "<<std::endl<<std::endl;
		std::cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<std::endl;
		exit(1);
	}

	TRint *theApp = 0;
	theApp = new TRint("ROOT example", 0, 0);

	TEveManager::Create();
	gGeoManager->SetVerboseLevel(0);
	TGLViewer* Ori = gEve->GetDefaultGLViewer();
	Ori->SetIgnoreSizesOnUpdate(kTRUE);

	make_gui();

	flagslcio = kFALSE;
	bool flaggeoroot = kFALSE;
	bool flaggeoxml = kFALSE;
	bool flagEventNumber = kFALSE;
	bool flagRunNumber = kFALSE; 

	if(argc == 2)
	{
		TString sa(argv[1]);    

		if(ExtendName(sa) == "steer")
		{
			char * steeringFileName=argv[1];
			SteerSetup(steeringFileName);
			flagslcio = kTRUE; 
			flagEventNumber = kTRUE; 
			flagRunNumber = kTRUE; 
			if(gdmlrootfileName.find("root") != string::npos) 
			{
				flaggeoroot = kTRUE;
			}
			else if(gearFileName.find("xml") != string::npos) 
			{	
				flaggeoxml = kTRUE;    
			}
		} 
		else if(ExtendName(sa) == "root")
		{
			std::cout<<"Display only the geometry according to root gdml file"<<std::endl;
			gdmlrootfileName = argv[1];   flaggeoroot = kTRUE; 
		} 
		else if(ExtendName(sa) == "gearxml")
		{
			std::cout<<"Display only the geometry according to gearxml file"<<std::endl;
			gearFileName = argv[1];   flaggeoxml = kTRUE;
		}
		else if(ExtendName(sa) == "slcio")
		{
			std::cout<<"Display the first event in LCIO file without geometry!"<<std::endl;
			lcioFileName = argv[1];   flagslcio = kTRUE; 
			runNumber = -1;
			event_id = 1;
		}

	}
	else if(argc == 3)
	{
		lcioFileName = argv[1];
		TString lciofile = lcioFileName;
		if(ExtendName(lciofile) == "slcio")
		{
			flagslcio = kTRUE; 
			runNumber = -1;

			TString sa(argv[2]);
			if(ExtendName(sa) == "root") 
			{gdmlrootfileName = argv[2];
				flaggeoroot = kTRUE;
				event_id = 1;
			}else if(ExtendName(sa) == "gearxml"){
				gearFileName = argv[2];   flaggeoxml = kTRUE;
				event_id = 1;
			}else{
				event_id = atoi( argv[2] );
				flagEventNumber = kTRUE;
			}
		}
		else
		{ 
			std::cout<<"Check your input, the first argument should be lcio file"<<std::endl;
		}

	}
	else if(argc == 4)
	{
		lcioFileName = argv[1];   
		flagslcio = kTRUE;
		TString sa(argv[2]);

		if(ExtendName(sa) == "root")
		{ 
			gdmlrootfileName = argv[2];   
			flaggeoroot = kTRUE;
		}
		else if(ExtendName(sa) == "gearxml")
		{
			gearFileName = argv[2];   
			flaggeoxml = kTRUE;
		} 
		event_id = atoi( argv[3] );   
		flagEventNumber = kTRUE;
		runNumber = -1;
	}
	else if(argc == 5)
	{
		lcioFileName = argv[1];	
		flagslcio = kTRUE; 
		TString sa(argv[2]);
		if(ExtendName(sa) == "root")
		{
			gdmlrootfileName = argv[2];   
			flaggeoroot = kTRUE;
		}
		else if(ExtendName(sa) == "gearxml")
		{
			gearFileName = argv[2];   
			flaggeoxml = kTRUE;
		}	
		runNumber = atoi( argv[3] );    
		flagRunNumber = kTRUE; 
		event_id = atoi( argv[4] );     
		flagEventNumber = kTRUE; 
	}

	std::cout<<" Flag Status Summary: FlagGeoROOT = "<<flaggeoroot<<"\t"<<"FlagGeoXML = "<<flaggeoxml<<"\t"<<"FlagSLCIO = "<<flagslcio<<std::endl;

	gMultiView = new MultiView;

	if(flagslcio)
	{
		lcReader = LCFactory::getInstance()->createLCReader(LCReader::directAccess) ;
		lcReader->open( lcioFileName.c_str() );
		load_event(event_id);
	}

	if(flaggeoroot)
	{
		BuildGeoGDMLRoot(gdmlrootfileName);
	}
	else if(flaggeoxml)
	{
		std::cout<<" Gear Geometry is no longer supporting since Druid 2.0! Please use gdml root file"<<std::endl;
	}

	gEve->FullRedraw3D(kTRUE);

	theApp->Run();

	if(flagslcio) 
	{
		lcReader->close();
	}

	return(0);

}


