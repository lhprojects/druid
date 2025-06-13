
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//    Build Simulated CaloHits and CaloHits into TEveBox,
//    // Specify the orientation, color and size according to different
//    // settings (detector concepts, subdetctors, information to emphasize...)
//    //
//															                 //
//    Date:   03 Dec 2010                                                    //
//    Author: Manqi Ruan (LLR)                                               //
//																			 //
//	  Last Modified, 17, Nov 2011, Cleanning
////
//    13, May 2011, adding Timing info support
//    //
//																			 //
///////////////////////////////////////////////////////////////////////////////

#include "EVENT/CalorimeterHit.h"
#include "EVENT/LCCollection.h"
#include "EVENT/LCEvent.h"
#include "EVENT/RawCalorimeterHit.h"
#include "EVENT/SimCalorimeterHit.h"
#include "GlobalDefs.hh"
#include "IMPL/LCTOOLS.h"
#include "TColor.h"
#include "TEveBox.h"
#include "TEveBoxSet.h"
#include "TEveRGBAPalette.h"
#include "TEveRGBAPaletteOverlay.h"
#include "TStyle.h"
#include "TVector3.h"
#include "UTIL/CellIDDecoder.h"
#include "lcio.h"

using namespace lcio;
using namespace EVENT;
using namespace std;

double cellEnergyThresh = 0.0;
int HitColourType = 0;
float HitColourFactor = 0.3;
float HitColourUnitFactor = 1000; ///< change the unit of over/under flow limit by into MeV/HitColourUnitFactor
float HitColourOverflowLimit = 10;
float HitColourUnderflowLimit = 0.;
bool HitSizeLog = kFALSE;
bool Flag_AttachTextToHit = kFALSE;
int GlobalRandomColorIndex = 0;
float cellfactor = 1.0;
float HitEnergyCut = 0.2;  // Mips
bool HiddenLowESimCell = kTRUE;

extern int flagdetectortype;
extern int event_id;
extern std::map<float, int>
    randomColor;  // used to give an random color to each MCParticle
extern std::map<MCParticle*, int>
    OriginColor;  // used to give ... to Origin Particle

extern TEveRGBAPalette* p;
// extern TEveRGBAPaletteOverlay* po;
extern TEveBoxSet* cal_shell;

namespace PrototypeTB {
// namespace PrototypeTB
const double SCECALCELLWIDTH = 0.5;         // cm
const double SCECALCELLLENGTH = 4.5;        // cm
const double SC_ECAL_CELL_THICKNESS = 0.2;  // cm
const float MIP_SC_ECAL = 0.305e-3;         // GeV
const float MIP_ADC_SC_ECAL = 450;

const double AHCAL_CELL_WIDTH = 4;        // cm
const double AHCAL_CELL_LENGTH = 4;       // cm
const double AHCAL_CELL_THICKNESS = 0.3;  // cm
const float MIP_AHCAL = 0.426e-3;         // GeV
const float MIP_ADC_AHCAL = 450;
};  // namespace PrototypeTB

const float MIPSIMECAL = 1.8e-4;
const float MIPSIMAHCAL = 1.0e-3;
const float MIPSIMDHCAL = 5.5e-7;
const float MIPSIMMUON = 1.0e-6;
const float MIPDIGIECAL = 0.007;
// const float MIP_SC_ECAL = 0.00167;//GeV
const float MIPDIGIAHCAL = 0.031;
// const float MIPDIGIDHCAL = 5.5e-5;	//not sure yet...depending on the
// digitization
const float MIPDIGIDHCAL = 1;
const float MIPDIGIMUON = 0.01;
// const float MIPSIMLCAL = 1.2e-4;	//similar for LHcal and LumiCal

// TEveBox* BoxPhi( TVector3 &HitPos, TVector3 &Scale, int Type, int
// SegOrStaveNumber, float HitEnergy );

TEveElementList* CaloHits(LCCollection* col, string name) {
  bool isSimHit = col->getTypeName() == LCIO::SIMCALORIMETERHIT;
  bool isRawHit = col->getTypeName() == LCIO::RAWCALORIMETERHIT;

  TEveElementList* cal = new TEveElementList;

  cal_shell->Reset(TEveBoxSet::kBT_FreeBox, kFALSE, 64);
  cal_shell->SetMainAlpha(1);
                      
  cal_shell->SetPalette(p);
  cal_shell->SetPickable(false);

  int icol(0);
  int mothercol(0);

  float HitX, HitY, HitZ, HitEn;
  char* HitRegion = 0;
  int MCPID, StaveNum, IndexI,
      IndexJ;  // IndexI and StaveNum, as well as Module Number and IndexK are
               // coded in ID0 of Hit...
  int MCHitType;  // EM, Had, Neutron, Proton;
  float MCenergy;
  int MotherPID;
  int LayerNum;
  float Mipcount;      // Meause Energy deposition in unit of MIP
  float MotherEnergy;  // The energy and PID of Origin
  float MCContTime;
  // OriginColor.clear();

  double MaxHitEnergy = 0;
  double MinHitEnergy = 1000000;  // GeV

  HitX = 0;
  HitY = 0;
  HitZ = 0;
  HitEn = 0;
  MCPID = 0;
  MCContTime = -10;
  MCHitType = 0;
  Mipcount = 0;
  LayerNum = 0;
  MCenergy = 0;
  StaveNum = 0;
  IndexI = 0;
  IndexJ = 0;
  MotherPID = 0;
  MotherEnergy = 0;

  int nHits = col->getNumberOfElements();

  cal->SetName(name.c_str());

  string SubDCollection(name, 0, 30);
  string SubDAHcal(name, 0, 13);
  string SubD(name, 0, 4);
  string SubDPos(name, 4, 3);
  string SubBarrel(name, 0, 10);

  gStyle->SetPalette(1, 0);
  TEveRGBAPalette* pal = new TEveRGBAPalette(0, 100);

  float HalfX = 0;
  float HalfY = 0;
  float HalfZ = 0;
  float s1, s2;
  float SF = isSimHit ? 1.0 : 1.2;  // Scale Factor
  float s1X, s1Y, s2X, s2Y;

  cout << "  Calo hits collection : " << name << ". Number of Hits: " << nHits
       << ", typeName:" << col->getTypeName() << endl;
  //cout << endl; Hao

  int MaxEnergyDepoID = 0;  // To record the MCPID of Simulated Calo Hits
  int MCPIDindex = 1;
  int Originindex = 0;
  int MotherEnergySignature = 0;

  // CellIDDecoder<SimCalorimeterHit> idDecoder(col);

  TEveElementList* cal2 = new TEveElementList;
  cal2->SetName("FiberHit");

  /*
  for(int i=0; i<nHits; i++)
{
          CalorimeterHit* hit = dynamic_cast<CalorimeterHit*>(
col->getElementAt( i ) ) ; CellIDDecoder<CalorimeterHit> idDecoder(col);
          HitX=hit->getPosition()[0]/10.0;
          HitY=hit->getPosition()[1]/10.0;
          HitZ=hit->getPosition()[2]/10.0;
          HitEn=hit->getEnergy();

          Mipcount=HitEn/MIP_SC_ECAL;
          int IndexHitColorLine = int(HitColourFactor*Mipcount*0.5 + 51);
          int IndexHitColorMain = int(HitColourFactor*Mipcount*0.5 + 51);

          if(IndexHitColorLine > 100) IndexHitColorLine = 100;
          if(IndexHitColorMain > 100) IndexHitColorMain = 100;

          if(name == "AllHits") HitEn = HitEn * 30000;	//Only used for all
cleaned hits check in arbor

          TVector3 HitPosition(HitX, HitY, HitZ);


          TEveBox * q = new TEveBox;

          TVector3 HitScale(0.5,0.5,0.2);

          {
                  if((int)(HitZ*100 -10)%199 ==0)
                  {

                          HitScale.SetXYZ(SCECALCELLLENGTH,SCECALCELLWIDTH,SC_ECAL_CELL_THICKNESS);
                          //HitPosition.SetY(0);
                  }
                  else
                  {
                          HitScale.SetXYZ(SCECALCELLWIDTH,SCECALCELLLENGTH,SC_ECAL_CELL_THICKNESS);
                          //HitPosition.SetX(0);
                  }
                  q = BoxPhi( HitPosition, HitScale, -1, 0, HitEn );
          }

                  q->SetPickable(kTRUE);		//Set boxset pickable by
mouse q->SetMainColor(19); q->SetLineColor(19);

                  cal->AddElement(q);
  }
  */
  if(isRawHit) {
  	std::cout << ", RawHit will be skiped" << std::endl;
  }
  for (int i = 0; i < nHits; i++) {
    //------------------ Set Hit Energy & CellIDs ------------------
    if (isRawHit && false) {
      RawCalorimeterHit* hit =
          dynamic_cast<RawCalorimeterHit*>(col->getElementAt(i));
      CellIDDecoder<RawCalorimeterHit> idDecoder(col);

      /*
       * !!! Should be varify!
       */
      if (idDecoder(hit)["hit"] == 0) {
        continue;
      }
      LayerNum = idDecoder(hit)["layer"];
      IndexI = idDecoder(hit)["chip"] -
               1;  ///< For test-beam prototype, IndexI denots chip index.
                   ///< WARNING: should be valiated
      IndexJ = idDecoder(hit)["channel"];  ///< ..., ... channel index.
      //< decode position with dedicate functions.
      std::cout << LayerNum << ", " << IndexI << ", " << IndexJ << std::endl;
      TVector3 pos;
      if (name == "ECALRawHits") {
        pos = TVector3(GetScEcalHitPos(LayerNum, IndexI, IndexJ));
      } else {
        pos = TVector3(GetAhcalHitPos(LayerNum, IndexI, IndexJ));
      }

      HitX = pos.X() / 10.;
      HitY = pos.Y() / 10.;
      HitZ = pos.Z() / 10.;
      HitEn = hit->getAmplitude();
      MCContTime = hit->getTimeStamp();
    } else if (isSimHit) {
      SimCalorimeterHit* hit =
          dynamic_cast<SimCalorimeterHit*>(col->getElementAt(i));
      HitX = hit->getPosition()[0] / 10.0;
      HitY = hit->getPosition()[1] / 10.0;
      HitZ = hit->getPosition()[2] / 10.0;
      HitEn = hit->getEnergy();

      CellIDDecoder<SimCalorimeterHit> idDecoder(col);
      //IndexI = idDecoder(hit)["I"];
      //IndexJ = idDecoder(hit)["J"];
      //if (SubD == "Ecal" || SubD == "Hcal") LayerNum = idDecoder(hit)["K-1"];
      IndexI = 0;
      IndexJ = 0;
      StaveNum = 0;

      // IndexI = (hit->getCellID0() & 0x00007FC0)>>6;
      if (flagdetectortype == 3 || flagdetectortype == 4)
        StaveNum = (IndexI % 128) / 8;  // SID ECAL Stave
      else {
	 StaveNum = IndexI;
        //StaveNum = idDecoder(
        //   hit)["S-1"];  // StaveNum=(hit->getCellID0() & 0x00000038)>>3;
      }
      // LayerNum=(hit->getCellID0() & 0x3F000000)>>24;

      if (hit->getNMCContributions() > 0 && name != "LumiCalCollection") {
        float Emax = -0.01;
        MaxEnergyDepoID = 0;
        for (int k = 0; k < hit->getNMCContributions(); k++) {
          MCParticle* hitMCPart =
              dynamic_cast<MCParticle*>(hit->getParticleCont(k));

          if (hit->getEnergyCont(k) > Emax) {
            MCPID = hitMCPart->getPDG();  // According to assigned Track index:
                                          // lots of mother contribution
            MCHitType = hit->getPDGCont(
                k);  // According to local deposited track PID: e-, e+, pi, K...
            MCenergy = hitMCPart->getEnergy();
            Emax = hit->getEnergyCont(k);
            MaxEnergyDepoID = k;
            MCContTime = hit->getTimeCont(k);

            if (randomColor.find(MCenergy) != randomColor.end()) {
              MCPIDindex = randomColor[MCenergy];
            } else {
              randomColor[MCenergy] = icol++;
              MCPIDindex = randomColor[MCenergy];
            }

            // MCPIDindex = MCPID;				//Color
            // according to the ID
          }
        }

        // Try to get Mother

        if (HitColourType == 2) {
          MCParticle* hitMCPart0 =
              dynamic_cast<MCParticle*>(hit->getParticleCont(MaxEnergyDepoID));

          for (int s = 0; s < 100; s++) {
            MCParticleVec mothers = hitMCPart0->getParents();
            if (mothers.size() > 0) {
              MCParticle* Mother = mothers[0];

              if (Mother->getPDG() == 92 || Mother->getParents().size() == 0 ||
                  Mother->getPDG() == 23 ||
                  Mother->getPDG() ==
                      25)  // || abs(Mother->getPDG())<6 ) 5, -5, 25;
              // if( (fabs(Mother->getPDG()) < 6 || Mother->getPDG() == 25) &&
              // Mother->getParents().size() == 0 )
              {
                MotherPID = Mother->getPDG();
                MotherEnergy = Mother->getEnergy();

                if (OriginColor.find(Mother) == OriginColor.end()) {
                  if (MotherPID == 23)
                    Originindex = 2;
                  else if (MotherPID == 25)
                    Originindex = 4;
                  else
                    Originindex = 9;
                  OriginColor[Mother] = Originindex;
                } else {
                  Originindex = OriginColor[Mother];
                }
                break;
              } else {
                hitMCPart0 = Mother;
              }
            }
          }
          // cout<<"Mother "<<MCPID<<"\t"<<MotherPID<<"\t"<<MotherEnergy<<endl;
        }
      }

    } else {
      CalorimeterHit* hit = dynamic_cast<CalorimeterHit*>(col->getElementAt(i));
      CellIDDecoder<CalorimeterHit> idDecoder(col);

      HitX = hit->getPosition()[0] / 10.0;
      HitY = hit->getPosition()[1] / 10.0;
      HitZ = hit->getPosition()[2] / 10.0;
      HitEn = hit->getEnergy();

      if (HitEn > MaxHitEnergy) MaxHitEnergy = HitEn;
      if (HitEn < MinHitEnergy) MinHitEnergy = HitEn;

      if (name == "AllHits")
        HitEn = HitEn * 30000;  // Only used for all cleaned hits check in arbor

      if (name == "ScEcalHit" or name == "AhcalHit") {
        LayerNum = int(hit->getCellID0() / 1e5);
        IndexI =
            int(hit->getCellID0() % 100000 /
                1e4);  ///< For test-beam prototype, IndexI denots chip index.
        IndexJ = int(hit->getCellID0() % 100);  ///< ..., ... channel index.
      } else if (SubD == "ECAL" || SubD == "HCAL") {
        //IndexI = idDecoder(hit)["I"];
        //IndexJ = idDecoder(hit)["J"];
        //LayerNum = idDecoder(hit)["K-1"]; Hao
	IndexI = 0;
	IndexJ = 0;
        LayerNum = 0;

        /*
        IndexI = (hit->getCellID0() & 0x00007FC0)>>6;
        if (flagdetectortype == 3 || flagdetectortype == 4)  StaveNum =
        (IndexI%128)/8;			//Sid ECAL StaveNum else
        StaveNum=(hit->getCellID0() & 0x00000038)>>3;
        LayerNum=(hit->getCellID0() & 0x3F000000)>>24;
        */
      }
    }

    TVector3 HitPosition(HitX, HitY, HitZ);

    if (SubD == "Muon" || SubD == "MUON") {
      SF = 10.0;
    } else {
      //      if (HitSizeLog == 0) {
      if (((SubD == "Hcal" || SubD == "HCAL") && flagdetectortype == 0) ||
          SubD == "AhcC") {
        SF = 3.0 * cellfactor;  // AHCAL Hits
      } else
        SF = cellfactor;
      //		  }
      //      else if (HitSizeLog == 1) {
      //		  if(SubD=="Hcal") SF = (8.0 + 0.1*log(HitEn));
      //		  else SF = (6.0 + 0.2*log(HitEn));
      //      }
    }

    int IndexHitColorLine = 0;
    int IndexHitColorMain = 0;

    float ecalcali = 1.;
    float hcalcali = 1.;
    float tcmtcali = 1.;
    Mipcount = 0;

    if (isRawHit) {
      ecalcali = PrototypeTB::MIP_ADC_SC_ECAL;
      hcalcali = PrototypeTB::MIP_ADC_AHCAL;
    } else if (isSimHit) {
      ecalcali = MIPSIMECAL;
      tcmtcali = MIPSIMMUON;
      if (flagdetectortype == 1 || flagdetectortype == 3) {
        hcalcali = MIPSIMDHCAL;
      }
      //			else if (flagdetectortype == 0 ||
      // flagdetectortype == 2 || flagdetectortype == 4) {hcalcali =
      // MIPSIMAHCAL;}	//ILD00(AHCAL) and CALICE TB with AHCAL Prototype
      else {
        hcalcali = MIPSIMAHCAL;
      }
    } else {
      if (flagdetectortype == 2 ||
          flagdetectortype == 10)  // Test beams, actually only
      {
        ecalcali = 1.;
        hcalcali = 1.;
        tcmtcali = 1.;
      }  // already calibrated TB Hits
      else {
        ecalcali = MIPDIGIECAL;
        tcmtcali = MIPDIGIMUON;
        if (flagdetectortype == 1 || flagdetectortype == 3) {
          hcalcali = MIPDIGIDHCAL;
        }  // Assume Cali Const for DHCAL/SDHCAL are the same
        //				else if ( flagdetectortype == 0 ||
        // flagdetectortype == 4) {hcalcali = MIPDIGIAHCAL;}
        else {
          hcalcali = MIPDIGIAHCAL;
        }
      }
    }

    //------------------ Set MIP Count ------------------
    if (name == "HCALRawHits") {
      Mipcount = HitEn / hcalcali;
    } else if (name == "ECALRawHits") {
      Mipcount = HitEn / ecalcali;
    } else if (name == "ScEcalHit") {
      Mipcount = HitEn / PrototypeTB::MIP_SC_ECAL;
    } else if (name == "AhcalHit") {
      Mipcount = HitEn / PrototypeTB::MIP_AHCAL;
    } else if (SubD == "Ecal" || SubD == "ECAL" || SubD == "ecal" ||
               SubD == "Prot")  // Prot rep sth as "ProtoSD03Collection", ECAL
                                // Prototype
    {
      Mipcount = HitEn / ecalcali;
    } else if (SubD == "Hcal" || SubD == "HCAL" || SubD == "hcal") {
      Mipcount = HitEn / hcalcali;
    } else  // muon, ecal...
    {
      Mipcount = HitEn / tcmtcali;
      ;
    }

    if (HitColourType == 1 || !isSimHit) {  // rep Hit Energy in Color
      // TColor Color;
      // int RGB = (int)((HitEn-MinHitEnergy)/(MaxHitEnergy-MinHitEnergy))*256;

      // Color.SetRGB(RGB,0,0); //it's red
      // IndexHitColorLine = Color.GetNumber(); IndexHitColorMain =
      // Color.GetNumber();

      // int color = (int)((HitEn-300)/(400-300))*256;
      // IndexHitColorLine = color; IndexHitColorMain = color;

      /*//-----------set hit energy in color----------------
      if(Mipcount < 0.2)
              {
                      IndexHitColorLine = 51; IndexHitColorMain = 51;//grey
              }
              else if(Mipcount<5 && Mipcount>=0.2)
              {
                      IndexHitColorLine = 4; IndexHitColorMain = 4;//blue
              }
              else if(Mipcount<10 && Mipcount>=5)
              {
                      IndexHitColorLine = 3; IndexHitColorMain = 3;//green
              }
              else if(Mipcount>=10)
              {
                      IndexHitColorLine = 2; IndexHitColorMain = 2;//red
              }
      //--------------------------------------------------
      */

      if (isSimHit && Mipcount < HitEnergyCut) {
        IndexHitColorLine = 10;
        IndexHitColorMain = 10;
      }  // White to indicate the Hitcolor for simulated Hits: with changing
         // bkgrd to white it will disappear.
      else {
        IndexHitColorLine = int(HitColourFactor * Mipcount * 0.5 + 51);
        IndexHitColorMain = int(HitColourFactor * Mipcount * 0.5 + 51);
      }

      if (IndexHitColorLine > 100) IndexHitColorLine = 100;
      if (IndexHitColorMain > 100) IndexHitColorMain = 100;

      if ((SubD == "Hcal" || SubD == "HCAL") && cellEnergyThresh != 0) {
        IndexHitColorMain = 10;  // to denote SemiDHCAL case

        if (Mipcount * cellEnergyThresh < 0.2) {
          IndexHitColorLine = 51;
          IndexHitColorMain = 51;
        } else if (Mipcount * cellEnergyThresh < 5 &&
                   Mipcount * cellEnergyThresh >= 0.2) {
          IndexHitColorLine = 4;
          IndexHitColorMain = 4;
        } else if (Mipcount * cellEnergyThresh < 10 &&
                   Mipcount * cellEnergyThresh >= 5) {
          IndexHitColorLine = 3;
          IndexHitColorMain = 3;
        } else if (Mipcount * cellEnergyThresh >= 10) {
          IndexHitColorLine = 2;
          IndexHitColorMain = 2;
        }
      } else if ((SubD == "SDHC"))  // From Trivent
      {
        if (HitEn == 1)  // Inverted ... 1 = 2
        {
          IndexHitColorLine = 3;
          IndexHitColorMain = 3;
        } else if (HitEn == 2) {
          IndexHitColorLine = 4;
          IndexHitColorMain = 4;
        } else if (HitEn == 3) {
          IndexHitColorLine = 2;
          IndexHitColorMain = 2;
        } else {
          IndexHitColorLine = 5;
          IndexHitColorMain = 5;
        }
      }

    } else if (HitColourType ==
               0) {  // rep direct deposition Particle Type with Color

      switch (MCPID) {
        case (12):  // Neutrinos: Normally doesn't creates hits...
        case (14):
        case (16):
        case (-12):
        case (-14):
        case (-16):
          IndexHitColorMain = 66;
          break;

        case (22):  // Gamma!
          IndexHitColorMain = 89;
          break;

        case (11):  // electrons & positrons
          IndexHitColorMain = 56;
          break;

        case (-11):
          IndexHitColorMain = 99;
          break;

        case (211):  // Pions
          IndexHitColorMain = 97;
          break;

        case (-211):
          IndexHitColorMain = 53;
          break;

        case (13):  // Muons
          IndexHitColorMain = 100;
          break;

        case (-13):
          IndexHitColorMain = 51;
          break;

        case (2212):  // Protons
          IndexHitColorMain = 96;
          break;

        case (-2212):
          IndexHitColorMain = 64;
          break;

        case (2122):  // Neutrons
          IndexHitColorMain = 85;
          break;

        case (130):  // Klong
          IndexHitColorMain = 80;
          break;

        case (321):  // Charged Kaon
          IndexHitColorMain = 98;
          break;

        case (-321):
          IndexHitColorMain = 54;
          break;

        default:  // Any other possibles
          IndexHitColorMain = 38;
          break;
      }
      IndexHitColorLine = IndexHitColorMain;
    } else if (HitColourType == 2) {
      // IndexHitColorMain = (Originindex*13 + GlobalRandomColorIndex*7)%50+51;
      IndexHitColorMain = Originindex;
    } else if (HitColourType == 3) {
      IndexHitColorMain = ((MCPIDindex % 2) * 25 + MCPIDindex * 5 +
                           GlobalRandomColorIndex * 13) %
                              50 +
                          51;
    } else if (HitColourType == 4)  // Uniformed Color for Simulated Hit: Blue
    {
      IndexHitColorMain = 5;
    } else if (HitColourType == 5) {
      if (MCHitType == 11 || MCHitType == -11 || MCHitType == 22) {
        IndexHitColorMain = 4;       // Blue for EM
      } else if (MCHitType == 2112)  // Neutron
      {
        IndexHitColorMain = 5;  // Yellow
      } else {
        IndexHitColorMain = 2;
      }
    } else if (HitColourType == 6) {
      if (MCContTime == -10) {
        std::cout << "Timing info not available in current data file"
                  << std::endl;
        IndexHitColorMain = 3;
      } else {
        if (MCContTime < 150)  // 150ns as electronic integration time...
        {
          IndexHitColorMain =
              int(0.5 * HitColourFactor * (MCContTime - 5.0) + 51);
        } else {
          IndexHitColorMain = 21;
        }
      }
      if (IndexHitColorMain > 100) IndexHitColorMain = 100;
      IndexHitColorLine = IndexHitColorMain;
    }

    if (Mipcount > HitEnergyCut || !HiddenLowESimCell || !isSimHit) {
      TEveBox* q = new TEveBox;

      TVector3 HitScale(SF, SF, 0.1 * SF);
      TVector3 TCMTScale(100.0, 2.0, 5.0);

      if (SubBarrel == "ECALBarrel" || SubBarrel == "EcalBarrel") {
        if (flagdetectortype < 2) {
          q = BoxPhi(HitPosition, HitScale, 0, StaveNum, HitEn);
        } else {
          q = BoxPhi(HitPosition, HitScale, 2, StaveNum, HitEn);
        }
      } else if (SubBarrel == "HCALBarrel" || SubBarrel == "HcalBarrel") {
        if (flagdetectortype == 1)  // a la Videau
        {
          q = BoxPhi(HitPosition, HitScale, 0, StaveNum, HitEn);
        } else if (flagdetectortype == 0)  // Tesla
        {
          q = BoxPhi(HitPosition, HitScale, 1, 8, HitEn);
        } else  // CLIC HCAL
        {
          q = BoxPhi(HitPosition, HitScale, 1, 12, HitEn);
        }
      } else if (SubDPos == "End") {
        q = BoxPhi(HitPosition, HitScale, -1, 0, HitEn);
      } else if (SubD == "MUON" || SubD == "Muon") {
        if (flagdetectortype < 2) {
          q = BoxPhi(HitPosition, HitScale, 1, 12, HitEn);

        } else {
          q = BoxPhi(HitPosition, HitScale, 1, 8, HitEn);
        }
      } else if (SubBarrel == "TcmtCalori") {
        q = BoxPhi(HitPosition, TCMTScale, 1, 8, HitEn);
      } else if (name == "ScEcalHit" or
                 name ==
                     "ECALRawHit")  // MyECALHits is for upward compatibility
      {
        // HitEn *= 1e-3; ///< hit energy from Jiaxuan's root files are under
        // the unit of MeV

        IndexHitColorMain = int(HitColourFactor * Mipcount * 0.5 + 51);
        IndexHitColorLine = int(HitColourFactor * Mipcount * 0.5 + 51);
        // if((int)(HitZ*100 -10)%199 ==0)
        if (LayerNum % 2 == 0) {
          HitScale.SetXYZ(PrototypeTB::SCECALCELLLENGTH,
                          PrototypeTB::SCECALCELLWIDTH,
                          PrototypeTB::SC_ECAL_CELL_THICKNESS);
          // HitPosition.SetY(0);
        } else {
          HitScale.SetXYZ(PrototypeTB::SCECALCELLWIDTH,
                          PrototypeTB::SCECALCELLLENGTH,
                          PrototypeTB::SC_ECAL_CELL_THICKNESS);
          // HitPosition.SetX(0);
        }

        q = BoxPhi(HitPosition, HitScale, -1, 0, HitEn);
      } else if (name == "AhcalHit" or name == "HCALRawHits") {
        // HitEn *= 1e-3; ///< hit energy from Jiaxuan's root files are under
        // the unit of MeV
        // std::cout << "HCAL Raw: " << HitX << ", " << HitY << ", " << HitZ <<
        // std::endl;

        IndexHitColorMain = int(HitColourFactor * Mipcount * 0.5 + 51);
        IndexHitColorLine = int(HitColourFactor * Mipcount * 0.5 + 51);
        // if((int)(HitZ*100 -10)%199 ==0)
        HitScale.SetXYZ(PrototypeTB::AHCAL_CELL_WIDTH,
                        PrototypeTB::AHCAL_CELL_LENGTH,
                        PrototypeTB::AHCAL_CELL_THICKNESS);

        q = BoxPhi(HitPosition, HitScale, -1, 0, HitEn);
        Float_t verts[24] = {
          q->GetVertex(0)[0], q->GetVertex(0)[1], q->GetVertex(0)[2],
          q->GetVertex(1)[0], q->GetVertex(1)[1], q->GetVertex(1)[2],
          q->GetVertex(2)[0], q->GetVertex(2)[1], q->GetVertex(2)[2],
          q->GetVertex(3)[0], q->GetVertex(3)[1], q->GetVertex(3)[2],
          q->GetVertex(4)[0], q->GetVertex(4)[1], q->GetVertex(4)[2],
          q->GetVertex(5)[0], q->GetVertex(5)[1], q->GetVertex(5)[2],
          q->GetVertex(6)[0], q->GetVertex(6)[1], q->GetVertex(6)[2],
          q->GetVertex(7)[0], q->GetVertex(7)[1], q->GetVertex(7)[2]
        };
        cal_shell->AddBox(verts);
        cal_shell->DigitValue(HitEn * 1e3 * HitColourUnitFactor);
      } else  // Including all the Endcap Hits
      {
        HitScale.SetXYZ(4.0, 4.0, 0.3);
        q = BoxPhi(HitPosition, HitScale, -1, 0, HitEn);
      }

      if (IndexHitColorLine == 0) IndexHitColorLine = IndexHitColorMain;

      q->SetPickable(kTRUE);  // Set boxset pickable by mouse
      q->SetMainColor(IndexHitColorMain);
      q->SetLineColor(0);
      q->SetMainAlpha(0.8);

      if (Flag_AttachTextToHit) {
        if (isSimHit) {
          q->SetTitle(Form(
              "Simulated Calo Hits %s\n"
              "EventNum=%d, SubDetector=%s\n"
              "Hit Energy = %.3e keV ~ %.3e Mip, Thresh = %3.e keV\n"
              "MCPID = %d, MCenergy = %.3f\n "
              "OriginPID = %d, OriginEnergy = %.3f\n"
              "PosX = %.3f mm, PosY = %.3f mm, PosZ = %.3f mm, StaveNum = %d, "
              "LayerNum = %d, IndexI = %d, IndexJ = %d, Originindex = %d, "
              "HitTime = %.3fns",
              name.c_str(), event_id, SubDCollection.c_str(), HitEn * 1000000,
              Mipcount, cellEnergyThresh * 1000000, MCPID, MCenergy, MotherPID,
              MotherEnergy, 10 * HitX, 10 * HitY, 10 * HitZ, StaveNum, LayerNum,
              IndexI, IndexJ, Originindex, MCContTime));
        } else {
          /*
             if(flagdetectortype != 2)
             q->SetTitle(Form( "Reconstructed Calo Hit %s \n"
             "EventNum=%d, SubDetector=%s\n"
             "Hit Energy = %.3e keV ~ %.3e Mip, StaveNum = %d\n"
             "PosX = %.3f mm, PosY = %.3f mm, PosZ = %.3f mm",
             name.c_str(), event_id, SubDCollection.c_str(), HitEn*1000000,
             Mipcount, StaveNum, 10*HitX, 10*HitY, 10*HitZ));

             else
             */
          q->SetTitle(
              Form("Reconstructed Calo Hit %s \n"
                   "EventNum=%d, SubDetector=%s, (Layer, Chip, "
                   "Channel)=(%d,%d,%d) \n"
                   "Hit Energy = %.3e GeV ~ %.3e Mip , StaveNum = %d\n"
                   "PosX = %.3f mm, PosY = %.3f mm, PosZ = %.3f mm",
                   name.c_str(), event_id, SubDCollection.c_str(), LayerNum,
                   IndexI, IndexJ, HitEn, Mipcount, StaveNum, 10 * HitX,
                   10 * HitY, 10 * HitZ));
        }
      }
      cal->AddElement(q);
    }
  }

  /*
     TEveBox * calicesdhcal = new TEveBox;
     TVector3 geoscale, geopos;
     geoscale.SetXYZ(100,100,100);
     geopos.SetXYZ(0, 0, 0);
     calicesdhcal = BoxPhi(geopos, geoscale, -1, 0, 20);
  //calicesdhcal->SetTransparency(70);

  cal->AddElement(calicesdhcal);
  */

  cout << "OriCoSize " << OriginColor.size() << ", " << mothercol << endl;

  return cal;
}

TEveBox* BoxPhi(TVector3& HitPos, TVector3& Scale, int Type,
                int SegOrStaveNumber, float HitEnergy) {
  // Type = -1; EndCap
  // Type = 0; Barrel, Based On Segment Number, ultilized for a la Videau Model
  // HCAL Type = 1; Barrel, Based On Phi Angle, ultilized for TESLA Model HCAL
  // Type = 2; Barrel, Based On Segment Number, ultilized for SID ECAL

  TEveBox* q = new TEveBox();
  q->SetName(Form("HitE = %.3e MeV", HitEnergy * 1000));

  gStyle->SetPalette(1, 0);
  q->SetMainTransparency(60);

  const float DegToRad = 1.74532925199432781e-02;
  const float Pi = 3.1415926535;

  float HitX = HitPos(0);
  float HitY = HitPos(1);
  float HitZ = HitPos(2);

  float s1X = 0;
  float s1Y = 0;
  float s1Z = 0;
  float s2X = 0;
  float s2Y = 0;
  float s2Z = 0;

  float phiAngle = 0;
  float SX = 0;
  float SY = 0;
  float SZ = 0;
  float SX1 = 0;
  float SY1 = 0;
  float SZ1 = 0;

  if (Type == -1)  // Based on EndCap
  {
    s1X = -1 * Scale(0) / 2.0;
    s2X = -1 * s1X;
    s1Y = -1 * Scale(1) / 2.0;
    s2Y = s1Y;
    s1Z = Scale(2) / 2.0;
    s2Z = s1Z;

    SX = fabs(0.5 * Scale(0));
    SY = fabs(0.5 * Scale(1));
    SZ = fabs(0.5 * Scale(2));
    SX1 = SX;
    SY1 = SY;
    SZ1 = SZ;
  } else {
    if (Type == 1)  // Based on Phi & SegmentNumbers(>=4 at least)
    {
      if (HitPos.Phi() > 0) {
        phiAngle = 2 * Pi / SegOrStaveNumber *
                   int(HitPos.Phi() * SegOrStaveNumber / 2 / Pi + 0.5);
      } else if (HitPos.Phi() <= 0) {
        phiAngle = 2 * Pi / SegOrStaveNumber *
                   int(HitPos.Phi() * SegOrStaveNumber / 2 / Pi - 0.5);
      }
    }
    if (Type == 0)  // Based on StaveNumber		Currently Only used for
                    // ILD00 a la Videau model
    {
      if (SegOrStaveNumber == 2 || SegOrStaveNumber == 6) {
        phiAngle = 0;
      }
      if (SegOrStaveNumber == 0 || SegOrStaveNumber == 4) {
        phiAngle = 90 * DegToRad;
      }
      if (SegOrStaveNumber == 3 || SegOrStaveNumber == 7) {
        phiAngle = 45 * DegToRad;
      }
      if (SegOrStaveNumber == 1 || SegOrStaveNumber == 5) {
        phiAngle = 135 * DegToRad;
      }
    }
    if (Type == 2)  // Based on StaveNumber		Used for SID ECAL
    {
      if (SegOrStaveNumber == 3 || SegOrStaveNumber == 9) {
        phiAngle = 0;
      }
      if (SegOrStaveNumber == 2 || SegOrStaveNumber == 8) {
        phiAngle = 30 * DegToRad;
      }
      if (SegOrStaveNumber == 1 || SegOrStaveNumber == 7) {
        phiAngle = 60 * DegToRad;
      }
      if (SegOrStaveNumber == 0 || SegOrStaveNumber == 6) {
        phiAngle = 90 * DegToRad;
      }
      if (SegOrStaveNumber == 5 || SegOrStaveNumber == 11) {
        phiAngle = 120 * DegToRad;
      }
      if (SegOrStaveNumber == 4 || SegOrStaveNumber == 10) {
        phiAngle = 150 * DegToRad;
      }
    }

    s1X = -0.5 * (Scale(0) * sin(phiAngle) + Scale(2) * cos(phiAngle));
    s1Y = 0.5 * (Scale(0) * cos(phiAngle) - Scale(2) * sin(phiAngle));
    s2X = -0.5 * (Scale(0) * sin(phiAngle) - Scale(2) * cos(phiAngle));
    s2Y = 0.5 * (Scale(0) * cos(phiAngle) + Scale(2) * sin(phiAngle));
    s1Z = Scale(1) / 2.0;
    s2Z = s1Z;

    SX = fabs(0.5 * (Scale(0) * sin(phiAngle) + Scale(2) * cos(phiAngle)));
    SX1 = fabs(0.5 * (Scale(0) * sin(phiAngle) - Scale(2) * cos(phiAngle)));
    SY = fabs(0.5 * (Scale(0) * cos(phiAngle) - Scale(2) * sin(phiAngle)));
    SY1 = fabs(0.5 * (Scale(0) * cos(phiAngle) + Scale(2) * sin(phiAngle)));
    SZ = fabs(0.5 * Scale(1));
    SZ1 = SZ;
  }

  q->SetVertex(5, HitX + s2X, HitY + s2Y, HitZ - s1Z);
  q->SetVertex(6, HitX + s2X, HitY + s2Y, HitZ + s1Z);

  q->SetVertex(4, HitX + s1X, HitY + s1Y, HitZ - s1Z);
  q->SetVertex(7, HitX + s1X, HitY + s1Y, HitZ + s1Z);

  q->SetVertex(3, HitX - s2X, HitY - s2Y, HitZ + s1Z);
  q->SetVertex(0, HitX - s2X, HitY - s2Y, HitZ - s1Z);

  q->SetVertex(2, HitX - s1X, HitY - s1Y, HitZ + s1Z);
  q->SetVertex(1, HitX - s1X, HitY - s1Y, HitZ - s1Z);

  return q;
}

// oooOO0OOooooooOO0OOooooooOO0OOooooooOO0OOooooooOO0OOooooooOO0OOooooooOO0OOooooooOO0OOooooooOO0OOooo
//************************************************************
//  Position Decoder for Test-Beam Prototype Sc-ECAL and AHCAL
//************************************************************

//------------------ For ScECAL ------------------
// Copy from Jiaxuan Wang @ USTC
// scintillator strips wrt. 6 SPIROC2E chips * 36 channels
TVector3 GetScEcalHitPos(int LayerIDs, int ChipIDs, int ChannelIDs) {
  // static const int layerNu = 32;
  static const int chipNu = 6;
  static const int chnNu = 36;
  int decodeID[chipNu][chnNu] = {
      0,   42,  1,   43,  2,   44,  3,   4,   5,   6,   7,   8,   9,   10,  11,
      12,  54,  13,  55,  14,  56,  15,  57,  16,  58,  17,  59,  18,  60,  19,
      61,  20,  62,  21,  22,  23,  24,  66,  25,  67,  26,  68,  27,  69,  28,
      70,  29,  71,  30,  72,  31,  73,  32,  74,  33,  75,  34,  76,  35,  77,
      36,  78,  37,  79,  38,  80,  39,  81,  40,  82,  41,  83,  149, 148, 147,
      96,  97,  98,  99,  100, 101, 102, 103, 104, 105, 106, 107, 63,  64,  65,
      108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122,
      123, 124, 125, 150, 192, 151, 193, 152, 194, 153, 195, 154, 196, 155, 197,
      156, 198, 157, 199, 158, 200, 159, 201, 160, 202, 161, 203, 162, 204, 163,
      205, 164, 206, 165, 207, 166, 208, 167, 209, 191, 190, 189, 188, 146, 187,
      145, 186, 144, 185, 143, 184, 142, 183, 141, 182, 140, 181, 139, 180, 138,
      179, 178, 177, 176, 175, 174, 173, 172, 171, 170, 128, 169, 127, 168, 126,
      137, 136, 135, 134, 133, 132, 131, 130, 129, 84,  85,  86,  87,  88,  89,
      90,  91,  92,  93,  94,  95,  45,  46,  47,  48,  49,  50,  51,  52,  53,
      210, 210, 210, 210, 210, 210};

  int ScintillatorIDs = decodeID[ChipIDs][ChannelIDs];
  double layerZ;

  const double _xInterval = 5.3;   // 300 um gap in width direction
  const double _yInterval = 45.4;  // 400 um gap in length direction
  const int rowNu = 42;
  const int columnNu = 5;
  int _yID = ScintillatorIDs / rowNu;
  int _xID = ScintillatorIDs % rowNu;
  TVector3 _position;
  double x0 = _xInterval * _xID - _xInterval * (rowNu - 1) / 2.;
  double y0 = _yInterval * _yID - _yInterval * (columnNu - 1) / 2.;

  // for prototype test
  if (LayerIDs % 2 == 0) {
    _position[0] = -y0;
    _position[1] = -x0;
  }
  if (LayerIDs % 2 == 1) {
    _position[0] = -x0;
    _position[1] = -y0;
  }
  if (LayerIDs % 2 == 0) layerZ = 1 + LayerIDs / 2 * 19.9;
  // else layerZ = 12.95+(LayerIDs-1)/2*19.9;
  else
    layerZ = 12.2 + (LayerIDs - 1) / 2 * 19.9;
  _position[2] = layerZ;
  // for grouped CR test
  // if(LayerIDs%4==0) {
  //     _position[0] = -y0;
  //     _position[1] = -x0;
  // }
  // if(LayerIDs%4==1) {
  //     _position[0] = x0;
  //     _position[1] = y0;
  // }
  // if(LayerIDs%4==2) {
  //     _position[0] = y0;
  //     _position[1] = x0;
  // }
  // if(LayerIDs%4==3) {
  //     _position[0] = -x0;
  //     _position[1] = -y0;
  // }
  // _position[2] = layerZ[LayerIDs];

  return _position;
}

//------------------ For AHCAL ------------------
// Copy from Yukun @ USTC
// scintillator strips wrt. 6 SPIROC2E chips * 36 channels
TVector3 GetAhcalHitPos(int layer_ID, int chip_ID, int channel_ID) {
  const int cell_SP = 16;
  const int chip_No = 9;
  const int channel_No = 36;
  const int Layer_No = 40;
  const double _Pos_X[channel_No] = {
      100.2411,  100.2411,  100.2411,  59.94146,  59.94146,  59.94146,
      19.64182,  19.64182,  19.64182,  19.64182,  59.94146,  100.2411,
      100.2411,  59.94146,  19.64182,  100.2411,  59.94146,  19.64182,
      -20.65782, -60.95746, -101.2571, -20.65782, -60.95746, -101.2571,
      -101.2571, -60.95746, -20.65782, -20.65782, -20.65782, -20.65782,
      -60.95746, -60.95746, -60.95746, -101.2571, -101.2571, -101.2571};
  const double _Pos_Y[channel_No] = {
      141.04874, 181.34838, 221.64802, 141.04874, 181.34838, 221.64802,
      141.04874, 181.34838, 221.64802, 261.94766, 261.94766, 261.94766,
      302.2473,  302.2473,  302.2473,  342.54694, 342.54694, 342.54694,
      342.54694, 342.54694, 342.54694, 302.2473,  302.2473,  302.2473,
      261.94766, 261.94766, 261.94766, 221.64802, 181.34838, 141.04874,
      221.64802, 181.34838, 141.04874, 221.64802, 181.34838, 141.04874};
  const double chip_dis_X = 239.3;
  const double chip_dis_Y = 241.8;
  const double HBU_X = 239.3;
  const double HBU_Y = 725.4;
  const double HBU_Z = 26;

  TVector3 pos;
  int HBU_ID = chip_ID / 3;
  chip_ID = chip_ID % 3;
  pos.SetX(_Pos_Y[channel_ID] - chip_ID * chip_dis_Y);
  pos.SetY(-(-_Pos_X[channel_ID] + (HBU_ID - 1) * HBU_X));
  pos.SetZ(layer_ID * HBU_Z);
  return pos;
}
