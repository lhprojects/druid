
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//    Build Simulated CaloHits into TEveBox visualization                   //
//    Handles SimCalorimeterHit collection processing and visualization      //
//                                                                           //
//    Date:   27 Nov 2024                                                    //
//    Extracted from CaloHits.cc                                             //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "EVENT/CalorimeterHit.h"
#include "EVENT/LCCollection.h"
#include "EVENT/SimCalorimeterHit.h"
#include "GlobalDefs.hh"
#include "TColor.h"
#include "TEveBox.h"
#include "TEveBoxSet.h"
#include "TEveRGBAPalette.h"
#include "TEveManager.h"
#include "TStyle.h"
#include "TVector3.h"
#include "UTIL/CellIDDecoder.h"
#include "lcio.h"
#include "TruthHelper.h"

using namespace lcio;
using namespace EVENT;
using namespace std;

// External declarations
extern double cellEnergyThresh;
extern float HitColourFactor;
extern bool Flag_AttachTextToHit;
extern int GlobalRandomColorIndex;
extern float cellfactor;
extern float HitEnergyCut;
extern bool HiddenLowESimCell;
extern int flagdetectortype;
extern std::map<float, int> randomColor;
extern std::map<MCParticle*, int> OriginColor;
extern TEveRGBAPalette* p;
extern TEveBoxSet* cal_shell;

extern const float MIPSIMECAL;
extern const float MIPSIMAHCAL;
extern const float MIPSIMDHCAL;

// Forward declaration
TEveBox *createSimCaloHitBox(SimCalorimeterHit *hit, float MIPS, int type, int caloType);

TEveElementList *createSimCaloHits(LCCollection *col, string name)
{
    bool isSimHit = true;

    int collIndex = gGUIManager.m_SimCaloHitCollNames.size();
    gGUIManager.m_SimCaloHitCollNames.push_back(name);

    // Clear SimCalorimeterHit boxes when processing SimCalorimeterHit collection

    TEveElementList *cal = new TEveElementList;

    cal_shell->Reset(TEveBoxSet::kBT_FreeBox, kFALSE, 64);
    cal_shell->SetMainAlpha(1);

    cal_shell->SetPalette(p);
    cal_shell->SetPickable(false);

    int icol(0);
    int mothercol(0);

    float HitX, HitY, HitZ, HitEn;
    char *HitRegion = 0;
    int MCPID, StaveNum, IndexI,
        IndexJ;    // IndexI and StaveNum, as well as Module Number and IndexK are
                   // coded in ID0 of Hit...
    int MCHitType; // EM, Had, Neutron, Proton;
    float MCenergy;
    int MotherPID;
    int LayerNum;
    float Mipcount;     // Meause Energy deposition in unit of MIP
    float MotherEnergy; // The energy and PID of Origin
    float MCContTime;
    // OriginColor.clear();

    double MaxHitEnergy = 0;
    double MinHitEnergy = 1000000; // GeV

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
    TEveRGBAPalette *pal = new TEveRGBAPalette(0, 100);

    float HalfX = 0;
    float HalfY = 0;
    float HalfZ = 0;
    float s1, s2;
    float SF = isSimHit ? 1.0 : 1.2; // Scale Factor
    float s1X, s1Y, s2X, s2Y;

    cout << "  Calo hits collection : " << name << ". Number of Hits: " << nHits
         << ", typeName:" << col->getTypeName() << endl;
    // cout << endl; Hao

    int MaxEnergyDepoID = 0; // To record the MCPID of Simulated Calo Hits
    int MCPIDindex = 1;
    int Originindex = 0;
    int MotherEnergySignature = 0;

    // CellIDDecoder<SimCalorimeterHit> idDecoder(col);

    for (int i = 0; i < nHits; i++)
    {
        //------------------ Set Hit Energy & CellIDs ------------------
        if (true)
        {
            SimCalorimeterHit *hit =
                dynamic_cast<SimCalorimeterHit *>(col->getElementAt(i));

            CellIDDecoder<SimCalorimeterHit> idDecoder(col);
            // IndexI = idDecoder(hit)["I"];
            // IndexJ = idDecoder(hit)["J"];
            // if (SubD == "Ecal" || SubD == "Hcal") LayerNum = idDecoder(hit)["K-1"];
            IndexI = 0;
            IndexJ = 0;


            if (hit->getNMCContributions() > 0 && name != "LumiCalCollection")
            {
                float Emax = -0.01;
                MaxEnergyDepoID = 0;
                for (int k = 0; k < hit->getNMCContributions(); k++)
                {
                    MCParticle *hitMCPart =
                        dynamic_cast<MCParticle *>(hit->getParticleCont(k));

                    if (hit->getEnergyCont(k) > Emax)
                    {
                        MCPID = hitMCPart->getPDG(); // According to assigned Track index:
                                                     // lots of mother contribution
                        MCHitType = hit->getPDGCont(
                            k); // According to local deposited track PID: e-, e+, pi, K...
                        MCenergy = hitMCPart->getEnergy();
                        Emax = hit->getEnergyCont(k);
                        MaxEnergyDepoID = k;
                        MCContTime = hit->getTimeCont(k);

                        if (randomColor.find(MCenergy) != randomColor.end())
                        {
                            MCPIDindex = randomColor[MCenergy];
                        }
                        else
                        {
                            randomColor[MCenergy] = icol++;
                            MCPIDindex = randomColor[MCenergy];
                        }

                        // MCPIDindex = MCPID;				//Color
                        // according to the ID
                    }
                }

                // Try to get Mother

                if (gGUIManager.HitColourType == 2)
                {
                    MCParticle *hitMCPart0 =
                        dynamic_cast<MCParticle *>(hit->getParticleCont(MaxEnergyDepoID));

                    for (int s = 0; s < 100; s++)
                    {
                        MCParticleVec mothers = hitMCPart0->getParents();
                        if (mothers.size() > 0)
                        {
                            MCParticle *Mother = mothers[0];

                            if (Mother->getPDG() == 92 || Mother->getParents().size() == 0 ||
                                Mother->getPDG() == 23 ||
                                Mother->getPDG() ==
                                    25) // || abs(Mother->getPDG())<6 ) 5, -5, 25;
                            // if( (fabs(Mother->getPDG()) < 6 || Mother->getPDG() == 25) &&
                            // Mother->getParents().size() == 0 )
                            {
                                MotherPID = Mother->getPDG();
                                MotherEnergy = Mother->getEnergy();

                                if (OriginColor.find(Mother) == OriginColor.end())
                                {
                                    if (MotherPID == 23)
                                        Originindex = 2;
                                    else if (MotherPID == 25)
                                        Originindex = 4;
                                    else
                                        Originindex = 9;
                                    OriginColor[Mother] = Originindex;
                                }
                                else
                                {
                                    Originindex = OriginColor[Mother];
                                }
                                break;
                            }
                            else
                            {
                                hitMCPart0 = Mother;
                            }
                        }
                    }
                    // cout<<"Mother "<<MCPID<<"\t"<<MotherPID<<"\t"<<MotherEnergy<<endl;
                }
            }
        }


        if (SubD == "Muon" || SubD == "MUON")
        {
            SF = 10.0;
        }
        else
        {
            //      if (HitSizeLog == 0) {
            if (((SubD == "Hcal" || SubD == "HCAL") && flagdetectortype == 0) ||
                SubD == "Ahcal")
            {
                SF = 3.0 * cellfactor; // AHCAL Hits
            }
            else
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

        bool isECAL = contains(name, "Ecal") || contains(name, "ECAL") ||
                        contains(name, "ScEcal") || contains(name, "SCECAL");
        bool isHCAL = contains(name, "Hcal") || contains(name, "HCAL") ||
                        contains(name, "Ahcal") || contains(name, "AHCAL");
        bool isBarrel = contains(name, "Barrel");
        bool isEndcap = contains(name, "Endcap");

    
        SimCalorimeterHit *simCaloHit = dynamic_cast<SimCalorimeterHit *>(col->getElementAt(i));
        float const hitEnergy = simCaloHit->getEnergy();
        float const MIP = isECAL ? MIPSIMECAL : (isHCAL ? MIPSIMAHCAL : MIPSIMDHCAL);
        Mipcount = hitEnergy / MIP;
        
        // Extract hit position and energy for later use
        HitX = simCaloHit->getPosition()[0] / 10.0;  // Convert to cm
        HitY = simCaloHit->getPosition()[1] / 10.0;
        HitZ = simCaloHit->getPosition()[2] / 10.0;
        HitEn = hitEnergy;
        
        //std::cout << "Hit " << i << ": Energy = " << hitEnergy << " GeV, Mipcount = " << Mipcount << std::endl;

        if (gGUIManager.HitColourType == 1)
        { 
            if (isSimHit && Mipcount < HitEnergyCut)
            {
                IndexHitColorLine = 10;
                IndexHitColorMain = 10;
            } // White to indicate the Hitcolor for simulated Hits: with changing
              // bkgrd to white it will disappear.
            else
            {
                IndexHitColorLine = int(HitColourFactor * Mipcount * 0.5 + 51);
                IndexHitColorMain = int(HitColourFactor * Mipcount * 0.5 + 51);
            }

            if (IndexHitColorLine > 100)
                IndexHitColorLine = 100;
            if (IndexHitColorMain > 100)
                IndexHitColorMain = 100;

            if ((SubD == "Hcal" || SubD == "HCAL") && cellEnergyThresh != 0)
            {
                IndexHitColorMain = 10; // to denote SemiDHCAL case

                if (Mipcount * cellEnergyThresh < 0.2)
                {
                    IndexHitColorLine = 51;
                    IndexHitColorMain = 51;
                }
                else if (Mipcount * cellEnergyThresh < 5 &&
                         Mipcount * cellEnergyThresh >= 0.2)
                {
                    IndexHitColorLine = 4;
                    IndexHitColorMain = 4;
                }
                else if (Mipcount * cellEnergyThresh < 10 &&
                         Mipcount * cellEnergyThresh >= 5)
                {
                    IndexHitColorLine = 3;
                    IndexHitColorMain = 3;
                }
                else if (Mipcount * cellEnergyThresh >= 10)
                {
                    IndexHitColorLine = 2;
                    IndexHitColorMain = 2;
                }
            }
            else if ((SubD == "SDHC")) // From Trivent
            {
                if (HitEn == 1) // Inverted ... 1 = 2
                {
                    IndexHitColorLine = 3;
                    IndexHitColorMain = 3;
                }
                else if (HitEn == 2)
                {
                    IndexHitColorLine = 4;
                    IndexHitColorMain = 4;
                }
                else if (HitEn == 3)
                {
                    IndexHitColorLine = 2;
                    IndexHitColorMain = 2;
                }
                else
                {
                    IndexHitColorLine = 5;
                    IndexHitColorMain = 5;
                }
            }
        }
        else if (gGUIManager.HitColourType == 0)
        { // rep direct deposition Particle Type with Color

            switch (MCPID)
            {
            case (12): // Neutrinos: Normally doesn't creates hits...
            case (14):
            case (16):
            case (-12):
            case (-14):
            case (-16):
                IndexHitColorMain = 66;
                break;

            case (22): // Gamma!
                IndexHitColorMain = 89;
                break;

            case (11): // electrons & positrons
                IndexHitColorMain = 56;
                break;

            case (-11):
                IndexHitColorMain = 99;
                break;

            case (211): // Pions
                IndexHitColorMain = 97;
                break;

            case (-211):
                IndexHitColorMain = 53;
                break;

            case (13): // Muons
                IndexHitColorMain = 100;
                break;

            case (-13):
                IndexHitColorMain = 51;
                break;

            case (2212): // Protons
                IndexHitColorMain = 96;
                break;

            case (-2212):
                IndexHitColorMain = 64;
                break;

            case (2122): // Neutrons
                IndexHitColorMain = 85;
                break;

            case (130): // Klong
                IndexHitColorMain = 80;
                break;

            case (321): // Charged Kaon
                IndexHitColorMain = 98;
                break;

            case (-321):
                IndexHitColorMain = 54;
                break;

            default: // Any other possibles
                IndexHitColorMain = 38;
                break;
            }
            IndexHitColorLine = IndexHitColorMain;
        }
        else if (gGUIManager.HitColourType == 2)
        {
            // IndexHitColorMain = (Originindex*13 + GlobalRandomColorIndex*7)%50+51;
            IndexHitColorMain = Originindex;
        }
        else if (gGUIManager.HitColourType == 3)
        {
            IndexHitColorMain = ((MCPIDindex % 2) * 25 + MCPIDindex * 5 +
                                 GlobalRandomColorIndex * 13) %
                                    50 +
                                51;
        }
        else if (gGUIManager.HitColourType == 4) // Uniformed Color for Simulated Hit: Blue
        {
            IndexHitColorMain = 5;
        }
        else if (gGUIManager.HitColourType == 5)
        {
            if (MCHitType == 11 || MCHitType == -11 || MCHitType == 22)
            {
                IndexHitColorMain = 4; // Blue for EM
            }
            else if (MCHitType == 2112) // Neutron
            {
                IndexHitColorMain = 5; // Yellow
            }
            else
            {
                IndexHitColorMain = 2;
            }
        }
        else if (gGUIManager.HitColourType == 6)
        {
            if (MCContTime == -10)
            {
                std::cout << "Timing info not available in current data file"
                          << std::endl;
                IndexHitColorMain = 3;
            }
            else
            {
                if (MCContTime < 150) // 150ns as electronic integration time...
                {
                    IndexHitColorMain =
                        int(0.5 * HitColourFactor * (MCContTime - 5.0) + 51);
                }
                else
                {
                    IndexHitColorMain = 21;
                }
            }
            if (IndexHitColorMain > 100)
                IndexHitColorMain = 100;
            IndexHitColorLine = IndexHitColorMain;
        }
        else if (gGUIManager.HitColourType == 7)
        {
            // Color scheme 7: Use MCParticle colors (handled by UpdateHitColorsToMatchMCParticles)
            // Set initial color to gray, will be updated by UpdateHitColorsToMatchMCParticles()
            IndexHitColorMain = kGray;
            IndexHitColorLine = IndexHitColorMain;
        }

        if (Mipcount > HitEnergyCut || !HiddenLowESimCell || !isSimHit)
        {
            TEveBox *q = new TEveBox;

            TVector3 HitScale(SF, SF, 0.1 * SF);
            TVector3 TCMTScale(100.0, 2.0, 5.0);



            if (IndexHitColorLine == 0)
                IndexHitColorLine = IndexHitColorMain;

            // type: 0 = barrel, 1 = endcap
            int hitType = isBarrel ? 0 : (isEndcap ? 1 : 0);
            // caloType: 0 = ECAL, 1 = HCAL
            int caloType = isECAL ? 0 : (isHCAL ? 1 : 0);
            
            q = createSimCaloHitBox(simCaloHit, Mipcount, hitType, caloType);

            q->SetPickable(kTRUE);
            q->SetMainColor(IndexHitColorMain);
            q->SetLineColor(IndexHitColorLine);
            q->SetMainAlpha(0.8);

            if (Flag_AttachTextToHit)
            {
                EVENT::MCParticle *mcp = gTruthHelper.GetMainMCP(simCaloHit);
                q->SetTitle(Form(
                    "SimCaloHit %s\n"
                    "EventNum=%d, SubDetector=%s\n"
                    "Hit Energy = %.3e keV ~ %.3e Mip, Thresh = %3.e keV\n"
                    "MCP %s, MCPID = %d, MCenergy = %.3f\n "
                    "OriginPID = %d, OriginEnergy = %.3f\n"
                    "PosX = %.3f mm, PosY = %.3f mm, PosZ = %.3f mm,\n"
                    "HitTime = %.3fns",
                    name.c_str(),
                    gDisplayState.getEventNumber(),
                    SubDCollection.c_str(),
                    HitEn * 1000000, Mipcount, cellEnergyThresh * 1000000,
                    gTruthHelper.GetStringID(mcp).c_str(), MCPID, MCenergy,
                    MotherPID,
                    MotherEnergy,
                    10 * HitX, 10 * HitY, 10 * HitZ, MCContTime));
            }

            // Store SimCalorimeterHit box for later updates
            SimCalorimeterHit *hit =
                dynamic_cast<SimCalorimeterHit *>(col->getElementAt(i));

            gGUIManager._SimCaloHitBoxes[q] = hit;
            gGUIManager._SimCaloHitType[q] = collIndex;

            cal->AddElement(q);
        }
    }

    cout << "OriCoSize " << OriginColor.size() << ", " << mothercol << endl;

    return cal;
}

// Function to update SimCalorimeterHit and CalorimeterHit colors to match MCParticle visibility and colors
void UpdateHitColorsToMatchMCParticles() {
    std::cout << "Updating SimCalorimeterHit and CalorimeterHit colors to match MCParticle visibility..." << std::endl;
    
    int simUpdatedCount = 0;
    int simGrayCount = 0;
    int caloUpdatedCount = 0;
    int caloGrayCount = 0;

    for(std::string name : gGUIManager.m_SimCaloHitCollNames) {
        std::cout << "Current SimCaloHit collection name: " << name << std::endl;
    }
    
    for(std::string name : gGUIManager.m_CaloHitCollNames) {
        std::cout << "Current CaloHit collection name: " << name << std::endl;
    }
    
    // Iterate through all stored SimCalorimeterHit boxes
    for (auto& hitBox : gGUIManager._SimCaloHitBoxes) {
        EVENT::SimCalorimeterHit* hit = hitBox.second;
        TEveBox* box = hitBox.first;
        
        //std::cout << "Processing hit: " << hit << ", box: " << box << std::endl;
        if (!hit || !box) continue;
        
        // Get the main MCParticle for this hit
        EVENT::MCParticle* mainMCP = gTruthHelper.GetMainMCP(hit);
        bool mcpIsVisible = false;
        int newColor = kGray;

        //std::cout << "Main MCP for hit " << hit << ": " << mainMCP << std::endl;
        std::string name = gGUIManager.m_SimCaloHitCollNames.at(gGUIManager._SimCaloHitType[box]);
        
        if (mainMCP && gGUIManager._MCPTracks.find(mainMCP) != gGUIManager._MCPTracks.end()) {
            TEveTrack* track = gGUIManager._MCPTracks[mainMCP];
            
            // Check if the track is visible
            bool trackVisible = track && track->GetRnrSelf();
            
            // Check if the particle's group is showing
            bool groupVisible = true;
            auto groupIt = gGUIManager._MCPGroups.find(mainMCP);
            if (groupIt != gGUIManager._MCPGroups.end()) {
                TEveElement* group = groupIt->second;
                if (group) {
                    groupVisible = group->GetRnrSelf() && group->GetRnrChildren();
                    if(false) {
                        std::cout << "group self visible: " << group->GetRnrSelf()
                                << ", children visible: " << group->GetRnrChildren() << std::endl;
                    }
                }
            }
            
            if (trackVisible && groupVisible) {
                mcpIsVisible = true;
                // Use the same color as the MCParticle track
                newColor = track->GetLineColor();                
                //std::cout << "MCP is visible. Using track color: " << newColor << std::endl;
            }
        }
        
        if (!mcpIsVisible) {
            newColor = kGray;
            simGrayCount++;
        } else {
            simUpdatedCount++;
        }
        
        // Update the box color
        box->SetMainColor(newColor);
        box->SetMainAlpha(0.8);
        box->SetLineColor(newColor);
    }
    
    // Iterate through all stored CalorimeterHit boxes
    for (auto& hitBox : gGUIManager._CaloHitBoxes) {
        EVENT::CalorimeterHit* hit = hitBox.second;
        TEveBox* box = hitBox.first;
        
        if (!hit || !box) continue;
        
        // Get the main MCParticle for this hit
        EVENT::MCParticle* mainMCP = gTruthHelper.GetMainMCP(hit);
        bool mcpIsVisible = false;
        int newColor = kGray;
        
        if (mainMCP && gGUIManager._MCPTracks.find(mainMCP) != gGUIManager._MCPTracks.end()) {
            TEveTrack* track = gGUIManager._MCPTracks[mainMCP];
            
            // Check if the track is visible
            bool trackVisible = track && track->GetRnrSelf();
            
            // Check if the particle's group is showing
            bool groupVisible = true;
            auto groupIt = gGUIManager._MCPGroups.find(mainMCP);
            if (groupIt != gGUIManager._MCPGroups.end()) {
                TEveElement* group = groupIt->second;
                if (group) {
                    groupVisible = group->GetRnrSelf() && group->GetRnrChildren();
                }
            }
            
            if (trackVisible && groupVisible) {
                mcpIsVisible = true;
                // Use the same color as the MCParticle track
                newColor = track->GetLineColor();
            }
        }
        
        if (!mcpIsVisible) {
            newColor = kGray;
            caloGrayCount++;
        } else {
            caloUpdatedCount++;
        }
        
        // Update the box color
        box->SetMainColor(newColor);
        box->SetMainAlpha(0.8);
        box->SetLineColor(newColor);
    }
    
    std::cout << "SimCalorimeterHit: Updated " << simUpdatedCount << " hits to match MCParticle colors, " 
              << simGrayCount << " hits set to gray" << std::endl;
    std::cout << "CalorimeterHit: Updated " << caloUpdatedCount << " hits to match MCParticle colors, " 
              << caloGrayCount << " hits set to gray" << std::endl;
    
    // Redraw without resetting camera view
    if (gEve) {
        gEve->Redraw3D(kFALSE);
    }
}
