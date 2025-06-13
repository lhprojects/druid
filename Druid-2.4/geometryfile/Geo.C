
{
//TGeoManager::Import("sidloi3.gdml");
//TGeoManager::Import("World.gdml");
//TGeoManager::Import("MokkaGDML/ILD_00Dhcal.gdml");
//TGeoManager::Import("tb_angela.gdml");
TGeoManager::Import("ShortTPC.gdml");
//TGeoManager::Import("clic_sid_cdr_b.gdml");
gGeoManager->GetTopVolume()->Draw("ogl");
TFile *f = new TFile("ShortTPC.root","recreate");
gGeoManager->Write();
f->Close();
}
