#include "Tools.cxx"

void JetRate(){
 
 TFile* file = new TFile("latino_stepB_latinosYieldSkim_MC_ggHww.root","READ");
 
 TTree* latino = (TTree*) file->Get("latino");

 std::vector<float> pt_edges;
 
//  pt_edges.push_back(5.0);
//  //  pt_edges.push_back(10.0);
//  pt_edges.push_back(15.0);
//  pt_edges.push_back(20.0);
//  pt_edges.push_back(25.0);
//  pt_edges.push_back(30.0);
//  pt_edges.push_back(35.0);
//  pt_edges.push_back(50.0);
//  pt_edges.push_back(100.0);
//  pt_edges.push_back(200.0);
//  pt_edges.push_back(300.0);
//  
 
 pt_edges.push_back(5.0);
 pt_edges.push_back(20.0);
 pt_edges.push_back(30.0);
 pt_edges.push_back(40.0);
 pt_edges.push_back(50.0);
 pt_edges.push_back(100.0);
 pt_edges.push_back(150.0);
 pt_edges.push_back(200.0);
 pt_edges.push_back(250.0);
 pt_edges.push_back(300.0);
 
 
  
 Float_t bins[1000];
 for (int ibin = 0; ibin < pt_edges.size(); ibin++) {
  bins[ibin] = pt_edges.at(ibin);
 }
 bins[pt_edges.size()] = pt_edges.at(pt_edges.size()-1) + 100;
 
 TH1F* histo_efficiency_standard = new TH1F("histo_efficiency_standard","efficiency standard", pt_edges.size(), bins);
 TH1F* histo_efficiency_standard_num = new TH1F("histo_efficiency_standard_num","efficiency standard", pt_edges.size(), bins);
 TH1F* histo_efficiency_standard_den = new TH1F("histo_efficiency_standard_den","efficiency standard", pt_edges.size(), bins);
 
 TH1F* histo_efficiency_puppi = new TH1F("histo_efficiency_puppi","efficiency puppi", pt_edges.size(), bins);
 TH1F* histo_efficiency_puppi_num = new TH1F("histo_efficiency_puppi_num","efficiency puppi", pt_edges.size(), bins);
 TH1F* histo_efficiency_puppi_den = new TH1F("histo_efficiency_puppi_den","efficiency puppi", pt_edges.size(), bins);
 
 TH1F* histo_efficiency_standard_over_puppi = new TH1F("histo_efficiency_standard_over_puppi","efficiency standard / puppi", pt_edges.size(), bins);
 
 
 TH1F* histo_fake_standard = new TH1F("histo_fake_standard","fake standard", pt_edges.size(), bins);
 TH1F* histo_fake_puppi = new TH1F("histo_fake_puppi","fake puppi", pt_edges.size(), bins);
 TH1F* histo_fake_standard_over_puppi = new TH1F("histo_fake_standard_over_puppi","fake standard / puppi", pt_edges.size(), bins);
 
 
 
 std::vector<float>* std_vector_jetGen_pt = 0; 
 std::vector<float>* std_vector_jetGen_eta = 0;
 std::vector<float>* std_vector_jetGen_phi = 0;

 std::vector<float>* std_vector_jet_pt = 0; 
 std::vector<float>* std_vector_jet_eta = 0;
 std::vector<float>* std_vector_jet_phi = 0;

 std::vector<float>* std_vector_puppijet_pt = 0; 
 std::vector<float>* std_vector_puppijet_eta = 0;
 std::vector<float>* std_vector_puppijet_phi = 0;

 
 latino->SetBranchAddress("std_vector_jetGen_pt" ,&std_vector_jetGen_pt);
 latino->SetBranchAddress("std_vector_jetGen_phi",&std_vector_jetGen_phi);
 latino->SetBranchAddress("std_vector_jetGen_eta",&std_vector_jetGen_eta);

 latino->SetBranchAddress("std_vector_jet_pt" ,&std_vector_jet_pt);
 latino->SetBranchAddress("std_vector_jet_phi",&std_vector_jet_phi);
 latino->SetBranchAddress("std_vector_jet_eta",&std_vector_jet_eta);

 latino->SetBranchAddress("std_vector_puppijet_pt" ,&std_vector_puppijet_pt);
 latino->SetBranchAddress("std_vector_puppijet_phi",&std_vector_puppijet_phi);
 latino->SetBranchAddress("std_vector_puppijet_eta",&std_vector_puppijet_eta);
 
 
 Long64_t nentries = latino->GetEntries();
  for (Long64_t i=0; i<nentries; i++) {
   latino->GetEntry(i);
   
   //---- match efficiency
   for (UInt_t j = 0; j < std_vector_jetGen_phi->size(); ++j) {
    float pt_gen = std_vector_jetGen_pt->at(j);
    if (pt_gen > pt_edges.at(pt_edges.size()-1)) { //---- to deal with overflow bin in "view" zone
     pt_gen = pt_edges.at(pt_edges.size()-1)+0.5;
    }
    
    histo_efficiency_standard_den->Fill(pt_gen);
    histo_efficiency_puppi_den->Fill(pt_gen);
    
    std::pair<int, float> closest_standard_jet = getClosestIndexAndDR(std_vector_jetGen_eta->at(j), std_vector_jetGen_phi->at(j),    *std_vector_jet_eta, *std_vector_jet_phi);
    if (closest_standard_jet.second < 0.4) {
     histo_efficiency_standard_num->Fill(pt_gen);
    }
    std::pair<int, float> closest_puppi_jet = getClosestIndexAndDR(std_vector_jetGen_eta->at(j), std_vector_jetGen_phi->at(j),    *std_vector_puppijet_eta, *std_vector_puppijet_phi);
    if (closest_puppi_jet.second < 0.4) {
     histo_efficiency_puppi_num->Fill(pt_gen);
    }
    
   }
   
   
   //---- fake
   for (UInt_t j = 0; j < std_vector_jet_phi->size(); ++j) {
    float pt_gen = std_vector_jet_pt->at(j);
    if (pt_gen > pt_edges.at(pt_edges.size()-1)) { //---- to deal with overflow bin in "view" zone
     pt_gen = pt_edges.at(pt_edges.size()-1)+0.5;
//      std::cout << " pt_gen = " << pt_gen << std::endl;
    }
    
    std::pair<int, float> closest_standard_jet = getClosestIndexAndDR(std_vector_jet_eta->at(j), std_vector_jet_phi->at(j),    *std_vector_jetGen_eta, *std_vector_jetGen_phi);
    if (closest_standard_jet.second > 0.4) {
     histo_fake_standard->Fill(pt_gen);
    }    
//     else {
//      std::cout << " closest_standard_jet.second = " << closest_standard_jet.second << std::endl;
//     }
   }
   
   
   for (UInt_t j = 0; j < std_vector_puppijet_phi->size(); ++j) {
    float pt_gen = std_vector_puppijet_pt->at(j);
    if (pt_gen > pt_edges.at(pt_edges.size()-1)) { //---- to deal with overflow bin in "view" zone
     pt_gen = pt_edges.at(pt_edges.size()-1)+0.5;
    }
    
    std::pair<int, float> closest_puppi_jet = getClosestIndexAndDR(std_vector_puppijet_eta->at(j), std_vector_puppijet_phi->at(j),    *std_vector_jetGen_eta, *std_vector_jetGen_phi);
    if (closest_puppi_jet.second > 0.4) {
     histo_fake_puppi->Fill(pt_gen);
    }    
   }
   
    
    
  }
 
 
 for (int ibin = 0; ibin <= pt_edges.size(); ibin++) {
  float num = histo_efficiency_standard_num->GetBinContent(ibin+1);
  float den = histo_efficiency_standard_den->GetBinContent(ibin+1);
  float eff = 0.;
  if (den != 0) eff = num / den;
  histo_efficiency_standard->SetBinContent(ibin+1, eff);
 }
 
 histo_efficiency_standard->SetLineColor(kBlue);
 histo_efficiency_standard->SetLineWidth(2);
 

 for (int ibin = 0; ibin <= pt_edges.size(); ibin++) {
  float num = histo_efficiency_puppi_num->GetBinContent(ibin+1);
  float den = histo_efficiency_puppi_den->GetBinContent(ibin+1);
  float eff = 0.;
  if (den != 0) eff = num / den;
  histo_efficiency_puppi->SetBinContent(ibin+1, eff);
 }
 
 histo_efficiency_puppi->SetLineColor(kRed);
 histo_efficiency_puppi->SetLineWidth(2);
 histo_efficiency_puppi->SetLineStyle(2);
 
 
 
 TCanvas* ccResult_efficiency = new TCanvas ("ccResult_efficiency","efficiency",800,800);
 ccResult_efficiency->Divide(1,2);
 
 ccResult_efficiency->cd(1);
 histo_efficiency_standard->Draw();
 histo_efficiency_puppi->Draw("same");
 histo_efficiency_standard->GetXaxis()->SetTitle("gen jet p_{T} [GeV]");
 histo_efficiency_standard->GetYaxis()->SetTitle("matching efficiency");
 histo_efficiency_standard->GetYaxis()->SetRangeUser(0.0, 1.0);
 gPad->SetGrid();
 gPad->BuildLegend();
 
 ccResult_efficiency->cd(2);
 
 for (int ibin = 0; ibin <= pt_edges.size(); ibin++) {
  float num = histo_efficiency_standard->GetBinContent(ibin+1);
  float den = histo_efficiency_puppi->GetBinContent(ibin+1);
  float eff = 0.;
  if (den != 0) eff = num / den;
  histo_efficiency_standard_over_puppi->SetBinContent(ibin+1, eff);
 }
 histo_efficiency_standard_over_puppi->SetLineColor(kGreen+2);
 histo_efficiency_standard_over_puppi->SetLineWidth(2);
 histo_efficiency_standard_over_puppi->Draw();
 histo_efficiency_standard_over_puppi->GetXaxis()->SetTitle("gen jet p_{T} [GeV]");
 histo_efficiency_standard_over_puppi->GetYaxis()->SetTitle("standard / puppi efficiency");
 histo_efficiency_standard_over_puppi->GetYaxis()->SetRangeUser(0.0, 2.0);
 
 gPad->SetGrid();
 
 
 
 
 
 
 TCanvas* ccResult_fake = new TCanvas ("ccResult_fake","fake",600,800);
 ccResult_fake->Divide(1,2);
 
 ccResult_fake->cd(1);
 histo_fake_standard->SetLineColor(kBlue);
 histo_fake_standard->SetLineWidth(2);

 histo_fake_puppi->SetLineColor(kRed);
 histo_fake_puppi->SetLineWidth(2);
 histo_fake_puppi->SetLineStyle(2);
 
 histo_fake_standard->Draw();
 histo_fake_puppi->Draw("same");
 histo_fake_standard->GetXaxis()->SetTitle("reco jet p_{T} [GeV]");
 histo_fake_standard->GetYaxis()->SetTitle("fake (a.u.)");
 gPad->SetGrid();
 gPad->SetLogy();
 gPad->BuildLegend();
 
 ccResult_fake->cd(2);
 
 for (int ibin = 0; ibin <= pt_edges.size(); ibin++) {
  float num = histo_fake_standard->GetBinContent(ibin+1);
  float den = histo_fake_puppi->GetBinContent(ibin+1);
  float eff = 0.;
  if (den != 0) eff = num / den;
  histo_fake_standard_over_puppi->SetBinContent(ibin+1, eff);
 }
 histo_fake_standard_over_puppi->SetLineColor(kGreen+2);
 histo_fake_standard_over_puppi->SetLineWidth(2);
 histo_fake_standard_over_puppi->Draw();
 histo_fake_standard_over_puppi->GetXaxis()->SetTitle("reco jet p_{T} [GeV]");
 histo_fake_standard_over_puppi->GetYaxis()->SetTitle("standard / puppi fake");
 histo_fake_standard_over_puppi->GetYaxis()->SetRangeUser(0.0, 10.0);
 
 gPad->SetGrid();
 
}


