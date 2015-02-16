#include "Tools.cxx"

void JetRateEta(){
 
 TFile* file = new TFile("latino_stepB_latinosYieldSkim_MC_ggHww.root","READ");
 
 TTree* latino = (TTree*) file->Get("latino");
 
 
 std::vector<float> eta_edges;
 
 eta_edges.push_back(0.0);
 eta_edges.push_back(1.3);
 eta_edges.push_back(2.5);
 eta_edges.push_back(3.0);
 eta_edges.push_back(4.5);
 
 
 
 
 std::vector<float> pt_edges;
 
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
 
 
 
 TH1F* histo_efficiency_standard[100];
 TH1F* histo_efficiency_standard_num[100];
 TH1F* histo_efficiency_standard_den[100];
 
 TH1F* histo_efficiency_puppi[100];
 TH1F* histo_efficiency_puppi_num[100];
 TH1F* histo_efficiency_puppi_den[100];
 
 TH1F* histo_efficiency_standard_over_puppi[100];
 
 TH1F* histo_fake_standard[100];
 TH1F* histo_fake_puppi[100];
 TH1F* histo_fake_standard_over_puppi[100];
 
 
 TString name;
 
 for (int ibineta = 0; ibineta < (eta_edges.size()-1); ibineta++) {
  
  name = Form ("histo_efficiency_standard_%d",ibineta);
  histo_efficiency_standard[ibineta] = new TH1F (name.Data(),"standard",pt_edges.size(), bins);   
  name = Form ("histo_efficiency_standard_num_%d",ibineta);
  histo_efficiency_standard_num[ibineta] = new TH1F (name.Data(),"standard",pt_edges.size(), bins);
  name = Form ("histo_efficiency_standard_den_%d",ibineta);
  histo_efficiency_standard_den[ibineta] = new TH1F (name.Data(),"standard",pt_edges.size(), bins);
  
  
  name = Form ("histo_efficiency_puppi_%d",ibineta);
  histo_efficiency_puppi[ibineta] = new TH1F (name.Data(),"puppi",pt_edges.size(), bins);
  name = Form ("histo_efficiency_puppi_num_%d",ibineta);
  histo_efficiency_puppi_num[ibineta] = new TH1F (name.Data(),"puppi",pt_edges.size(), bins);
  name = Form ("histo_efficiency_puppi_den_%d",ibineta);
  histo_efficiency_puppi_den[ibineta] = new TH1F (name.Data(),"puppi",pt_edges.size(), bins);
  
  
  name = Form ("histo_efficiency_standard_over_puppi_%d",ibineta);
  histo_efficiency_standard_over_puppi[ibineta] = new TH1F (name.Data(),"efficiency standard / puppi",pt_edges.size(), bins);
  
  
  name = Form ("histo_fake_puppi_%d",ibineta);
  histo_fake_puppi[ibineta] = new TH1F (name.Data(),"puppi",pt_edges.size(), bins);
  name = Form ("histo_fake_standard_%d",ibineta);
  histo_fake_standard[ibineta] = new TH1F (name.Data(),"standard",pt_edges.size(), bins);
  name = Form ("histo_fake_standard_over_puppi_%d",ibineta);
  histo_fake_standard_over_puppi[ibineta] = new TH1F (name.Data(),"fake standard / puppi",pt_edges.size(), bins);
  
 }
 
 
 
 
 
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
   
   
   float eta_gen = fabs(std_vector_jetGen_eta->at(j));
   int iGenEtaBin = -1;
   for (int ibineta = 0; ibineta < (eta_edges.size()-1); ibineta++) {
    if (ibineta != (eta_edges.size()-1)) {
     if (eta_gen >= eta_edges.at(ibineta) && eta_gen < eta_edges.at(ibineta+1)) {
      iGenEtaBin = ibineta;
      break;
     }
    }
   }
   
   
   if (iGenEtaBin != -1) {
    
    histo_efficiency_standard_den[iGenEtaBin]->Fill(pt_gen);
    histo_efficiency_puppi_den[iGenEtaBin]->Fill(pt_gen);
    
    std::pair<int, float> closest_standard_jet = getClosestIndexAndDR(std_vector_jetGen_eta->at(j), std_vector_jetGen_phi->at(j),    *std_vector_jet_eta, *std_vector_jet_phi);
    if (closest_standard_jet.second < 0.4) {
     histo_efficiency_standard_num[iGenEtaBin]->Fill(pt_gen);
    }
    std::pair<int, float> closest_puppi_jet = getClosestIndexAndDR(std_vector_jetGen_eta->at(j), std_vector_jetGen_phi->at(j),    *std_vector_puppijet_eta, *std_vector_puppijet_phi);
    if (closest_puppi_jet.second < 0.4) {
     histo_efficiency_puppi_num[iGenEtaBin]->Fill(pt_gen);
    }
   }
   
  }
  
  
  //---- fake
  for (UInt_t j = 0; j < std_vector_jet_phi->size(); ++j) {
   float pt_gen = std_vector_jet_pt->at(j);
   if (pt_gen > pt_edges.at(pt_edges.size()-1)) { //---- to deal with overflow bin in "view" zone
    pt_gen = pt_edges.at(pt_edges.size()-1)+0.5;
    //      std::cout << " pt_gen = " << pt_gen << std::endl;
   }
   
   
   float eta_jet = fabs(std_vector_jet_eta->at(j));
   int iRecoEtaBin = -1;
   for (int ibineta = 0; ibineta < (eta_edges.size()-1); ibineta++) {
    if (ibineta != (eta_edges.size()-1)) {
     if (eta_jet >= eta_edges.at(ibineta) && eta_jet < eta_edges.at(ibineta+1)) {
      iRecoEtaBin = ibineta;
      break;
     }
    }
   }
   
   if (iRecoEtaBin != -1) {
    
    std::pair<int, float> closest_standard_jet = getClosestIndexAndDR(std_vector_jet_eta->at(j), std_vector_jet_phi->at(j),    *std_vector_jetGen_eta, *std_vector_jetGen_phi);
    if (closest_standard_jet.second > 0.4) {
     histo_fake_standard[iRecoEtaBin]->Fill(pt_gen);
    }    
   }
  }
  
  for (UInt_t j = 0; j < std_vector_puppijet_phi->size(); ++j) {
   float pt_gen = std_vector_puppijet_pt->at(j);
   if (pt_gen > pt_edges.at(pt_edges.size()-1)) { //---- to deal with overflow bin in "view" zone
    pt_gen = pt_edges.at(pt_edges.size()-1)+0.5;
   }
   
   float eta_jet = fabs(std_vector_puppijet_eta->at(j));
   int iRecoEtaBin = -1;
   for (int ibineta = 0; ibineta < (eta_edges.size()-1); ibineta++) {
    if (ibineta != (eta_edges.size()-1)) {
     if (eta_jet >= eta_edges.at(ibineta) && eta_jet < eta_edges.at(ibineta+1)) {
      iRecoEtaBin = ibineta;
      break;
     }
    }
   }
   
   if (iRecoEtaBin != -1) {
    
    std::pair<int, float> closest_puppi_jet = getClosestIndexAndDR(std_vector_puppijet_eta->at(j), std_vector_puppijet_phi->at(j),    *std_vector_jetGen_eta, *std_vector_jetGen_phi);
    if (closest_puppi_jet.second > 0.4) {
     histo_fake_puppi[iRecoEtaBin]->Fill(pt_gen);
    }
   }
  }
  
  
 }
 
 
 for (int ibineta = 0; ibineta < (eta_edges.size()-1); ibineta++) {
  
  for (int ibin = 0; ibin <= pt_edges.size(); ibin++) {
   float num = histo_efficiency_standard_num[ibineta]->GetBinContent(ibin+1);
   float den = histo_efficiency_standard_den[ibineta]->GetBinContent(ibin+1);
   float eff = 0.;
   if (den != 0) eff = num / den;
   histo_efficiency_standard[ibineta]->SetBinContent(ibin+1, eff);
  }
  
  histo_efficiency_standard[ibineta]->SetLineColor(kBlue);
  histo_efficiency_standard[ibineta]->SetLineWidth(2);
  
  
  for (int ibin = 0; ibin <= pt_edges.size(); ibin++) {
   float num = histo_efficiency_puppi_num[ibineta]->GetBinContent(ibin+1);
   float den = histo_efficiency_puppi_den[ibineta]->GetBinContent(ibin+1);
   float eff = 0.;
   if (den != 0) eff = num / den;
   histo_efficiency_puppi[ibineta]->SetBinContent(ibin+1, eff);
  }
  
  histo_efficiency_puppi[ibineta]->SetLineColor(kRed);
  histo_efficiency_puppi[ibineta]->SetLineWidth(2);
  histo_efficiency_puppi[ibineta]->SetLineStyle(2);
 }
 
 
 
 TCanvas* ccResult_efficiency[100];
 TCanvas* ccResult_fake[100];
 
 for (int ibineta = 0; ibineta < (eta_edges.size()-1); ibineta++) {
  name = Form ("ccResult_efficiency_%d",ibineta);
  TString cut;
  cut = Form ("#eta [%f, %f]", eta_edges.at(ibineta), eta_edges.at(ibineta+1) );
  
  ccResult_efficiency[ibineta] = new TCanvas (name.Data(),cut.Data(),800,800);
  ccResult_efficiency[ibineta]->Divide(1,2);
  
  ccResult_efficiency[ibineta]->cd(1);
  histo_efficiency_standard[ibineta]->Draw();
  histo_efficiency_puppi[ibineta]->Draw("same");
  histo_efficiency_standard[ibineta]->GetXaxis()->SetTitle("gen jet p_{T} [GeV]");
  histo_efficiency_standard[ibineta]->GetYaxis()->SetTitle("matching efficiency");
  histo_efficiency_standard[ibineta]->GetYaxis()->SetRangeUser(0.0, 1.0);
  gPad->SetGrid();
  gPad->BuildLegend();
  
  ccResult_efficiency[ibineta]->cd(2);
  
  for (int ibin = 0; ibin <= pt_edges.size(); ibin++) {
   float num = histo_efficiency_standard[ibineta]->GetBinContent(ibin+1);
   float den = histo_efficiency_puppi[ibineta]->GetBinContent(ibin+1);
   float eff = 0.;
   if (den != 0) eff = num / den;
   histo_efficiency_standard_over_puppi[ibineta]->SetBinContent(ibin+1, eff);
  }
  histo_efficiency_standard_over_puppi[ibineta]->SetLineColor(kGreen+2);
  histo_efficiency_standard_over_puppi[ibineta]->SetLineWidth(2);
  histo_efficiency_standard_over_puppi[ibineta]->Draw();
  histo_efficiency_standard_over_puppi[ibineta]->GetXaxis()->SetTitle("reco jet p_{T} [GeV]");
  histo_efficiency_standard_over_puppi[ibineta]->GetYaxis()->SetTitle("puppi / standard efficiency");
  histo_efficiency_standard_over_puppi[ibineta]->GetYaxis()->SetRangeUser(0.0, 2.0);
  
  gPad->SetGrid();
  
  
  
  
  
  name = Form ("ccResult_fake_%d",ibineta);
  ccResult_fake[ibineta] = new TCanvas (name.Data(),cut.Data(),600,800);
  ccResult_fake[ibineta]->Divide(1,2);
  
  
  ccResult_fake[ibineta]->cd(1);
  histo_fake_standard[ibineta]->SetLineColor(kBlue);
  histo_fake_standard[ibineta]->SetLineWidth(2);
  
  histo_fake_puppi[ibineta]->SetLineColor(kRed);
  histo_fake_puppi[ibineta]->SetLineWidth(2);
  histo_fake_puppi[ibineta]->SetLineStyle(2);
  
  histo_fake_standard[ibineta]->Draw();
  histo_fake_puppi[ibineta]->Draw("same");
  histo_fake_standard[ibineta]->GetXaxis()->SetTitle("reco jet p_{T} [GeV]");
  histo_fake_standard[ibineta]->GetYaxis()->SetTitle("fake (a.u.)");
  gPad->SetGrid();
  gPad->SetLogy();
  gPad->BuildLegend();
  
  ccResult_fake[ibineta]->cd(2);
  
  for (int ibin = 0; ibin <= pt_edges.size(); ibin++) {
   float num = histo_fake_standard[ibineta]->GetBinContent(ibin+1);
   float den = histo_fake_puppi[ibineta]->GetBinContent(ibin+1);
   float eff = 0.;
   if (den != 0) eff = num / den;
   histo_fake_standard_over_puppi[ibineta]->SetBinContent(ibin+1, eff);
  }
  histo_fake_standard_over_puppi[ibineta]->SetLineColor(kGreen+2);
  histo_fake_standard_over_puppi[ibineta]->SetLineWidth(2);
  histo_fake_standard_over_puppi[ibineta]->Draw();
  histo_fake_standard_over_puppi[ibineta]->GetXaxis()->SetTitle("reco jet p_{T} [GeV]");
  histo_fake_standard_over_puppi[ibineta]->GetYaxis()->SetTitle("puppi / standard fake");
  histo_fake_standard_over_puppi[ibineta]->GetYaxis()->SetRangeUser(0.0, 10.0);
  
  gPad->SetGrid();
  
 }
 
 
}




