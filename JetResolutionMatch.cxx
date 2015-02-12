#include "Tools.cxx"

void JetResolutionMatch(){

 TFile* file = new TFile("latino_stepB_latinosYieldSkim_MC_ggHww.root","READ");
 
 TTree* latino = (TTree*) file->Get("latino");
 
 std::vector<float> pt_edges;
 
 pt_edges.push_back(5.0);
 pt_edges.push_back(20.0);
 pt_edges.push_back(30.0);
 pt_edges.push_back(40.0);
 pt_edges.push_back(50.0);
 pt_edges.push_back(100.0);
 pt_edges.push_back(200.0);
 pt_edges.push_back(300.0);
 
 TH1F* histo_standard[100];
 TH1F* histo[100];
 TString name;
 
 TF1* fit_histo_standard[100];
 TF1* fit_histo[100];
 
 
 
 Float_t bins[1000];
 for (int ibin = 0; ibin < pt_edges.size(); ibin++) {
  
  TString cut;
  if (ibin != (pt_edges.size()-1)) cut = Form ("pt [%f, %f]", pt_edges.at(ibin), pt_edges.at(ibin+1) );
  else                             cut = Form ("pt [%f, -]",  pt_edges.at(ibin) ); 

  name = Form ("histo_%d",ibin);
  histo[ibin] = new TH1F (name.Data(),cut.Data(),100,0,3);
  
  name = Form ("histo_standard_%d",ibin);
  histo_standard[ibin] = new TH1F (name.Data(),cut.Data(),100,0,3);
  
  
  name = Form ("fit_histo_%d",ibin);
  if      (ibin <= 1) fit_histo[ibin] = new TF1 (name.Data(),"gaus(0)+pol2(3)",0.2,3.0);
  else if (ibin == 3) fit_histo[ibin] = new TF1 (name.Data(),"gaus(0)+pol2(3)",0.3,1.4);
  else                fit_histo[ibin] = new TF1 (name.Data(),"gaus(0)+pol2(3)",0.2,1.5);
  fit_histo[ibin]->SetParameter(1,1.0);
  fit_histo[ibin]->SetParameter(2,0.5);
  
  name = Form ("fit_histo_standard_%d",ibin);
  if      (ibin <= 1) fit_histo_standard[ibin] = new TF1 (name.Data(),"gaus+pol2(3)",0.2,3.0);
  else if (ibin == 2) fit_histo_standard[ibin] = new TF1 (name.Data(),"gaus(0)+pol2(3)",0.5,1.7);
  else if (ibin == 3) fit_histo_standard[ibin] = new TF1 (name.Data(),"gaus(0)+pol2(3)",0.3,1.4);
  else                fit_histo_standard[ibin] = new TF1 (name.Data(),"gaus+pol2(3)",0.2,1.5);
  fit_histo_standard[ibin]->SetParameter(1,1.0);
  fit_histo_standard[ibin]->SetParameter(2,0.5);
  
  
  fit_histo_standard[ibin]->SetLineColor(kBlue);
  fit_histo[ibin]->SetLineColor(kRed);
  
  bins[ibin] = pt_edges.at(ibin);
 }
 
 bins[pt_edges.size()] = pt_edges.at(pt_edges.size()-1) + 100;
 
 
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
 
 
 
 TCanvas* cc_standard = new TCanvas ("cc_standard","standard",800,800);
 cc_standard->Divide(4,4);
 
 TCanvas* cc = new TCanvas ("cc","puppi",800,800);
 cc->Divide(4,4);
 
 
 
 Long64_t nentries = latino->GetEntries();
 for (Long64_t i=0; i<nentries; i++) {
  latino->GetEntry(i);
  
  //---- match efficiency
  for (UInt_t j = 0; j < std_vector_jetGen_phi->size(); ++j) {
   
   float pt_gen = std_vector_jetGen_pt->at(j);
   int iGenPtBin = 0;
   for (int ibin = 0; ibin < pt_edges.size(); ibin++) {
    if (ibin != (pt_edges.size()-1)) {
     if (pt_gen < pt_edges.at(ibin+1)) {
      iGenPtBin = ibin;
      break;
     }
    }
    else {
     if (pt_gen > pt_edges.at(ibin)) {
      iGenPtBin = ibin;
      break;
     }
    } 
   }
   
   
   if (pt_gen > pt_edges.at(pt_edges.size()-1)) { //---- to deal with overflow bin in "view" zone
    pt_gen = pt_edges.at(pt_edges.size()-1)+0.5;
   }
      
   std::pair<int, float> closest_standard_jet = getClosestIndexAndDR(std_vector_jetGen_eta->at(j), std_vector_jetGen_phi->at(j),    *std_vector_jet_eta, *std_vector_jet_phi);
   if (closest_standard_jet.second < 0.4) {
    histo_standard[iGenPtBin]->Fill(std_vector_jet_pt->at(closest_standard_jet.first) / std_vector_jetGen_pt->at(j));
   }
   std::pair<int, float> closest_puppi_jet = getClosestIndexAndDR(std_vector_jetGen_eta->at(j), std_vector_jetGen_phi->at(j),    *std_vector_puppijet_eta, *std_vector_puppijet_phi);
   if (closest_puppi_jet.second < 0.4) {
    histo[iGenPtBin]->Fill(std_vector_puppijet_pt->at(closest_puppi_jet.first) / std_vector_jetGen_pt->at(j));
   }
   
  } 
  
 }
 
 
 
 
 for (int ibin = 0; ibin < pt_edges.size(); ibin++) {
  
  cc-> cd(ibin+1);
  histo[ibin]->Draw();  
  name = Form ("fit_histo_%d",ibin);
  histo[ibin]->Fit(name.Data(),"RMQ");
  histo[ibin]->Fit(name.Data(),"RMQ");
  
  cc_standard-> cd(ibin+1);
  histo_standard[ibin]->Draw();  
  name = Form ("fit_histo_standard_%d",ibin);
  histo_standard[ibin]->Fit(name.Data(),"RMQ");
  histo_standard[ibin]->Fit(name.Data(),"RMQ");
  
 }
 
 
 
 
 
 
 
 
 
 TGraphErrors* gr_response = new TGraphErrors();     gr_response->SetTitle("puppi");
 TGraphErrors* gr_response_fit = new TGraphErrors(); gr_response_fit->SetTitle("puppi");
 for (int ibin = 0; ibin < pt_edges.size(); ibin++) {
  float x;
  float delta;
  if (ibin!= (pt_edges.size()-1)) {
   x     =  (pt_edges.at(ibin) + pt_edges.at(ibin+1)) / 2.;
   delta = -(pt_edges.at(ibin) - pt_edges.at(ibin+1)) / 2.;
  }
  else {
   x = (pt_edges.at(ibin)) ;
   delta = 0;
  }
  float mean = histo[ibin]->GetMean();
  float RMS  = histo[ibin]->GetRMS();
  gr_response->SetPoint      (ibin, x, mean);
  gr_response->SetPointError (ibin, delta, RMS);
  //   std::cout << " delta = " << delta << std::endl;
  
  gr_response_fit->SetPoint      (ibin, x, fit_histo[ibin]->GetParameter(1));
  gr_response_fit->SetPointError (ibin, delta, fabs(fit_histo[ibin]->GetParameter(2)));
 }
 
 gr_response->SetFillColor(0);
 gr_response->SetMarkerSize(1);
 gr_response->SetMarkerStyle(22);
 gr_response->SetMarkerColor(kRed);
 gr_response->SetLineColor(kRed);
 
 gr_response_fit->SetFillColor(0);
 gr_response_fit->SetMarkerSize(1);
 gr_response_fit->SetMarkerStyle(24);
 gr_response_fit->SetMarkerColor(kRed);
 gr_response_fit->SetLineColor(kRed);
 
 
 
 
 TGraphErrors* gr_response_standard = new TGraphErrors();     gr_response_standard->SetTitle("standard");
 TGraphErrors* gr_response_standard_fit = new TGraphErrors(); gr_response_standard_fit->SetTitle("standard");
 for (int ibin = 0; ibin < pt_edges.size(); ibin++) {
  float x;
  float delta;
  if (ibin!= (pt_edges.size()-1)) {
   x     =  (pt_edges.at(ibin) + pt_edges.at(ibin+1)) / 2.;
   delta = -(pt_edges.at(ibin) - pt_edges.at(ibin+1)) / 2.;
  }
  else {
   x = (pt_edges.at(ibin)) ;
   delta = 0;
  }
  float mean = histo_standard[ibin]->GetMean();
  float RMS  = histo_standard[ibin]->GetRMS();
  gr_response_standard->SetPoint      (ibin, x, mean);
  gr_response_standard->SetPointError (ibin, delta, RMS);
  //   std::cout << " delta = " << delta << std::endl;
  
  gr_response_standard_fit->SetPoint      (ibin, x, fit_histo_standard[ibin]->GetParameter(1));
  gr_response_standard_fit->SetPointError (ibin, delta, fabs(fit_histo_standard[ibin]->GetParameter(2)));
 }
 
 gr_response_standard->SetFillColor(0);
 gr_response_standard->SetMarkerSize(1);
 gr_response_standard->SetMarkerStyle(22);
 gr_response_standard->SetMarkerColor(kBlue);
 gr_response_standard->SetLineColor(kBlue);
 
 gr_response_standard_fit->SetFillColor(0);
 gr_response_standard_fit->SetMarkerSize(1);
 gr_response_standard_fit->SetMarkerStyle(24);
 gr_response_standard_fit->SetMarkerColor(kBlue);
 gr_response_standard_fit->SetLineColor(kBlue);
 
 
 
 
 TCanvas* ccResult = new TCanvas ("ccResult","Response",800,800);
 gr_response->Draw("AP");
 gr_response->GetXaxis()->SetTitle("gen jet p_{T} [GeV]");
 gr_response->GetYaxis()->SetTitle("jet p_{T} / gen jet p_{T}");
 gr_response_standard->Draw("P");
 ccResult->SetGrid();
 ccResult->BuildLegend();
 
 
 TCanvas* ccResult_fit = new TCanvas ("ccResult_fit","Response fit",800,800);
 gr_response_fit->Draw("AP");
 gr_response_fit->GetXaxis()->SetTitle("gen jet p_{T} [GeV]");
 gr_response_fit->GetYaxis()->SetTitle("jet p_{T} / gen jet p_{T}");
 gr_response_standard_fit->Draw("P");
 ccResult_fit->SetGrid();
 ccResult_fit->SetLogx();
 ccResult_fit->BuildLegend();
 
 
 
 
 
 
 
 //---- resolution
 
 TGraphErrors* gr_resolution = new TGraphErrors();     gr_resolution->SetTitle("puppi");
 TGraphErrors* gr_resolution_fit = new TGraphErrors(); gr_resolution_fit->SetTitle("puppi");
 for (int ibin = 0; ibin < pt_edges.size(); ibin++) {
  float x;
  float delta;
  if (ibin!= (pt_edges.size()-1)) {
   x     =  (pt_edges.at(ibin) + pt_edges.at(ibin+1)) / 2.;
   delta = -(pt_edges.at(ibin) - pt_edges.at(ibin+1)) / 2.;
  }
  else {
   x = (pt_edges.at(ibin)) ;
   delta = 0;
  }
  float mean = histo[ibin]->GetMean();
  float RMS  = histo[ibin]->GetRMS();
  gr_resolution->SetPoint      (ibin, x, RMS/mean);
  gr_resolution->SetPointError (ibin, delta, 0);
  //   std::cout << " delta = " << delta << std::endl;
  
  gr_resolution_fit->SetPoint      (ibin, x, fabs(fit_histo[ibin]->GetParameter(2)) / fit_histo[ibin]->GetParameter(1));
  gr_resolution_fit->SetPointError (ibin, delta, 0.);
 }
 
 gr_resolution->SetFillColor(0);
 gr_resolution->SetMarkerSize(1);
 gr_resolution->SetMarkerStyle(22);
 gr_resolution->SetMarkerColor(kRed);
 gr_resolution->SetLineColor(kRed);
 
 gr_resolution_fit->SetFillColor(0);
 gr_resolution_fit->SetMarkerSize(1);
 gr_resolution_fit->SetMarkerStyle(24);
 gr_resolution_fit->SetMarkerColor(kRed);
 gr_resolution_fit->SetLineColor(kRed);
 
 
 
 
 TGraphErrors* gr_resolution_standard = new TGraphErrors();      gr_resolution_standard->SetTitle("standard");
 TGraphErrors* gr_resolution_standard_fit = new TGraphErrors();  gr_resolution_standard_fit->SetTitle("standard");
 for (int ibin = 0; ibin < pt_edges.size(); ibin++) {
  float x;
  float delta;
  if (ibin!= (pt_edges.size()-1)) {
   x     =  (pt_edges.at(ibin) + pt_edges.at(ibin+1)) / 2.;
   delta = -(pt_edges.at(ibin) - pt_edges.at(ibin+1)) / 2.;
  }
  else {
   x = (pt_edges.at(ibin)) ;
   delta = 0;
  }
  float mean = histo_standard[ibin]->GetMean();
  float RMS  = histo_standard[ibin]->GetRMS();
  gr_resolution_standard->SetPoint      (ibin, x, RMS/mean);
  gr_resolution_standard->SetPointError (ibin, delta, 0);
  //   std::cout << " delta = " << delta << std::endl;
  
  gr_resolution_standard_fit->SetPoint      (ibin, x, fabs(fit_histo_standard[ibin]->GetParameter(2)) / fit_histo_standard[ibin]->GetParameter(1));
  gr_resolution_standard_fit->SetPointError (ibin, delta, 0.);
 }
 
 gr_resolution_standard->SetFillColor(0);
 gr_resolution_standard->SetMarkerSize(1);
 gr_resolution_standard->SetMarkerStyle(22);
 gr_resolution_standard->SetMarkerColor(kBlue);
 gr_resolution_standard->SetLineColor(kBlue);
 
 gr_resolution_standard_fit->SetFillColor(0);
 gr_resolution_standard_fit->SetMarkerSize(1);
 gr_resolution_standard_fit->SetMarkerStyle(24);
 gr_resolution_standard_fit->SetMarkerColor(kBlue);
 gr_resolution_standard_fit->SetLineColor(kBlue);
 
 
 
 TCanvas* ccResult_resolution = new TCanvas ("ccResult_resolution","Resolution",800,800);
 gr_resolution->Draw("AP");
 gr_resolution->GetXaxis()->SetTitle("gen jet p_{T} [GeV]");
 gr_resolution->GetYaxis()->SetTitle("resolution jet p_{T} / gen jet p_{T}");
 gr_resolution_standard->Draw("P");
 ccResult_resolution->SetGrid();
 ccResult_resolution->BuildLegend();
 
 
 TCanvas* ccResult_resolution_fit = new TCanvas ("ccResult_resolution_fit","Resolution",800,800);
 gr_resolution_fit->Draw("AP");
 gr_resolution_fit->GetXaxis()->SetTitle("gen jet p_{T} [GeV]");
 gr_resolution_fit->GetYaxis()->SetTitle("resolution fit jet p_{T} / gen jet p_{T}");
 gr_resolution_standard_fit->Draw("P");
 ccResult_resolution_fit->SetGrid();
 ccResult_resolution_fit->SetLogx();
 ccResult_resolution_fit->BuildLegend(); 
 
 
 
 
 
 
 
 
 
}