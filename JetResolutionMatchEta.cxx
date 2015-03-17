#include "Tools.cxx"

void JetResolutionMatchEta(){

 TFile* file = new TFile("latino_stepB_latinosYieldSkim_MC_ggHww.root","READ");
 
 TTree* latino = (TTree*) file->Get("latino");
 
 std::vector<float> eta_edges;
 
 eta_edges.push_back(0.0);
 eta_edges.push_back(1.3);
 eta_edges.push_back(2.5);
 eta_edges.push_back(3.0);
 eta_edges.push_back(4.5);
 
 
 
 
 std::vector<float> pt_edges;


 pt_edges.push_back(10.0);
 pt_edges.push_back(30.0);
 pt_edges.push_back(50.0);
 pt_edges.push_back(75.0);
 pt_edges.push_back(100.0);
 
 
 TH1F* histo_standard[100][100];
 TH1F* histo[100][100];
 TString name;
 
 TF1* fit_histo_standard[100][100];
 TF1* fit_histo[100][100];
 
 
 TCanvas* cc_standard[100];
 TCanvas* cc[100];
 
 TCanvas* cc_summary_standard;
 TCanvas* cc_summary;
 
 
 
 Float_t bins[1000];
 
  
 for (int ibinpt = 0; ibinpt < pt_edges.size(); ibinpt++) {
  TString cut;
  if (ibinpt != (pt_edges.size()-1)) cut = Form ("pt [%f, %f]", pt_edges.at(ibinpt), pt_edges.at(ibinpt+1) );
  else                             cut = Form ("pt [%f, -]",  pt_edges.at(ibinpt) ); 
  
  for (int ibineta = 0; ibineta < (eta_edges.size()-1); ibineta++) {
   
     
   name = Form ("histo_%d_%d",ibineta,ibinpt);
   histo[ibineta][ibinpt] = new TH1F (name.Data(),cut.Data(),100,0,3);
   
   name = Form ("histo_standard_%d_%d",ibineta,ibinpt);
   histo_standard[ibineta][ibinpt] = new TH1F (name.Data(),cut.Data(),100,0,3);
   
   
   name = Form ("fit_histo_%d_%d",ibineta,ibinpt);
   if      (ibinpt == 2) {
    fit_histo[ibineta][ibinpt] = new TF1 (name.Data(),"gaus(0)+pol2(3)",0.3,1.2);
    fit_histo[ibineta][ibinpt]->SetParameter(3,0.1);
   }
   else if (ibinpt == 1) fit_histo[ibineta][ibinpt] = new TF1 (name.Data(),"gaus(0)+pol2(3)",0.3,1.5);
   else             fit_histo[ibineta][ibinpt] = new TF1 (name.Data(),"gaus(0)+pol2(3)",0.3,2.0);
   fit_histo[ibineta][ibinpt]->SetParameter(1,1.0);
   fit_histo[ibineta][ibinpt]->SetParameter(2,0.5);

   
   name = Form ("fit_histo_standard_%d_%d",ibineta,ibinpt);
   if (ibinpt == 2) fit_histo_standard[ibineta][ibinpt] = new TF1 (name.Data(),"gaus(0)+pol2(3)",0.5,1.5);
   else             fit_histo_standard[ibineta][ibinpt] = new TF1 (name.Data(),"gaus+pol2(3)",0.3,2.0);
   fit_histo_standard[ibineta][ibinpt]->SetParameter(1,1.0);
   fit_histo_standard[ibineta][ibinpt]->SetParameter(2,0.5);
   
   
   fit_histo_standard[ibineta][ibinpt]->SetLineColor(kBlue);
   fit_histo[ibineta][ibinpt]->SetLineColor(kRed);   
  }
  bins[ibinpt] = pt_edges.at(ibinpt);
 }
 
 
 for (int ibineta = 0; ibineta < (eta_edges.size()-1); ibineta++) {
  name = Form ("cc_standard_%d",ibineta);
  cc_standard[ibineta] = new TCanvas (name.Data(),name.Data(),800,800);
  cc_standard[ibineta]->Divide(3,3);
 
  name = Form ("cc_puppi_%d",ibineta);
  cc[ibineta] = new TCanvas (name.Data(),name.Data(),800,800);
  cc[ibineta]->Divide(3,3);
 }
 
 cc_summary_standard = new TCanvas ("cc_summary_standard","standard",800,800);
 cc_summary          = new TCanvas ("cc_summary",         "puppi",   800,800);
 
 
 
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
 
 
 
 Long64_t nentries = latino->GetEntries();
 for (Long64_t i=0; i<nentries; i++) {
  latino->GetEntry(i);
  
  //---- match efficiency
  for (UInt_t j = 0; j < std_vector_jetGen_phi->size(); ++j) {
   
   float pt_gen = std_vector_jetGen_pt->at(j);
   int iGenPtBin = -1;
   for (int ibinpt = 0; ibinpt < pt_edges.size(); ibinpt++) {
    if (ibinpt != (pt_edges.size()-1)) {
     if (pt_gen >= pt_edges.at(ibinpt) && pt_gen < pt_edges.at(ibinpt+1)) {
      iGenPtBin = ibinpt;
      break;
     }
    }
    else {
     if (pt_gen >= pt_edges.at(ibinpt)) {
      iGenPtBin = ibinpt;
     }
    } 
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
   
   //     if (eta_gen > 1.5) std::cout << " eta_gen = " << eta_gen << std::endl;
   
   
   if (iGenPtBin != -1 && iGenEtaBin != -1) {
    std::pair<int, float> closest_standard_jet = getClosestIndexAndDR(std_vector_jetGen_eta->at(j), std_vector_jetGen_phi->at(j),    *std_vector_jet_eta, *std_vector_jet_phi);
    if (closest_standard_jet.second < 0.4) {
     histo_standard[iGenEtaBin][iGenPtBin]->Fill(std_vector_jet_pt->at(closest_standard_jet.first) / std_vector_jetGen_pt->at(j));
    }
    std::pair<int, float> closest_puppi_jet = getClosestIndexAndDR(std_vector_jetGen_eta->at(j), std_vector_jetGen_phi->at(j),    *std_vector_puppijet_eta, *std_vector_puppijet_phi);
    if (closest_puppi_jet.second < 0.4) {
     histo[iGenEtaBin][iGenPtBin]->Fill(std_vector_puppijet_pt->at(closest_puppi_jet.first) / std_vector_jetGen_pt->at(j));
    }
   }
   
  } 
  
 }
 
 
 
 for (int ibineta = 0; ibineta < (eta_edges.size()-1); ibineta++) {
  for (int ibinpt = 0; ibinpt < pt_edges.size(); ibinpt++) {
   
   cc[ibineta]-> cd(ibinpt+1);
   histo[ibineta][ibinpt]->Draw();  
   name = Form ("fit_histo_%d_%d",ibineta,ibinpt);
   histo[ibineta][ibinpt]->Fit(name.Data(),"RMQ");
   histo[ibineta][ibinpt]->Fit(name.Data(),"RMQ");
   histo[ibineta][ibinpt]->Fit(name.Data(),"RMQ");
   
   cc_standard[ibineta]-> cd(ibinpt+1);
   histo_standard[ibineta][ibinpt]->Draw();  
   name = Form ("fit_histo_standard_%d_%d",ibineta,ibinpt);
   histo_standard[ibineta][ibinpt]->Fit(name.Data(),"RMQ");
   histo_standard[ibineta][ibinpt]->Fit(name.Data(),"RMQ");
   histo_standard[ibineta][ibinpt]->Fit(name.Data(),"RMQ");
   
  }
 }
 
 
 
 
 TGraphErrors* gr_response[100];
 TGraphErrors* gr_response_fit[100];
 
 for (int ibineta = 0; ibineta < (eta_edges.size()-1); ibineta++) { 
  gr_response[ibineta] = new TGraphErrors();     gr_response[ibineta]->SetTitle("puppi");
  gr_response_fit[ibineta] = new TGraphErrors(); gr_response_fit[ibineta]->SetTitle("puppi");
  
  
  for (int ibinpt = 0; ibinpt < pt_edges.size(); ibinpt++) {
   float x;
   float delta;
   if (ibinpt!= (pt_edges.size()-1)) {
    x     =  (pt_edges.at(ibinpt) + pt_edges.at(ibinpt+1)) / 2.;
    delta = -(pt_edges.at(ibinpt) - pt_edges.at(ibinpt+1)) / 2.;
   }
   else {
    x = (pt_edges.at(ibinpt)) ;
    delta = 0;
   }
   float mean = histo[ibineta][ibinpt]->GetMean();
   float RMS  = histo[ibineta][ibinpt]->GetRMS();
   gr_response[ibineta]->SetPoint      (ibinpt, x, mean);
   gr_response[ibineta]->SetPointError (ibinpt, delta, RMS);
   //   std::cout << " delta = " << delta << std::endl;
   
   gr_response_fit[ibineta]->SetPoint      (ibinpt, x, fit_histo[ibineta][ibinpt]->GetParameter(1));
   gr_response_fit[ibineta]->SetPointError (ibinpt, delta, fabs(fit_histo[ibineta][ibinpt]->GetParameter(2)));
  }
  
  gr_response[ibineta]->SetFillColor(0);
  gr_response[ibineta]->SetMarkerSize(1);
  gr_response[ibineta]->SetMarkerStyle(22);
  gr_response[ibineta]->SetMarkerColor(kRed);
  gr_response[ibineta]->SetLineColor(kRed);
  
  gr_response_fit[ibineta]->SetFillColor(0);
  gr_response_fit[ibineta]->SetMarkerSize(1);
  gr_response_fit[ibineta]->SetMarkerStyle(24);
  gr_response_fit[ibineta]->SetMarkerColor(kRed);
  gr_response_fit[ibineta]->SetLineColor(kRed);
 }
 

 TGraphErrors* gr_response_standard[100];
 TGraphErrors* gr_response_standard_fit[100];


 for (int ibineta = 0; ibineta < (eta_edges.size()-1); ibineta++) { 
  
  gr_response_standard[ibineta] = new TGraphErrors();     gr_response_standard[ibineta]->SetTitle("standard");
  gr_response_standard_fit[ibineta] = new TGraphErrors(); gr_response_standard_fit[ibineta]->SetTitle("standard");
  for (int ibinpt = 0; ibinpt < pt_edges.size(); ibinpt++) {
   float x;
   float delta;
   if (ibinpt!= (pt_edges.size()-1)) {
    x     =  (pt_edges.at(ibinpt) + pt_edges.at(ibinpt+1)) / 2.;
    delta = -(pt_edges.at(ibinpt) - pt_edges.at(ibinpt+1)) / 2.;
   }
   else {
    x = (pt_edges.at(ibinpt)) ;
    delta = 0;
   }
   float mean = histo_standard[ibineta][ibinpt]->GetMean();
   float RMS  = histo_standard[ibineta][ibinpt]->GetRMS();
   gr_response_standard[ibineta]->SetPoint      (ibinpt, x, mean);
   gr_response_standard[ibineta]->SetPointError (ibinpt, delta, RMS);
   //   std::cout << " delta = " << delta << std::endl;
   
   gr_response_standard_fit[ibineta]->SetPoint      (ibinpt, x, fit_histo_standard[ibineta][ibinpt]->GetParameter(1));
   gr_response_standard_fit[ibineta]->SetPointError (ibinpt, delta, fabs(fit_histo_standard[ibineta][ibinpt]->GetParameter(2)));
  }
  
  gr_response_standard[ibineta]->SetFillColor(0);
  gr_response_standard[ibineta]->SetMarkerSize(1);
  gr_response_standard[ibineta]->SetMarkerStyle(22);
  gr_response_standard[ibineta]->SetMarkerColor(kBlue);
  gr_response_standard[ibineta]->SetLineColor(kBlue);
  
  gr_response_standard_fit[ibineta]->SetFillColor(0);
  gr_response_standard_fit[ibineta]->SetMarkerSize(1);
  gr_response_standard_fit[ibineta]->SetMarkerStyle(24);
  gr_response_standard_fit[ibineta]->SetMarkerColor(kBlue);
  gr_response_standard_fit[ibineta]->SetLineColor(kBlue);
 }
 
 

 
 TCanvas* ccResult = new TCanvas ("ccResult","Response",400,800);
 ccResult->Divide (1,(eta_edges.size()-1));
 
 TCanvas* ccResult_fit = new TCanvas ("ccResult_fit","Response fit",400,800);
 ccResult_fit->Divide (1,(eta_edges.size()-1));
 
 for (int ibineta = 0; ibineta < (eta_edges.size()-1); ibineta++) { 
  
  ccResult->cd(ibineta+1);
  gr_response[ibineta]->Draw("AP");
  gr_response[ibineta]->GetXaxis()->SetTitle("gen jet p_{T} [GeV]");
  gr_response[ibineta]->GetYaxis()->SetTitle("jet p_{T} / gen jet p_{T}");
  gr_response_standard[ibineta]->Draw("P");
  gPad->SetGrid();
  gPad->BuildLegend();
  
  
  ccResult_fit->cd(ibineta+1);
  gr_response_fit[ibineta]->Draw("AP");
  gr_response_fit[ibineta]->GetXaxis()->SetTitle("gen jet p_{T} [GeV]");
  gr_response_fit[ibineta]->GetYaxis()->SetTitle("jet p_{T} / gen jet p_{T}");
  gr_response_standard_fit[ibineta]->Draw("P");
  gPad->SetGrid();
//   gPad->SetLogx();
  gPad->BuildLegend();
  
 }
 
 
 //---- resolution
 
 TGraphErrors* gr_resolution[100];
 TGraphErrors* gr_resolution_fit[100];
 
 for (int ibineta = 0; ibineta < (eta_edges.size()-1); ibineta++) { 
  gr_resolution[ibineta] = new TGraphErrors();     gr_resolution[ibineta]->SetTitle("puppi");
  gr_resolution_fit[ibineta] = new TGraphErrors(); gr_resolution_fit[ibineta]->SetTitle("puppi");
  
  
  for (int ibinpt = 0; ibinpt < pt_edges.size(); ibinpt++) {
   float x;
   float delta;
   if (ibinpt!= (pt_edges.size()-1)) {
    x     =  (pt_edges.at(ibinpt) + pt_edges.at(ibinpt+1)) / 2.;
    delta = -(pt_edges.at(ibinpt) - pt_edges.at(ibinpt+1)) / 2.;
   }
   else {
    x = (pt_edges.at(ibinpt)) ;
    delta = 0;
   }
   float mean = histo[ibineta][ibinpt]->GetMean();
   float RMS  = histo[ibineta][ibinpt]->GetRMS();
   gr_resolution[ibineta]->SetPoint      (ibinpt, x, RMS/mean);
   gr_resolution[ibineta]->SetPointError (ibinpt, delta, 0);
   //   std::cout << " delta = " << delta << std::endl;
   
   gr_resolution_fit[ibineta]->SetPoint      (ibinpt, x, fabs(fit_histo[ibineta][ibinpt]->GetParameter(2)) / fit_histo[ibineta][ibinpt]->GetParameter(1));
   gr_resolution_fit[ibineta]->SetPointError (ibinpt, delta, 0);
  }
  
  gr_resolution[ibineta]->SetFillColor(0);
  gr_resolution[ibineta]->SetMarkerSize(1);
  gr_resolution[ibineta]->SetMarkerStyle(22+ibineta);
  gr_resolution[ibineta]->SetMarkerColor(kRed);
  gr_resolution[ibineta]->SetLineColor(kRed);
  
  gr_resolution_fit[ibineta]->SetFillColor(0);
  gr_resolution_fit[ibineta]->SetMarkerSize(1);
  gr_resolution_fit[ibineta]->SetMarkerStyle(24+ibineta);
  gr_resolution_fit[ibineta]->SetMarkerColor(kRed);
  gr_resolution_fit[ibineta]->SetLineColor(kRed);
 }
 
 
 TGraphErrors* gr_resolution_standard[100];
 TGraphErrors* gr_resolution_standard_fit[100];
 
 
 for (int ibineta = 0; ibineta < (eta_edges.size()-1); ibineta++) { 
  
  gr_resolution_standard[ibineta] = new TGraphErrors();     gr_resolution_standard[ibineta]->SetTitle("standard");
  gr_resolution_standard_fit[ibineta] = new TGraphErrors(); gr_resolution_standard_fit[ibineta]->SetTitle("standard");
  for (int ibinpt = 0; ibinpt < pt_edges.size(); ibinpt++) {
   float x;
   float delta;
   if (ibinpt!= (pt_edges.size()-1)) {
    x     =  (pt_edges.at(ibinpt) + pt_edges.at(ibinpt+1)) / 2.;
    delta = -(pt_edges.at(ibinpt) - pt_edges.at(ibinpt+1)) / 2.;
   }
   else {
    x = (pt_edges.at(ibinpt)) ;
    delta = 0;
   }
   float mean = histo_standard[ibineta][ibinpt]->GetMean();
   float RMS  = histo_standard[ibineta][ibinpt]->GetRMS();
   gr_resolution_standard[ibineta]->SetPoint      (ibinpt, x, RMS/mean);
   gr_resolution_standard[ibineta]->SetPointError (ibinpt, delta, 0);
   //   std::cout << " delta = " << delta << std::endl;
   
   gr_resolution_standard_fit[ibineta]->SetPoint      (ibinpt, x, fabs(fit_histo_standard[ibineta][ibinpt]->GetParameter(2)) / fit_histo_standard[ibineta][ibinpt]->GetParameter(1));
   gr_resolution_standard_fit[ibineta]->SetPointError (ibinpt, delta, 0);
  }
  
  gr_resolution_standard[ibineta]->SetFillColor(0);
  gr_resolution_standard[ibineta]->SetMarkerSize(1);
  gr_resolution_standard[ibineta]->SetMarkerStyle(22+ibineta);
  gr_resolution_standard[ibineta]->SetMarkerColor(kBlue);
  gr_resolution_standard[ibineta]->SetLineColor(kBlue);
  
  gr_resolution_standard_fit[ibineta]->SetFillColor(0);
  gr_resolution_standard_fit[ibineta]->SetMarkerSize(1);
  gr_resolution_standard_fit[ibineta]->SetMarkerStyle(24+ibineta);
  gr_resolution_standard_fit[ibineta]->SetMarkerColor(kBlue);
  gr_resolution_standard_fit[ibineta]->SetLineColor(kBlue);
 }
 
 

 
 TCanvas* ccResult_resolution = new TCanvas ("ccResult_resolution","Resolution",400,800);
 ccResult_resolution->Divide (1,(eta_edges.size()-1));
 
 TCanvas* ccResult_resolution_fit = new TCanvas ("ccResult_resolution_fit","Resolution fit",400,800);
 ccResult_resolution_fit->Divide (1,(eta_edges.size()-1));
 
 for (int ibineta = 0; ibineta < (eta_edges.size()-1); ibineta++) { 
  
  ccResult_resolution->cd(ibineta+1);
  gr_resolution[ibineta]->Draw("AP");
  gr_resolution[ibineta]->GetXaxis()->SetTitle("gen jet p_{T} [GeV]");
  gr_resolution[ibineta]->GetYaxis()->SetTitle("resolution jet p_{T} / gen jet p_{T}");
  gr_resolution_standard[ibineta]->Draw("P");
  gPad->SetGrid();
  gPad->BuildLegend();
  
  
  ccResult_resolution_fit->cd(ibineta+1);
  gr_resolution_fit[ibineta]->Draw("AP");
  gr_resolution_fit[ibineta]->GetXaxis()->SetTitle("gen jet p_{T} [GeV]");
  gr_resolution_fit[ibineta]->GetYaxis()->SetTitle("resolution fit jet p_{T} / gen jet p_{T}");
  gr_resolution_standard_fit[ibineta]->Draw("P");
  gPad->SetGrid();
//   gPad->SetLogx();
  gPad->BuildLegend();
  
 }
 
 
 TLegend* leg_summary_standard = new TLegend(0.1,0.7,0.48,0.9);
 
 cc_summary_standard->cd();
 for (int ibineta = 0; ibineta < (eta_edges.size()-1); ibineta++) {
  gr_resolution_standard_fit[ibineta]->SetMarkerColor(kBlue+ibineta); 
  gr_resolution_standard_fit[ibineta]->SetLineColor(kBlue+ibineta); 
  if (ibineta == 0) gr_resolution_standard_fit[ibineta]->Draw("AP");
  else              gr_resolution_standard_fit[ibineta]->Draw("P");
  gr_resolution_standard_fit[ibineta]->GetYaxis()->SetRangeUser(0.0,0.6);
  gr_resolution_standard_fit[ibineta]->GetXaxis()->SetTitle("gen jet p_{T} [GeV]");
  gr_resolution_standard_fit[ibineta]->GetYaxis()->SetTitle("resolution jet p_{T} / gen jet p_{T}");
  name = Form ("#eta [%f, %f]", eta_edges.at(ibineta), eta_edges.at(ibineta+1) );
  leg_summary_standard->AddEntry(gr_resolution_standard_fit[ibineta],name.Data(),"lep");  
 }
 leg_summary_standard->Draw(); 
 cc_summary_standard->SetGrid();


 TLegend* leg_summary = new TLegend(0.5,0.7,0.9,0.9);
 
 cc_summary->cd();
 for (int ibineta = 0; ibineta < (eta_edges.size()-1); ibineta++) { 
  gr_resolution_fit[ibineta]->SetMarkerColor(kRed+ibineta);
  gr_resolution_fit[ibineta]->SetLineColor(kRed+ibineta);
  if (ibineta == 0) gr_resolution_fit[ibineta]->Draw("AP");
  else              gr_resolution_fit[ibineta]->Draw("P");
  gr_resolution_fit[ibineta]->GetYaxis()->SetRangeUser(0.0,0.6);
  gr_resolution_fit[ibineta]->GetXaxis()->SetTitle("gen jet p_{T} [GeV]");
  gr_resolution_fit[ibineta]->GetYaxis()->SetTitle("resolution jet p_{T} / gen jet p_{T}");
  name = Form ("#eta [%f, %f]", eta_edges.at(ibineta), eta_edges.at(ibineta+1) );
  leg_summary->AddEntry(gr_resolution_fit[ibineta],name.Data(),"lep");
  
 }
 leg_summary->Draw(); 
 cc_summary->SetGrid();
 
 
 
 cc_summary_standard->SaveAs("summary_standard.png");
 cc_summary         ->SaveAs("summary_puppi.png");
 
 
 for (int ibineta = 0; ibineta < (eta_edges.size()-1); ibineta++) {
  name = Form ("puppi_eta_%f_%f.png", eta_edges.at(ibineta), eta_edges.at(ibineta+1) ); 
  cc[ibineta]-> SaveAs(name.Data());
  name = Form ("eta_%f_%f.png", eta_edges.at(ibineta), eta_edges.at(ibineta+1) );
  cc_standard[ibineta]-> SaveAs(name.Data());  
 }
 
}