{
 
 TFile* file = new TFile("latino_stepB_latinosYieldSkim_MC_ggHww.root","READ");
 
 TTree* latino = (TTree*) file->Get("latino");

 std::vector<float> pt_edges;
 
 pt_edges.push_back(0.0);
 pt_edges.push_back(10.0);
 pt_edges.push_back(15.0);
 pt_edges.push_back(20.0);
 pt_edges.push_back(25.0);
 pt_edges.push_back(30.0);
 pt_edges.push_back(35.0);
 pt_edges.push_back(50.0);
 pt_edges.push_back(100.0);
 pt_edges.push_back(200.0);
 pt_edges.push_back(300.0);
 
 TCanvas* cc = new TCanvas ("cc","",800,800);
 cc->Divide(4,4);
 
 TH1F* histo[100];
 
 for (int ibin = 0; ibin < pt_edges.size(); ibin++) {
  TString cut;
  if (ibin!= (pt_edges.size()-1)) {
   cut  = Form ("std_vector_jetGen_pt[0]>=0 && std_vector_puppijet_pt[0]>0 && std_vector_jetGen_pt[0]>%f && std_vector_jetGen_pt[0]<=%f", pt_edges.at(ibin), pt_edges.at(ibin+1));
  }
  else {
   cut  = Form ("std_vector_jetGen_pt[0]>=0 && std_vector_puppijet_pt[0]>0 && std_vector_jetGen_pt[0]>%f", pt_edges.at(ibin));
  }
  
  TString name = Form ("histo_%d",ibin);
  TString nameToDraw = Form ("std_vector_puppijet_pt[0] / std_vector_jetGen_pt[0] >> histo_%d",ibin);
  histo[ibin] = new TH1F (name.Data(),cut.Data(),100,0,3);
  latino->Draw(nameToDraw.Data(), cut.Data(), "goff");
  cc-> cd(ibin+1);
  histo[ibin]->Draw();  
 }
 
 
 TGraphErrors* gr_response = new TGraphErrors();
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
 }
 
 gr_response->SetMarkerSize(1);
 gr_response->SetMarkerStyle(22);
 gr_response->SetMarkerColor(kRed);
 gr_response->SetLineColor(kRed);
 
 TCanvas* ccResult = new TCanvas ("ccResult","",800,800);
 gr_response->Draw("AP");
 gr_response->GetXaxis()->SetTitle("gen jet p_{T}");
 gr_response->GetYaxis()->SetTitle("jet p_{T} / gen jet p_{T}");
 ccResult->SetGrid();
 
 
 //  latino->Draw("std_vector_jetGen_pt[0]", "std_vector_jetGen_pt[0]>=0");
//  latino->Draw("std_vector_puppijet_pt[0] / std_vector_jetGen_pt[0]", "std_vector_jetGen_pt[0]>=0 && std_vector_puppijet_pt[0]>0");

 
}