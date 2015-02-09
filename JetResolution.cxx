{
 
 TFile* file = new TFile("latino_stepB_latinosYieldSkim_MC_ggHww.root","READ");
 
 TTree* latino = (TTree*) file->Get("latino");
 
//  latino->Draw("std_vector_jetGen_pt[0]", "std_vector_jetGen_pt[0]>=0");
 latino->Draw("std_vector_puppijet_pt[0] / std_vector_jetGen_pt[0]", "std_vector_jetGen_pt[0]>=0 && std_vector_puppijet_pt[0]>0");

 
}