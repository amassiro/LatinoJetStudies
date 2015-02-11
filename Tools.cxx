// 
// Tools
// 

float DPhi (float phi1, float phi2) {
 float delta = fabs(phi1-phi2);
 if (delta > 3.14159265) {
  return (delta-3.14159265);
 }
 else {
  return delta;
 }
}
 
 
float DR (float eta1, float eta2, float phi1, float phi2) {
 float deta = eta1 - eta2;
 float dphi = DPhi(phi1, phi2);
 return sqrt(deta*deta + dphi*dphi);
}



int getClosest(float eta, float phi, std::vector<float>& v_eta, std::vector<float>& v_phi) {
 int iObj = -1;
 if (v_eta.size() != v_phi.size()) return iObj;
 
 float mindr = 100000;
 for (UInt_t j = 0; j < v_eta.size(); ++j) {
  float temp_dr = DR(eta, v_eta.at(j), phi, v_phi.at(j));
  if (temp_dr < mindr) {
   mindr = temp_dr;
   iObj = j;
  }
 }
 
 return iObj;
}



std::pair<int, float> getClosestIndexAndDR(float eta, float phi, std::vector<float>& v_eta, std::vector<float>& v_phi) {
 int iObj = -1;
 float mindr = 100000;
 
 std::pair <int, float> result;
 result.first = iObj;
 result.second = mindr;
 
 if (v_eta.size() != v_phi.size()) return result;
 
 for (UInt_t j = 0; j < v_eta.size(); ++j) {
  float temp_dr = DR(eta, v_eta.at(j), phi, v_phi.at(j));
  if (temp_dr < mindr) {
   mindr = temp_dr;
   iObj = j;
  }
 }
 
 result.first = iObj;
 result.second = mindr;
 
 return result;
}
