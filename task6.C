#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <TTree.h>
#include <TFile.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TVector.h>
#include <TLorentzVector.h>

using namespace std;

int task6(){
    double Ecm = 1020;
    double MKs = 497.6;
    double MPi = 139.57;
    double pKs = sqrt(Ecm*Ecm/4-MKs*MKs);
    double V = pKs/sqrt(MKs*MKs+pKs*pKs);
    double pPi = sqrt(Ecm*Ecm/16-MPi*MPi);
    TLorentzVector Ks[100000];
    TLorentzVector PiPos[100000];
    TLorentzVector PiNeg[100000];
  
    TF1 u("u","sin(x)*sin(x)", 0, 3.1415);
    TCanvas *s = new TCanvas();
    TH1D* ThetaKs = new TH1D("h","h",1000,0,1);
    TH1D* PhiKs = new TH1D("h","h",1000,0,1);
    double phi, theta, phiPi, thetaPi;

    for(int i = 0; i < 100000; i++){
      phi = 2*3.1415*gRandom->Rndm();
      theta = u.GetRandom();
      PhiKs -> Fill(phi);
      ThetaKs->Fill(theta);
      Ks[i] = TLorentzVector(pKs*sin(theta)*cos(phi), pKs*sin(phi)*sin(theta),pKs*cos(theta),Ecm/2);

      phiPi = 2*3.1415*gRandom->Rndm();
      thetaPi= 2*3.1415*gRandom->Rndm();

      PiPos[i] = TLorentzVector(pPi*sin(theta-thetaPi)*cos(phiPi), pPi*sin(phiPi)*sin(theta-thetaPi),pPi*cos(theta-thetaPi),Ecm/4); 
      PiNeg[i] = TLorentzVector(pPi*sin(theta+thetaPi)*cos(phiPi), pPi*sin(phiPi)*sin(theta+thetaPi),pPi*cos(theta+thetaPi),Ecm/4);

    }
    ThetaKs->Draw();
    s->SaveAs("ThetaKs.png");
    PhiKs -> Draw();
    s->SaveAs("PhiKs.png");
    
  return 0;
}
