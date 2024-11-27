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

int task8(){
    TFile file("newroot.root","UPDATE");
    TTree *tree = (TTree *)file.Get("MyTree");
    Int_t nph;
    Float_t GoodPh[2];
    
    tree->SetBranchAddress("nph", &nph);
    Float_t eph[nph], thetaph[nph], phiph[nph];
    tree->SetBranchAddress("eph", eph);
    tree->SetBranchAddress("thetaph", thetaph);
    tree->SetBranchAddress("phiph", phiph);
 
    TCanvas *s = new TCanvas();
    TH1D* hM2ph= new TH1D("heph","Energy of photons",100,0.1,0.2);
    TH1D* hangle= new TH1D("hangle","Angle of photons",100,-4.,1.);

    Int_t ncandidate = 0;
    double e2ph, MinvPh, angle;

    Long64_t nentries = tree->GetEntries();
    for (Long64_t ientry = 0; ientry < nentries; ++ientry) {
        
        tree->GetEntry(ientry);
        
        for(int i =0; i<nph-1; ++i){
		 
		    for(int j=i+1; j<nph; ++j){
                TLorentzVector ph1, ph2, ph3, ph4;
                ph1 = TLorentzVector(eph[i]*sin(thetaph[i])*cos(phiph[i]),
                                     eph[i]*sin(phiph[i])*sin(thetaph[i]),
                                     eph[i]*cos(thetaph[i]),eph[i]);
                ph2 = TLorentzVector(eph[j]*sin(thetaph[j])*cos(phiph[j]),
                                     eph[j]*sin(phiph[j])*sin(thetaph[j]),
                                     eph[j]*cos(thetaph[j]),eph[j]);

                //ph1.SetE(phen[i]); ph1.SetPhi(phphi[i]); ph1.SetTheta(phth[i]);                                                                                                 
                //ph2.SetE(phen[j]); ph2.SetPhi(phphi[j]); ph2.SetTheta(phth[j]); 
                //TLorentzVector g1(eph), g2();
                //g1.SetPhi(phi); g1.SetTheta(th);
                //(g1+g2).M();  

                e2ph = (ph1+ph2).E(); 
                MinvPh = (ph1+ph2).M(); 
                angle = 1 - MinvPh/(2*eph[i]*eph[j]);
                hangle -> Fill(angle);

                if (MinvPh<0.2 && MinvPh>0.1) {
                    if (ncandidate < 2) GoodPh[ncandidate] = MinvPh;
                    ncandidate++;
                }

            }
        }  
        if (ncandidate == 2){
            hM2ph->Fill(GoodPh[0]);
            hM2ph->Fill(GoodPh[1]);
        }
        ncandidate = 0;  
    }

    hM2ph->GetYaxis()->SetTitle("Number of particles");
    hM2ph->GetXaxis()->SetTitle("Energy");
    //hM2ph->Draw();
    //s->SaveAs("hM2ph.png");

    hangle->GetYaxis()->SetTitle("Number of particles");
    hangle->GetXaxis()->SetTitle("Angle");
    //hangle->Draw();
    //s->SaveAs("hangle.png");

    file.Write();
    file.Save();

    return 0;
}
