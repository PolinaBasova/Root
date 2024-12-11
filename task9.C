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
#include <TMath.h>
#include <TGraphErrors.h>
#include <malloc.h>

using namespace std;

int task9(){
    TFile file("newroot.root","UPDATE");
    TTree *tree = (TTree *)file.Get("MyTree");

    Int_t nph;
    
    tree->SetBranchAddress("nph", &nph);
    Float_t eph[nph], thetaph[nph], phiph[nph];
    TLorentzVector GoodPh[2];
    tree->SetBranchAddress("eph", eph);
    tree->SetBranchAddress("thetaph", thetaph);
    tree->SetBranchAddress("phiph", phiph);
 
    TDirectory* subdir = file.mkdir("subdir","subdir title");

    //TH1D* hangle = new TH1D("hangle","Angle of photons",100,-4.,1.);

    Int_t ncandidate = 0;
    double e2ph, MinvPh, angle;

    int n_pair = 0;
    double* theta = NULL;
    double* phi = NULL;

    Long64_t nentries = tree->GetEntries();
    for (Long64_t ientry = 0; ientry < nentries; ++ientry) {
        
        tree->GetEntry(ientry);
        
        for(int i =0; i<nph-1; ++i){
		 
		    for(int j=i+1; j<nph; ++j){
                TLorentzVector ph1, ph2;
                ph1 = TLorentzVector(eph[i]*sin(thetaph[i])*cos(phiph[i]),
                                     eph[i]*sin(phiph[i])*sin(thetaph[i]),
                                     eph[i]*cos(thetaph[i]),eph[i]);
                ph2 = TLorentzVector(eph[j]*sin(thetaph[j])*cos(phiph[j]),
                                     eph[j]*sin(phiph[j])*sin(thetaph[j]),
                                     eph[j]*cos(thetaph[j]),eph[j]);  

                e2ph = (ph1+ph2).E(); 
                MinvPh = (ph1+ph2).M(); 
                angle = 1 - MinvPh/(2*eph[i]*eph[j]);
                //hangle -> Fill(angle);

                if (MinvPh<0.2 && MinvPh>0.1) {
                    if (ncandidate < 2) {
                        GoodPh[ncandidate] = ph1+ph2; 
                        }
                    ncandidate++;
                }

            }
        }  
        if (ncandidate == 2){
            theta = (double*)realloc(theta, (n_pair + 1) * sizeof(double));
            theta[n_pair] = GoodPh[0].Theta();
            phi = (double*)realloc(phi, (n_pair + 1) * sizeof(double));
            phi[n_pair] = GoodPh[0].Phi();
            n_pair++;
            theta = (double*)realloc(theta, (n_pair + 1) * sizeof(double));
            theta[n_pair] = GoodPh[1].Theta();
            phi = (double*)realloc(phi, (n_pair + 1) * sizeof(double));
            phi[n_pair] = GoodPh[1].Phi();
            n_pair++;
        }
        ncandidate = 0;  
    }
    
    TCanvas *s2 = new TCanvas("s2","s2",500,300);
    TLegend *legend2 = new TLegend(0.7,0.3,0.9,0.5);
    double x[18];
    double ex[18];
    double y[18];
    double ey[18];
    for (int j = 0; j<18; ++j){
        y[j]=0;
    }

    for (int i = 0; i<18; ++i){
        x[i]=5+i*10;
        ex[i]=0;
        ey[i]=TMath::RMS(n_pair-1, theta);
        for (int j = 0; j<n_pair-1; ++j){
            if ( (x[i]-5)*3.14/180 < theta[j] < (x[i]+5)*3.14/180 ) {y[i]+=1;}
        }
    }
    
    auto theta_gr = new TGraphErrors(18, x, y, ex, ey);
    theta_gr->SetTitle("TGraphErrors of theta");
    theta_gr->GetYaxis()->SetTitle("Number of particles");
    theta_gr->GetXaxis()->SetTitle("Theta");
    theta_gr->SetMarkerColor(2);
    theta_gr->SetMarkerStyle(21);
    theta_gr->Draw();
    legend2->AddEntry(theta_gr, "Histogram \"Theta\"", "l");
    legend2->AddEntry((TObject*)0, "", "");
    legend2->AddEntry((TObject*)0, "Points", "lep");
    legend2->Draw();

    s2->Write();

    TCanvas *s = new TCanvas("s","s",500,300);
    TLegend *legend = new TLegend(0.7,0.3,0.9,0.5);

    double x2[35];
    double ex2[35];
    double y2[35];
    double ey2[35];
    x2[0]=-175;
    for (int j = 1; j<35; ++j){
        x2[j]=x2[j-1]+10;
        y2[j]=0;
    }

    for (int j = 0; j<35; ++j){
        for (int i = 0; i<n_pair-1; ++i){
            if ( (x2[j]-5)*3.14/180 < phi[i] < (x2[j]+5)*3.14/180 ) {y2[j]+=1;}
        }
        ex2[j]=0;   
        ey2[j]=TMath::RMS(n_pair-1, phi); 
    }
    auto phi_gr = new TGraphErrors(35, x2, y2, ex2, ey2);
    phi_gr->SetTitle("TGraphErrors of phi");
    phi_gr->GetYaxis()->SetTitle("Number of particles");
    phi_gr->GetXaxis()->SetTitle("Phi");
    phi_gr->SetMarkerStyle(21);
    phi_gr->SetMarkerColor(2);
    phi_gr->Draw();
    legend->AddEntry(phi_gr, "Histogram \"Phi\"", "l");
    legend->AddEntry((TObject*)0, "", "");
    legend->AddEntry((TObject*)0, "Points", "lep");
    legend->Draw();

    s->Write();

    file.Write();
    file.Save();

    free(theta);
    free(phi);

    return 0;
}
