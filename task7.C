#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <TTree.h>
#include <TFile.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

using namespace std;

int task7(){
    TFile file("m3pimc.root","read");
    TTree *tree = (TTree *)file.Get("h10");
    
    Int_t nph;
    UChar_t Isrfilter;
    Float_t chi2_3p;
    tree->SetBranchAddress("nph", &nph);
    tree->SetBranchAddress("Isrfilter", &Isrfilter);
    tree->SetBranchAddress("chi2_3p", &chi2_3p);

    TFile* newroot = new TFile("newroot.root", "RECREATE");
    TTree* MyTree = new TTree("MyTree", "MyTree");

    MyTree->Branch("nph", &nph);

    file.cd();

    Float_t eph[nph];
    Float_t thetaph[nph];
    Float_t phiph[nph];  
    tree->SetBranchAddress("eph", eph);
    tree->SetBranchAddress("thetaph", thetaph);
    tree->SetBranchAddress("phiph", phiph);

    newroot->cd();

    MyTree->Branch("eph", eph, "eph[nph]/F");
    MyTree->Branch("thetaph", thetaph, "thetaph[nph]/F");
    MyTree->Branch("phiph", phiph, "phiph[nph]/F");
    
    TCanvas *s = new TCanvas();
    TH1D* heph= new TH1D("heph","Energy of photons",100,0.,10.);

    Long64_t nentries = tree->GetEntries();
    for (Long64_t ientry = 0; ientry < nentries; ++ientry) {
        
        tree->GetEntry(ientry);
        if (Isrfilter == 1 && chi2_3p < 30) {
            MyTree->Fill();
            for (Long64_t i = 0; i < nph; ++i){
                heph->Fill(eph[i]);
            }
        }
    }
    heph->GetYaxis()->SetTitle("Nunber of entries");
    heph->GetXaxis()->SetTitle("Energy of photons");
    //heph->SetLogy();
    s->SetLogy();
    s->Draw();

    TF1 *func = new TF1("func","500*[0]/x",0., 10.);
    heph->TH1D::Fit(func);
    
    newroot->Write();
    newroot->Save();

    cout<< "Maximum energy of photon: "<< MyTree->TTree::GetMaximum("eph") << endl;
    cout<< "Minimum energy of photon: "<< MyTree->TTree::GetMinimum("eph") << endl;

    //file.Print();
    //newroot->Print();
    cout<<"New file entries:"<< newroot-> GetSize() << endl;
    cout<<"Old file entries:"<< file.GetSize() << endl;

    newroot->Close();
    file.Close();

    //logarifm scale 
    return 0;
}
