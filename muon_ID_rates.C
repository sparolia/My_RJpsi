#include <iostream>
#include <ostream>
#include <cmath>
#include <TMath.h>
#include <TRandom2.h>
#include <TH1.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <fstream>
#include <string>
#include <TFile.h>
#include <TArrayD.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <math.h>
#include <TLegend.h>



using namespace std;

void muon_ID_rates(){

  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain
  
  //Opening File
  TFile *f1 = TFile::Open("RootupleBcTo3Mu_pionChannel.root"); //opening the normalization root file for nn output
  TTree *tree = (TTree*)f1->Get("rootuple/ntuple"); //reading the tree from the ntuple
  //166209 events
  
  //Declaring leaf Variables
  //double Bc_mu_pt=0,Bc_mu_eta=0;
  float pion_pt_gen=0,pion_eta_gen=0,pion_phi_gen=0;
   
  vector<short> *truthMatchMu;
  vector<short> *isMuSoft;
  vector<short> *isMuTight;
  vector<short> *isMuLoose;
  vector<short> *isMuPF;
  vector<short> *isMuTracker;
  vector<short> *isMuGlobal;
  vector<short> *isMuMedium;
  vector<double> *Bc_mu_pt;
  vector<double> *Bc_mu_eta;
  vector<short> *triggerMatchDimuon0;
  vector<short> *triggerMatchJpsiTk;
  
  //TLorentzVector *gen_mu_p4;
  TLorentzVector *gen_pion_p4;
   
  //List of Branches
  TBranch *b_truthMatchMu;
  TBranch *b_isMuSoft;
  TBranch *b_isMuTight;
  TBranch *b_isMuLoose;
  TBranch *b_isMuPF;
  TBranch *b_isMuTracker;
  TBranch *b_isMuGlobal;
  TBranch *b_isMuMedium;
  TBranch *b_triggerMatchDimuon0;
  TBranch *b_triggerMatchJpsiTk;

  TBranch *b_Bc_mu_pt;
  TBranch *b_Bc_mu_eta;  
  //TBranch *b_gen_mu_p4;
  TBranch *b_gen_pion_p4;

  //Initializing values of the variable pointers
  truthMatchMu = 0;
  isMuSoft = 0;
  isMuTight = 0;
  isMuLoose = 0;
  isMuPF = 0;
  isMuTracker = 0;
  isMuGlobal = 0;
  isMuMedium = 0;
  Bc_mu_pt = 0;
  Bc_mu_eta = 0;
  triggerMatchDimuon0 = 0;
  triggerMatchJpsiTk = 0;
  //gen_mu_p4 = 0;
  gen_pion_p4 = 0;
 
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  
  //reading variables from the normalization tree for gen level
  //fChain->SetBranchAddress("gen_mu_p4",&gen_mu_p4,&b_gen_mu_p4);
  fChain->SetBranchAddress("gen_pion_p4",&gen_pion_p4,&b_gen_pion_p4); 

  tree->SetBranchAddress("Bc_mu_pt",&Bc_mu_pt,&b_Bc_mu_pt);
  tree->SetBranchAddress("Bc_mu_eta",&Bc_mu_eta,&b_Bc_mu_eta);
  
  fChain->SetBranchAddress("truthMatchMu",&truthMatchMu,&b_truthMatchMu);
  fChain->SetBranchAddress("isMuSoft",&isMuSoft,&b_isMuSoft);
  fChain->SetBranchAddress("isMuTight",&isMuTight,&b_isMuTight);
  fChain->SetBranchAddress("isMuLoose",&isMuLoose,&b_isMuLoose);
  fChain->SetBranchAddress("isMuPF",&isMuPF,&b_isMuPF);
  fChain->SetBranchAddress("isMuGlobal",&isMuGlobal,&b_isMuGlobal);
  fChain->SetBranchAddress("isMuTracker",&isMuTracker,&b_isMuTracker);
  fChain->SetBranchAddress("isMuMedium",&isMuMedium,&b_isMuMedium);
  fChain->SetBranchAddress("triggerMatchJpsiTk",&triggerMatchJpsiTk,&b_triggerMatchJpsiTk);
  fChain->SetBranchAddress("triggerMatchDimuon0",&triggerMatchDimuon0,&b_triggerMatchDimuon0);

  //Defining Canvas, removing statistics and setting the legend
  auto c1 = new TCanvas("c", "c", 500,500);
  //gStyle->SetOptStat(0);
  auto* legend = new TLegend(0.0, 0.7, 0.2, 0.9);


  TH2D* hmu_true_pion_pt_eta = new TH2D("hmu_true_pion_pt_eta","True Pt Vs Eta of Pion",60,-3,3,40,0,20);
  TH2D* hmu_recco_mu_pt_eta = new TH2D("hmu_recco_mu_pt_eta","Recco Pt Vs Eta of unpaired muon",60,-3,3,40,0,20);
  TH2D* hmu_mu_pt_eta_Soft = new TH2D("hmu_mu_pt_eta_Soft","Pt Vs Eta of Soft unpaired muon",60,-3,3,40,0,20);
  TH2D* hmu_mu_pt_eta_Tight = new TH2D("hmu_mu_pt_eta_Tight","Pt Vs Eta of Tight unpaired muon",60,-3,3,40,0,20);
  TH2D* hmu_mu_pt_eta_Loose = new TH2D("hmu_mu_pt_eta_Loose","Pt Vs Eta of Loose unpaired muon",60,-3,3,40,0,20);
  TH2D* hmu_mu_pt_eta_Global = new TH2D("hmu_mu_pt_eta_Global","Pt Vs Eta of Global unpaired muon",60,-3,3,40,0,20);
  TH2D* hmu_mu_pt_eta_PF = new TH2D("hmu_mu_pt_eta_PF","Pt Vs Eta of PF unpaired muon",60,-3,3,40,0,20);
  TH2D* hmu_mu_pt_eta_Tracker = new TH2D("hmu_mu_pt_eta_Tracker","Pt Vs Eta of Tracker unpaired muon",60,-3,3,40,0,20);
  TH2D* hmu_mu_pt_eta_Medium = new TH2D("hmu_mu_pt_eta_Medium","Pt Vs Eta of Medium unpaired muon",60,-3,3,40,0,20);
   
  int fake_num = 0,fake_den_acp = 0,fake_den_full=0;
  //count_true_pion=0;
  int MuSoft_num=0,MuTight_num=0,MuLoose_num=0,MuGlobal_num=0,MuPF_num=0,MuTracker_num=0,MuMedium_num=0;
    
  for(int i=0; i<tree->GetEntries();i++){
    tree->GetEntry(i);

    if((triggerMatchJpsiTk->at(0)) and (triggerMatchDimuon0->at(0))){
      fake_den_full +=1;
    
    pion_pt_gen = gen_pion_p4->Pt();
    pion_eta_gen = gen_pion_p4->Eta();
    pion_phi_gen = gen_pion_p4->Phi();

    int j_match = -1;
    for(int j=0;j < Bc_mu_pt->size();j++){
      
      if (truthMatchMu->at(j))
	j_match = j; 
    }

    if(j_match == -1)
      continue;
    fake_den_acp +=1;

    float recco_mu_pt = Bc_mu_pt->at(j_match);
    float recco_mu_eta = Bc_mu_eta->at(j_match);

    //cout << recco_mu_pt<<endl;
    if ((abs(recco_mu_eta) < 1.2 and recco_mu_pt > 3.0) or (1.2 < abs(recco_mu_eta) < 2.4 and recco_mu_pt > 0.8)){//acceptance cuts
      	
      //if (truthMatchMu->at(0)){
	fake_num +=1;
	
	  if(isMuSoft->at(j_match)){
	    MuSoft_num +=1;
	    hmu_mu_pt_eta_Soft->Fill(recco_mu_eta,recco_mu_pt);
	  }
	  if(isMuTight->at(j_match)){
	    MuTight_num +=1;
	    hmu_mu_pt_eta_Tight->Fill(recco_mu_eta,recco_mu_pt);
	  }
	  if(isMuLoose->at(j_match)){
	    MuLoose_num +=1;
	    hmu_mu_pt_eta_Loose->Fill(recco_mu_eta,recco_mu_pt);
	  }
	  if(isMuGlobal->at(j_match)){
	    MuGlobal_num +=1;
	    hmu_mu_pt_eta_Global->Fill(recco_mu_eta,recco_mu_pt);
	  }
	  if(isMuPF->at(j_match)){
	    MuPF_num +=1;
	    hmu_mu_pt_eta_PF->Fill(recco_mu_eta,recco_mu_pt);
	  }
	  if(isMuTracker->at(j_match)){
	    MuTracker_num +=1;
	    hmu_mu_pt_eta_Tracker->Fill(recco_mu_eta,recco_mu_pt);
	  }
	  if(isMuMedium->at(j_match)){
	    MuMedium_num +=1;
	    hmu_mu_pt_eta_Medium->Fill(recco_mu_eta,recco_mu_pt);
	  }
	  // }

      hmu_recco_mu_pt_eta->Fill(recco_mu_eta,recco_mu_pt);
      hmu_true_pion_pt_eta->Fill(pion_eta_gen,pion_pt_gen);
    }
    
     }
  }
  TH2D* hmu_fake_eff = (TH2D*)hmu_recco_mu_pt_eta->Clone();
  TH2D* hmu_mu_eff_Soft = (TH2D*)hmu_mu_pt_eta_Soft->Clone();
  TH2D* hmu_mu_eff_Tight = (TH2D*)hmu_mu_pt_eta_Tight->Clone();
  TH2D* hmu_mu_eff_Loose = (TH2D*)hmu_mu_pt_eta_Loose->Clone();
  TH2D* hmu_mu_eff_Global = (TH2D*)hmu_mu_pt_eta_Global->Clone();
  TH2D* hmu_mu_eff_PF = (TH2D*)hmu_mu_pt_eta_PF->Clone();
  TH2D* hmu_mu_eff_Tracker = (TH2D*)hmu_mu_pt_eta_Tracker->Clone();
  TH2D* hmu_mu_eff_Medium = (TH2D*)hmu_mu_pt_eta_Medium->Clone();

  hmu_fake_eff->Divide(hmu_true_pion_pt_eta);
  hmu_mu_eff_Soft->Divide(hmu_true_pion_pt_eta);
  hmu_mu_eff_Tight->Divide(hmu_true_pion_pt_eta);
  hmu_mu_eff_Loose->Divide(hmu_true_pion_pt_eta);
  hmu_mu_eff_Global->Divide(hmu_true_pion_pt_eta);
  hmu_mu_eff_PF->Divide(hmu_true_pion_pt_eta);
  hmu_mu_eff_Tracker->Divide(hmu_true_pion_pt_eta);
  hmu_mu_eff_Medium->Divide(hmu_true_pion_pt_eta);

  float eff_MuSoft = MuSoft_num/(fake_num*1.);
  float eff_MuTight = MuTight_num/(fake_num*1.);
  float eff_MuLoose = MuLoose_num/(fake_num*1.);
  float eff_MuGlobal = MuGlobal_num/(fake_num*1.);
  float eff_MuPF = MuPF_num/(fake_num*1.);
  float eff_MuTracker = MuTracker_num/(fake_num*1.);
  float eff_MuMedium = MuMedium_num/(fake_num*1.);

  float fake_rate = fake_num/(fake_den_full*1.);

  //cout<< MuSoft_num << endl;
  //cout<< MuTight_num << endl;
  
  
  cout<<"Eff_Soft = "<< eff_MuSoft << endl;
  cout<<"Eff_Medium = "<< eff_MuMedium << endl;
  cout<<"Eff_Tight = "<< eff_MuTight << endl;
  cout<<"Eff_Loose = "<< eff_MuLoose << endl;
  cout<<"Eff_PF = "<< eff_MuPF << endl;
  cout<<"Eff_Global = "<< eff_MuGlobal << endl;
  cout<<"Eff_Tracker = "<< eff_MuTracker << endl;
 
  cout<<"fake rate = "<<fake_rate<<endl;
  cout<<"Number of pions misidentified as muons = "<<fake_num<<endl;
  cout<<"Total number of events crossing the acceptance cuts = "<<fake_den_acp<<endl;
  cout<<"Total number of events = "<<fake_den_full<<endl;

  cout<<eff_MuSoft << endl;
  cout<<eff_MuMedium << endl;
  cout<<eff_MuTight << endl;
  cout<<eff_MuLoose << endl;
  cout<<eff_MuPF << endl;
  cout<<eff_MuGlobal << endl;
  cout<<eff_MuTracker << endl;
 
  cout<<fake_rate<<endl;
  cout<<fake_num<<endl;
  cout<<fake_den_acp<<endl;
  cout<<fake_den_full<<endl;

  /**
  //Pt_eta distribution for recco muons
  hmu_recco_mu_pt_eta->SetLineColor(kRed);
  hmu_recco_mu_pt_eta->GetXaxis()->SetTitle("Eta");
  hmu_recco_mu_pt_eta->GetYaxis()->SetTitle("Pt");
  hmu_recco_mu_pt_eta->Draw("colz");
  c1->SaveAs("recco_mu_pt_eta.png");

  //Pt_eta distribution for true pions
  hmu_true_pion_pt_eta->SetLineColor(kRed);
  hmu_true_pion_pt_eta->GetXaxis()->SetTitle("Eta");
  hmu_true_pion_pt_eta->GetYaxis()->SetTitle("Pt");
  hmu_true_pion_pt_eta->Draw("colz");
  c1->SaveAs("True_pion_pt_eta.png");

  //Efficiency for Soft muons
  hmu_mu_pt_eta_Soft->SetLineColor(kRed);
  hmu_mu_pt_eta_Soft->GetXaxis()->SetTitle("Efficiency");
  hmu_mu_pt_eta_Soft->GetYaxis()->SetTitle("a.u.");
  hmu_mu_pt_eta_Soft->Draw("colz");
  c1->SaveAs("mu_soft.png");

  //Efficiency for Soft muons
  hmu_mu_eff_Soft->SetLineColor(kRed);
  hmu_mu_eff_Soft->GetXaxis()->SetTitle("Efficiency");
  hmu_mu_eff_Soft->GetYaxis()->SetTitle("a.u.");
  hmu_mu_eff_Soft->Draw("colz");
  c1->SaveAs("mu_eff_soft.png");
  //Efficiency for Tight muons
  hmu_mu_eff_Tight->SetLineColor(kRed);
  hmu_mu_eff_Tight->GetXaxis()->SetTitle("Efficiency");
  hmu_mu_eff_Tight->GetYaxis()->SetTitle("a.u.");
  hmu_mu_eff_Tight->Draw("colz");
  c1->SaveAs("mu_eff_tight.png");
  //Efficiency for Loose muons
  hmu_mu_eff_Loose->SetLineColor(kRed);
  hmu_mu_eff_Loose->GetXaxis()->SetTitle("Efficiency");
  hmu_mu_eff_Loose->GetYaxis()->SetTitle("a.u.");
  hmu_mu_eff_Loose->Draw("colz");
  c1->SaveAs("mu_eff_loose.png");
  //Efficiency for Global muons
  hmu_mu_eff_Global->SetLineColor(kRed);
  hmu_mu_eff_Global->GetXaxis()->SetTitle("Efficiency");
  hmu_mu_eff_Global->GetYaxis()->SetTitle("a.u.");
  hmu_mu_eff_Global->Draw("colz");
  c1->SaveAs("mu_eff_global.png");
  //Efficiency for PF muons
  hmu_mu_eff_PF->SetLineColor(kRed);
  hmu_mu_eff_PF->GetXaxis()->SetTitle("Efficiency");
  hmu_mu_eff_PF->GetYaxis()->SetTitle("a.u.");
  hmu_mu_eff_PF->Draw("colz");
  c1->SaveAs("mu_eff_PF.png");
  //Efficiency for Tracker muons
  hmu_mu_eff_Tracker->SetLineColor(kRed);
  hmu_mu_eff_Tracker->GetXaxis()->SetTitle("Efficiency");
  hmu_mu_eff_Tracker->GetYaxis()->SetTitle("a.u.");
  hmu_mu_eff_Tracker->Draw("colz");
  c1->SaveAs("mu_eff_tracker.png");
  //Efficiency for Medium muons
  hmu_mu_eff_Medium->SetLineColor(kRed);
  hmu_mu_eff_Medium->GetXaxis()->SetTitle("Efficiency");
  hmu_mu_eff_Medium->GetYaxis()->SetTitle("a.u.");
  hmu_mu_eff_Medium->Draw("colz");
  c1->SaveAs("mu_eff_medium.png");
  **/
  f1->Close();
}
      
