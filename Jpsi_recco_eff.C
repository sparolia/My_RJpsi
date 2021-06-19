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
#include <vector>

using namespace std;

void Jpsi_recco_eff(){

  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain
  
  //Opening File
  TFile *f1 = TFile::Open("RootupleBcTo3Mu_muonChannel.root"); //opening the normalization root file for nn output
  TTree *tree = (TTree*)f1->Get("rootuple/ntuple"); //reading the tree from the ntuple
  //166209 events
  
  //Declaring leaf Variables
   
  vector<short> *truthMatchMu1;
  vector<short> *isMu1Soft;
  vector<short> *isMu1Tight;
  vector<short> *isMu1Loose;
  vector<short> *isMu1PF;
  vector<short> *isMu1Tracker;
  vector<short> *isMu1Global;
  vector<short> *isMu1Medium;
  vector<double> *Bc_jpsi_mu1_pt;
  vector<double> *Bc_jpsi_mu1_eta;

  vector<short> *truthMatchMu2;
  vector<short> *isMu2Soft;
  vector<short> *isMu2Tight;
  vector<short> *isMu2Loose;
  vector<short> *isMu2PF;
  vector<short> *isMu2Tracker;
  vector<short> *isMu2Global;
  vector<short> *isMu2Medium;
  vector<double> *Bc_jpsi_mu2_pt;
  vector<double> *Bc_jpsi_mu2_eta;

  vector<double> *Bc_jpsi_mass;
  vector<short> *triggerMatchDimuon0;
  vector<short> *triggerMatchJpsiTk;
  
   
  //List of Branches
  TBranch *b_truthMatchMu1;
  TBranch *b_isMu1Soft;
  TBranch *b_isMu1Tight;
  TBranch *b_isMu1Loose;
  TBranch *b_isMu1PF;
  TBranch *b_isMu1Tracker;
  TBranch *b_isMu1Global;
  TBranch *b_isMu1Medium;
  TBranch *b_Bc_jpsi_mu1_pt;
  TBranch *b_Bc_jpsi_mu1_eta;

  TBranch *b_truthMatchMu2;
  TBranch *b_isMu2Soft;
  TBranch *b_isMu2Tight;
  TBranch *b_isMu2Loose;
  TBranch *b_isMu2PF;
  TBranch *b_isMu2Tracker;
  TBranch *b_isMu2Global;
  TBranch *b_isMu2Medium;
  TBranch *b_Bc_jpsi_mu2_pt;
  TBranch *b_Bc_jpsi_mu2_eta;

  TBranch *b_Bc_jpsi_mass;
  TBranch *b_triggerMatchDimuon0;
  TBranch *b_triggerMatchJpsiTk;
  

  //Initializing values of the variable pointers
  truthMatchMu1 = 0;
  isMu1Soft = 0;
  isMu1Tight = 0;
  isMu1Loose = 0;
  isMu1PF = 0;
  isMu1Tracker = 0;
  isMu1Global = 0;
  isMu1Medium = 0;
  Bc_jpsi_mu1_pt = 0;
  Bc_jpsi_mu1_eta = 0;

  truthMatchMu2 = 0;
  isMu2Soft = 0;
  isMu2Tight = 0;
  isMu2Loose = 0;
  isMu2PF = 0;
  isMu2Tracker = 0;
  isMu2Global = 0;
  isMu2Medium = 0;
  Bc_jpsi_mu2_pt = 0;
  Bc_jpsi_mu2_eta = 0;

  Bc_jpsi_mass = 0;
  triggerMatchDimuon0 = 0;
  triggerMatchJpsiTk = 0;
  
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  
  //reading variables from the normalization tree for gen level
  
  fChain->SetBranchAddress("truthMatchMu1",&truthMatchMu1,&b_truthMatchMu1);
  fChain->SetBranchAddress("isMu1Soft",&isMu1Soft,&b_isMu1Soft);
  fChain->SetBranchAddress("isMu1Tight",&isMu1Tight,&b_isMu1Tight);
  fChain->SetBranchAddress("isMu1Loose",&isMu1Loose,&b_isMu1Loose);
  fChain->SetBranchAddress("isMu1PF",&isMu1PF,&b_isMu1PF);
  fChain->SetBranchAddress("isMu1Global",&isMu1Global,&b_isMu1Global);
  fChain->SetBranchAddress("isMu1Medium",&isMu1Medium,&b_isMu1Medium);
  fChain->SetBranchAddress("isMu1Tracker",&isMu1Tracker,&b_isMu1Tracker);
  fChain->SetBranchAddress("Bc_jpsi_mu1_pt",&Bc_jpsi_mu1_pt,&b_Bc_jpsi_mu1_pt);
  fChain->SetBranchAddress("Bc_jpsi_mu1_eta",&Bc_jpsi_mu1_eta,&b_Bc_jpsi_mu1_eta);
  
  fChain->SetBranchAddress("truthMatchMu2",&truthMatchMu2,&b_truthMatchMu2);
  fChain->SetBranchAddress("isMu2Soft",&isMu2Soft,&b_isMu2Soft);
  fChain->SetBranchAddress("isMu2Tight",&isMu2Tight,&b_isMu2Tight);
  fChain->SetBranchAddress("isMu2Loose",&isMu2Loose,&b_isMu2Loose);
  fChain->SetBranchAddress("isMu2PF",&isMu2PF,&b_isMu2PF);
  fChain->SetBranchAddress("isMu2Global",&isMu2Global,&b_isMu2Global);
  fChain->SetBranchAddress("isMu2Medium",&isMu2Medium,&b_isMu2Medium);
  fChain->SetBranchAddress("isMu2Tracker",&isMu2Tracker,&b_isMu2Tracker);
  fChain->SetBranchAddress("Bc_jpsi_mu2_pt",&Bc_jpsi_mu2_pt,&b_Bc_jpsi_mu2_pt);
  fChain->SetBranchAddress("Bc_jpsi_mu2_eta",&Bc_jpsi_mu2_eta,&b_Bc_jpsi_mu2_eta);
  
  fChain->SetBranchAddress("triggerMatchJpsiTk",&triggerMatchJpsiTk,&b_triggerMatchJpsiTk);
  fChain->SetBranchAddress("triggerMatchDimuon0",&triggerMatchDimuon0,&b_triggerMatchDimuon0);
  fChain->SetBranchAddress("Bc_jpsi_mass",&Bc_jpsi_mass,&b_Bc_jpsi_mass);
  

  //Defining Canvas, removing statistics and setting the legend
  auto c1 = new TCanvas("c", "c", 500,500);
  //gStyle->SetOptStat(0);
  auto* legend = new TLegend(0.0, 0.7, 0.2, 0.9);

  TH2D* hmu_recco_mu_pt_eta = new TH2D("hmu_recco_mu_pt_eta","Recco Pt Vs Eta of unpaired muon",60,-3,3,40,0,20);
  TH2D* hmu_recco_mu_high_pt_eta = new TH2D("hmu_recco_mu_high_pt_eta","Recco Pt Vs Eta of unpaired high Pt muon",60,-3,3,40,0,20);
  TH2D* hmu_recco_mu_low_pt_eta = new TH2D("hmu_recco_mu_low_pt_eta","Recco Pt Vs Eta of unpaired low Pt muon",60,-3,3,40,0,20);
  TH1F* hmu_jpsi_events_ntrig = new TH1F("hmu_jpsi_events_ntrig","Number of reconstructed Jpsi per event without trigger", 6,0,6);
  TH1F* hmu_jpsi_events_trig = new TH1F("hmu_jpsi_events_trig","Number of reconstructed Jpsi per event with trigger", 6,0,6);
  
  int tot = 0, tot_accp = 0, tot_trig = 0, tot_truth = 0;
  //count_true_pion=0;
  int Soft_num=0, Tight_num=0, Loose_num=0, Global_num=0, PF_num=0, Tracker_num=0, Medium_num=0;
    
  for(int i=0; i<tree->GetEntries();i++){
    tree->GetEntry(i);
    
    //std::vector<double> *jpsi = Bc_jpsi_mass;
    //set<double> s( Bc_jpsi_mass->begin(), Bc_jpsi_mass->end() );
    //Bc_jpsi_mass->assign( s.begin(), s.end() );
    int jpsi_events = 0;
    //cout<<Bc_jpsi_mass->size()<<"*"<<endl;
    //std::sort(Bc_jpsi_mass->begin(), Bc_jpsi_mass->end());
    Bc_jpsi_mass->erase(std::unique(Bc_jpsi_mass->begin(), Bc_jpsi_mass->end()), Bc_jpsi_mass->end());
    //cout<<Bc_jpsi_mass->size()<<endl;
    
    for(int k =0; k<Bc_jpsi_mass->size(); k++){
      if(2.9 < Bc_jpsi_mass->at(k) < 3.1)
	jpsi_events +=1;
    }
    
    
    int j_match = -1;
    for(int j=0; j<Bc_jpsi_mass->size();j++){
      if (truthMatchMu1->at(j) and truthMatchMu2->at(j))
	j_match = j; 
    }

    if(j_match == -1)
      continue;
    tot +=1;
     
    float recco_mu1_pt = Bc_jpsi_mu1_pt->at(j_match);
    float recco_mu1_eta = Bc_jpsi_mu1_eta->at(j_match);

    float recco_mu2_pt = Bc_jpsi_mu2_pt->at(j_match);
    float recco_mu2_eta = Bc_jpsi_mu2_eta->at(j_match);

    if (((abs(recco_mu1_eta) < 1.2 and recco_mu1_pt > 3.0) or (1.2 < abs(recco_mu1_eta) < 2.4 and recco_mu1_pt > 0.8)) and ((abs(recco_mu2_eta) < 1.2 and recco_mu2_pt > 3.0) or (1.2 < abs(recco_mu2_eta) < 2.4 and recco_mu2_pt > 0.8))) {//acceptance cuts
      tot_accp +=1;
      
      if((triggerMatchJpsiTk->at(0)) and !(triggerMatchDimuon0->at(0))){
	tot_trig +=1;
    
	if(isMu1Soft->at(j_match) and isMu2Soft->at(j_match))
	  //if(isMu1Soft->at(j_match)){
	  Soft_num +=1;
	
	if(isMu1Tight->at(j_match) and isMu2Tight->at(j_match))
	  Tight_num +=1;
	
	if(isMu1Loose->at(j_match) and isMu2Loose->at(j_match))
	  Loose_num +=1;

	if(isMu1Global->at(j_match) and isMu2Global->at(j_match))
	  Global_num +=1;

	if(isMu1PF->at(j_match) and isMu2PF->at(j_match))
	  PF_num +=1;

	if(isMu1Tracker->at(j_match) and isMu2Tracker->at(j_match))
	  Tracker_num +=1;

	if(isMu1Medium->at(j_match) and isMu2Medium->at(j_match))
	  Medium_num +=1;

	hmu_jpsi_events_trig->Fill(jpsi_events);
    
      }

      hmu_recco_mu_pt_eta->Fill(recco_mu1_eta,recco_mu1_pt);
      hmu_recco_mu_pt_eta->Fill(recco_mu2_eta,recco_mu2_pt);
      if(recco_mu1_pt > recco_mu2_pt){
	hmu_recco_mu_high_pt_eta->Fill(recco_mu1_eta,recco_mu1_pt);
	hmu_recco_mu_low_pt_eta->Fill(recco_mu2_eta,recco_mu2_pt);
      }
      else{
	hmu_recco_mu_high_pt_eta->Fill(recco_mu2_eta,recco_mu2_pt);
	hmu_recco_mu_low_pt_eta->Fill(recco_mu1_eta,recco_mu1_pt);
      }
      hmu_jpsi_events_ntrig->Fill(jpsi_events);
      
    }
  }
  
  float eff_Soft = Soft_num/(tot_trig*1.);
  float eff_Tight = Tight_num/(tot_trig*1.);
  float eff_Loose = Loose_num/(tot_trig*1.);
  float eff_Global = Global_num/(tot_trig*1.);
  float eff_PF = PF_num/(tot_trig*1.);
  float eff_Tracker = Tracker_num/(tot_trig*1.);
  float eff_Medium = Medium_num/(tot_trig*1.);
  
  cout<< tot <<endl;
  cout<< tot_accp <<endl;
  cout<< tot_trig <<endl;
  cout<< eff_Soft << endl;
  cout<< eff_Medium << endl;
  cout<< eff_Tight << endl;
  cout<< eff_Loose << endl;
  cout<< eff_PF << endl;
  cout<< eff_Global << endl;
  cout<< eff_Tracker << endl;

  //Pt_eta distribution for recco muons
  hmu_recco_mu_pt_eta->SetLineColor(kRed);
  hmu_recco_mu_pt_eta->GetXaxis()->SetTitle("Eta");
  hmu_recco_mu_pt_eta->GetYaxis()->SetTitle("Pt");
  hmu_recco_mu_pt_eta->Draw("colz");
  c1->SaveAs("recco_mu_pt_eta.png");

  //Pt_eta distribution for high Pt recco muons
  hmu_recco_mu_high_pt_eta->SetLineColor(kRed);
  hmu_recco_mu_high_pt_eta->GetXaxis()->SetTitle("Eta");
  hmu_recco_mu_high_pt_eta->GetYaxis()->SetTitle("Pt");
  hmu_recco_mu_high_pt_eta->Draw("colz");
  c1->SaveAs("recco_mu_high_pt_eta.png");

  //Pt_eta distribution for low Pt recco muons
  hmu_recco_mu_low_pt_eta->SetLineColor(kRed);
  hmu_recco_mu_low_pt_eta->GetXaxis()->SetTitle("Eta");
  hmu_recco_mu_low_pt_eta->GetYaxis()->SetTitle("Pt");
  hmu_recco_mu_low_pt_eta->Draw("colz");
  c1->SaveAs("recco_mu_low_pt_eta.png");
  
  //Number of recco Jpsi per event
  hmu_jpsi_events_ntrig->SetLineColor(kRed);
  hmu_jpsi_events_ntrig->GetXaxis()->SetTitle("Events");
  hmu_jpsi_events_ntrig->GetYaxis()->SetTitle("Frequency");
  hmu_jpsi_events_ntrig->Draw();
  c1->SaveAs("Jpsi_frequency_ntrig.png");

  //Number of recco Jpsi per event
  hmu_jpsi_events_trig->SetLineColor(kRed);
  hmu_jpsi_events_trig->GetXaxis()->SetTitle("Events");
  hmu_jpsi_events_trig->GetYaxis()->SetTitle("Frequency");
  hmu_jpsi_events_trig->Draw();
  c1->SaveAs("Jpsi_frequency_trig.png");
  
  f1->Close();
}
      
