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

void muon_ID_rates_sig_norm(){

  //Opening File
  TFile *f1 = TFile::Open("results-tau_channel-nodes_100_100_100_60-vfrac_0p3-dropoutrate_0p1-bsize_200.root"); //opening the signal root file for nn output
  TTree *tree = (TTree*)f1->Get("tree"); //reading the tree from the ntuple
  //166209 events
  
  //Declaring leaf Variables
  float Bc_mu_pt=0,Bc_mu_eta=0;
  float gen_mu_pt=0,gen_mu_eta=0;
    
  float truthMatchMu=0;
  float isMuSoft=0,isMuTight=0,isMuLoose=0;
  float isMuPF=0,isMuTracker=0,isMuGlobal=0;
  float triggerMatchDimuon0=0,triggerMatchJpsiTk=0;

  
  //reading variables from the normalization tree for gen level
  tree->SetBranchAddress("gen_mu_pt",&gen_mu_pt);
  tree->SetBranchAddress("gen_mu_Eta",&gen_mu_eta);
  tree->SetBranchAddress("Bc_mu_pt",&Bc_mu_pt);
  tree->SetBranchAddress("Bc_mu_eta",&Bc_mu_eta);
  
  tree->SetBranchAddress("truthMatchMu",&truthMatchMu);
  tree->SetBranchAddress("isMuSoft",&isMuSoft);
  tree->SetBranchAddress("isMuTight",&isMuTight);
  tree->SetBranchAddress("isMuLoose",&isMuLoose);
  tree->SetBranchAddress("isMuPF",&isMuPF);
  tree->SetBranchAddress("isMuGlobal",&isMuGlobal);
  tree->SetBranchAddress("isMuTracker",&isMuTracker);
  tree->SetBranchAddress("triggerMatchJpsiTk",&triggerMatchJpsiTk);
  tree->SetBranchAddress("triggerMatchDimuon0",&triggerMatchDimuon0);

  //Defining Canvas, removing statistics and setting the legend
  auto c1 = new TCanvas("c", "c", 600,600);
  //gStyle->SetOptStat(0);
  auto* legend = new TLegend(0.0, 0.7, 0.2, 0.9);


  TH2D* hmu_true_muon_pt_eta = new TH2D("hmu_true_muon_pt_eta","True Pt Vs Eta of Muon",60,-3,3,40,0,20);
  TH2D* hmu_recco_mu_pt_eta = new TH2D("hmu_recco_mu_pt_eta","Recco Pt Vs Eta of unpaired muon",60,-3,3,40,0,20);
  TH2D* hmu_mu_pt_eta_Soft = new TH2D("hmu_mu_pt_eta_Soft","Pt Vs Eta of Soft unpaired muon",60,-3,3,40,0,20);
  TH2D* hmu_mu_pt_eta_Tight = new TH2D("hmu_mu_pt_eta_Tight","Pt Vs Eta of Tight unpaired muon",60,-3,3,40,0,20);
  TH2D* hmu_mu_pt_eta_Loose = new TH2D("hmu_mu_pt_eta_Loose","Pt Vs Eta of Loose unpaired muon",60,-3,3,40,0,20);
  TH2D* hmu_mu_pt_eta_Global = new TH2D("hmu_mu_pt_eta_Global","Pt Vs Eta of Global unpaired muon",60,-3,3,40,0,20);
  TH2D* hmu_mu_pt_eta_PF = new TH2D("hmu_mu_pt_eta_PF","Pt Vs Eta of PF unpaired muon",60,-3,3,40,0,20);
  TH2D* hmu_mu_pt_eta_Tracker = new TH2D("hmu_mu_pt_eta_Tracker","Pt Vs Eta of Tracker unpaired muon",60,-3,3,40,0,20);
   
  float fake_num = 0,fake_den_full=0,fake_den_acp = 0;
  //count_true_pion=0;
  int MuSoft_num=0,MuTight_num=0,MuLoose_num=0,MuGlobal_num=0,MuPF_num=0,MuTracker_num=0;
    
  for(int i=0; i<tree->GetEntries();i++){
    tree->GetEntry(i);

    if(triggerMatchJpsiTk or triggerMatchDimuon0){
    fake_den_full +=1;
    
    float recco_mu_pt = Bc_mu_pt;
    float recco_mu_eta = Bc_mu_eta;

    //cout << recco_mu_pt<<endl;
    if ((abs(recco_mu_eta) < 1.2 and recco_mu_pt > 3.0) or (1.2 < abs(recco_mu_eta) < 2.4 and recco_mu_pt > 0.8)){//acceptance cuts
      fake_den_acp +=1;
      if (truthMatchMu){
	fake_num +=1;
	//cout<<isMuSoft->at(0)<<endl;
	if(isMuSoft){
	    MuSoft_num +=1;
	    //cout<<recco_mu_pt<<endl;
	    hmu_mu_pt_eta_Soft->Fill(recco_mu_eta,recco_mu_pt);
	  }
	  if(isMuTight){
	    MuTight_num +=1;
	    hmu_mu_pt_eta_Tight->Fill(recco_mu_eta,recco_mu_pt);
	  }
	  if(isMuLoose){
	    MuLoose_num +=1;
	    hmu_mu_pt_eta_Loose->Fill(recco_mu_eta,recco_mu_pt);
	  }
	  if(isMuGlobal){
	    MuGlobal_num +=1;
	    hmu_mu_pt_eta_Global->Fill(recco_mu_eta,recco_mu_pt);
	  }
	  if(isMuPF){
	    MuPF_num +=1;
	    hmu_mu_pt_eta_PF->Fill(recco_mu_eta,recco_mu_pt);
	  }
	  if(isMuTracker){
	    MuTracker_num +=1;
	    hmu_mu_pt_eta_Tracker->Fill(recco_mu_eta,recco_mu_pt);
	  }
      }

      hmu_recco_mu_pt_eta->Fill(recco_mu_eta,recco_mu_pt);
      hmu_true_muon_pt_eta->Fill(gen_mu_eta,gen_mu_pt);
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

  hmu_fake_eff->Divide(hmu_true_muon_pt_eta);
  hmu_mu_eff_Soft->Divide(hmu_true_muon_pt_eta);
  hmu_mu_eff_Tight->Divide(hmu_true_muon_pt_eta);
  hmu_mu_eff_Loose->Divide(hmu_true_muon_pt_eta);
  hmu_mu_eff_Global->Divide(hmu_true_muon_pt_eta);
  hmu_mu_eff_PF->Divide(hmu_true_muon_pt_eta);
  hmu_mu_eff_Tracker->Divide(hmu_true_muon_pt_eta);

  float eff_MuSoft = MuSoft_num/fake_num;
  float eff_MuTight = MuTight_num/fake_num;
  float eff_MuLoose = MuLoose_num/fake_num;
  float eff_MuGlobal = MuGlobal_num/fake_num;
  float eff_MuPF = MuPF_num/fake_num;
  float eff_MuTracker = MuTracker_num/fake_num;

  float fake_rate = fake_num/fake_den_acp;

  cout<< eff_MuSoft << endl;
  cout<< eff_MuTight << endl;
  cout<< eff_MuLoose << endl;
  cout<< eff_MuPF << endl;
  cout<< eff_MuGlobal << endl;
  cout<< eff_MuTracker << endl;
  cout<<fake_rate<<endl;
  cout<<fake_num<<endl;
  cout<<fake_den_acp<<endl;
  cout<<fake_den_full<<endl;

  //Pt_eta distribution for recco muons
  hmu_recco_mu_pt_eta->SetLineColor(kRed);
  hmu_recco_mu_pt_eta->GetXaxis()->SetTitle("Eta");
  hmu_recco_mu_pt_eta->GetYaxis()->SetTitle("Pt");
  hmu_recco_mu_pt_eta->Draw("colz");
  c1->SaveAs("recco_mu_pt_eta.png");

  //Pt_eta distribution for true pions
  hmu_true_muon_pt_eta->SetLineColor(kRed);
  hmu_true_muon_pt_eta->GetXaxis()->SetTitle("Eta");
  hmu_true_muon_pt_eta->GetYaxis()->SetTitle("Pt");
  hmu_true_muon_pt_eta->Draw("colz");
  c1->SaveAs("True_muon_pt_eta.png");

  
  //Efficeincy for Soft muons
  hmu_mu_eff_Soft->SetLineColor(kRed);
  hmu_mu_eff_Soft->GetXaxis()->SetTitle("Efficiency");
  hmu_mu_eff_Soft->GetYaxis()->SetTitle("a.u.");
  hmu_mu_eff_Soft->Draw("colz");
  c1->SaveAs("mu_eff_soft.png");
  //Efficeincy for Tight muons
  hmu_mu_eff_Tight->SetLineColor(kRed);
  hmu_mu_eff_Tight->GetXaxis()->SetTitle("Efficiency");
  hmu_mu_eff_Tight->GetYaxis()->SetTitle("a.u.");
  hmu_mu_eff_Tight->Draw("colz");
  c1->SaveAs("mu_eff_tight.png");
  //Efficeincy for Loose muons
  hmu_mu_eff_Loose->SetLineColor(kRed);
  hmu_mu_eff_Loose->GetXaxis()->SetTitle("Efficiency");
  hmu_mu_eff_Loose->GetYaxis()->SetTitle("a.u.");
  hmu_mu_eff_Loose->Draw("colz");
  c1->SaveAs("mu_eff_loose.png");
  //Efficeincy for Global muons
  hmu_mu_eff_Global->SetLineColor(kRed);
  hmu_mu_eff_Global->GetXaxis()->SetTitle("Efficiency");
  hmu_mu_eff_Global->GetYaxis()->SetTitle("a.u.");
  hmu_mu_eff_Global->Draw("colz");
  c1->SaveAs("mu_eff_global.png");
  //Efficeincy for PF muons
  hmu_mu_eff_PF->SetLineColor(kRed);
  hmu_mu_eff_PF->GetXaxis()->SetTitle("Efficiency");
  hmu_mu_eff_PF->GetYaxis()->SetTitle("a.u.");
  hmu_mu_eff_PF->Draw("colz");
  c1->SaveAs("mu_eff_PF.png");
  //Efficeincy for Tracker muons
  hmu_mu_eff_Tracker->SetLineColor(kRed);
  hmu_mu_eff_Tracker->GetXaxis()->SetTitle("Efficiency");
  hmu_mu_eff_Tracker->GetYaxis()->SetTitle("a.u.");
  hmu_mu_eff_Tracker->Draw("colz");
  c1->SaveAs("mu_eff_tracker.png");

  f1->Close();
}
      
