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

void Macro_nn(){

  //declaring all the required variable for normalization channel
  //float bc_pt_norm_nn=0;
  float bc_eta_norm_nn=0,bc_phi_norm_nn=0,bc_mass_norm_nn=0,bc_true_pt_norm_nn=0;
  float mu_pt_norm_nn=0,mu_eta_norm_nn=0,mu_phi_norm_nn=0;
  float mu1_pt_norm_nn=0,mu1_eta_norm_nn=0,mu1_phi_norm_nn=0;
  float mu2_pt_norm_nn=0,mu2_eta_norm_nn=0,mu2_phi_norm_nn=0;
  float jpsi_pt_norm_nn=0,jpsi_eta_norm_nn=0,jpsi_phi_norm_nn=0,jpsi_mass_norm_nn=0;
  float Trig_dimuon_norm_nn=0,Trig_JpsiTk_norm_nn=0;
  float munu_pt_norm_nn=0,munu_eta_norm_nn=0,munu_phi_norm_nn=0;
  double bc_pt_norm_nn;

  //declaring all the required variable for signal channel
  //float bc_pt_sig_nn=0;
  float bc_eta_sig_nn=0,bc_phi_sig_nn=0,bc_mass_sig_nn=0,bc_true_pt_sig_nn=0;
  float mu_pt_sig_nn=0,mu_eta_sig_nn=0,mu_phi_sig_nn=0;
  float mu1_pt_sig_nn=0,mu1_eta_sig_nn=0,mu1_phi_sig_nn=0;
  float mu2_pt_sig_nn=0,mu2_eta_sig_nn=0,mu2_phi_sig_nn=0;
  float jpsi_pt_sig_nn=0,jpsi_eta_sig_nn=0,jpsi_phi_sig_nn=0,jpsi_mass_sig_nn=0;
  float Trig_dimuon_sig_nn=0,Trig_JpsiTk_sig_nn=0;
  float taunu_pt_sig_nn=0,taunu_eta_sig_nn=0,taunu_phi_sig_nn=0;
  float taunu2_pt_sig_nn=0,taunu2_eta_sig_nn=0,taunu2_phi_sig_nn=0;
  float munu_pt_sig_nn=0,munu_eta_sig_nn=0,munu_phi_sig_nn=0;
  double bc_pt_sig_nn;
  
  TFile *f1 = TFile::Open("results-muonChannel-nodes_80_60_40_20.root"); //opening the normalization root file
  TTree *tree_norm_nn = (TTree*)f1->Get("tree"); //reading the tree from the ntuple
  TFile *f2 = TFile::Open("results-tauChannel-nodes_80_60_40_20.root"); //opening the signal root file
  TTree *tree_sig_nn = (TTree*)f2->Get("tree"); //reading the tree from the ntuple


  //reading variables from the normalization tree
  //tree_norm_nn->SetBranchAddress("Bc_pt",&bc_pt_norm_nn); 
  tree_norm_nn_nn->SetBranchAddress("Bc_eta",&bc_eta_norm_nn); 
  tree_norm_nn->SetBranchAddress("Bc_phi",&bc_phi_norm_nn); 
  tree_norm_nn->SetBranchAddress("Bc_mass",&bc_mass_norm_nn); 

  tree_norm_nn->SetBranchAddress("Bc_mu_pt",&mu_pt_norm_nn); 
  tree_norm_nn->SetBranchAddress("Bc_mu_eta",&mu_eta_norm_nn); 
  tree_norm_nn->SetBranchAddress("Bc_mu_phi",&mu_phi_norm_nn);

  tree_norm_nn->SetBranchAddress("Bc_jpsi_mu1_pt",&mu1_pt_norm_nn); 
  tree_norm_nn->SetBranchAddress("Bc_jpsi_mu1_eta",&mu1_eta_norm_nn); 
  tree_norm_nn->SetBranchAddress("Bc_jpsi_mu1_phi",&mu1_phi_norm_nn);

  tree_norm_nn->SetBranchAddress("Bc_jpsi_mu2_pt",&mu2_pt_norm_nn); 
  tree_norm_nn->SetBranchAddress("Bc_jpsi_mu2_eta",&mu2_eta_norm_nn); 
  tree_norm_nn->SetBranchAddress("Bc_jpsi_mu2_phi",&mu2_phi_norm_nn);

  tree_norm_nn->SetBranchAddress("Bc_jpsi_pt",&jpsi_pt_norm_nn); 
  tree_norm_nn->SetBranchAddress("Bc_jpsi_eta",&jpsi_eta_norm_nn); 
  tree_norm_nn->SetBranchAddress("Bc_jpsi_phi",&jpsi_phi_norm_nn);
  tree_norm_nn->SetBranchAddress("Bc_jpsi_mass",&jpsi_mass_norm_nn);

  tree_norm_nn->SetBranchAddress("triggerMatchDimuon0",&Trig_dimuon_norm_nn);
  tree_norm_nn->SetBranchAddress("triggerMatchJpsiTk",&Trig_JpsiTk_norm_nn);

  //tree_norm_nn->SetBranchAddress("nn_predicted_pt",&bc_pt_predicted_norm_nn);
  tree_norm_nn->SetBranchAddress("bc_pt_predicted",&bc_pt_norm_nn);
  tree_norm_nn->SetBranchAddress("gen_b_pt",&bc_true_pt_norm_nn);
 

  //reading variables from the signal tree
  //tree_sig_nn->SetBranchAddress("Bc_pt",&bc_pt_sig_nn); 
  tree_sig_nn->SetBranchAddress("Bc_eta",&bc_eta_sig_nn); 
  tree_sig_nn->SetBranchAddress("Bc_phi",&bc_phi_sig_nn); 
  tree_sig_nn->SetBranchAddress("Bc_mass",&bc_mass_sig_nn);

  tree_sig_nn->SetBranchAddress("Bc_mu_pt",&mu_pt_sig_nn); 
  tree_sig_nn->SetBranchAddress("Bc_mu_eta",&mu_eta_sig_nn); 
  tree_sig_nn->SetBranchAddress("Bc_mu_phi",&mu_phi_sig_nn);

  tree_sig_nn->SetBranchAddress("Bc_jpsi_mu1_pt",&mu1_pt_sig_nn); 
  tree_sig_nn->SetBranchAddress("Bc_jpsi_mu1_eta",&mu1_eta_sig_nn); 
  tree_sig_nn->SetBranchAddress("Bc_jpsi_mu1_phi",&mu1_phi_sig_nn);

  tree_sig_nn->SetBranchAddress("Bc_jpsi_mu2_pt",&mu2_pt_sig_nn); 
  tree_sig_nn->SetBranchAddress("Bc_jpsi_mu2_eta",&mu2_eta_sig_nn); 
  tree_sig_nn->SetBranchAddress("Bc_jpsi_mu2_phi",&mu2_phi_sig_nn);

  tree_sig_nn->SetBranchAddress("Bc_jpsi_pt",&jpsi_pt_sig_nn); 
  tree_sig_nn->SetBranchAddress("Bc_jpsi_eta",&jpsi_eta_sig_nn); 
  tree_sig_nn->SetBranchAddress("Bc_jpsi_phi",&jpsi_phi_sig_nn); 
  tree_sig_nn->SetBranchAddress("Bc_jpsi_mass",&jpsi_mass_sig_nn);

  tree_sig_nn->SetBranchAddress("triggerMatchDimuon0",&Trig_dimuon_sig_nn);
  tree_sig_nn->SetBranchAddress("triggerMatchJpsiTk",&Trig_JpsiTk_sig_nn);
  
  //tree_sig_nn->SetBranchAddress("nn_predicted_pt",&bc_pt_predicted_sig_nn);
  tree_sig_nn->SetBranchAddress("bc_pt_predicted",&bc_pt_sig_nn);
  tree_sig_nn->SetBranchAddress("gen_b_pt",&bc_true_pt_sig_nn);
  

  //Defining Canvas, removing statistics and setting the legend
  auto c1 = new TCanvas("c", "c", 500,500);
  //gStyle->SetOptStat(0);
  auto* legend = new TLegend(0.7, 0.3, 0.9, 0.5);
   

  //booking the histograms
  TH1D* hmu_m2_miss_norm_nn = new TH1D("hmu_m2_miss_norm_nn","Missing mass squared",20,-10,10); 
  TH1D* hmu_m2_miss_sig_nn = new TH1D("hmu_m2_miss_sig_nn","Missing mass squared",20,-10,10);
  TH1D* hmu_Q2_norm_nn = new TH1D("hmu_Q2_norm_nn","Squared four momentum transfered to lepton system",50,-10,20);
  TH1D* hmu_Q2_sig_nn = new TH1D("hmu_Q2_sig_nn","Squared four momentum transfered to lepton system",50,-10,20);
  TH1D* hmu_Pt_miss_norm_nn = new TH1D("hmu_Pt_miss_norm_nn","Missing transverse momentum",120,-10,50);
  TH1D* hmu_Pt_miss_sig_nn = new TH1D("hmu_Pt_miss_sig_nn","Missing transverse momentum",120,-10,50);
  TH1D* hmu_Pt_var_norm_nn = new TH1D("hmu_Pt_var_norm_nn","Tranverse variable momentum",100,-30,50);
  TH1D* hmu_Pt_var_sig_nn = new TH1D("hmu_Pt_var_sig_nn","Tranverse variable momentum",100,-30,50);
  TH1D* hmu_del_R_norm_nn = new TH1D("hmu_del_R_norm_nn","Delta R",100,0,2);
  TH1D* hmu_del_R_sig_nn = new TH1D("hmu_del_R_sig_nn","Delta R",100,0,2);
  TH1D* hmu_mu_pt_norm_nn = new TH1D("hmu_mu_pt_norm_nn","Unpaired muon tranverse momentum",200,0,50);
  TH1D* hmu_mu_pt_sig_nn = new TH1D("hmu_mu_pt_sig_nn","Unpaired muon tranverse momentum",200,0,50);
  TH1D* hmu_mu_px_norm_nn = new TH1D("hmu_mu_px_norm_nn","Unpaired muon momentum in x-direction",200,0,50);
  TH1D* hmu_mu_px_sig_nn = new TH1D("hmu_mu_px_sig_nn","Unpaired muon momentum in x-direction",200,0,50);
  TH1D* hmu_mu_py_norm_nn = new TH1D("hmu_mu_py_norm_nn","Unpaired muon momentum in y-direction",200,0,50);
  TH1D* hmu_mu_py_sig_nn = new TH1D("hmu_mu_py_sig_nn","Unpaired muon momentum in y-direction",200,0,50);
  TH1D* hmu_mu_pz_norm_nn = new TH1D("hmu_mu_pz_norm_nn","Unpaired muon momentum in z-direction",200,0,50);
  TH1D* hmu_mu_pz_sig_nn = new TH1D("hmu_mu_pz_sig_nn","Unpaired muon momentum in z-direction",200,0,50);
  TH1D* hmu_mu_eta_norm_nn = new TH1D("hmu_mu_eta_norm_nn","Unpaired muon psuedorapidity",90,-3,3);
  TH1D* hmu_mu_eta_sig_nn = new TH1D("hmu_mu_eta_sig_nn","Unpaired muon psuedorapidity",90,-3,3);
  TH1D* hmu_mu_phi_norm_nn = new TH1D("hmu_mu_phi_norm_nn","Unpaired muon phi",100,-5,5);
  TH1D* hmu_mu_phi_sig_nn = new TH1D("hmu_mu_phi_sig_nn","Unpaired muon phi",100,-5,5);
  TH1D* hmu_mu_en_bc_frame_norm_nn = new TH1D("hmu_mu_en_bc_frame_norm_nn","Energy of Unpaired muon in Bc rest frame",200,0,20);
  TH1D* hmu_mu_en_bc_frame_sig_nn = new TH1D("hmu_mu_en_bc_frame_sig_nn","Energy of Unpaired muon in Bc rest frame",200,0,20);
  TH1D* hmu_mu_en_jpsi_frame_norm_nn = new TH1D("hmu_mu_en_jpsi_frame_norm_nn","Energy of Unpaired muon in Jpsi rest frame",200,0,20);
  TH1D* hmu_mu_en_jpsi_frame_sig_nn = new TH1D("hmu_mu_en_jpsi_frame_sig_nn","Energy of Unpaired muon in Jpsi rest frame",200,0,20);
  TH1D* hmu_bc_pt_norm_nn = new TH1D("hmu_bc_pt_norm_nn","Transverse momentum of B_{c}",100,0,100);
  TH1D* hmu_bc_pt_sig_nn = new TH1D("hmu_bc_pt_sig_nn","Transverse momentum of B_{c}",100,0,100);
  TH1D* hmu_bc_pt_ratio_norm_nn = new TH1D("hmu_bc_pt_ratio_norm_nn","Ratio of true and predicted_nn transverse momentum of B_{c}",100,0,3);
  TH1D* hmu_bc_pt_ratio_sig_nn = new TH1D("hmu_bc_pt_ratio_sig_nn","Ratio of true and predicted_nn transverse momentum of B_{c}",100,0,3);
  auto hprof_norm_nn  = new TProfile("hprof_norm_nn","Profile of Transverse momentum of B_{c}",100,0,100,0,1000);
  auto hprof_sig_nn  = new TProfile("hprof_sig_nn","Profile of Transverse momentum of B_{c}",100,0,100,0,1000);
  TH2D* hmu_bc_pt_scatter_norm_nn = new TH2D("hmu_bc_pt_scatter_norm_nn","Scatter plot of transverse momentum of B_{c}",100,0,100,100,0,100);
  TH2D* hmu_bc_pt_scatter_sig_nn = new TH2D("hmu_bc_pt_scatter_sig_nn","Scatter plot of transverse momentum of B_{c}",100,0,100,100,0,100);
  TH1D* hmu_nn_bc_pt_norm_nn = new TH1D("hmu_nn_bc_pt_norm_nn","Predicted transverse momentum of B_{c}",100,0,100);
  TH1D* hmu_nn_bc_pt_sig_nn = new TH1D("hmu_nn_bc_pt_sig_nn","Predicted transverse momentum of B_{c}",100,0,100);
  
 


  //ParticleMass muonMass = 0.10565837;
  float mu_mass = 0.113428;
  float bcPdgMass = 6.2756;

  //loop over events for norm_nnalized channel
  for(int i=0; i<tree_norm_nn->GetEntries();i++){
    tree_norm_nn->GetEntry(i);
    
    //TLorentzVector class for calculating four momenta
    TLorentzVector momntm_bc_norm_nn;
    momntm_bc_norm_nn.SetPtEtaPhiM(bc_pt_norm_nn, bc_eta_norm_nn,bc_phi_norm_nn,bcPdgMass);
    TLorentzVector momntm_mu1_norm_nn;
    momntm_mu1_norm_nn.SetPtEtaPhiM(mu1_pt_norm_nn, mu1_eta_norm_nn,mu1_phi_norm_nn,mu_mass);
    TLorentzVector momntm_mu2_norm_nn;
    momntm_mu2_norm_nn.SetPtEtaPhiM(mu2_pt_norm_nn, mu2_eta_norm_nn,mu2_phi_norm_nn,mu_mass);
    TLorentzVector momntm_mu_norm_nn;
    momntm_mu_norm_nn.SetPtEtaPhiM(mu_pt_norm_nn, mu_eta_norm_nn,mu_phi_norm_nn,mu_mass);
    TLorentzVector momntm_jpsi_norm_nn;
    momntm_jpsi_norm_nn.SetPtEtaPhiM(jpsi_pt_norm_nn, jpsi_eta_norm_nn,jpsi_phi_norm_nn,jpsi_mass_norm_nn);

    //calculating the kinematic variables
    float m2_miss_norm_nn =(momntm_bc_norm_nn - momntm_mu1_norm_nn - momntm_mu2_norm_nn - momntm_mu_norm_nn).M2();
    float Q2_norm_nn = (momntm_bc_norm_nn - momntm_mu1_norm_nn - momntm_mu2_norm_nn).M2();
    float Pt_miss_norm_nn = (bc_pt_norm_nn - mu1_pt_norm_nn - mu2_pt_norm_nn - mu_pt_norm_nn);
    float Pt_var_norm_nn = (jpsi_pt_norm_nn - mu_pt_norm_nn );
    float del_R_norm_nn = momntm_mu1_norm_nn.DeltaR(momntm_mu2_norm_nn);
    float bc_pt_ratio_norm_nn = bc_pt_norm_nn/bc_true_pt_norm_nn;
    
    TVector3 X_bc_norm_nn = -momntm_bc_corr_norm_nn.BoostVector();
    TVector3 X_jpsi_norm_nn = -momntm_jpsi_norm_nn.BoostVector();
    TLorentzVector momntm_mu_bc_frame_norm_nn = momntm_mu_norm_nn;
    momntm_mu_bc_frame_norm_nn.Boost(X_bc_norm_nn);
    TLorentzVector momntm_mu_jpsi_frame_norm_nn = momntm_mu_norm_nn;
    momntm_mu_jpsi_frame_norm_nn.Boost(X_jpsi_norm_nn);

    float mu_en_bc_frame_norm_nn = momntm_mu_bc_frame_norm_nn.E();
    float mu_en_jpsi_frame_norm_nn = momntm_mu_jpsi_frame_norm_nn.E();

    float mu_px_norm_nn = momntm_mu_norm_nn.Px();
    float mu_py_norm_nn = momntm_mu_norm_nn.Py();
    float mu_pz_norm_nn = momntm_mu_norm_nn.Pz();
        
    //Filling the histograms

    //if ((abs(mu_eta_norm_nn) < 1.3 and mu_pt_norm_nn > 3.3) or (1.3 < abs(mu_eta_norm_nn) < 2.2 and mu_pt_norm_nn > 2.9) or (2.2 < abs(mu_eta_norm_nn) < 2.4 and mu_pt_norm_nn > 2.4)){ //acceptance cuts
      //if ((Trig_dimuon_norm_nn == false) and (Trig_JpsiTk_norm_nn == true)){ //trigger flags
	hmu_m2_miss_norm_nn->Fill(m2_miss_norm_nn);
	hmu_Q2_norm_nn->Fill(Q2_norm_nn);
	hmu_Pt_miss_norm_nn->Fill(Pt_miss_norm_nn);
	hmu_Pt_var_norm_nn->Fill(Pt_var_norm_nn);
	hmu_mu_pt_norm_nn->Fill(mu_pt_norm_nn);
	hmu_mu_px_norm_nn->Fill(mu_px_norm_nn);
	hmu_mu_py_norm_nn->Fill(mu_py_norm_nn);
	hmu_mu_pz_norm_nn->Fill(mu_pz_norm_nn);
	hmu_mu_eta_norm_nn->Fill(mu_eta_norm_nn);
	hmu_mu_phi_norm_nn->Fill(mu_phi_norm_nn);
	hmu_del_R_norm_nn->Fill(del_R_norm_nn);
	hmu_mu_en_bc_frame_norm_nn->Fill(mu_en_bc_frame_norm_nn);
	hmu_mu_en_jpsi_frame_norm_nn->Fill(mu_en_jpsi_frame_norm_nn);
	hmu_bc_pt_norm_nn->Fill(bc_pt_norm_nn);
	hmu_bc_pt_ratio_norm_nn->Fill(bc_pt_ratio_norm_nn);
	hprof_norm_nn->Fill(bc_pt_norm_nn,bc_pt_ratio_norm_nn,1);
	hmu_bc_pt_scatter_norm_nn->Fill(bc_true_pt_norm_nn,bc_pt_norm_nn);
	//}
	//}
  }

  //loop over events for sig_nnnal channel
  for(int i=0; i<tree_sig_nn->GetEntries();i++){
    tree_sig_nn->GetEntry(i);

    //TLorentzVector class for calculating four momenta
    TLorentzVector momntm_bc_sig_nn;
    momntm_bc_sig_nn.SetPtEtaPhiM(bc_pt_sig_nn, bc_eta_sig_nn,bc_phi_sig_nn,bcPdgMass);
    TLorentzVector momntm_mu1_sig_nn;
    momntm_mu1_sig_nn.SetPtEtaPhiM(mu1_pt_sig_nn, mu1_eta_sig_nn,mu1_phi_sig_nn,mu_mass);
    TLorentzVector momntm_mu2_sig_nn;
    momntm_mu2_sig_nn.SetPtEtaPhiM(mu2_pt_sig_nn, mu2_eta_sig_nn,mu2_phi_sig_nn,mu_mass);
    TLorentzVector momntm_mu_sig_nn;
    momntm_mu_sig_nn.SetPtEtaPhiM(mu_pt_sig_nn, mu_eta_sig_nn,mu_phi_sig_nn,mu_mass);
    TLorentzVector momntm_jpsi_sig_nn;
    momntm_jpsi_sig_nn.SetPtEtaPhiM(jpsi_pt_sig_nn, jpsi_eta_sig_nn,jpsi_phi_sig_nn,jpsi_mass_sig_nn);

    //calculating the kinematic variables
    float m2_miss_sig_nn = (momntm_bc_sig_nn - momntm_mu1_sig_nn - momntm_mu2_sig_nn - momntm_mu_sig_nn).M2();
    float Q2_sig_nn = (momntm_bc_sig_nn - momntm_mu1_sig_nn - momntm_mu2_sig_nn).M2();
    float Pt_miss_sig_nn = (bc_pt_sig_nn-mu1_pt_sig_nn-mu2_pt_sig_nn-mu_pt_sig_nn);
    float Pt_var_sig_nn = (jpsi_pt_sig_nn - mu_pt_sig_nn );
    float del_R_sig_nn = momntm_mu1_sig_nn.DeltaR(momntm_mu2_sig_nn);
   
    // float bc_pt_ratio_sig_nn = bc_pt_corr_sig_nn/bc_pt_sig_nn;

    TVector3 X_bc_sig_nn = -momntm_bc_corr_sig_nn.BoostVector();
    TVector3 X_jpsi_sig_nn = -momntm_jpsi_sig_nn.BoostVector();
    TLorentzVector momntm_mu_bc_frame_sig_nn = momntm_mu_sig_nn;
    momntm_mu_bc_frame_sig_nn.Boost(X_bc_sig_nn);
    TLorentzVector momntm_mu_jpsi_frame_sig_nn = momntm_mu_sig_nn;
    momntm_mu_jpsi_frame_sig_nn.Boost(X_jpsi_sig_nn);
    
    float mu_en_bc_frame_sig_nn = momntm_mu_bc_frame_sig_nn.E();
    float mu_en_jpsi_frame_sig_nn = momntm_mu_jpsi_frame_sig_nn.E();

    float mu_px_sig_nn = momntm_mu_sig_nn.Px();
    float mu_py_sig_nn = momntm_mu_sig_nn.Py();
    float mu_pz_sig_nn = momntm_mu_sig_nn.Pz();
    
    //Filling the histograms
    //if ((abs(mu_eta_sig_nn) < 1.3 and mu_pt_sig_nn > 3.3) or (1.3 < abs(mu_eta_sig_nn) < 2.2 and mu_pt_sig_nn > 2.9) or (2.2 < abs(mu_eta_sig_nn) < 2.4 and mu_pt_sig_nn > 2.4)){ //acceptance cuts
      //if ((Trig_dimuon_sig_nn == false) and (Trig_JpsiTk_sig_nn == true)){ //trigger flags
	hmu_m2_miss_sig_nn->Fill(m2_miss_sig_nn);
	hmu_Q2_sig_nn->Fill(Q2_sig_nn);
	hmu_Pt_miss_sig_nn->Fill(Pt_miss_sig_nn);
	hmu_Pt_var_sig_nn->Fill(Pt_var_sig_nn);
	hmu_mu_pt_sig_nn->Fill(mu_pt_sig_nn);
	hmu_mu_px_sig_nn->Fill(mu_px_sig_nn);
	hmu_mu_py_sig_nn->Fill(mu_py_sig_nn);
	hmu_mu_pz_sig_nn->Fill(mu_pz_sig_nn);
	hmu_mu_eta_sig_nn->Fill(mu_eta_sig_nn);
	hmu_mu_phi_sig_nn->Fill(mu_phi_sig_nn);
	hmu_del_R_sig_nn->Fill(del_R_sig_nn);
	hmu_mu_en_bc_frame_sig_nn->Fill(mu_en_bc_frame_sig_nn);
	hmu_mu_en_jpsi_frame_sig_nn->Fill(mu_en_jpsi_frame_sig_nn);
	hmu_bc_pt_sig_nn->Fill(bc_pt_sig_nn);
	hprof_sig_nn->Fill(bc_pt_sig_nn,bc_pt_ratio_sig_nn,1);
	hmu_bc_pt_scatter_sig_nn->Fill(bc_true_pt_sig_nn,bc_pt_sig_nn);
	//}
	//}
  }  
  
  //Drawing Histograms

  
  //Missing mass Squared
  //float norm_nn = hmu_m2_miss_norm_nn->GetEntries();
  hmu_m2_miss_norm_nn->SetLineColor(kRed);
  hmu_m2_miss_norm_nn->GetXaxis()->SetTitle("m^{2}_{miss}[GeV^{2}]");
  hmu_m2_miss_norm_nn->GetYaxis()->SetTitle("a.u.");
  hmu_m2_miss_norm_nn->DrawNorm_Nnalized();
  hmu_m2_miss_sig_nn->SetLineColor(kBlue);
  hmu_m2_miss_sig_nn->DrawNorm_Nnalized("same");
  legend->AddEntry(hmu_m2_miss_norm_nn, "#mu channel", "l");
  legend->AddEntry(hmu_m2_miss_sig_nn, "#tau channel", "l");
  legend->Draw();
  c1->SaveAs("M2_miss.png");

  
  //Squared four momentum transfered to lepton system
  hmu_Q2_norm_nn->SetLineColor(kRed);
  hmu_Q2_norm_nn->GetXaxis()->SetTitle("Q^{2}[GeV^{2}]");
  hmu_Q2_norm_nn->GetYaxis()->SetTitle("a.u.");
  hmu_Q2_norm_nn->DrawNorm_Nnalized();
  hmu_Q2_sig_nn->SetLineColor(kBlue);
  hmu_Q2_sig_nn->DrawNorm_Nnalized("same");
  legend->Draw();
  c1->SaveAs("Q2.png");
  

  
  //Missing transverse momentum
  hmu_Pt_miss_norm_nn->SetLineColor(kRed);
  hmu_Pt_miss_norm_nn->GetXaxis()->SetTitle("P_{T}^{miss}[GeV]");
  hmu_Pt_miss_norm_nn->GetYaxis()->SetTitle("a.u.");
  hmu_Pt_miss_norm_nn->DrawNorm_Nnalized();
  hmu_Pt_miss_sig_nn->SetLineColor(kBlue);
  hmu_Pt_miss_sig_nn->DrawNorm_Nnalized("same");
  legend->Draw();
  c1->SaveAs("Pt_miss.png");
  

  //Tranverse variable momentum
  hmu_Pt_var_norm_nn->SetLineColor(kRed);
  hmu_Pt_var_norm_nn->GetXaxis()->SetTitle("P_{T}^{var}[GeV]");
  hmu_Pt_var_norm_nn->GetYaxis()->SetTitle("a.u.");
  hmu_Pt_var_norm_nn->DrawNorm_Nnalized();
  hmu_Pt_var_sig_nn->SetLineColor(kBlue);
  hmu_Pt_var_sig_nn->DrawNorm_Nnalized("same");
  legend->Draw();
  c1->SaveAs("Pt_var.png");
  
  //Delta R
  hmu_del_R_norm_nn->SetLineColor(kRed);
  hmu_del_R_norm_nn->GetXaxis()->SetTitle("#Delta R(#mu_{1},#mu_{2})");
  hmu_del_R_norm_nn->GetYaxis()->SetTitle("a.u.");
  hmu_del_R_norm_nn->DrawNorm_Nnalized();
  hmu_del_R_sig_nn->SetLineColor(kBlue);
  hmu_del_R_sig_nn->DrawNorm_Nnalized("same");
  legend->Draw();
  c1->SaveAs("del_R.png");


  //Unpaired muon tranverse momentum
  hmu_mu_pt_norm_nn->SetLineColor(kRed);
  hmu_mu_pt_norm_nn->GetXaxis()->SetTitle("P_{T}^{#mu}[GeV]");
  hmu_mu_pt_norm_nn->GetYaxis()->SetTitle("a.u.");
  hmu_mu_pt_norm_nn->DrawNorm_Nnalized();
  hmu_mu_pt_sig_nn->SetLineColor(kBlue);
  hmu_mu_pt_sig_nn->DrawNorm_Nnalized("same");
  legend->Draw();
  c1->SaveAs("mu_pt.png");


  //Unpaired muon momentum in x-direction
  hmu_mu_px_norm_nn->SetLineColor(kRed);
  hmu_mu_px_norm_nn->GetXaxis()->SetTitle("P_{x}^{#mu}[GeV]");
  hmu_mu_px_norm_nn->GetYaxis()->SetTitle("a.u.");
  hmu_mu_px_norm_nn->DrawNorm_Nnalized();
  hmu_mu_px_sig_nn->SetLineColor(kBlue);
  hmu_mu_px_sig_nn->DrawNorm_Nnalized("same");
  legend->Draw();
  c1->SaveAs("mu_px.png");


  //Unpaired muon momentum in y-direction
  hmu_mu_py_norm_nn->SetLineColor(kRed);
  hmu_mu_py_norm_nn->GetXaxis()->SetTitle("P_{y}^{#mu}[GeV]");
  hmu_mu_py_norm_nn->GetYaxis()->SetTitle("a.u.");
  hmu_mu_py_norm_nn->DrawNorm_Nnalized();
  hmu_mu_py_sig_nn->SetLineColor(kBlue);
  hmu_mu_py_sig_nn->DrawNorm_Nnalized("same");
  legend->Draw();
  c1->SaveAs("mu_py.png");


  //Unpaired muon momentum in z-direction
  hmu_mu_pz_norm_nn->SetLineColor(kRed);
  hmu_mu_pz_norm_nn->GetXaxis()->SetTitle("P_{z}^{#mu}[GeV]");
  hmu_mu_pz_norm_nn->GetYaxis()->SetTitle("a.u.");
  hmu_mu_pz_norm_nn->DrawNorm_Nnalized();
  hmu_mu_pz_sig_nn->SetLineColor(kBlue);
  hmu_mu_pz_sig_nn->DrawNorm_Nnalized("same");
  legend->Draw();
  c1->SaveAs("mu_pz.png");

  
  //Unpaired muon psuedorapidity
  hmu_mu_eta_norm_nn->SetLineColor(kRed);
  hmu_mu_eta_norm_nn->GetXaxis()->SetTitle("#eta^{#mu}");
  hmu_mu_eta_norm_nn->GetYaxis()->SetTitle("a.u.");
  hmu_mu_eta_norm_nn->DrawNorm_Nnalized();
  hmu_mu_eta_sig_nn->SetLineColor(kBlue);
  hmu_mu_eta_sig_nn->DrawNorm_Nnalized("same");
  legend->Draw();
  c1->SaveAs("mu_eta.png");

  
  //Unpaired muon phi
  hmu_mu_phi_norm_nn->SetLineColor(kRed);
  hmu_mu_phi_norm_nn->GetXaxis()->SetTitle("#phi^{#mu}");
  hmu_mu_phi_norm_nn->GetYaxis()->SetTitle("a.u.");
  hmu_mu_phi_norm_nn->DrawNorm_Nnalized();
  hmu_mu_phi_sig_nn->SetLineColor(kBlue);
  hmu_mu_phi_sig_nn->DrawNorm_Nnalized("same");
  legend->Draw();
  c1->SaveAs("mu_phi.png");

  
  //Energy of unpaired muon in Bc rest frame
  hmu_mu_en_bc_frame_norm_nn->SetLineColor(kRed);
  hmu_mu_en_bc_frame_norm_nn->GetXaxis()->SetTitle("E(#mu) in rest frame of B_{c} (GeV)");
  hmu_mu_en_bc_frame_norm_nn->GetYaxis()->SetTitle("a.u.");
  hmu_mu_en_bc_frame_norm_nn->DrawNorm_Nnalized();
  hmu_mu_en_bc_frame_sig_nn->SetLineColor(kBlue);
  hmu_mu_en_bc_frame_sig_nn->DrawNorm_Nnalized("same");
  legend->Draw();
  c1->SaveAs("E_bc.png");

  //Energy of unpaired muon in Jpsi rest frame
  hmu_mu_en_jpsi_frame_norm_nn->SetLineColor(kRed);
  hmu_mu_en_jpsi_frame_norm_nn->GetXaxis()->SetTitle("E(#mu) in rest frame of J/(#Psi) (GeV)");
  hmu_mu_en_jpsi_frame_norm_nn->GetYaxis()->SetTitle("a.u.");
  hmu_mu_en_jpsi_frame_norm_nn->DrawNorm_Nnalized();
  hmu_mu_en_jpsi_frame_sig_nn->SetLineColor(kBlue);
  hmu_mu_en_jpsi_frame_sig_nn->DrawNorm_Nnalized("same");
  legend->Draw();
  c1->SaveAs("E_jpsi.png");\

  
  //Transverse momentum of Bc
  hmu_bc_pt_norm_nn->SetLineColor(kRed);
  hmu_bc_pt_norm_nn->GetXaxis()->SetTitle("P_{T}^{Bc} (GeV)");
  hmu_bc_pt_norm_nn->GetYaxis()->SetTitle("a.u.");
  hmu_bc_pt_norm_nn->DrawNorm_Nnalized();
  hmu_bc_pt_sig_nn->SetLineColor(kBlue);
  hmu_bc_pt_sig_nn->DrawNorm_Nnalized("same");
  legend->Draw();
  c1->SaveAs("Pt_Bc.png");

  
  //Ratio of true and corrected transverse momentum of Bc
  hmu_bc_pt_ratio_norm_nn->SetLineColor(kRed);
  hmu_bc_pt_ratio_norm_nn->GetXaxis()->SetTitle("Ratio of true and corrected P_{T}^{Bc} (GeV)");
  hmu_bc_pt_ratio_norm_nn->GetYaxis()->SetTitle("a.u.");
  hmu_bc_pt_ratio_norm_nn->DrawNorm_Nnalized();
  hmu_bc_pt_ratio_sig_nn->SetLineColor(kBlue);
  hmu_bc_pt_ratio_sig_nn->DrawNorm_Nnalized("same");
  legend->Draw();
  c1->SaveAs("Ratio_pt_bc.png");
  

  //Profile histogram of transverse momentum of Bc
  hprof_norm_nn->SetLineColor(kRed);
  hprof_norm_nn->GetXaxis()->SetTitle("P_{T}^{Bc} (GeV)");
  hprof_norm_nn->GetYaxis()->SetTitle("Ratio of reconstructed and corrected P_{T}^{Bc}");
  hprof_norm_nn->Draw();
  hprof_sig_nn->SetLineColor(kBlue);
  hprof_sig_nn->Draw("same");
  legend->Draw();
  c1->SaveAs("Profile_pt_bc.png");
  
  
  //Comparison of Reconstructed Bc Pt, Jona recipe and prediction from neural network for norm_nnalization channel
  hmu_bc_pt_norm_nn->SetLineColor(kRed);
  hmu_bc_pt_norm_nn->GetXaxis()->SetTitle("P_{T}^{Bc} (GeV)");
  hmu_bc_pt_norm_nn->GetYaxis()->SetTitle("a.u.");
  hmu_bc_pt_norm_nn->DrawNorm_Nnalized();
  hmu_nn_bc_pt_norm_nn->SetLineColor(kGreen);
  hmu_nn_bc_pt_norm_nn->DrawNorm_Nnalized("same");
  legend->AddEntry(hmu_nn_bc_pt_norm_nn, "Predicted", "l");
  legend->Draw();
  c1->SaveAs("Compare_recco_pt_bc_norm_nn.png");

  //Comparison of Reconstructed Bc Pt, Jona recipe and prediction from neural network for sig_nnnal channel
  //Comparison of Bc Pt and Reconstructed Bc Pt
  hmu_bc_pt_sig_nn->SetLineColor(kBlue);
  hmu_bc_pt_sig_nn->GetXaxis()->SetTitle("P_{T}^{Bc} (GeV)");
  hmu_bc_pt_sig_nn->GetYaxis()->SetTitle("a.u.");
  hmu_bc_pt_sig_nn->DrawNorm_Nnalized();
  hmu_nn_bc_pt_sig_nn->SetLineColor(kGreen);
  hmu_nn_bc_pt_sig_nn->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("Compare_recco_pt_bc_sig_nn.png");


  //Closing the files
  f1->Close();
  f2->Close();
  
}
