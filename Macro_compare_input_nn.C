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

void Macro_compare_input_nn(){

  //declaring all the required variable for signal channel
  double bc_pt_sig=0,bc_eta_sig=0,bc_phi_sig=0,bc_mass_sig=0;
  float mu_pt_sig=0,mu_eta_sig=0,mu_phi_sig=0;
  float mu1_pt_sig=0,mu1_eta_sig=0,mu1_phi_sig=0;
  float mu2_pt_sig=0,mu2_eta_sig=0,mu2_phi_sig=0;
  float jpsi_pt_sig=0,jpsi_eta_sig=0,jpsi_phi_sig=0,jpsi_mass_sig=0;
  float Trig_dimuon_sig=0,Trig_JpsiTk_sig=0;
  float taunu_pt_sig=0,taunu_eta_sig=0,taunu_phi_sig=0;
  float taunu2_pt_sig=0,taunu2_eta_sig=0,taunu2_phi_sig=0;
  float munu_pt_sig=0,munu_eta_sig=0,munu_phi_sig=0;

  /**
  //declaring all the required variable for signal 2017 channel
  float bc_pt_sig_17=0,bc_eta_sig_17=0,bc_phi_sig_17=0,bc_mass_sig_17=0;
  float mu_pt_sig_17=0,mu_eta_sig_17=0,mu_phi_sig_17=0;
  float mu1_pt_sig_17=0,mu1_eta_sig_17=0,mu1_phi_sig_17=0;
  float mu2_pt_sig_17=0,mu2_eta_sig_17=0,mu2_phi_sig_17=0;
  float jpsi_pt_sig_17=0,jpsi_eta_sig_17=0,jpsi_phi_sig_17=0,jpsi_mass_sig_17=0;
  float Trig_dimuon_sig_17=0,Trig_JpsiTk_sig_17=0;
  float taunu_pt_sig_17=0,taunu_eta_sig_17=0,taunu_phi_sig_17=0;
  float taunu2_pt_sig_17=0,taunu2_eta_sig_17=0,taunu2_phi_sig_17=0;
  float munu_pt_sig_17=0,munu_eta_sig_17=0,munu_phi_sig_17=0;
  **/

  //declaring all the required variable for signal 2018 channel
  float bc_pt_sig_18=0,bc_eta_sig_18=0,bc_phi_sig_18=0,bc_mass_sig_18=0;
  float mu_pt_sig_18=0,mu_eta_sig_18=0,mu_phi_sig_18=0;
  float mu1_pt_sig_18=0,mu1_eta_sig_18=0,mu1_phi_sig_18=0;
  float mu2_pt_sig_18=0,mu2_eta_sig_18=0,mu2_phi_sig_18=0;
  float jpsi_pt_sig_18=0,jpsi_eta_sig_18=0,jpsi_phi_sig_18=0,jpsi_mass_sig_18=0;
  float Trig_dimuon_sig_18=0,Trig_JpsiTk_sig_18=0;
  float taunu_pt_sig_18=0,taunu_eta_sig_18=0,taunu_phi_sig_18=0;
  float taunu2_pt_sig_18=0,taunu2_eta_sig_18=0,taunu2_phi_sig_18=0;
  float munu_pt_sig_18=0,munu_eta_sig_18=0,munu_phi_sig_18=0;
  
  TFile *f1 = TFile::Open("RootupleBcTo3Mu_tauChannel_old.root"); //opening the signal root file of old data
  //TNtuple *f1_ntuple = (TNtuple*)f1->Get("rootuple");
  TTree *tree_sig = (TTree*)f1->Get("rootuple/ntuple"); //reading the tree from the ntuple
  /**
  TFile *f2 = TFile::Open("bc_jpsi_tau_nu_gen_2017.root"); //opening the signal root file of 2017 data
  TTree *tree_sig_17 = (TTree*)f2->Get("tree"); //reading the tree from the ntuple
  **/
  TFile *f3 = TFile::Open("RootupleBcTo3Mu_tauChannel_2018.root"); //opening the signal root file of 2018 data
  //TNtuple *f3_ntuple = (TNtuple*)f3->Get("rootuple");
  TTree *tree_sig_18 = (TTree*)f3->Get("rootuple/ntuple"); //reading the tree from the ntuple



  //reading variables from the signal tree
  tree_sig->SetBranchAddress("Bc_pt",&bc_pt_sig); 
  tree_sig->SetBranchAddress("Bc_eta",&bc_eta_sig); 
  tree_sig->SetBranchAddress("Bc_phi",&bc_phi_sig); 
  tree_sig->SetBranchAddress("Bc_mass",&bc_mass_sig);

  tree_sig->SetBranchAddress("Bc_mu_pt",&mu_pt_sig); 
  tree_sig->SetBranchAddress("Bc_mu_eta",&mu_eta_sig); 
  tree_sig->SetBranchAddress("mu_phi",&mu_phi_sig);

  tree_sig->SetBranchAddress("Bc_jpsi_mu1_pt",&mu1_pt_sig); 
  tree_sig->SetBranchAddress("Bc_jpsi_mu1_eta",&mu1_eta_sig); 
  tree_sig->SetBranchAddress("Bc_jpsi_mu1_phi",&mu1_phi_sig);

  tree_sig->SetBranchAddress("Bc_jpsi_mu2_pt",&mu2_pt_sig); 
  tree_sig->SetBranchAddress("Bc_jpsi_mu2_eta",&mu2_eta_sig); 
  tree_sig->SetBranchAddress("Bc_jpsi_mu2_phi",&mu2_phi_sig);

  tree_sig->SetBranchAddress("Bc_jpsi_pt",&jpsi_pt_sig); 
  tree_sig->SetBranchAddress("Bc_jpsi_eta",&jpsi_eta_sig); 
  tree_sig->SetBranchAddress("Bc_jpsi_phi",&jpsi_phi_sig); 
  tree_sig->SetBranchAddress("Bc_jpsi_mass",&jpsi_mass_sig);

  tree_sig->SetBranchAddress("triggerMatchDimuon0",&Trig_dimuon_sig);
  tree_sig->SetBranchAddress("triggerMatchJpsiTk",&Trig_JpsiTk_sig);

  /**
  //reading variables from the signal 2017 tree
  tree_sig_17->SetBranchAddress("bc_pt",&bc_pt_sig_17); 
  tree_sig_17->SetBranchAddress("bc_eta",&bc_eta_sig_17); 
  tree_sig_17->SetBranchAddress("bc_phi",&bc_phi_sig_17); 
  tree_sig_17->SetBranchAddress("bc_mass",&bc_mass_sig_17); 

  tree_sig_17->SetBranchAddress("mu_pt",&mu_pt_sig_17); 
  tree_sig_17->SetBranchAddress("mu_eta",&mu_eta_sig_17); 
  tree_sig_17->SetBranchAddress("mu_phi",&mu_phi_sig_17);

  tree_sig_17->SetBranchAddress("mu1_pt",&mu1_pt_sig_17); 
  tree_sig_17->SetBranchAddress("mu1_eta",&mu1_eta_sig_17); 
  tree_sig_17->SetBranchAddress("mu1_phi",&mu1_phi_sig_17);

  tree_sig_17->SetBranchAddress("mu2_pt",&mu2_pt_sig_17); 
  tree_sig_17->SetBranchAddress("mu2_eta",&mu2_eta_sig_17); 
  tree_sig_17->SetBranchAddress("mu2_phi",&mu2_phi_sig_17);

  tree_sig_17->SetBranchAddress("jpsi_pt",&jpsi_pt_sig_17); 
  tree_sig_17->SetBranchAddress("jpsi_eta",&jpsi_eta_sig_17); 
  tree_sig_17->SetBranchAddress("jpsi_phi",&jpsi_phi_sig_17);
  tree_sig_17->SetBranchAddress("jpsi_mass",&jpsi_mass_sig_17);

  tree_sig_17->SetBranchAddress("trigDimuon0",&Trig_dimuon_sig_17);
  tree_sig_17->SetBranchAddress("trigJpsiTk",&Trig_JpsiTk_sig_17);

  tree_sig_17->SetBranchAddress("munu_pt",&munu_pt_sig_17); 
  tree_sig_17->SetBranchAddress("munu_eta",&munu_eta_sig_17); 
  tree_sig_17->SetBranchAddress("munu_phi",&munu_phi_sig_17);
  **/

  //reading variables from the signal 2018 tree
  tree_sig_18->SetBranchAddress("Bc_pt",&bc_pt_sig_18); 
  tree_sig_18->SetBranchAddress("Bc_eta",&bc_eta_sig_18); 
  tree_sig_18->SetBranchAddress("Bc_phi",&bc_phi_sig_18); 
  tree_sig_18->SetBranchAddress("Bc_mass",&bc_mass_sig_18);

  tree_sig_18->SetBranchAddress("Bc_mu_pt",&mu_pt_sig_18); 
  tree_sig_18->SetBranchAddress("Bc_mu_eta",&mu_eta_sig_18); 
  tree_sig_18->SetBranchAddress("Bc_mu_phi",&mu_phi_sig_18);

  tree_sig_18->SetBranchAddress("Bc_jpsi_mu1_pt",&mu1_pt_sig_18); 
  tree_sig_18->SetBranchAddress("Bc_jpsi_mu1_eta",&mu1_eta_sig_18); 
  tree_sig_18->SetBranchAddress("Bc_jpsi_mu1_phi",&mu1_phi_sig_18);

  tree_sig_18->SetBranchAddress("Bc_jpsi_mu2_pt",&mu2_pt_sig_18); 
  tree_sig_18->SetBranchAddress("Bc_jpsi_mu2_eta",&mu2_eta_sig_18); 
  tree_sig_18->SetBranchAddress("Bc_jpsi_mu2_phi",&mu2_phi_sig_18);

  tree_sig_18->SetBranchAddress("Bc_jpsi_pt",&jpsi_pt_sig_18); 
  tree_sig_18->SetBranchAddress("Bc_jpsi_eta",&jpsi_eta_sig_18); 
  tree_sig_18->SetBranchAddress("Bc_jpsi_phi",&jpsi_phi_sig_18); 
  tree_sig_18->SetBranchAddress("Bc_jpsi_mass",&jpsi_mass_sig_18);

  tree_sig_18->SetBranchAddress("triggerMatchDimuon0",&Trig_dimuon_sig_18);
  tree_sig_18->SetBranchAddress("triggerMatchJpsiTk",&Trig_JpsiTk_sig_18);

  //Defining Canvas, removing statistics and setting the legend
  auto c1 = new TCanvas("c", "c", 600,600);
  //gStyle->SetOptStat(0);
  auto* legend = new TLegend(0.7, 0.3, 0.9, 0.5);
   

  //booking the histograms
  TH1D* hmu_m2_miss_sig_17 = new TH1D("hmu_m2_miss_sig_17","Missing mass squared",80,0,8);
  TH1D* hmu_m2_miss_sig_18 = new TH1D("hmu_m2_miss_sig_18","Missing mass squared",80,0,8); 
  TH1D* hmu_m2_miss_sig = new TH1D("hmu_m2_miss_sig","Missing mass squared",80,0,8);
  TH1D* hmu_Q2_sig_17 = new TH1D("hmu_Q2_sig_17","Squared four momentum transfered to lepton system",200,0,20);
  TH1D* hmu_Q2_sig_18 = new TH1D("hmu_Q2_sig_18","Squared four momentum transfered to lepton system",200,0,20);
  TH1D* hmu_Q2_sig = new TH1D("hmu_Q2_sig","Squared four momentum transfered to lepton system",200,0,20);
  TH1D* hmu_Pt_miss_sig_17 = new TH1D("hmu_Pt_miss_sig_17","Missing transverse momentum",200,-5,25);
  TH1D* hmu_Pt_miss_sig_18 = new TH1D("hmu_Pt_miss_sig_18","Missing transverse momentum",200,-5,25);
  TH1D* hmu_Pt_miss_sig = new TH1D("hmu_Pt_miss_sig","Missing transverse momentum",200,-5,25);
  TH1D* hmu_m2_miss_corr_sig_17 = new TH1D("hmu_m2_miss_corr_sig_17","Missing mass squared corrected",150,0,15);
  TH1D* hmu_m2_miss_corr_sig_18 = new TH1D("hmu_m2_miss_corr_sig_18","Missing mass squared corrected",150,0,15); 
  TH1D* hmu_m2_miss_corr_sig = new TH1D("hmu_m2_miss_corr_sig","Missing mass squared corrected",150,0,15);
  TH1D* hmu_Q2_corr_sig_17 = new TH1D("hmu_Q2_corr_sig_17","Squared four momentum transfered to lepton system corrected",200,0,20);
  TH1D* hmu_Q2_corr_sig_18 = new TH1D("hmu_Q2_corr_sig_18","Squared four momentum transfered to lepton system corrected",200,0,20);
  TH1D* hmu_Q2_corr_sig = new TH1D("hmu_Q2_corr_sig","Squared four momentum transfered to lepton system corrected",200,0,20);
  TH1D* hmu_Pt_miss_corr_sig_17 = new TH1D("hmu_Pt_miss_corr_sig_17","Missing transverse momentum corrected",500,0,50);
  TH1D* hmu_Pt_miss_corr_sig_18 = new TH1D("hmu_Pt_miss_corr_sig_18","Missing transverse momentum corrected",500,0,50);
  TH1D* hmu_Pt_miss_corr_sig = new TH1D("hmu_Pt_miss_corr_sig","Missing transverse momentum corrected",500,0,50);
  TH1D* hmu_Pt_var_sig_17 = new TH1D("hmu_Pt_var_sig_17","Tranverse variable momentum",600,-30,30);
  TH1D* hmu_Pt_var_sig_18 = new TH1D("hmu_Pt_var_sig_18","Tranverse variable momentum",600,-30,30);
  TH1D* hmu_Pt_var_sig = new TH1D("hmu_Pt_var_sig","Tranverse variable momentum",600,-30,30);
  TH1D* hmu_del_R_sig_17 = new TH1D("hmu_del_R_sig_17","Delta R",200,0,2);
  TH1D* hmu_del_R_sig_18 = new TH1D("hmu_del_R_sig_18","Delta R",200,0,2);
  TH1D* hmu_del_R_sig = new TH1D("hmu_del_R_sig","Delta R",200,0,2);
  TH1D* hmu_mu_pt_sig_17 = new TH1D("hmu_mu_pt_sig_17","Unpaired muon tranverse momentum",200,0,20);
  TH1D* hmu_mu_pt_sig_18 = new TH1D("hmu_mu_pt_sig_18","Unpaired muon tranverse momentum",200,0,20);
  TH1D* hmu_mu_pt_sig = new TH1D("hmu_mu_pt_sig","Unpaired muon tranverse momentum",200,0,20);
  TH1D* hmu_mu_px_sig_17 = new TH1D("hmu_mu_px_sig_17","Unpaired muon momentum in x-direction",200,0,50);
  TH1D* hmu_mu_px_sig_18 = new TH1D("hmu_mu_px_sig_18","Unpaired muon momentum in x-direction",200,0,50);
  TH1D* hmu_mu_px_sig = new TH1D("hmu_mu_px_sig","Unpaired muon momentum in x-direction",200,0,50);
  TH1D* hmu_mu_py_sig_17 = new TH1D("hmu_mu_py_sig_17","Unpaired muon momentum in y-direction",200,0,50);
  TH1D* hmu_mu_py_sig_18 = new TH1D("hmu_mu_py_sig_18","Unpaired muon momentum in y-direction",200,0,50);
  TH1D* hmu_mu_py_sig = new TH1D("hmu_mu_py_sig","Unpaired muon momentum in y-direction",200,0,50);
  TH1D* hmu_mu_pz_sig_17 = new TH1D("hmu_mu_pz_sig_17","Unpaired muon momentum in z-direction",200,0,50);
  TH1D* hmu_mu_pz_sig_18 = new TH1D("hmu_mu_pz_sig_18","Unpaired muon momentum in z-direction",200,0,50);
  TH1D* hmu_mu_pz_sig = new TH1D("hmu_mu_pz_sig","Unpaired muon momentum in z-direction",200,0,50);
  TH1D* hmu_mu_eta_sig_17 = new TH1D("hmu_mu_eta_sig_17","Unpaired muon psuedorapidity",60,-3,3);
  TH1D* hmu_mu_eta_sig_18 = new TH1D("hmu_mu_eta_sig_18","Unpaired muon psuedorapidity",60,-3,3);
  TH1D* hmu_mu_eta_sig = new TH1D("hmu_mu_eta_sig","Unpaired muon psuedorapidity",60,-3,3);
  TH1D* hmu_mu_phi_sig_17 = new TH1D("hmu_mu_phi_sig_17","Unpaired muon phi",100,-5,5);
  TH1D* hmu_mu_phi_sig_18 = new TH1D("hmu_mu_phi_sig_18","Unpaired muon phi",100,-5,5);
  TH1D* hmu_mu_phi_sig = new TH1D("hmu_mu_phi_sig","Unpaired muon phi",100,-5,5);
  TH1D* hmu_mu_en_bc_frame_sig_17 = new TH1D("hmu_mu_en_bc_frame_sig_17","Energy of Unpaired muon in Bc rest frame",40,0,4);
  TH1D* hmu_mu_en_bc_frame_sig_18 = new TH1D("hmu_mu_en_bc_frame_sig_18","Energy of Unpaired muon in Bc rest frame",40,0,4);
  TH1D* hmu_mu_en_bc_frame_sig = new TH1D("hmu_mu_en_bc_frame_sig","Energy of Unpaired muon in Bc rest frame",40,0,4);
  TH1D* hmu_mu_en_jpsi_frame_sig_17 = new TH1D("hmu_mu_en_jpsi_frame_sig_17","Energy of Unpaired muon in Jpsi rest frame",50,0,5);
  TH1D* hmu_mu_en_jpsi_frame_sig_18 = new TH1D("hmu_mu_en_jpsi_frame_sig_18","Energy of Unpaired muon in Jpsi rest frame",50,0,5);
  TH1D* hmu_mu_en_jpsi_frame_sig = new TH1D("hmu_mu_en_jpsi_frame_sig","Energy of Unpaired muon in Jpsi rest frame",50,0,5);
  TH1D* hmu_bc_pt_sig_17 = new TH1D("hmu_bc_pt_sig_17","Transverse momentum of B_{c}",100,0,100);
  TH1D* hmu_bc_pt_sig_18 = new TH1D("hmu_bc_pt_sig_18","Transverse momentum of B_{c}",100,0,100);
  TH1D* hmu_bc_pt_sig = new TH1D("hmu_bc_pt_sig","Transverse momentum of B_{c}",100,0,100);
  TH1D* hmu_bc_pt_ratio_sig_17 = new TH1D("hmu_bc_pt_ratio_sig_17","Ratio of true and corrected transverse momentum of B_{c}",100,0,2);
  TH1D* hmu_bc_pt_ratio_sig_18 = new TH1D("hmu_bc_pt_ratio_sig_18","Ratio of true and corrected transverse momentum of B_{c}",100,0,2);
  TH1D* hmu_bc_pt_ratio_sig = new TH1D("hmu_bc_pt_ratio_sig","Ratio of true and corrected transverse momentum of B_{c}",100,0,2);
  TH1D* hmu_bc_pt_corr_sig_17 = new TH1D("hmu_bc_pt_corr_sig_17","Corrected transverse momentum of B_{c}",100,0,100);
  TH1D* hmu_bc_pt_corr_sig_18 = new TH1D("hmu_bc_pt_corr_sig_18","Corrected transverse momentum of B_{c}",100,0,100);
  TH1D* hmu_bc_pt_corr_sig = new TH1D("hmu_bc_pt_corr_sig","Corrected transverse momentum of B_{c}",100,0,100);
  auto hprof_sig_17  = new TProfile("hprof_sig_17","Profile of Transverse momentum of B_{c}",120,0,120,0,1000);
  auto hprof_sig_18  = new TProfile("hprof_sig_18","Profile of Transverse momentum of B_{c}",120,0,120,0,1000);
  auto hprof_sig  = new TProfile("hprof_sig","Profile of Transverse momentum of B_{c}",120,0,120,0,1000);
  //TH2D* hmu_bc_pt_scatter_sig_17 = new TH2D("hmu_bc_pt_scatter_sig_17","Scatter plot of transverse momentum of B_{c}",100,-5,5,100,0,20);
  TH2D* hmu_bc_pt_scatter_sig_17 = new TH2D("hmu_bc_pt_scatter_sig_17","Scatter plot of transverse momentum of B_{c}",100,0,100,100,0,100);
  TH2D* hmu_bc_pt_scatter_sig_18 = new TH2D("hmu_bc_pt_scatter_sig_18","Scatter plot of transverse momentum of B_{c}",100,0,100,100,0,100);
  TH2D* hmu_bc_pt_scatter_sig = new TH2D("hmu_bc_pt_scatter_sig","Scatter plot of transverse momentum of B_{c}",100,0,100,100,0,100);
  TH2D* hmu_bc_pt_eta_scatter_sig_17 = new TH2D("hmu_bc_pt_eta_scatter_sig_17","Scatter plot of transverse momentum of B_{c} with eta",100,-5,5,100,0,60);
  TH2D* hmu_bc_pt_eta_scatter_sig_18 = new TH2D("hmu_bc_pt_eta_scatter_sig_18","Scatter plot of transverse momentum of B_{c} with eta",100,-5,5,100,0,60);
  TH2D* hmu_bc_pt_eta_scatter_sig = new TH2D("hmu_bc_pt_eta_scatter_sig","Scatter plot of transverse momentum of B_{c} with eta",100,-5,5,100,0,60);
  TH2D* hmu_mu_pt_eta_scatter_sig_17 = new TH2D("hmu_mu_pt_eta_scatter_sig_17","Scatter plot of transverse momentum of unpaired muon",100,-5,5,100,0,40);
  TH2D* hmu_mu_pt_eta_scatter_sig_18 = new TH2D("hmu_mu_pt_eta_scatter_sig_18","Scatter plot of transverse momentum of unpaired muon",100,-5,5,100,0,40);
  TH2D* hmu_mu_pt_eta_scatter_sig = new TH2D("hmu_mu_pt_eta_scatter_sig","Scatter plot of transverse momentum of unpaired muon",100,-5,5,100,0,40);
  TH2D* hmu_jpsi_pt_eta_scatter_sig_17 = new TH2D("hmu_jpsi_pt_eta_scatter_sig_17","Scatter plot of transverse momentum of J/#Psi with eta",100,-5,5,100,0,20);
  TH2D* hmu_jpsi_pt_eta_scatter_sig_18 = new TH2D("hmu_jpsi_pt_eta_scatter_sig_18","Scatter plot of transverse momentum of J/#Psi with eta",100,-5,5,100,0,20);
  TH2D* hmu_jpsi_pt_eta_scatter_sig = new TH2D("hmu_jpsi_pt_eta_scatter_sig","Scatter plot of transverse momentum of J/#Psi with eta",100,-5,5,100,0,20);
  TH1D* hmu_recco_bc_pt_sig_17 = new TH1D("hmu_recco_bc_pt_sig_17","Reconstructed transverse momentum of B_{c}",100,0,100);
  TH1D* hmu_recco_bc_pt_sig_18 = new TH1D("hmu_recco_bc_pt_sig_18","Reconstructed transverse momentum of B_{c}",100,0,100);
  TH1D* hmu_recco_bc_pt_sig = new TH1D("hmu_recco_bc_pt_sig","Reconstructed transverse momentum of B_{c}",100,0,100);
  TH1D* hmu_jpsi_pt_sig_17 = new TH1D("hmu_jpsi_pt_sig_17","Transverse momentum of J/#Psi",100,0,100);
  TH1D* hmu_jpsi_pt_sig_18 = new TH1D("hmu_jpsi_pt_sig_18","Transverse momentum of J/#Psi",100,0,100);
  TH1D* hmu_jpsi_pt_sig = new TH1D("hmu_bc_jpsi_sig","Transverse momentum of J/#Psi",100,0,100);
  
 


  //ParticleMass muonMass = 0.10565837;
  float mu_mass = 0.113428;

  /**
  //loop over events for signal 2017 channel
  for(int i=0; i<tree_sig_17->GetEntries();i++){
    tree_sig_17->GetEntry(i);
    
    //TLorentzVector class for calculating four momenta
    TLorentzVector momntm_bc_sig_17;
    momntm_bc_sig_17.SetPtEtaPhiM(bc_pt_sig_17, bc_eta_sig_17,bc_phi_sig_17,bc_mass_sig_17);
    TLorentzVector momntm_mu1_sig_17;
    momntm_mu1_sig_17.SetPtEtaPhiM(mu1_pt_sig_17, mu1_eta_sig_17,mu1_phi_sig_17,mu_mass);
    TLorentzVector momntm_mu2_sig_17;
    momntm_mu2_sig_17.SetPtEtaPhiM(mu2_pt_sig_17, mu2_eta_sig_17,mu2_phi_sig_17,mu_mass);
    TLorentzVector momntm_mu_sig_17;
    momntm_mu_sig_17.SetPtEtaPhiM(mu_pt_sig_17, mu_eta_sig_17,mu_phi_sig_17,mu_mass);
    TLorentzVector momntm_jpsi_sig_17;
    momntm_jpsi_sig_17.SetPtEtaPhiM(jpsi_pt_sig_17, jpsi_eta_sig_17,jpsi_phi_sig_17,jpsi_mass_sig_17);

    //calculating the kinematic variables
    float m2_miss_sig_17 =(momntm_bc_sig_17 - momntm_mu1_sig_17 - momntm_mu2_sig_17 - momntm_mu_sig_17).M2();
    float Q2_sig_17 = (momntm_bc_sig_17 - momntm_mu1_sig_17 - momntm_mu2_sig_17).M2();
    float Pt_miss_sig_17 = (bc_pt_sig_17 - mu1_pt_sig_17 - mu2_pt_sig_17 - mu_pt_sig_17);
    float Pt_var_sig_17 = (jpsi_pt_sig_17 - mu_pt_sig_17 );
    float del_R_sig_17 = momntm_mu1_sig_17.DeltaR(momntm_mu2_sig_17);
    TLorentzVector muon_system_sig_17 = (momntm_mu1_sig_17+momntm_mu2_sig_17+momntm_mu_sig_17);
    float bc_pt_corr_sig_17 = (bc_mass_sig_17 * muon_system_sig_17.Pt())/ (muon_system_sig_17.M()); //corrected Pt according to recipe in Jona thesis
    float bc_pt_ratio_sig_17 = bc_pt_corr_sig_17/bc_pt_sig_17;
    
    TLorentzVector momntm_bc_corr_sig_17;
    momntm_bc_corr_sig_17.SetPtEtaPhiM(bc_pt_corr_sig_17, muon_system_sig_17.Eta(), muon_system_sig_17.Phi(),momntm_bc_sig_17.M());
    float Q2_corr_sig_17 = (momntm_bc_corr_sig_17 - momntm_mu1_sig_17 - momntm_mu2_sig_17).M2();
    float m2_miss_corr_sig_17 = (momntm_bc_corr_sig_17 - momntm_mu1_sig_17 - momntm_mu2_sig_17 - momntm_mu_sig_17).M2();
    float Pt_miss_corr_sig_17 = (bc_pt_corr_sig_17 - mu1_pt_sig_17 - mu2_pt_sig_17 - mu_pt_sig_17);

    TVector3 X_bc_sig_17 = -momntm_bc_corr_sig_17.BoostVector();
    TVector3 X_jpsi_sig_17 = -momntm_jpsi_sig_17.BoostVector();
    TLorentzVector momntm_mu_bc_frame_sig_17 = momntm_mu_sig_17;
    momntm_mu_bc_frame_sig_17.Boost(X_bc_sig_17);
    TLorentzVector momntm_mu_jpsi_frame_sig_17 = momntm_mu_sig_17;
    momntm_mu_jpsi_frame_sig_17.Boost(X_jpsi_sig_17);

    float mu_en_bc_frame_sig_17 = momntm_mu_bc_frame_sig_17.E();
    float mu_en_jpsi_frame_sig_17 = momntm_mu_jpsi_frame_sig_17.E();

    float recco_bc_pt_sig_17 = mu1_pt_sig_17 + mu2_pt_sig_17 + mu_pt_sig_17 + munu_pt_sig_17;
    float mu_px_sig_17 = momntm_mu_sig_17.Px();
    float mu_py_sig_17 = momntm_mu_sig_17.Py();
    float mu_pz_sig_17 = momntm_mu_sig_17.Pz();
        
    //Filling the histograms

    //if ((abs(mu_eta_sig_17) < 1.3 and mu_pt_sig_17 > 3.3) or (1.3 < abs(mu_eta_sig_17) < 2.2 and mu_pt_sig_17 > 2.9) or (2.2 < abs(mu_eta_sig_17) < 2.4 and mu_pt_sig_17 > 2.4)){ //acceptance cuts
      //if ((Trig_dimuon_sig_17 == false) and (Trig_JpsiTk_sig_17 == true)){ //trigger flags
	hmu_m2_miss_sig_17->Fill(m2_miss_sig_17);
	hmu_Q2_sig_17->Fill(Q2_sig_17);
	hmu_Pt_miss_sig_17->Fill(Pt_miss_sig_17);
	hmu_m2_miss_corr_sig_17->Fill(m2_miss_corr_sig_17);
	hmu_Q2_corr_sig_17->Fill(Q2_corr_sig_17);
	hmu_Pt_miss_corr_sig_17->Fill(Pt_miss_corr_sig_17);
	hmu_Pt_var_sig_17->Fill(Pt_var_sig_17);
	hmu_mu_pt_sig_17->Fill(mu_pt_sig_17);
	hmu_mu_px_sig_17->Fill(mu_px_sig_17);
	hmu_mu_py_sig_17->Fill(mu_py_sig_17);
	hmu_mu_pz_sig_17->Fill(mu_pz_sig_17);
	hmu_mu_eta_sig_17->Fill(mu_eta_sig_17);
	hmu_mu_phi_sig_17->Fill(mu_phi_sig_17);
	hmu_del_R_sig_17->Fill(del_R_sig_17);
	hmu_mu_en_bc_frame_sig_17->Fill(mu_en_bc_frame_sig_17);
	hmu_mu_en_jpsi_frame_sig_17->Fill(mu_en_jpsi_frame_sig_17);
	hmu_bc_pt_sig_17->Fill(bc_pt_sig_17);
	hmu_bc_pt_corr_sig_17->Fill(bc_pt_corr_sig_17);
	hmu_bc_pt_ratio_sig_17->Fill(bc_pt_ratio_sig_17);
	hprof_sig_17->Fill(bc_pt_sig_17,bc_pt_ratio_sig_17,1);
	hmu_bc_pt_scatter_sig_17->Fill(bc_pt_sig_17,bc_pt_corr_sig_17);
	hmu_bc_pt_eta_scatter_sig_17->Fill(bc_eta_sig_17,bc_pt_sig_17);
	hmu_mu_pt_eta_scatter_sig_17->Fill(mu_eta_sig_17,mu_pt_sig_17);
	hmu_jpsi_pt_eta_scatter_sig_17->Fill(jpsi_eta_sig_17,jpsi_pt_sig_17);
	hmu_recco_bc_pt_sig_17->Fill(recco_bc_pt_sig_17);
	hmu_jpsi_pt_sig_17->Fill(jpsi_pt_sig_17);
	// }
	//}
  }

  **/
  //loop over events for signal 2018 channel
  for(int i=0; i<tree_sig_18->GetEntries();i++){
    tree_sig_18->GetEntry(i);
    
    //TLorentzVector class for calculating four momenta
    TLorentzVector momntm_bc_sig_18;
    momntm_bc_sig_18.SetPtEtaPhiM(bc_pt_sig_18, bc_eta_sig_18,bc_phi_sig_18,bc_mass_sig_18);
    TLorentzVector momntm_mu1_sig_18;
    momntm_mu1_sig_18.SetPtEtaPhiM(mu1_pt_sig_18, mu1_eta_sig_18,mu1_phi_sig_18,mu_mass);
    TLorentzVector momntm_mu2_sig_18;
    momntm_mu2_sig_18.SetPtEtaPhiM(mu2_pt_sig_18, mu2_eta_sig_18,mu2_phi_sig_18,mu_mass);
    TLorentzVector momntm_mu_sig_18;
    momntm_mu_sig_18.SetPtEtaPhiM(mu_pt_sig_18, mu_eta_sig_18,mu_phi_sig_18,mu_mass);
    TLorentzVector momntm_jpsi_sig_18;
    momntm_jpsi_sig_18.SetPtEtaPhiM(jpsi_pt_sig_18, jpsi_eta_sig_18,jpsi_phi_sig_18,jpsi_mass_sig_18);

    //calculating the kinematic variables
    float m2_miss_sig_18 =(momntm_bc_sig_18 - momntm_mu1_sig_18 - momntm_mu2_sig_18 - momntm_mu_sig_18).M2();
    float Q2_sig_18 = (momntm_bc_sig_18 - momntm_mu1_sig_18 - momntm_mu2_sig_18).M2();
    float Pt_miss_sig_18 = (bc_pt_sig_18 - mu1_pt_sig_18 - mu2_pt_sig_18 - mu_pt_sig_18);
    float Pt_var_sig_18 = (jpsi_pt_sig_18 - mu_pt_sig_18 );
    float del_R_sig_18 = momntm_mu1_sig_18.DeltaR(momntm_mu2_sig_18);
    TLorentzVector muon_system_sig_18 = (momntm_mu1_sig_18+momntm_mu2_sig_18+momntm_mu_sig_18);
    float bc_pt_corr_sig_18 = (bc_mass_sig_18 * muon_system_sig_18.Pt())/ (muon_system_sig_18.M()); //corrected Pt according to recipe in Jona thesis
    float bc_pt_ratio_sig_18 = bc_pt_corr_sig_18/bc_pt_sig_18;
    
    TLorentzVector momntm_bc_corr_sig_18;
    momntm_bc_corr_sig_18.SetPtEtaPhiM(bc_pt_corr_sig_18, muon_system_sig_18.Eta(), muon_system_sig_18.Phi(),momntm_bc_sig_18.M());
    float Q2_corr_sig_18 = (momntm_bc_corr_sig_18 - momntm_mu1_sig_18 - momntm_mu2_sig_18).M2();
    float m2_miss_corr_sig_18 = (momntm_bc_corr_sig_18 - momntm_mu1_sig_18 - momntm_mu2_sig_18 - momntm_mu_sig_18).M2();
    float Pt_miss_corr_sig_18 = (bc_pt_corr_sig_18 - mu1_pt_sig_18 - mu2_pt_sig_18 - mu_pt_sig_18);

    TVector3 X_bc_sig_18 = -momntm_bc_corr_sig_18.BoostVector();
    TVector3 X_jpsi_sig_18 = -momntm_jpsi_sig_18.BoostVector();
    TLorentzVector momntm_mu_bc_frame_sig_18 = momntm_mu_sig_18;
    momntm_mu_bc_frame_sig_18.Boost(X_bc_sig_18);
    TLorentzVector momntm_mu_jpsi_frame_sig_18 = momntm_mu_sig_18;
    momntm_mu_jpsi_frame_sig_18.Boost(X_jpsi_sig_18);

    float mu_en_bc_frame_sig_18 = momntm_mu_bc_frame_sig_18.E();
    float mu_en_jpsi_frame_sig_18 = momntm_mu_jpsi_frame_sig_18.E();

    float recco_bc_pt_sig_18 = mu1_pt_sig_18 + mu2_pt_sig_18 + mu_pt_sig_18 + munu_pt_sig_18;
    float mu_px_sig_18 = momntm_mu_sig_18.Px();
    float mu_py_sig_18 = momntm_mu_sig_18.Py();
    float mu_pz_sig_18 = momntm_mu_sig_18.Pz();
        
    //Filling the histograms

    //if ((abs(mu_eta_sig_18) < 1.3 and mu_pt_sig_18 > 3.3) or (1.3 < abs(mu_eta_sig_18) < 2.2 and mu_pt_sig_18 > 2.9) or (2.2 < abs(mu_eta_sig_18) < 2.4 and mu_pt_sig_18 > 2.4)){ //acceptance cuts
    // if ((Trig_dimuon_sig_18 == false) and (Trig_JpsiTk_sig_18 == true)){ //trigger flags
	hmu_m2_miss_sig_18->Fill(m2_miss_sig_18);
	hmu_Q2_sig_18->Fill(Q2_sig_18);
	hmu_Pt_miss_sig_18->Fill(Pt_miss_sig_18);
	hmu_m2_miss_corr_sig_18->Fill(m2_miss_corr_sig_18);
	hmu_Q2_corr_sig_18->Fill(Q2_corr_sig_18);
	hmu_Pt_miss_corr_sig_18->Fill(Pt_miss_corr_sig_18);
	hmu_Pt_var_sig_18->Fill(Pt_var_sig_18);
	hmu_mu_pt_sig_18->Fill(mu_pt_sig_18);
	hmu_mu_px_sig_18->Fill(mu_px_sig_18);
	hmu_mu_py_sig_18->Fill(mu_py_sig_18);
	hmu_mu_pz_sig_18->Fill(mu_pz_sig_18);
	hmu_mu_eta_sig_18->Fill(mu_eta_sig_18);
	hmu_mu_phi_sig_18->Fill(mu_phi_sig_18);
	hmu_del_R_sig_18->Fill(del_R_sig_18);
	hmu_mu_en_bc_frame_sig_18->Fill(mu_en_bc_frame_sig_18);
	hmu_mu_en_jpsi_frame_sig_18->Fill(mu_en_jpsi_frame_sig_18);
	hmu_bc_pt_sig_18->Fill(bc_pt_sig_18);
	hmu_bc_pt_corr_sig_18->Fill(bc_pt_corr_sig_18);
	hmu_bc_pt_ratio_sig_18->Fill(bc_pt_ratio_sig_18);
	hprof_sig_18->Fill(bc_pt_sig_18,bc_pt_ratio_sig_18,1);
	hmu_bc_pt_scatter_sig_18->Fill(bc_pt_sig_18,bc_pt_corr_sig_18);
	hmu_bc_pt_eta_scatter_sig_18->Fill(bc_eta_sig_18,bc_pt_sig_18);
	hmu_mu_pt_eta_scatter_sig_18->Fill(mu_eta_sig_18,mu_pt_sig_18);
	hmu_jpsi_pt_eta_scatter_sig_18->Fill(jpsi_eta_sig_18,jpsi_pt_sig_18);
	hmu_recco_bc_pt_sig_18->Fill(recco_bc_pt_sig_18);
	hmu_jpsi_pt_sig_18->Fill(jpsi_pt_sig_18);
	//}
	//}
  }

  //loop over events for signal channel
  for(int i=0; i<tree_sig->GetEntries();i++){
    tree_sig->GetEntry(i);

    //TLorentzVector class for calculating four momenta
    TLorentzVector momntm_bc_sig;
    momntm_bc_sig.SetPtEtaPhiM(bc_pt_sig, bc_eta_sig,bc_phi_sig,bc_mass_sig);
    TLorentzVector momntm_mu1_sig;
    momntm_mu1_sig.SetPtEtaPhiM(mu1_pt_sig, mu1_eta_sig,mu1_phi_sig,mu_mass);
    TLorentzVector momntm_mu2_sig;
    momntm_mu2_sig.SetPtEtaPhiM(mu2_pt_sig, mu2_eta_sig,mu2_phi_sig,mu_mass);
    TLorentzVector momntm_mu_sig;
    momntm_mu_sig.SetPtEtaPhiM(mu_pt_sig, mu_eta_sig,mu_phi_sig,mu_mass);
    TLorentzVector momntm_jpsi_sig;
    momntm_jpsi_sig.SetPtEtaPhiM(jpsi_pt_sig, jpsi_eta_sig,jpsi_phi_sig,jpsi_mass_sig);

    //calculating the kinematic variables
    float m2_miss_sig = (momntm_bc_sig - momntm_mu1_sig - momntm_mu2_sig - momntm_mu_sig).M2();
    float Q2_sig = (momntm_bc_sig - momntm_mu1_sig - momntm_mu2_sig).M2();
    float Pt_miss_sig = (bc_pt_sig-mu1_pt_sig-mu2_pt_sig-mu_pt_sig);
    float Pt_var_sig = (jpsi_pt_sig - mu_pt_sig );
    float del_R_sig = momntm_mu1_sig.DeltaR(momntm_mu2_sig);
    TLorentzVector muon_system_sig = (momntm_mu1_sig+momntm_mu2_sig+momntm_mu_sig);
    float bc_pt_corr_sig = (bc_mass_sig * muon_system_sig.Pt())/ (muon_system_sig.M());
    float bc_pt_ratio_sig = bc_pt_corr_sig/bc_pt_sig;

    TLorentzVector momntm_bc_corr_sig;
    momntm_bc_corr_sig.SetPtEtaPhiM(bc_pt_corr_sig, muon_system_sig.Eta(), muon_system_sig.Phi(), momntm_bc_sig.M());
    float Q2_corr_sig = (momntm_bc_corr_sig - momntm_mu1_sig - momntm_mu2_sig).M2();
    float m2_miss_corr_sig = (momntm_bc_corr_sig - momntm_mu1_sig - momntm_mu2_sig - momntm_mu_sig).M2();
    float Pt_miss_corr_sig = (bc_pt_corr_sig - mu1_pt_sig - mu2_pt_sig - mu_pt_sig);
   
    TVector3 X_bc_sig = -momntm_bc_corr_sig.BoostVector();
    TVector3 X_jpsi_sig = -momntm_jpsi_sig.BoostVector();
    TLorentzVector momntm_mu_bc_frame_sig = momntm_mu_sig;
    momntm_mu_bc_frame_sig.Boost(X_bc_sig);
    TLorentzVector momntm_mu_jpsi_frame_sig = momntm_mu_sig;
    momntm_mu_jpsi_frame_sig.Boost(X_jpsi_sig);
    
    float mu_en_bc_frame_sig = momntm_mu_bc_frame_sig.E();
    float mu_en_jpsi_frame_sig = momntm_mu_jpsi_frame_sig.E();

    float recco_bc_pt_sig = mu1_pt_sig + mu2_pt_sig + mu_pt_sig + taunu_pt_sig + taunu2_pt_sig + munu_pt_sig;
    float mu_px_sig = momntm_mu_sig.Px();
    float mu_py_sig = momntm_mu_sig.Py();
    float mu_pz_sig = momntm_mu_sig.Pz();
    
    //Filling the histograms
    //if ((abs(mu_eta_sig) < 1.3 and mu_pt_sig > 3.3) or (1.3 < abs(mu_eta_sig) < 2.2 and mu_pt_sig > 2.9) or (2.2 < abs(mu_eta_sig) < 2.4 and mu_pt_sig > 2.4)){ //acceptance cuts
    //if ((Trig_dimuon_sig == false) and (Trig_JpsiTk_sig == true)){ //trigger flags
    if ((-2.7 < mu1_eta_sig < 2.7) and (-2.7 < mu2_eta_sig < 2.7) and (mu1_pt_sig > 3.0) and (mu2_pt_sig > 3.0)){
	hmu_m2_miss_sig->Fill(m2_miss_sig);
	hmu_Q2_sig->Fill(Q2_sig);
	hmu_Pt_miss_sig->Fill(Pt_miss_sig);
	hmu_m2_miss_corr_sig->Fill(m2_miss_corr_sig);
	hmu_Q2_corr_sig->Fill(Q2_corr_sig);
	hmu_Pt_miss_corr_sig->Fill(Pt_miss_corr_sig);
	hmu_Pt_var_sig->Fill(Pt_var_sig);
	hmu_mu_pt_sig->Fill(mu_pt_sig);
	hmu_mu_px_sig->Fill(mu_px_sig);
	hmu_mu_py_sig->Fill(mu_py_sig);
	hmu_mu_pz_sig->Fill(mu_pz_sig);
	hmu_mu_eta_sig->Fill(mu_eta_sig);
	hmu_mu_phi_sig->Fill(mu_phi_sig);
	hmu_del_R_sig->Fill(del_R_sig);
	hmu_mu_en_bc_frame_sig->Fill(mu_en_bc_frame_sig);
	hmu_mu_en_jpsi_frame_sig->Fill(mu_en_jpsi_frame_sig);
	hmu_bc_pt_sig->Fill(bc_pt_sig);
	hmu_bc_pt_corr_sig->Fill(bc_pt_corr_sig);
	hmu_bc_pt_ratio_sig->Fill(bc_pt_ratio_sig);
	hprof_sig->Fill(bc_pt_sig,bc_pt_ratio_sig,1);
	hmu_bc_pt_scatter_sig->Fill(bc_pt_sig,bc_pt_corr_sig);
	hmu_bc_pt_eta_scatter_sig->Fill(bc_eta_sig,bc_pt_sig);
	hmu_mu_pt_eta_scatter_sig->Fill(mu_eta_sig,mu_pt_sig);
	hmu_jpsi_pt_eta_scatter_sig->Fill(jpsi_eta_sig,jpsi_pt_sig);
	hmu_recco_bc_pt_sig->Fill(recco_bc_pt_sig);
	hmu_jpsi_pt_sig->Fill(jpsi_pt_sig);
    }
	//}
	//    }
  }  
  
  //Drawing Histograms

  /**
  //Missing mass Squared
  hmu_m2_miss_sig_18->SetLineColor(kRed);
  hmu_m2_miss_sig_18->GetXaxis()->SetTitle("m^{2}_{miss}[GeV^{2}]");
  hmu_m2_miss_sig_18->GetYaxis()->SetTitle("a.u.");
  hmu_m2_miss_sig_18->DrawNormalized();
  hmu_m2_miss_sig_17->SetLineColor(kBlue);
  hmu_m2_miss_sig_17->DrawNormalized("same");
  hmu_m2_miss_sig->SetLineColor(kGreen);
  hmu_m2_miss_sig->DrawNormalized("same");
  legend->AddEntry(hmu_m2_miss_sig_18, "2018 data", "l");
  legend->AddEntry(hmu_m2_miss_sig_17, "2017 data", "l");
  legend->AddEntry(hmu_m2_miss_sig, "old data", "l");
  legend->Draw();
  c1->SaveAs("M2_miss.png");

  
  //Squared four momentum transfered to lepton system
  hmu_Q2_sig_18->SetLineColor(kRed);
  hmu_Q2_sig_18->GetXaxis()->SetTitle("Q^{2}[GeV^{2}]");
  hmu_Q2_sig_18->GetYaxis()->SetTitle("a.u.");
  hmu_Q2_sig_18->DrawNormalized();
  hmu_Q2_sig_17->SetLineColor(kBlue);
  hmu_Q2_sig_17->DrawNormalized("same");
  hmu_Q2_sig->SetLineColor(kGreen);
  hmu_Q2_sig->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("Q2.png");
  


  //Missing transverse momentum
  hmu_Pt_miss_sig_18->SetLineColor(kRed);
  hmu_Pt_miss_sig_18->GetXaxis()->SetTitle("P_{T}^{miss}[GeV]");
  hmu_Pt_miss_sig_18->GetYaxis()->SetTitle("a.u.");
  hmu_Pt_miss_sig_18->DrawNormalized();
  hmu_Pt_miss_sig_17->SetLineColor(kBlue);
  hmu_Pt_miss_sig_17->DrawNormalized("same");
  hmu_Pt_miss_sig->SetLineColor(kGreen);
  hmu_Pt_miss_sig->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("Pt_miss.png");


   //Missing mass Squared corrected
  hmu_m2_miss_corr_sig_18->SetLineColor(kRed);
  hmu_m2_miss_corr_sig_18->GetXaxis()->SetTitle("m^{2}_{miss}[GeV^{2}]");
  hmu_m2_miss_corr_sig_18->GetYaxis()->SetTitle("a.u.");
  hmu_m2_miss_corr_sig_18->DrawNormalized();
  hmu_m2_miss_corr_sig_17->SetLineColor(kBlue);
  hmu_m2_miss_corr_sig_17->DrawNormalized("same");
  hmu_m2_miss_corr_sig->SetLineColor(kGreen);
  hmu_m2_miss_corr_sig->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("M2_miss_corr.png");

  
  //Squared four momentum transfered to lepton system corrected
  hmu_Q2_corr_sig_18->SetLineColor(kRed);
  hmu_Q2_corr_sig_18->GetXaxis()->SetTitle("Q^{2}[GeV^{2}]");
  hmu_Q2_corr_sig_18->GetYaxis()->SetTitle("a.u.");
  hmu_Q2_corr_sig_18->DrawNormalized();
  hmu_Q2_corr_sig_17->SetLineColor(kBlue);
  hmu_Q2_corr_sig_17->DrawNormalized("same");
  hmu_Q2_corr_sig->SetLineColor(kGreen);
  hmu_Q2_corr_sig->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("Q2_corr.png");
  


  //Missing transverse momentum corrected
  hmu_Pt_miss_corr_sig_18->SetLineColor(kRed);
  hmu_Pt_miss_corr_sig_18->GetXaxis()->SetTitle("P_{T}^{miss}[GeV]");
  hmu_Pt_miss_corr_sig_18->GetYaxis()->SetTitle("a.u.");
  hmu_Pt_miss_corr_sig_18->DrawNormalized();
  hmu_Pt_miss_corr_sig_17->SetLineColor(kBlue);
  hmu_Pt_miss_corr_sig_17->DrawNormalized("same");
  hmu_Pt_miss_corr_sig->SetLineColor(kGreen);
  hmu_Pt_miss_corr_sig->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("Pt_miss_corr.png");
  

  //Tranverse variable momentum
  hmu_Pt_var_sig_18->SetLineColor(kRed);
  hmu_Pt_var_sig_18->GetXaxis()->SetTitle("P_{T}^{var}[GeV]");
  hmu_Pt_var_sig_18->GetYaxis()->SetTitle("a.u.");
  hmu_Pt_var_sig_18->DrawNormalized();
  hmu_Pt_var_sig_17->SetLineColor(kBlue);
  hmu_Pt_var_sig_17->DrawNormalized("same");
  hmu_Pt_var_sig->SetLineColor(kGreen);
  hmu_Pt_var_sig->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("Pt_var.png");
  
  //Delta R
  hmu_del_R_sig_18->SetLineColor(kRed);
  hmu_del_R_sig_18->GetXaxis()->SetTitle("#Delta R(#mu_{1},#mu_{2})");
  hmu_del_R_sig_18->GetYaxis()->SetTitle("a.u.");
  hmu_del_R_sig_18->DrawNormalized();
  hmu_del_R_sig_17->SetLineColor(kBlue);
  hmu_del_R_sig_17->DrawNormalized("same");
  hmu_del_R_sig->SetLineColor(kGreen);
  hmu_del_R_sig->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("del_R.png");

  
  //Unpaired muon tranverse momentum
  hmu_mu_pt_sig_18->SetLineColor(kRed);
  hmu_mu_pt_sig_18->GetXaxis()->SetTitle("P_{T}^{#mu}[GeV]");
  hmu_mu_pt_sig_18->GetYaxis()->SetTitle("a.u.");
  hmu_mu_pt_sig_18->DrawNormalized();
  hmu_mu_pt_sig_17->SetLineColor(kBlue);
  hmu_mu_pt_sig_17->DrawNormalized("same");
  hmu_mu_pt_sig->SetLineColor(kGreen);
  hmu_mu_pt_sig->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("mu_pt.png");

  
  //Unpaired muon momentum in x-direction
  hmu_mu_px_sig_18->SetLineColor(kRed);
  hmu_mu_px_sig_18->GetXaxis()->SetTitle("P_{x}^{#mu}[GeV]");
  hmu_mu_px_sig_18->GetYaxis()->SetTitle("a.u.");
  hmu_mu_px_sig_18->DrawNormalized();
  hmu_mu_px_sig->SetLineColor(kBlue);
  hmu_mu_px_sig->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("mu_px.png");


  //Unpaired muon momentum in y-direction
  hmu_mu_py_sig_18->SetLineColor(kRed);
  hmu_mu_py_sig_18->GetXaxis()->SetTitle("P_{y}^{#mu}[GeV]");
  hmu_mu_py_sig_18->GetYaxis()->SetTitle("a.u.");
  hmu_mu_py_sig_18->DrawSig_18alized();
  hmu_mu_py_sig->SetLineColor(kBlue);
  hmu_mu_py_sig->DrawSig_18alized("same");
  legend->Draw();
  c1->SaveAs("mu_py.png");


  //Unpaired muon momentum in z-direction
  hmu_mu_pz_sig_18->SetLineColor(kRed);
  hmu_mu_pz_sig_18->GetXaxis()->SetTitle("P_{z}^{#mu}[GeV]");
  hmu_mu_pz_sig_18->GetYaxis()->SetTitle("a.u.");
  hmu_mu_pz_sig_18->DrawNormalized();
  hmu_mu_pz_sig->SetLineColor(kBlue);
  hmu_mu_pz_sig->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("mu_pz.png");

  
  //Unpaired muon psuedorapidity
  hmu_mu_eta_sig_18->SetLineColor(kRed);
  hmu_mu_eta_sig_18->GetXaxis()->SetTitle("#eta^{#mu}");
  hmu_mu_eta_sig_18->GetYaxis()->SetTitle("a.u.");
  hmu_mu_eta_sig_18->DrawNormalized();
  hmu_mu_eta_sig_17->SetLineColor(kBlue);
  hmu_mu_eta_sig_17->DrawNormalized("same");
  hmu_mu_eta_sig->SetLineColor(kGreen);
  hmu_mu_eta_sig->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("mu_eta.png");

  
  //Unpaired muon phi
  hmu_mu_phi_sig_18->SetLineColor(kRed);
  hmu_mu_phi_sig_18->GetXaxis()->SetTitle("#phi^{#mu}");
  hmu_mu_phi_sig_18->GetYaxis()->SetTitle("a.u.");
  hmu_mu_phi_sig_18->DrawNormalized();
  hmu_mu_phi_sig_17->SetLineColor(kBlue);
  hmu_mu_phi_sig_17->DrawNormalized("same");
  hmu_mu_phi_sig->SetLineColor(kGreen);
  hmu_mu_phi_sig->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("mu_phi.png");
  
  
  //Energy of unpaired muon in Bc rest frame
  hmu_mu_en_bc_frame_sig_18->SetLineColor(kRed);
  hmu_mu_en_bc_frame_sig_18->GetXaxis()->SetTitle("E(#mu) in rest frame of B_{c} (GeV)");
  hmu_mu_en_bc_frame_sig_18->GetYaxis()->SetTitle("a.u.");
  hmu_mu_en_bc_frame_sig_18->DrawNormalized();
  hmu_mu_en_bc_frame_sig_17->SetLineColor(kBlue);
  hmu_mu_en_bc_frame_sig_17->DrawNormalized("same");
  hmu_mu_en_bc_frame_sig->SetLineColor(kGreen);
  hmu_mu_en_bc_frame_sig->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("E_bc.png");

  //Energy of unpaired muon in Jpsi rest frame
  hmu_mu_en_jpsi_frame_sig_18->SetLineColor(kRed);
  hmu_mu_en_jpsi_frame_sig_18->GetXaxis()->SetTitle("E(#mu) in rest frame of J/(#Psi) (GeV)");
  hmu_mu_en_jpsi_frame_sig_18->GetYaxis()->SetTitle("a.u.");
  hmu_mu_en_jpsi_frame_sig_18->DrawNormalized();
  hmu_mu_en_jpsi_frame_sig_17->SetLineColor(kBlue);
  hmu_mu_en_jpsi_frame_sig_17->DrawNormalized("same");
  hmu_mu_en_jpsi_frame_sig->SetLineColor(kGreen);
  hmu_mu_en_jpsi_frame_sig->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("E_jpsi.png");
**/
  //Transverse momentum of Bc
  hmu_jpsi_pt_sig_18->SetLineColor(kRed);
  hmu_jpsi_pt_sig_18->GetXaxis()->SetTitle("P_{T}^{J/#Psi} (GeV)");
  hmu_jpsi_pt_sig_18->GetYaxis()->SetTitle("a.u.");
  hmu_jpsi_pt_sig_18->DrawNormalized();
  hmu_jpsi_pt_sig->SetLineColor(kGreen);
  hmu_jpsi_pt_sig->DrawNormalized("same");
  legend->AddEntry(hmu_m2_miss_sig_18, "2018 data", "l");
  legend->AddEntry(hmu_m2_miss_sig, "old data", "l");
  legend->Draw();
  c1->SaveAs("Pt_jpsi.png");

  
  //Transverse momentum of Bc
  hmu_bc_pt_sig_18->SetLineColor(kRed);
  hmu_bc_pt_sig_18->GetXaxis()->SetTitle("P_{T}^{Bc} (GeV)");
  hmu_bc_pt_sig_18->GetYaxis()->SetTitle("a.u.");
  hmu_bc_pt_sig_18->DrawNormalized();
  hmu_bc_pt_sig->SetLineColor(kGreen);
  hmu_bc_pt_sig->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("Pt_Bc.png");

  //Unpaired muon tranverse momentum
  hmu_mu_pt_sig_18->SetLineColor(kRed);
  hmu_mu_pt_sig_18->GetXaxis()->SetTitle("P_{T}^{#mu}[GeV]");
  hmu_mu_pt_sig_18->GetYaxis()->SetTitle("a.u.");
  hmu_mu_pt_sig_18->DrawNormalized();
  hmu_mu_pt_sig->SetLineColor(kGreen);
  hmu_mu_pt_sig->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("mu_pt.png");

  
  //Corrected transverse momentum of Bc
  hmu_bc_pt_corr_sig_18->SetLineColor(kRed);
  hmu_bc_pt_corr_sig_18->GetXaxis()->SetTitle("Corrected P_{T}^{Bc} (GeV)");
  hmu_bc_pt_corr_sig_18->GetYaxis()->SetTitle("a.u.");
  hmu_bc_pt_corr_sig_18->DrawNormalized();
  hmu_bc_pt_corr_sig->SetLineColor(kGreen);
  hmu_bc_pt_corr_sig->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("Corr_pt_bc.png");

  //Ratio of true and corrected transverse momentum of Bc
  hmu_bc_pt_ratio_sig_18->SetLineColor(kRed);
  hmu_bc_pt_ratio_sig_18->GetXaxis()->SetTitle("Ratio of true and corrected P_{T}^{Bc} (GeV)");
  hmu_bc_pt_ratio_sig_18->GetYaxis()->SetTitle("a.u.");
  hmu_bc_pt_ratio_sig_18->DrawNormalized();
  hmu_bc_pt_ratio_sig->SetLineColor(kGreen);
  hmu_bc_pt_ratio_sig->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("Ratio_pt_bc.png");
  

  //Scatter plot of Transverse momentum of Bc
  hmu_bc_pt_scatter_sig_18->SetMarkerColor(kRed);
  hmu_bc_pt_scatter_sig_18->GetXaxis()->SetTitle("P_{T}^{B_{c}} (GeV)");
  hmu_bc_pt_scatter_sig_18->GetYaxis()->SetTitle("corrected P_{T}^{B_{c}} (GeV)");
  hmu_bc_pt_scatter_sig_18->DrawNormalized();
  hmu_bc_pt_scatter_sig->SetMarkerColor(kGreen);
  hmu_bc_pt_scatter_sig->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("scatter_Pt_Bc.png");


  //Scatter plot of Transverse momentum of Bc with eta for old data
  hmu_bc_pt_eta_scatter_sig->SetMarkerColor(kGreen);
  hmu_bc_pt_eta_scatter_sig->GetXaxis()->SetTitle("#eta^{B_{c}}");
  hmu_bc_pt_eta_scatter_sig->GetYaxis()->SetTitle("P_{T}^{B_{C}} (GeV)");
  hmu_bc_pt_eta_scatter_sig->Draw();
  legend->Draw();
  c1->SaveAs("scatter_Pt_eta_Bc_old.png");


  //Scatter plot of Transverse momentum of Bc with eta for 2018
  hmu_bc_pt_eta_scatter_sig_18->SetMarkerColor(kRed);
  hmu_bc_pt_eta_scatter_sig_18->GetXaxis()->SetTitle("#eta^{B_{c}}");
  hmu_bc_pt_eta_scatter_sig_18->GetYaxis()->SetTitle("P_{T}^{B_{c}} (GeV)");
  hmu_bc_pt_eta_scatter_sig_18->Draw();
  c1->SaveAs("scatter_Pt_eta_Bc_2018.png");

  //Scatter plot of Transverse momentum of mu with eta for old data
  hmu_mu_pt_eta_scatter_sig->SetMarkerColor(kGreen);
  hmu_mu_pt_eta_scatter_sig->GetXaxis()->SetTitle("#eta^{#mu}");
  hmu_mu_pt_eta_scatter_sig->GetYaxis()->SetTitle("P_{T}^{#mu} (GeV)");
  hmu_mu_pt_eta_scatter_sig->Draw();
  legend->Draw();
  c1->SaveAs("scatter_Pt_eta_mu_old.png");


  //Scatter plot of Transverse momentum of mu with eta for 2018
  hmu_mu_pt_eta_scatter_sig_18->SetMarkerColor(kRed);
  hmu_mu_pt_eta_scatter_sig_18->GetXaxis()->SetTitle("#eta^{#mu}");
  hmu_mu_pt_eta_scatter_sig_18->GetYaxis()->SetTitle("P_{T}^{#mu} (GeV)");
  hmu_mu_pt_eta_scatter_sig_18->Draw();
  c1->SaveAs("scatter_Pt_eta_mu_2018.png");

  //Scatter plot of Transverse momentum of jpsi with eta for old data
  hmu_jpsi_pt_eta_scatter_sig->SetMarkerColor(kGreen);
  hmu_jpsi_pt_eta_scatter_sig->GetXaxis()->SetTitle("#eta^{J/#Psi}");
  hmu_jpsi_pt_eta_scatter_sig->GetYaxis()->SetTitle("P_{T}^{J/#Psi} (GeV)");
  hmu_jpsi_pt_eta_scatter_sig->Draw();
  legend->Draw();
  c1->SaveAs("scatter_Pt_eta_jpsi_old.png");


  //Scatter plot of Transverse momentum of jpsi with eta for 2018
  hmu_jpsi_pt_eta_scatter_sig_18->SetMarkerColor(kRed);
  hmu_jpsi_pt_eta_scatter_sig_18->GetXaxis()->SetTitle("#eta^{J/#Psi}");
  hmu_jpsi_pt_eta_scatter_sig_18->GetYaxis()->SetTitle("P_{T}^{J/#Psi} (GeV)");
  hmu_jpsi_pt_eta_scatter_sig_18->Draw();
  c1->SaveAs("scatter_Pt_eta_jpsi_2018.png");


  //Profile histogram of transverse momentum of Bc
  hprof_sig_18->SetLineColor(kRed);
  hprof_sig_18->GetXaxis()->SetTitle("P_{T}^{Bc} (GeV)");
  hprof_sig_18->GetYaxis()->SetTitle("Ratio of true and corrected P_{T}^{Bc}");
  hprof_sig_18->Draw();
  hprof_sig->SetLineColor(kGreen);
  hprof_sig->Draw("same");
  legend->Draw();
  c1->SaveAs("Profile_pt_bc.png");
  
  //Closing the files
  f1->Close();
  f3->Close();
  
}
