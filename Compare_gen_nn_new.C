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

void Compare_gen_nn_new(){

  //declaring all the required variable for signal channel for generator level 
  float bc_pt_sig_gen=0,bc_eta_sig_gen=0,bc_phi_sig_gen=0,bc_mass_sig_gen=0;
  float mu_pt_sig_gen=0,mu_eta_sig_gen=0,mu_phi_sig_gen=0;
  float mu1_pt_sig_gen=0,mu1_eta_sig_gen=0,mu1_phi_sig_gen=0;
  float mu2_pt_sig_gen=0,mu2_eta_sig_gen=0,mu2_phi_sig_gen=0;
  float jpsi_pt_sig_gen=0,jpsi_eta_sig_gen=0,jpsi_phi_sig_gen=0,jpsi_E_sig_gen=0;
  float Trig_dimuon_sig_gen=0,Trig_JpsiTk_sig_gen=0;
  float taunu_pt_sig_gen=0,taunu_eta_sig_gen=0,taunu_phi_sig_gen=0;
  float taunu2_pt_sig_gen=0,taunu2_eta_sig_gen=0,taunu2_phi_sig_gen=0;
  float munu_pt_sig_gen=0,munu_eta_sig_gen=0,munu_phi_sig_gen=0;

  //declaring all the required variable for normalization channel for generator level
  float bc_pt_norm_gen=0,bc_eta_norm_gen=0,bc_phi_norm_gen=0,bc_mass_norm_gen=0;
  float mu_pt_norm_gen=0,mu_eta_norm_gen=0,mu_phi_norm_gen=0;
  float mu1_pt_norm_gen=0,mu1_eta_norm_gen=0,mu1_phi_norm_gen=0;
  float mu2_pt_norm_gen=0,mu2_eta_norm_gen=0,mu2_phi_norm_gen=0;
  float jpsi_pt_norm_gen=0,jpsi_eta_norm_gen=0,jpsi_phi_norm_gen=0,jpsi_E_norm_gen=0;
  float Trig_dimuon_norm_gen=0,Trig_JpsiTk_norm_gen=0;
  float munu_pt_norm_gen=0,munu_eta_norm_gen=0,munu_phi_norm_gen=0;


  //declaring all the required variable for normalization channel of nn output
  //float bc_pt_norm_nn=0;
  float bc_eta_norm_nn=0,bc_phi_norm_nn=0,bc_mass_norm_nn=0;
  float mu_pt_norm_nn=0,mu_eta_norm_nn=0,mu_phi_norm_nn=0;
  float mu1_pt_norm_nn=0,mu1_eta_norm_nn=0,mu1_phi_norm_nn=0;
  float mu2_pt_norm_nn=0,mu2_eta_norm_nn=0,mu2_phi_norm_nn=0;
  float jpsi_pt_norm_nn=0,jpsi_eta_norm_nn=0,jpsi_phi_norm_nn=0,jpsi_mass_norm_nn=0;
  float Trig_dimuon_norm_nn=0,Trig_JpsiTk_norm_nn=0;
  float munu_pt_norm_nn=0,munu_eta_norm_nn=0,munu_phi_norm_nn=0;
  float bc_pt_norm_nn,bc_px_norm_nn,bc_py_norm_nn,bc_pz_norm_nn;
  float bc_pt_norm_recco;

  //declaring all the required variable for signal channel of nn output
  //float bc_pt_sig_nn=0;
  float bc_eta_sig_nn=0,bc_phi_sig_nn=0,bc_mass_sig_nn=0;
  float mu_pt_sig_nn=0,mu_eta_sig_nn=0,mu_phi_sig_nn=0;
  float mu1_pt_sig_nn=0,mu1_eta_sig_nn=0,mu1_phi_sig_nn=0;
  float mu2_pt_sig_nn=0,mu2_eta_sig_nn=0,mu2_phi_sig_nn=0;
  float jpsi_pt_sig_nn=0,jpsi_eta_sig_nn=0,jpsi_phi_sig_nn=0,jpsi_mass_sig_nn=0;
  float Trig_dimuon_sig_nn=0,Trig_JpsiTk_sig_nn=0;
  float taunu_pt_sig_nn=0,taunu_eta_sig_nn=0,taunu_phi_sig_nn=0;
  float taunu2_pt_sig_nn=0,taunu2_eta_sig_nn=0,taunu2_phi_sig_nn=0;
  float munu_pt_sig_nn=0,munu_eta_sig_nn=0,munu_phi_sig_nn=0;
  float bc_pt_sig_nn,bc_px_sig_nn,bc_py_sig_nn,bc_pz_sig_nn;
  float bc_pt_sig_recco;

 
  /**
  TFile *f1 = TFile::Open("bc_jpsi_mu_nu_gen_v2.root"); //opening the normalization root file for gen level
  TTree *tree_norm_gen = (TTree*)f1->Get("tree"); //reading the tree from the ntuple
  TFile *f2 = TFile::Open("bc_jpsi_tau_nu_gen_v2.root"); //opening the signal root file for gen level
  TTree *tree_sig_gen = (TTree*)f2->Get("tree"); //reading the tree from the ntuple
  **/

  /**
  TFile *f3 = TFile::Open("results-muon_channel-nodes_100_100_100_60-vfrac_0p4-dropoutrate_0p1-bsize_125.root"); //opening the normalization root file for nn output
  TTree *tree_norm_nn = (TTree*)f3->Get("tree"); //reading the tree from the ntuple
  TFile *f4 = TFile::Open("results-tau_channel-nodes_100_100_100_60-vfrac_0p4-dropoutrate_0p1-bsize_125.root"); //opening the signal root file for nn output
  TTree *tree_sig_nn = (TTree*)f4->Get("tree"); //reading the tree from the ntuple
  **/
  
  TFile *f3 = TFile::Open("results-muon_channel-nodes_100_100_100_60-vfrac_0p3-dropoutrate_0p1-bsize_200.root"); //opening the normalization root file for nn output
  TTree *tree_norm_nn = (TTree*)f3->Get("tree"); //reading the tree from the ntuple
  TFile *f4 = TFile::Open("results-tau_channel-nodes_100_100_100_60-vfrac_0p3-dropoutrate_0p1-bsize_200.root"); //opening the signal root file for nn output
  TTree *tree_sig_nn = (TTree*)f4->Get("tree"); //reading the tree from the ntuple
  

  //reading variables from the normalization tree for gen level
  tree_norm_nn->SetBranchAddress("gen_b_pt",&bc_pt_norm_gen); 
  tree_norm_nn->SetBranchAddress("gen_b_Eta",&bc_eta_norm_gen); 
  tree_norm_nn->SetBranchAddress("gen_b_Phi",&bc_phi_norm_gen); 
  //tree_norm_nn->SetBranchAddress("bc_mass",&bc_mass_norm_gen); 

  tree_norm_nn->SetBranchAddress("gen_mu_pt",&mu_pt_norm_gen); 
  tree_norm_nn->SetBranchAddress("gen_mu_Eta",&mu_eta_norm_gen); 
  tree_norm_nn->SetBranchAddress("gen_mu_Phi",&mu_phi_norm_gen);

  tree_norm_nn->SetBranchAddress("gen_jpsi_mu1_pt",&mu1_pt_norm_gen); 
  tree_norm_nn->SetBranchAddress("gen_jpsi_mu1_Eta",&mu1_eta_norm_gen); 
  tree_norm_nn->SetBranchAddress("gen_jpsi_mu1_Phi",&mu1_phi_norm_gen);

  tree_norm_nn->SetBranchAddress("gen_jpsi_mu2_pt",&mu2_pt_norm_gen); 
  tree_norm_nn->SetBranchAddress("gen_jpsi_mu2_Eta",&mu2_eta_norm_gen); 
  tree_norm_nn->SetBranchAddress("gen_jpsi_mu2_Phi",&mu2_phi_norm_gen);

  tree_norm_nn->SetBranchAddress("gen_jpsi_pt",&jpsi_pt_norm_gen); 
  tree_norm_nn->SetBranchAddress("gen_jpsi_Eta",&jpsi_eta_norm_gen); 
  tree_norm_nn->SetBranchAddress("gen_jpsi_Phi",&jpsi_phi_norm_gen);
  tree_norm_nn->SetBranchAddress("gen_jpsi_E",&jpsi_E_norm_gen);
  
  
  //reading variables from the signal tree for gen level
  tree_sig_nn->SetBranchAddress("gen_b_pt",&bc_pt_sig_gen); 
  tree_sig_nn->SetBranchAddress("gen_b_Eta",&bc_eta_sig_gen); 
  tree_sig_nn->SetBranchAddress("gen_b_Phi",&bc_phi_sig_gen); 
  //tree_sig_nn->SetBranchAddress("gen_b_mass",&bc_mass_sig_gen);

  tree_sig_nn->SetBranchAddress("gen_mu_pt",&mu_pt_sig_gen); 
  tree_sig_nn->SetBranchAddress("gen_mu_Eta",&mu_eta_sig_gen); 
  tree_sig_nn->SetBranchAddress("gen_mu_Phi",&mu_phi_sig_gen);

  tree_sig_nn->SetBranchAddress("gen_jpsi_mu1_pt",&mu1_pt_sig_gen); 
  tree_sig_nn->SetBranchAddress("gen_jpsi_mu1_Eta",&mu1_eta_sig_gen); 
  tree_sig_nn->SetBranchAddress("gen_jpsi_mu1_Phi",&mu1_phi_sig_gen);

  tree_sig_nn->SetBranchAddress("gen_jpsi_mu2_pt",&mu2_pt_sig_gen); 
  tree_sig_nn->SetBranchAddress("gen_jpsi_mu2_Eta",&mu2_eta_sig_gen); 
  tree_sig_nn->SetBranchAddress("gen_jpsi_mu2_Phi",&mu2_phi_sig_gen);

  tree_sig_nn->SetBranchAddress("gen_jpsi_pt",&jpsi_pt_sig_gen); 
  tree_sig_nn->SetBranchAddress("gen_jpsi_Eta",&jpsi_eta_sig_gen); 
  tree_sig_nn->SetBranchAddress("gen_jpsi_Phi",&jpsi_phi_sig_gen);
  tree_sig_nn->SetBranchAddress("gen_jpsi_E",&jpsi_E_sig_gen);
  

  //reading variables from the normalization tree for nn output
  //tree_norm_nn->SetBranchAddress("Bc_pt",&bc_pt_norm_nn); 
  tree_norm_nn->SetBranchAddress("Bc_eta",&bc_eta_norm_nn); 
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

  //tree_norm_nn->SetBranchAddress("bc_pt_predicted",&bc_pt_norm_nn);
  tree_norm_nn->SetBranchAddress("bc_px_predicted",&bc_px_norm_nn);
  tree_norm_nn->SetBranchAddress("bc_py_predicted",&bc_py_norm_nn);
  tree_norm_nn->SetBranchAddress("bc_pz_predicted",&bc_pz_norm_nn);
  tree_norm_nn->SetBranchAddress("Bc_pt",&bc_pt_norm_recco);
  //tree_norm_nn->SetBranchAddress("gen_b_pt",&bc_true_pt_norm_nn);
 

  //reading variables from the signal tree for nn output
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
  
  //tree_sig_nn->SetBranchAddress("bc_pt_predicted",&bc_pt_sig_nn);
  tree_sig_nn->SetBranchAddress("bc_px_predicted",&bc_px_sig_nn);
  tree_sig_nn->SetBranchAddress("bc_py_predicted",&bc_py_sig_nn);
  tree_sig_nn->SetBranchAddress("bc_pz_predicted",&bc_pz_sig_nn);
  tree_sig_nn->SetBranchAddress("Bc_pt",&bc_pt_sig_recco);
  //tree_sig_nn->SetBranchAddress("gen_b_pt",&bc_true_pt_sig_nn);
  

  //Defining Canvas, removing statistics and setting the legend
  auto c1 = new TCanvas("c", "c", 500,500);
  //gStyle->SetOptStat(0);
  auto* legend = new TLegend(0.7, 0.3, 0.9, 0.5);

  //booking the histograms for gen level
  TH1D* hmu_m2_miss_norm_gen = new TH1D("hmu_m2_miss_norm_gen","Missing mass squared",20,0,10); 
  TH1D* hmu_m2_miss_sig_gen = new TH1D("hmu_m2_miss_sig_gen","Missing mass squared",20,0,10);
  TH1D* hmu_Q2_norm_gen = new TH1D("hmu_Q2_norm_gen","Squared four momentum transfered to lepton system",20,0,10);
  TH1D* hmu_Q2_sig_gen = new TH1D("hmu_Q2_sig_gen","Squared four momentum transfered to lepton system",20,0,10);
  TH1D* hmu_Pt_miss_norm_gen = new TH1D("hmu_Pt_miss_norm_gen","Missing transverse momentum",25,-5,20);
  TH1D* hmu_Pt_miss_sig_gen = new TH1D("hmu_Pt_miss_sig_gen","Missing transverse momentum",25,-5,20);
  TH1D* hmu_m2_miss_corr_norm_gen = new TH1D("hmu_m2_miss_corr_norm_gen","Missing mass squared corrected",20,0,10); 
  TH1D* hmu_m2_miss_corr_sig_gen = new TH1D("hmu_m2_miss_corr_sig_gen","Missing mass squared corrected",20,0,10);
  TH1D* hmu_Q2_corr_norm_gen = new TH1D("hmu_Q2_corr_norm_gen","Squared four momentum transfered to lepton system corrected",20,0,10);
  TH1D* hmu_Q2_corr_sig_gen = new TH1D("hmu_Q2_corr_sig_gen","Squared four momentum transfered to lepton system corrected",20,0,10);
  TH1D* hmu_Pt_miss_corr_norm_gen = new TH1D("hmu_Pt_miss_corr_norm_gen","Missing transverse momentum corrected",25,-5,20);
  TH1D* hmu_Pt_miss_corr_sig_gen = new TH1D("hmu_Pt_miss_corr_sig_gen","Missing transverse momentum corrected",25,-5,20);
  TH1D* hmu_Pt_var_norm_gen = new TH1D("hmu_Pt_var_norm_gen","Tranverse variable momentum",60,-30,30);
  TH1D* hmu_Pt_var_sig_gen = new TH1D("hmu_Pt_var_sig_gen","Tranverse variable momentum",60,-30,30);
  TH1D* hmu_del_R_norm_gen = new TH1D("hmu_del_R_norm_gen","Delta R",50,0,3);
  TH1D* hmu_del_R_sig_gen = new TH1D("hmu_del_R_sig_gen","Delta R",50,0,3);
  TH1D* hmu_mu_pt_norm_gen = new TH1D("hmu_mu_pt_norm_gen","Unpaired muon tranverse momentum",60,0,20);
  TH1D* hmu_mu_pt_sig_gen = new TH1D("hmu_mu_pt_sig_gen","Unpaired muon tranverse momentum",60,0,20);
  TH1D* hmu_mu_eta_norm_gen = new TH1D("hmu_mu_eta_norm_gen","Unpaired muon psuedorapidity",60,-3,3);
  TH1D* hmu_mu_eta_sig_gen = new TH1D("hmu_mu_eta_sig_gen","Unpaired muon psuedorapidity",60,-3,3);
  TH1D* hmu_mu_phi_norm_gen = new TH1D("hmu_mu_phi_norm_gen","Unpaired muon phi",50,-5,5);
  TH1D* hmu_mu_phi_sig_gen = new TH1D("hmu_mu_phi_sig_gen","Unpaired muon phi",50,-5,5);
  TH1D* hmu_mu_en_bc_frame_norm_gen = new TH1D("hmu_mu_en_bc_frame_norm_gen","Energy of Unpaired muon in Bc rest frame",50,0,5);
  TH1D* hmu_mu_en_bc_frame_sig_gen = new TH1D("hmu_mu_en_bc_frame_sig_gen","Energy of Unpaired muon in Bc rest frame",50,0,5);
  TH1D* hmu_mu_en_jpsi_frame_norm_gen = new TH1D("hmu_mu_en_jpsi_frame_norm_gen","Energy of Unpaired muon in Jpsi rest frame",70,0,7);
  TH1D* hmu_mu_en_jpsi_frame_sig_gen = new TH1D("hmu_mu_en_jpsi_frame_sig_gen","Energy of Unpaired muon in Jpsi rest frame",70,0,7);
  TH1D* hmu_mu_en_bc_frame_corr_norm_gen = new TH1D("hmu_mu_en_bc_frame_corr_norm_gen","Energy of Unpaired muon in Bc rest frame",50,0,5);
  TH1D* hmu_mu_en_bc_frame_corr_sig_gen = new TH1D("hmu_mu_en_bc_frame_corr_sig_gen","Energy of Unpaired muon in Bc rest frame",50,0,5);
  TH1D* hmu_bc_pt_norm_gen = new TH1D("hmu_bc_pt_norm_gen","Transverse momentum of B_{c}",100,0,100);
  TH1D* hmu_bc_pt_sig_gen = new TH1D("hmu_bc_pt_sig_gen","Transverse momentum of B_{c}",100,0,100);
  TH1D* hmu_bc_pt_ratio_norm_gen = new TH1D("hmu_bc_pt_ratio_norm_gen","Ratio of corrected and true transverse momentum of B_{c}",40,0,2);
  TH1D* hmu_bc_pt_ratio_sig_gen = new TH1D("hmu_bc_pt_ratio_sig_gen","Ratio of corrected and true transverse momentum of B_{c}",40,0,2);
  TH1D* hmu_bc_pt_corr_norm_gen = new TH1D("hmu_bc_pt_corr_norm_gen","Corrected transverse momentum of B_{c}",100,0,100);
  TH1D* hmu_bc_pt_corr_sig_gen = new TH1D("hmu_bc_pt_corr_sig_gen","Corrected transverse momentum of B_{c}",100,0,100);
  auto hprof_norm_gen  = new TProfile("hprof_norm_gen","Profile of Transverse momentum of B_{c}",100,0,100,0,100);
  auto hprof_sig_gen  = new TProfile("hprof_sig_gen","Profile of Transverse momentum of B_{c}",100,0,100,0,100);
  TH2D* hmu_bc_pt_scatter_norm_gen = new TH2D("hmu_bc_pt_scatter_norm_gen","Scatter plot of transverse momentum of B_{c}",100,0,100,100,0,100);
  TH2D* hmu_bc_pt_scatter_sig_gen = new TH2D("hmu_bc_pt_scatter_sig_gen","Scatter plot of transverse momentum of B_{c}",100,0,100,100,0,100);
  TH1D* hmu_m2_miss_jpsi_norm_gen = new TH1D("hmu_m2_miss_jpsi_norm_gen","Missing mass squared for Jpsi",50,0,10); 
  TH1D* hmu_m2_miss_jpsi_sig_gen = new TH1D("hmu_m2_miss_jpsi_sig_gen","Missing mass squared for Jpsi",50,0,10);
  TH1D* hmu_bc_pt_norm_reso_corr = new TH1D("hmu_bc_pt_norm_reso_corr","B_{c} Pt resolution for corr",50,-1,1); 
  TH1D* hmu_bc_pt_sig_reso_corr = new TH1D("hmu_bc_pt_sig_reso_corr","B_{c} Pt resolution for corr",50,-1,1);
  TH1D* hmu_bc_px_norm_reso_corr = new TH1D("hmu_bc_px_norm_reso_corr","B_{c} Px resolution for corr",50,-1,1); 
  TH1D* hmu_bc_px_sig_reso_corr = new TH1D("hmu_bc_px_sig_reso_corr","B_{c} Px resolution for corr",50,-1,1);
  TH1D* hmu_bc_py_norm_reso_corr = new TH1D("hmu_bc_py_norm_reso_corr","B_{c} Py resolution for corr",50,-1,1); 
  TH1D* hmu_bc_py_sig_reso_corr = new TH1D("hmu_bc_py_sig_reso_corr","B_{c} Py resolution for corr",50,-1,1);
  TH1D* hmu_bc_pz_norm_reso_corr = new TH1D("hmu_bc_pz_norm_reso_corr","B_{c} Pz resolution for corr",50,-1,1); 
  TH1D* hmu_bc_pz_sig_reso_corr = new TH1D("hmu_bc_pz_sig_reso_corr","B_{c} Pz resolution for corr",50,-1,1);
  TH1D* hmu_m2_miss_norm_reso_corr = new TH1D("hmu_m2_miss_norm_reso_corr","Missing mass squared resolution for corr",50,-10,10); 
  TH1D* hmu_m2_miss_sig_reso_corr = new TH1D("hmu_m2_miss_sig_reso_corr","Missing mass squared resolution for corr",50,-10,10);
  TH1D* hmu_Q2_norm_reso_corr = new TH1D("hmu_Q2_norm_reso_corr","Squared four momentum transfered to lepton system resolution for corr",100,-1,5);
  TH1D* hmu_Q2_sig_reso_corr = new TH1D("hmu_Q2_sig_reso_corr","Squared four momentum transfered to lepton system resolution for corr",100,-1,2);
  TH1D* hmu_Pt_miss_norm_reso_corr = new TH1D("hmu_Pt_miss_norm_reso_corr","Missing transverse momentum resolution for corr",50,-2,8);
  TH1D* hmu_Pt_miss_sig_reso_corr = new TH1D("hmu_Pt_miss_sig_reso_corr","Missing transverse momentum resolution for corr",50,-2,8);
  auto hprof_bc_pt_norm_reso_corr  = new TProfile("hprof_bc_pt_norm_reso_corr","Profile of resolution of Bc Pt ",100,0,100,-100,100);
  auto hprof_bc_pt_sig_reso_corr  = new TProfile("hprof_bc_pt_sig_reso_corr","Profile of resolution of Bc Pt",100,0,100,-100,100);
  auto hprof_bc_px_norm_reso_corr  = new TProfile("hprof_bc_px_norm_reso_corr","Profile of resolution of Bc Px ",100,0,100,-100,100,"S");
  auto hprof_bc_px_sig_reso_corr  = new TProfile("hprof_bc_px_sig_reso_corr","Profile of resolution of Bc Px",100,0,100,-100,100,"S");
  auto hprof_bc_py_norm_reso_corr  = new TProfile("hprof_bc_py_norm_reso_corr","Profile of resolution of Bc Py ",100,0,100,-100,100);
  auto hprof_bc_py_sig_reso_corr  = new TProfile("hprof_bc_py_sig_reso_corr","Profile of resolution of Bc Py",100,0,100,-100,100);
  auto hprof_bc_pz_norm_reso_corr  = new TProfile("hprof_bc_pz_norm_reso_corr","Profile of resolution of Bc Pz ",100,0,100,-100,100);
  auto hprof_bc_pz_sig_reso_corr  = new TProfile("hprof_bc_pz_sig_reso_corr","Profile of resolution of Bc Pz",100,0,100,-100,100);
  auto hprof_m2_miss_norm_reso_corr  = new TProfile("hprof_m2_miss_norm_reso_corr","Profile of resolution of missing mass squared ",100,0,10,-100,100,"S");
  auto hprof_m2_miss_sig_reso_corr  = new TProfile("hprof_m2_miss_sig_reso_corr","Profile of resolution of missing mass squared",100,0,10,-100,100,"S");
  auto hprof_Q2_norm_reso_corr  = new TProfile("hprof_Q2_norm_reso_corr","Profile of resolution of Q2 ",100,0,10,-100,100,"S");
  auto hprof_Q2_sig_reso_corr  = new TProfile("hprof_Q2_sig_reso_corr","Profile of resolution of Q2",100,0,10,-100,100,"S");
  auto hprof_Pt_miss_norm_reso_corr  = new TProfile("hprof_Pt_miss_norm_reso_corr","Profile of resolution of Pt_miss ",100,-1,20,-100,100,"S");
  auto hprof_Pt_miss_sig_reso_corr  = new TProfile("hprof_Pt_miss_sig_reso_corr","Profile of resolution of Pt_miss",100,-1,20,-100,100,"S");
 
  
  //booking the histograms for nn level
  TH1D* hmu_m2_miss_norm_nn = new TH1D("hmu_m2_miss_norm_nn","Missing mass squared",20,0,10); 
  TH1D* hmu_m2_miss_sig_nn = new TH1D("hmu_m2_miss_sig_nn","Missing mass squared",20,0,10);
  TH1D* hmu_Q2_norm_nn = new TH1D("hmu_Q2_norm_nn","Squared four momentum transfered to lepton system",20,0,10);
  TH1D* hmu_Q2_sig_nn = new TH1D("hmu_Q2_sig_nn","Squared four momentum transfered to lepton system",20,0,10);
  TH1D* hmu_Pt_miss_norm_nn = new TH1D("hmu_Pt_miss_norm_nn","Missing transverse momentum",25,-5,20);
  TH1D* hmu_Pt_miss_sig_nn = new TH1D("hmu_Pt_miss_sig_nn","Missing transverse momentum",25,-5,20);
  TH1D* hmu_Pt_var_norm_nn = new TH1D("hmu_Pt_var_norm_nn","Tranverse variable momentum",60,-30,30);
  TH1D* hmu_Pt_var_sig_nn = new TH1D("hmu_Pt_var_sig_nn","Tranverse variable momentum",60,-30,30);
  TH1D* hmu_del_R_norm_nn = new TH1D("hmu_del_R_norm_nn","Delta R",60,0,3);
  TH1D* hmu_del_R_sig_nn = new TH1D("hmu_del_R_sig_nn","Delta R",60,0,3);
  TH1D* hmu_mu_pt_norm_nn = new TH1D("hmu_mu_pt_norm_nn","Unpaired muon tranverse momentum",60,0,20);
  TH1D* hmu_mu_pt_sig_nn = new TH1D("hmu_mu_pt_sig_nn","Unpaired muon tranverse momentum",60,0,20);
  TH1D* hmu_mu_eta_norm_nn = new TH1D("hmu_mu_eta_norm_nn","Unpaired muon psuedorapidity",60,-3,3);
  TH1D* hmu_mu_eta_sig_nn = new TH1D("hmu_mu_eta_sig_nn","Unpaired muon psuedorapidity",60,-3,3);
  TH1D* hmu_mu_phi_norm_nn = new TH1D("hmu_mu_phi_norm_nn","Unpaired muon phi",50,-5,5);
  TH1D* hmu_mu_phi_sig_nn = new TH1D("hmu_mu_phi_sig_nn","Unpaired muon phi",50,-5,5);
  TH1D* hmu_mu_en_bc_frame_norm_nn = new TH1D("hmu_mu_en_bc_frame_norm_nn","Energy of Unpaired muon in Bc rest frame",50,0,5);
  TH1D* hmu_mu_en_bc_frame_sig_nn = new TH1D("hmu_mu_en_bc_frame_sig_nn","Energy of Unpaired muon in Bc rest frame",50,0,5);
  TH1D* hmu_mu_en_jpsi_frame_norm_nn = new TH1D("hmu_mu_en_jpsi_frame_norm_nn","Energy of Unpaired muon in Jpsi rest frame",70,0,7);
  TH1D* hmu_mu_en_jpsi_frame_sig_nn = new TH1D("hmu_mu_en_jpsi_frame_sig_nn","Energy of Unpaired muon in Jpsi rest frame",70,0,7);
  TH1D* hmu_bc_pt_norm_nn = new TH1D("hmu_bc_pt_norm_nn","Transverse momentum of B_{c}",100,0,100);
  TH1D* hmu_bc_pt_sig_nn = new TH1D("hmu_bc_pt_sig_nn","Transverse momentum of B_{c}",100,0,100);
  TH1D* hmu_bc_pt_ratio_norm_nn = new TH1D("hmu_bc_pt_ratio_norm_nn","Ratio of predicted_nn and true transverse momentum of B_{c}",40,0,2);
  TH1D* hmu_bc_pt_ratio_sig_nn = new TH1D("hmu_bc_pt_ratio_sig_nn","Ratio of predicted_nn and true transverse momentum of B_{c}",40,0,2);
  auto hprof_norm_nn  = new TProfile("hprof_norm_nn","Profile of Transverse momentum of B_{c}",100,0,100,0,100);
  auto hprof_sig_nn  = new TProfile("hprof_sig_nn","Profile of Transverse momentum of B_{c}",100,0,100,0,100);
  TH2D* hmu_bc_pt_scatter_norm_nn = new TH2D("hmu_bc_pt_scatter_norm_nn","Scatter plot of transverse momentum of B_{c}",100,0,100,100,0,100);
  TH2D* hmu_bc_pt_scatter_sig_nn = new TH2D("hmu_bc_pt_scatter_sig_nn","Scatter plot of transverse momentum of B_{c}",100,0,100,100,0,100);
  TH1D* hmu_m2_miss_jpsi_norm_nn = new TH1D("hmu_m2_miss_jpsi_norm_nn","Missing mass Jpsi squared",50,0,10); 
  TH1D* hmu_m2_miss_jpsi_sig_nn = new TH1D("hmu_m2_miss_jpsi_sig_nn","Missing mass Jpsi squared",50,0,10);
  TH1D* hmu_bc_pt_norm_recco = new TH1D("hmu_bc_pt_norm_recco","Transverse momentum of B_{c}",100,0,100);
  TH1D* hmu_bc_pt_sig_recco = new TH1D("hmu_bc_pt_sig_recco","Transverse momentum of B_{c}",100,0,100);
  TH1D* hmu_m2_miss_norm_recco = new TH1D("hmu_m2_miss_norm_recco","Missing mass squared",100,0,10); 
  TH1D* hmu_m2_miss_sig_recco = new TH1D("hmu_m2_miss_sig_recco","Missing mass squared",100,0,10);
  TH1D* hmu_Q2_norm_recco = new TH1D("hmu_Q2_norm_recco","Squared four momentum transfered to lepton system",200,-5,15);
  TH1D* hmu_Q2_sig_recco = new TH1D("hmu_Q2_sig_recco","Squared four momentum transfered to lepton system",200,-5,15);
  TH1D* hmu_Pt_miss_norm_recco = new TH1D("hmu_Pt_miss_norm_recco","Missing transverse momentum",250,-5,20);
  TH1D* hmu_Pt_miss_sig_recco = new TH1D("hmu_Pt_miss_sig_recco","Missing transverse momentum",250,-5,20);
  TH1D* hmu_ratio_m2_miss_nn = new TH1D("hmu_ratio_m2_miss_nn","Ratio of signal and total Missing mass square",250,-5,20);
  TH1D* hmu_bc_pt_norm_reso_nn = new TH1D("hmu_bc_pt_norm_reso_nn","B_{c} Pt resolution for nn",50,-1,1); 
  TH1D* hmu_bc_pt_sig_reso_nn = new TH1D("hmu_bc_pt_sig_reso_nn","B_{c} Pt resolution for nn",50,-1,1);
  TH1D* hmu_bc_px_norm_reso_nn = new TH1D("hmu_bc_px_norm_reso_nn","B_{c} Px resolution for nn",50,-1,1); 
  TH1D* hmu_bc_px_sig_reso_nn = new TH1D("hmu_bc_px_sig_reso_nn","B_{c} Px resolution for nn",50,-1,1);
  TH1D* hmu_bc_py_norm_reso_nn = new TH1D("hmu_bc_py_norm_reso_nn","B_{c} Py resolution for nn",50,-1,1); 
  TH1D* hmu_bc_py_sig_reso_nn = new TH1D("hmu_bc_py_sig_reso_nn","B_{c} Py resolution for nn",50,-1,1);
  TH1D* hmu_bc_pz_norm_reso_nn = new TH1D("hmu_bc_pz_norm_reso_nn","B_{c} Pz resolution for nn",50,-1,1); 
  TH1D* hmu_bc_pz_sig_reso_nn = new TH1D("hmu_bc_pz_sig_reso_nn","B_{c} Pz resolution for nn",50,-1,1);
  TH1D* hmu_m2_miss_norm_reso_nn = new TH1D("hmu_m2_miss_norm_reso_nn","Missing mass squared resolution for nn",50,-10,10); 
  TH1D* hmu_m2_miss_sig_reso_nn = new TH1D("hmu_m2_miss_sig_reso_nn","Missing mass squared resolution for nn",50,-10,10);
  TH1D* hmu_Q2_norm_reso_nn = new TH1D("hmu_Q2_norm_reso_nn","Squared four momentum transfered to lepton system resolution for nn",100,-1,5);
  TH1D* hmu_Q2_sig_reso_nn = new TH1D("hmu_Q2_sig_reso_nn","Squared four momentum transfered to lepton system resolution for nn",100,-1,2);
  TH1D* hmu_Pt_miss_norm_reso_nn = new TH1D("hmu_Pt_miss_norm_reso_nn","Missing transverse momentum resolution for nn",50,-2,8);
  TH1D* hmu_Pt_miss_sig_reso_nn = new TH1D("hmu_Pt_miss_sig_reso_nn","Missing transverse momentum resolution for nn",50,-2,8);
  auto hprof_bc_pt_norm_reso_nn  = new TProfile("hprof_bc_pt_norm_reso_nn","Profile of resolution of Bc Pt ",100,0,100,-100,100);
  auto hprof_bc_pt_sig_reso_nn  = new TProfile("hprof_bc_pt_sig_reso_nn","Profile of resolution of Bc Pt",100,0,100,-100,100);
  auto hprof_bc_px_norm_reso_nn  = new TProfile("hprof_bc_px_norm_reso_nn","Profile of resolution of Bc Px ",100,0,100,-100,100,"S");
  auto hprof_bc_px_sig_reso_nn  = new TProfile("hprof_bc_px_sig_reso_nn","Profile of resolution of Bc Px",100,0,100,-100,100,"S");
  auto hprof_bc_py_norm_reso_nn  = new TProfile("hprof_bc_py_norm_reso_nn","Profile of resolution of Bc Py ",100,0,100,-100,100);
  auto hprof_bc_py_sig_reso_nn  = new TProfile("hprof_bc_py_sig_reso_nn","Profile of resolution of Bc Py",100,0,100,-100,100);
  auto hprof_bc_pz_norm_reso_nn  = new TProfile("hprof_bc_pz_norm_reso_nn","Profile of resolution of Bc Pz ",100,0,100,-100,100);
  auto hprof_bc_pz_sig_reso_nn  = new TProfile("hprof_bc_pz_sig_reso_nn","Profile of resolution of Bc Pz",100,0,100,-100,100);
  auto hprof_m2_miss_norm_reso_nn  = new TProfile("hprof_m2_miss_norm_reso_nn","Profile of resolution of missing mass squared ",100,0,10,-100,100,"S");
  auto hprof_m2_miss_sig_reso_nn  = new TProfile("hprof_m2_miss_sig_reso_nn","Profile of resolution of missing mass squared",100,0,10,-100,100,"S");
  auto hprof_Q2_norm_reso_nn  = new TProfile("hprof_Q2_norm_reso_nn","Profile of resolution of Q2 ",100,0,10,-100,100,"s");
  auto hprof_Q2_sig_reso_nn  = new TProfile("hprof_Q2_sig_reso_nn","Profile of resolution of Q2",100,0,10,-100,100,"s");
  auto hprof_Pt_miss_norm_reso_nn  = new TProfile("hprof_Pt_miss_norm_reso_nn","Profile of resolution of Pt_miss ",100,-1,20,-100,100,"S");
  auto hprof_Pt_miss_sig_reso_nn  = new TProfile("hprof_Pt_miss_reso_sig_nn","Profile of resolution of Pt_miss",100,-1,20,-100,100,"S");
  TH2D* hmu_m2_miss_scatter_norm_nn = new TH2D("hmu_m2_miss_scatter_norm_nn","Scatter plot of m2 miss",100,-5,10,100,-5,10);
  TH2D* hmu_m2_miss_scatter_sig_nn = new TH2D("hmu_m2_miss_scatter_sig_nn","Scatter plot of m2 miss",100,-5,10,100,-5,10);
  TH2D* hmu_Q2_scatter_norm_nn = new TH2D("hmu_Q2_scatter_norm_nn","Scatter plot of Q2",100,-5,20,100,-5,20);
  TH2D* hmu_Q2_scatter_sig_nn = new TH2D("hmu_Q2_scatter_sig_nn","Scatter plot of Q2",100,-5,20,100,-5,20);
  TH2D* hmu_Pt_miss_scatter_norm_nn = new TH2D("hmu_Pt_miss_scatter_norm_nn","Scatter plot of Pt miss",100,-5,20,100,-5,20);
  TH2D* hmu_Pt_miss_scatter_sig_nn = new TH2D("hmu_Pt_miss_scatter_sig_nn","Scatter plot of Pt miss",100,-5,20,100,-5,20);

  TH1D* hmu_w_ratio_Q2_reso_nn = new TH1D("hmu_w_ratio_Q2_reso_nn","Resolution plot of w_ratio_Q2 for nn",50,-1,4);
  TH1D* hmu_w_ratio_Q2_reso_corr = new TH1D("hmu_w_ratio_Q2_reso_corr","Resolution plot of w_ratio_Q2 for Jona",50,-1,4);
  TH2D* hmu_w_ratio_Q2_scatter_nn = new TH2D("hmu_w_ratio_Q2_scatter_nn","Scatter plot of w_ratio_Q2 for nn",100,0,10,100,0,10);
  TH2D* hmu_w_ratio_Q2_scatter_corr = new TH2D("hmu_w_ratio_Q2_scatter_corr","Scatter plot of w_ratio_Q2 for Jona",100,0,10,100,0,10);

  
  //ParticleMass muonMass = 0.10565837;
  float mu_mass = 0.113428;
  float bcPdgMass = 6.2756;
  float tauPdgMass = 1.77682;

  
 for(int i=0; i<tree_norm_nn->GetEntries();i++){
    tree_norm_nn->GetEntry(i);

    //Generator level
    //TLorentzVector class for calculating four momenta
    TLorentzVector momntm_bc_norm_gen;
    momntm_bc_norm_gen.SetPtEtaPhiM(bc_pt_norm_gen, bc_eta_norm_gen,bc_phi_norm_gen,bcPdgMass);
    TLorentzVector momntm_mu1_norm_gen;
    momntm_mu1_norm_gen.SetPtEtaPhiM(mu1_pt_norm_gen, mu1_eta_norm_gen,mu1_phi_norm_gen,mu_mass);
    TLorentzVector momntm_mu2_norm_gen;
    momntm_mu2_norm_gen.SetPtEtaPhiM(mu2_pt_norm_gen, mu2_eta_norm_gen,mu2_phi_norm_gen,mu_mass);
    TLorentzVector momntm_mu_norm_gen;
    momntm_mu_norm_gen.SetPtEtaPhiM(mu_pt_norm_gen, mu_eta_norm_gen,mu_phi_norm_gen,mu_mass);
    TLorentzVector momntm_jpsi_norm_gen;
    momntm_jpsi_norm_gen.SetPtEtaPhiE(jpsi_pt_norm_gen, jpsi_eta_norm_gen,jpsi_phi_norm_gen,jpsi_E_norm_gen);

    float bc_px_norm_gen = momntm_bc_norm_gen.Px();
    float bc_py_norm_gen = momntm_bc_norm_gen.Py();
    float bc_pz_norm_gen = momntm_bc_norm_gen.Pz();
    
    //calculating the kinematic variables
    float m2_miss_norm_gen =(momntm_bc_norm_gen - momntm_mu1_norm_gen - momntm_mu2_norm_gen - momntm_mu_norm_gen).M2();
    float Q2_norm_gen = (momntm_bc_norm_gen - momntm_mu1_norm_gen - momntm_mu2_norm_gen).M2();
    float Pt_miss_norm_gen = (bc_pt_norm_gen - mu1_pt_norm_gen - mu2_pt_norm_gen - mu_pt_norm_gen);
    float Pt_var_norm_gen = (mu1_pt_norm_gen + mu2_pt_norm_gen - mu_pt_norm_gen );
    
    float del_R_norm_gen = momntm_mu1_norm_gen.DeltaR(momntm_mu2_norm_gen);
    float m2_miss_jpsi_norm_gen =(momntm_jpsi_norm_gen - momntm_mu1_norm_gen - momntm_mu2_norm_gen).M2();


    TVector3 X_bc_norm_gen = -momntm_bc_norm_gen.BoostVector();
    TVector3 X_jpsi_norm_gen = -momntm_jpsi_norm_gen.BoostVector();
    TLorentzVector momntm_mu_bc_frame_norm_gen = momntm_mu_norm_gen;
    momntm_mu_bc_frame_norm_gen.Boost(X_bc_norm_gen);
    TLorentzVector momntm_mu_jpsi_frame_norm_gen = momntm_mu_norm_gen;
    momntm_mu_jpsi_frame_norm_gen.Boost(X_jpsi_norm_gen);

    float mu_en_bc_frame_norm_gen = momntm_mu_bc_frame_norm_gen.E();
    float mu_en_jpsi_frame_norm_gen = momntm_mu_jpsi_frame_norm_gen.E();

    //Neural network
    //TLorentzVector class for calculating four momenta
    float bc_e_norm_nn = sqrt(pow(bcPdgMass,2)+pow(bc_px_norm_nn,2)+pow(bc_py_norm_nn,2)+pow(bc_pz_norm_nn,2));
    TLorentzVector momntm_bc_norm_nn;
    //momntm_bc_norm_nn.SetPtEtaPhiM(bc_pt_norm_nn, bc_eta_norm_nn, bc_phi_norm_nn,bcPdgMass);
    momntm_bc_norm_nn.SetPxPyPzE(bc_px_norm_nn, bc_py_norm_nn,bc_pz_norm_nn, bc_e_norm_nn);
    TLorentzVector momntm_mu1_norm_nn;
    momntm_mu1_norm_nn.SetPtEtaPhiM(mu1_pt_norm_nn, mu1_eta_norm_nn,mu1_phi_norm_nn,mu_mass);
    TLorentzVector momntm_mu2_norm_nn;
    momntm_mu2_norm_nn.SetPtEtaPhiM(mu2_pt_norm_nn, mu2_eta_norm_nn,mu2_phi_norm_nn,mu_mass);
    TLorentzVector momntm_mu_norm_nn;
    momntm_mu_norm_nn.SetPtEtaPhiM(mu_pt_norm_nn, mu_eta_norm_nn,mu_phi_norm_nn,mu_mass);
    TLorentzVector momntm_jpsi_norm_nn;
    momntm_jpsi_norm_nn.SetPtEtaPhiM(jpsi_pt_norm_nn, jpsi_eta_norm_nn,jpsi_phi_norm_nn,jpsi_mass_norm_nn);

    float bc_pt_norm_nn = momntm_bc_norm_nn.Pt();
    //calculating the kinematic variables
    float m2_miss_norm_nn =(momntm_bc_norm_nn - momntm_mu1_norm_nn - momntm_mu2_norm_nn - momntm_mu_norm_nn).M2();
    float Q2_norm_nn = (momntm_bc_norm_nn - momntm_mu1_norm_nn - momntm_mu2_norm_nn).M2();
    float Pt_miss_norm_nn = (bc_pt_norm_nn - mu1_pt_norm_nn - mu2_pt_norm_nn - mu_pt_norm_nn);
    float Pt_var_norm_nn = (mu1_pt_norm_nn + mu2_pt_norm_nn - mu_pt_norm_nn );
    float del_R_norm_nn = momntm_mu1_norm_nn.DeltaR(momntm_mu2_norm_nn);

    float bc_pt_ratio_norm_nn = bc_pt_norm_nn/bc_pt_norm_gen;
    float m2_miss_jpsi_norm_nn =(momntm_jpsi_norm_nn - momntm_mu1_norm_nn - momntm_mu2_norm_nn).M2();
    
    TVector3 X_bc_norm_nn = -momntm_bc_norm_nn.BoostVector();
    TVector3 X_jpsi_norm_nn = -momntm_jpsi_norm_nn.BoostVector();
    TLorentzVector momntm_mu_bc_frame_norm_nn = momntm_mu_norm_nn;
    momntm_mu_bc_frame_norm_nn.Boost(X_bc_norm_nn);
    TLorentzVector momntm_mu_jpsi_frame_norm_nn = momntm_mu_norm_nn;
    momntm_mu_jpsi_frame_norm_nn.Boost(X_jpsi_norm_nn);

    float mu_en_bc_frame_norm_nn = momntm_mu_bc_frame_norm_nn.E();
    float mu_en_jpsi_frame_norm_nn = momntm_mu_jpsi_frame_norm_nn.E();

    //reconstructed
    TLorentzVector momntm_bc_norm_recco;
    momntm_bc_norm_recco.SetPtEtaPhiM(bc_pt_norm_recco, bc_eta_norm_nn,bc_phi_norm_nn,bcPdgMass);
    float m2_miss_norm_recco =(momntm_bc_norm_recco - momntm_mu1_norm_nn - momntm_mu2_norm_nn - momntm_mu_norm_nn).M2();
    float Q2_norm_recco = (momntm_bc_norm_recco - momntm_mu1_norm_nn - momntm_mu2_norm_nn).M2();
    float Pt_miss_norm_recco = (bc_pt_norm_recco - mu1_pt_norm_nn - mu2_pt_norm_nn - mu_pt_norm_nn);
    

    //jona recipe
    float bc_pt_corr_norm_gen = (bcPdgMass*bc_pt_norm_recco)/(bc_mass_norm_nn); //corrected Pt according to recipe in Jona thesis
    float bc_pt_ratio_norm_gen = bc_pt_corr_norm_gen/bc_pt_norm_gen;
    
    TLorentzVector momntm_bc_corr_norm_gen;
    momntm_bc_corr_norm_gen.SetPtEtaPhiM(bc_pt_corr_norm_gen, bc_eta_norm_nn, bc_phi_norm_nn,bcPdgMass);
    float Q2_corr_norm_gen = (momntm_bc_corr_norm_gen - momntm_mu1_norm_nn - momntm_mu2_norm_nn).M2();
    float m2_miss_corr_norm_gen = (momntm_bc_corr_norm_gen - momntm_mu1_norm_nn - momntm_mu2_norm_nn - momntm_mu_norm_nn).M2();
    float Pt_miss_corr_norm_gen = (bc_pt_corr_norm_gen - mu1_pt_norm_nn - mu2_pt_norm_nn - mu_pt_norm_nn);

    TVector3 X_bc_corr_norm_gen = -momntm_bc_corr_norm_gen.BoostVector();
    TLorentzVector momntm_mu_bc_frame_corr_norm_gen = momntm_mu_norm_nn;
    momntm_mu_bc_frame_corr_norm_gen.Boost(X_bc_norm_gen);
  
    float mu_en_bc_frame_corr_norm_gen = momntm_mu_bc_frame_corr_norm_gen.E();

    float bc_px_corr_norm_gen =  momntm_bc_corr_norm_gen.Px();
    float bc_py_corr_norm_gen =  momntm_bc_corr_norm_gen.Py();
    float bc_pz_corr_norm_gen =  momntm_bc_corr_norm_gen.Pz();

    //resolution variables

    float bc_pt_norm_reso_nn = (bc_pt_norm_nn - bc_pt_norm_gen)/bc_pt_norm_gen;
    float bc_px_norm_reso_nn = (bc_px_norm_nn - bc_px_norm_gen)/bc_px_norm_gen;
    float bc_py_norm_reso_nn = (bc_py_norm_nn - bc_py_norm_gen)/bc_py_norm_gen;
    float bc_pz_norm_reso_nn = (bc_pz_norm_nn - bc_pz_norm_gen)/bc_pz_norm_gen;
    float m2_miss_norm_reso_nn = (m2_miss_norm_nn - m2_miss_norm_gen);
    float Q2_norm_reso_nn = (Q2_norm_nn - Q2_norm_gen)/Q2_norm_gen;
    float Pt_miss_norm_reso_nn = (Pt_miss_norm_nn - Pt_miss_norm_gen)/Pt_miss_norm_gen;

    float bc_pt_norm_reso_corr = (bc_pt_corr_norm_gen - bc_pt_norm_gen)/bc_pt_norm_gen;
    float bc_px_norm_reso_corr = (bc_px_corr_norm_gen - bc_px_norm_gen)/bc_px_norm_gen;
    float bc_py_norm_reso_corr = (bc_py_corr_norm_gen - bc_py_norm_gen)/bc_py_norm_gen;
    float bc_pz_norm_reso_corr = (bc_pz_corr_norm_gen - bc_pz_norm_gen)/bc_pz_norm_gen;
    float m2_miss_norm_reso_corr = (m2_miss_corr_norm_gen - m2_miss_norm_gen);
    float Q2_norm_reso_corr = (Q2_corr_norm_gen - Q2_norm_gen)/Q2_norm_gen;
    float Pt_miss_norm_reso_corr = (Pt_miss_corr_norm_gen - Pt_miss_norm_gen)/Pt_miss_norm_gen;

    //weight ratio
    float wmu_Q2_gen = pow((1-(pow(mu_mass,2)/Q2_norm_gen)),2)*(1+(pow(mu_mass,2)/(2*Q2_norm_gen)));
    float wtau_Q2_gen = pow((1-(pow(tauPdgMass,2)/Q2_norm_gen)),2)*(1+(pow(tauPdgMass,2)/(2*Q2_norm_gen)));
    float w_ratio_Q2_gen = wtau_Q2_gen/wmu_Q2_gen;

    float wmu_Q2_corr = pow((1-(pow(mu_mass,2)/Q2_corr_norm_gen)),2)*(1+(pow(mu_mass,2)/(2*Q2_corr_norm_gen)));
    float wtau_Q2_corr = pow((1-(pow(tauPdgMass,2)/Q2_corr_norm_gen)),2)*(1+(pow(tauPdgMass,2)/(2*Q2_corr_norm_gen)));
    float w_ratio_Q2_corr = wtau_Q2_corr/wmu_Q2_corr;

    float wmu_Q2_nn = pow((1-(pow(mu_mass,2)/Q2_norm_nn)),2)*(1+(pow(mu_mass,2)/(2*Q2_norm_nn)));
    float wtau_Q2_nn = pow((1-(pow(tauPdgMass,2)/Q2_norm_nn)),2)*(1+(pow(tauPdgMass,2)/(2*Q2_norm_nn)));
    float w_ratio_Q2_nn = wtau_Q2_nn/wmu_Q2_nn;

    float w_ratio_Q2_reso_nn = (w_ratio_Q2_nn - w_ratio_Q2_gen)/w_ratio_Q2_gen;
    float w_ratio_Q2_reso_corr = (w_ratio_Q2_corr - w_ratio_Q2_gen)/w_ratio_Q2_gen;


    //Filling the histograms

    //if ((abs(mu_eta_norm_gen) < 1.3 and mu_pt_norm_gen > 3.3) or (1.3 < abs(mu_eta_norm_gen) < 2.2 and mu_pt_norm_gen > 2.9) or (2.2 < abs(mu_eta_norm_gen) < 2.4 and mu_pt_norm_gen > 2.4)){ //acceptance cuts
      //if ((Trig_dimuon_norm_nn == true) and (Trig_JpsiTk_norm_nn == true)){ //trigger flags
        
    	hmu_m2_miss_norm_gen->Fill(m2_miss_norm_gen);
	hmu_Q2_norm_gen->Fill(Q2_norm_gen);
	hmu_Pt_miss_norm_gen->Fill(Pt_miss_norm_gen);
	hmu_m2_miss_corr_norm_gen->Fill(m2_miss_corr_norm_gen);
	hmu_Q2_corr_norm_gen->Fill(Q2_corr_norm_gen);
	hmu_Pt_miss_corr_norm_gen->Fill(Pt_miss_corr_norm_gen);
	hmu_Pt_var_norm_gen->Fill(Pt_var_norm_gen);
	hmu_mu_pt_norm_gen->Fill(mu_pt_norm_gen);
	hmu_mu_eta_norm_gen->Fill(mu_eta_norm_gen);
	hmu_mu_phi_norm_gen->Fill(mu_phi_norm_gen);
	hmu_del_R_norm_gen->Fill(del_R_norm_gen);
	hmu_mu_en_bc_frame_norm_gen->Fill(mu_en_bc_frame_norm_gen);
	hmu_mu_en_jpsi_frame_norm_gen->Fill(mu_en_jpsi_frame_norm_gen);
	hmu_mu_en_bc_frame_corr_norm_gen->Fill(mu_en_bc_frame_corr_norm_gen);
	hmu_bc_pt_norm_gen->Fill(bc_pt_norm_gen);
	hmu_bc_pt_corr_norm_gen->Fill(bc_pt_corr_norm_gen);
	hmu_bc_pt_ratio_norm_gen->Fill(bc_pt_ratio_norm_gen);
	hprof_norm_gen->Fill(bc_pt_norm_gen,bc_pt_ratio_norm_gen,1);
	hmu_bc_pt_scatter_norm_gen->Fill(bc_pt_norm_gen,bc_pt_corr_norm_gen);
	hmu_m2_miss_jpsi_norm_gen->Fill(m2_miss_jpsi_norm_gen);


	hmu_m2_miss_norm_nn->Fill(m2_miss_norm_nn);
	hmu_Q2_norm_nn->Fill(Q2_norm_nn);
	hmu_Pt_miss_norm_nn->Fill(Pt_miss_norm_nn);
	hmu_Pt_var_norm_nn->Fill(Pt_var_norm_nn);
	hmu_mu_pt_norm_nn->Fill(mu_pt_norm_nn);
	hmu_mu_eta_norm_nn->Fill(mu_eta_norm_nn);
	hmu_mu_phi_norm_nn->Fill(mu_phi_norm_nn);
	hmu_del_R_norm_nn->Fill(del_R_norm_nn);
	hmu_mu_en_bc_frame_norm_nn->Fill(mu_en_bc_frame_norm_nn);
	hmu_mu_en_jpsi_frame_norm_nn->Fill(mu_en_jpsi_frame_norm_nn);
	hmu_bc_pt_norm_nn->Fill(bc_pt_norm_nn);
	hmu_bc_pt_norm_recco->Fill(bc_pt_norm_recco);
	hmu_bc_pt_ratio_norm_nn->Fill(bc_pt_ratio_norm_nn);
	hprof_norm_nn->Fill(bc_pt_norm_gen,bc_pt_ratio_norm_nn,1);
	hmu_bc_pt_scatter_norm_nn->Fill(bc_pt_norm_gen,bc_pt_norm_nn);
	hmu_m2_miss_jpsi_norm_nn->Fill(m2_miss_jpsi_norm_nn);
	hmu_m2_miss_norm_recco->Fill(m2_miss_norm_recco);
	hmu_Q2_norm_recco->Fill(Q2_norm_recco);
	hmu_Pt_miss_norm_recco->Fill(Pt_miss_norm_recco);

	
	hmu_bc_pt_norm_reso_nn->Fill(bc_pt_norm_reso_nn);
	hmu_bc_px_norm_reso_nn->Fill(bc_px_norm_reso_nn);
	hmu_bc_py_norm_reso_nn->Fill(bc_py_norm_reso_nn);
	hmu_bc_pz_norm_reso_nn->Fill(bc_pz_norm_reso_nn);
	hmu_m2_miss_norm_reso_nn->Fill(m2_miss_norm_reso_nn);
	hmu_Q2_norm_reso_nn->Fill(Q2_norm_reso_nn);
	hmu_Pt_miss_norm_reso_nn->Fill(Pt_miss_norm_reso_nn);
	hprof_bc_pt_norm_reso_nn->Fill(bc_pt_norm_gen,bc_pt_norm_reso_nn,1);
	hprof_bc_px_norm_reso_nn->Fill(bc_px_norm_gen,bc_px_norm_reso_nn,1);
	hprof_bc_py_norm_reso_nn->Fill(bc_py_norm_gen,bc_py_norm_reso_nn,1);
	hprof_bc_pz_norm_reso_nn->Fill(bc_pz_norm_gen,bc_pz_norm_reso_nn,1);
	hprof_m2_miss_norm_reso_nn->Fill(m2_miss_norm_gen,m2_miss_norm_reso_nn,1);
	hprof_Q2_norm_reso_nn->Fill(Q2_norm_gen,Q2_norm_reso_nn,1);
	hprof_Pt_miss_norm_reso_nn->Fill(Pt_miss_norm_gen,Pt_miss_norm_reso_nn,1);

	hmu_bc_pt_norm_reso_corr->Fill(bc_pt_norm_reso_corr);
	hmu_bc_px_norm_reso_corr->Fill(bc_px_norm_reso_corr);
	hmu_bc_py_norm_reso_corr->Fill(bc_py_norm_reso_corr);
	hmu_bc_pz_norm_reso_corr->Fill(bc_pz_norm_reso_corr);
	hmu_m2_miss_norm_reso_corr->Fill(m2_miss_norm_reso_corr);
	hmu_Q2_norm_reso_corr->Fill(Q2_norm_reso_corr);
	hmu_Pt_miss_norm_reso_corr->Fill(Pt_miss_norm_reso_corr);
	hprof_bc_pt_norm_reso_corr->Fill(bc_pt_norm_gen,bc_pt_norm_reso_corr,1);
	hprof_bc_px_norm_reso_corr->Fill(bc_px_norm_gen,bc_px_norm_reso_corr,1);
	hprof_bc_py_norm_reso_corr->Fill(bc_py_norm_gen,bc_py_norm_reso_corr,1);
	hprof_bc_pz_norm_reso_corr->Fill(bc_pz_norm_gen,bc_pz_norm_reso_corr,1);
	hprof_m2_miss_norm_reso_corr->Fill(m2_miss_norm_gen,m2_miss_norm_reso_corr,1);
	hprof_Q2_norm_reso_corr->Fill(Q2_norm_gen,Q2_norm_reso_corr,1);
	hprof_Pt_miss_norm_reso_corr->Fill(Pt_miss_norm_gen,Pt_miss_norm_reso_corr,1);

	hmu_m2_miss_scatter_norm_nn->Fill(m2_miss_norm_gen,m2_miss_norm_nn);
	hmu_Q2_scatter_norm_nn->Fill(Q2_norm_gen,Q2_norm_nn);
	hmu_Pt_miss_scatter_norm_nn->Fill(Pt_miss_norm_gen,Pt_miss_norm_nn);

	//if(Q2_norm_gen != 0){
	  if(Q2_norm_gen > 3){
	    hmu_w_ratio_Q2_reso_nn->Fill(w_ratio_Q2_reso_nn);
	    hmu_w_ratio_Q2_reso_corr->Fill(w_ratio_Q2_reso_corr);
	    hmu_w_ratio_Q2_scatter_nn->Fill(w_ratio_Q2_gen,w_ratio_Q2_nn);
	    hmu_w_ratio_Q2_scatter_corr->Fill(w_ratio_Q2_gen,w_ratio_Q2_corr);
	    //}
	}
	//}
	//}
 }


  //loop over events for signal channel for nn output
  for(int i=0; i<tree_sig_nn->GetEntries();i++){
    tree_sig_nn->GetEntry(i);


    //Gen level
    //TLorentzVector class for calculating four momenta
    TLorentzVector momntm_bc_sig_gen;
    momntm_bc_sig_gen.SetPtEtaPhiM(bc_pt_sig_gen, bc_eta_sig_gen,bc_phi_sig_gen,bcPdgMass);
    TLorentzVector momntm_mu1_sig_gen;
    momntm_mu1_sig_gen.SetPtEtaPhiM(mu1_pt_sig_gen, mu1_eta_sig_gen,mu1_phi_sig_gen,mu_mass);
    TLorentzVector momntm_mu2_sig_gen;
    momntm_mu2_sig_gen.SetPtEtaPhiM(mu2_pt_sig_gen, mu2_eta_sig_gen,mu2_phi_sig_gen,mu_mass);
    TLorentzVector momntm_mu_sig_gen;
    momntm_mu_sig_gen.SetPtEtaPhiM(mu_pt_sig_gen, mu_eta_sig_gen,mu_phi_sig_gen,mu_mass);
    TLorentzVector momntm_jpsi_sig_gen;
    momntm_jpsi_sig_gen.SetPtEtaPhiE(jpsi_pt_sig_gen, jpsi_eta_sig_gen,jpsi_phi_sig_gen,jpsi_E_sig_gen);
   

    float bc_px_sig_gen = momntm_bc_sig_gen.Px();
    float bc_py_sig_gen = momntm_bc_sig_gen.Py();
    float bc_pz_sig_gen = momntm_bc_sig_gen.Pz();
    
    //calculating the kinematic variables
    float m2_miss_sig_gen = (momntm_bc_sig_gen - momntm_mu1_sig_gen - momntm_mu2_sig_gen - momntm_mu_sig_gen).M2();
    float Q2_sig_gen = (momntm_bc_sig_gen - momntm_mu1_sig_gen - momntm_mu2_sig_gen).M2();
    float Pt_miss_sig_gen = (bc_pt_sig_gen-mu1_pt_sig_gen-mu2_pt_sig_gen-mu_pt_sig_gen);
    float Pt_var_sig_gen = (jpsi_pt_sig_gen - mu_pt_sig_gen );
   
    float del_R_sig_gen = momntm_mu1_sig_gen.DeltaR(momntm_mu2_sig_gen);
    float m2_miss_jpsi_sig_gen = (momntm_jpsi_sig_gen - momntm_mu1_sig_gen - momntm_mu2_sig_gen).M2();

    TVector3 X_bc_sig_gen = -momntm_bc_sig_gen.BoostVector();
    TVector3 X_jpsi_sig_gen = -momntm_jpsi_sig_gen.BoostVector();
    TLorentzVector momntm_mu_bc_frame_sig_gen = momntm_mu_sig_gen;
    momntm_mu_bc_frame_sig_gen.Boost(X_bc_sig_gen);
    TLorentzVector momntm_mu_jpsi_frame_sig_gen = momntm_mu_sig_gen;
    momntm_mu_jpsi_frame_sig_gen.Boost(X_jpsi_sig_gen);
    
    float mu_en_bc_frame_sig_gen = momntm_mu_bc_frame_sig_gen.E();
    float mu_en_jpsi_frame_sig_gen = momntm_mu_jpsi_frame_sig_gen.E();
   

    //neural network
    //TLorentzVector class for calculating four momenta
    float bc_e_sig_nn = sqrt(pow(bcPdgMass,2)+pow(bc_px_sig_nn,2)+pow(bc_py_sig_nn,2)+pow(bc_pz_sig_nn,2));
    TLorentzVector momntm_bc_sig_nn;
    //momntm_bc_sig_nn.SetPtEtaPhiM(bc_pt_sig_nn, bc_eta_sig_nn, bc_phi_sig_nn,bcPdgMass);
    momntm_bc_sig_nn.SetPxPyPzE(bc_px_sig_nn, bc_py_sig_nn, bc_pz_sig_nn, bc_e_sig_nn);
    TLorentzVector momntm_mu1_sig_nn;
    momntm_mu1_sig_nn.SetPtEtaPhiM(mu1_pt_sig_nn, mu1_eta_sig_nn ,mu1_phi_sig_nn,mu_mass);
    TLorentzVector momntm_mu2_sig_nn;
    momntm_mu2_sig_nn.SetPtEtaPhiM(mu2_pt_sig_nn, mu2_eta_sig_nn, mu2_phi_sig_nn,mu_mass);
    TLorentzVector momntm_mu_sig_nn;
    momntm_mu_sig_nn.SetPtEtaPhiM(mu_pt_sig_nn, mu_eta_sig_nn, mu_phi_sig_nn,mu_mass);
    TLorentzVector momntm_jpsi_sig_nn;
    momntm_jpsi_sig_nn.SetPtEtaPhiM(jpsi_pt_sig_nn, jpsi_eta_sig_nn,jpsi_phi_sig_nn,jpsi_mass_sig_nn);

    bc_pt_sig_nn =  momntm_bc_sig_nn.Pt();
    
    //calculating the kinematic variables
    float m2_miss_sig_nn = (momntm_bc_sig_nn - momntm_mu1_sig_nn - momntm_mu2_sig_nn - momntm_mu_sig_nn).M2();
    float Q2_sig_nn = (momntm_bc_sig_nn - momntm_mu1_sig_nn - momntm_mu2_sig_nn).M2();
    float Pt_miss_sig_nn = (bc_pt_sig_nn-mu1_pt_sig_nn-mu2_pt_sig_nn-mu_pt_sig_nn);
    float Pt_var_sig_nn = (mu1_pt_sig_nn +mu2_pt_sig_nn- mu_pt_sig_nn );
    float del_R_sig_nn = momntm_mu1_sig_nn.DeltaR(momntm_mu2_sig_nn);
   
    float bc_pt_ratio_sig_nn = bc_pt_sig_nn/bc_pt_sig_gen;
    float m2_miss_jpsi_sig_nn = (momntm_jpsi_sig_nn - momntm_mu1_sig_nn - momntm_mu2_sig_nn).M2();

    TVector3 X_bc_sig_nn = -momntm_bc_sig_nn.BoostVector();
    TVector3 X_jpsi_sig_nn = -momntm_jpsi_sig_nn.BoostVector();
    TLorentzVector momntm_mu_bc_frame_sig_nn = momntm_mu_sig_nn;
    momntm_mu_bc_frame_sig_nn.Boost(X_bc_sig_nn);
    TLorentzVector momntm_mu_jpsi_frame_sig_nn = momntm_mu_sig_nn;
    momntm_mu_jpsi_frame_sig_nn.Boost(X_jpsi_sig_nn);
    
    float mu_en_bc_frame_sig_nn = momntm_mu_bc_frame_sig_nn.E();
    float mu_en_jpsi_frame_sig_nn = momntm_mu_jpsi_frame_sig_nn.E();

    //Reconstructed
    TLorentzVector momntm_bc_sig_recco;
    momntm_bc_sig_recco.SetPtEtaPhiM(bc_pt_sig_recco, bc_eta_sig_nn, bc_phi_sig_nn,bcPdgMass);
    float m2_miss_sig_recco = (momntm_bc_sig_recco - momntm_mu1_sig_nn - momntm_mu2_sig_nn - momntm_mu_sig_nn).M2();
    float Q2_sig_recco = (momntm_bc_sig_recco - momntm_mu1_sig_nn - momntm_mu2_sig_nn).M2();
    float Pt_miss_sig_recco = (bc_pt_sig_recco-mu1_pt_sig_nn-mu2_pt_sig_nn-mu_pt_sig_nn);

    //jona recipe
    float bc_pt_corr_sig_gen = ((bcPdgMass*bc_pt_sig_recco)/(bc_mass_sig_nn));
    float bc_pt_ratio_sig_gen = bc_pt_corr_sig_gen/bc_pt_sig_gen;

    TLorentzVector momntm_bc_corr_sig_gen;
    momntm_bc_corr_sig_gen.SetPtEtaPhiM(bc_pt_corr_sig_gen, bc_eta_sig_nn, bc_phi_sig_nn, bcPdgMass);
    float Q2_corr_sig_gen = (momntm_bc_corr_sig_gen - momntm_mu1_sig_nn - momntm_mu2_sig_nn).M2();
    float m2_miss_corr_sig_gen = (momntm_bc_corr_sig_gen - momntm_mu1_sig_nn - momntm_mu2_sig_nn - momntm_mu_sig_nn).M2();
    float Pt_miss_corr_sig_gen = (bc_pt_corr_sig_gen - mu1_pt_sig_nn - mu2_pt_sig_nn - mu_pt_sig_nn);
   
    TVector3 X_bc_corr_sig_gen = -momntm_bc_corr_sig_gen.BoostVector();
    TLorentzVector momntm_mu_bc_frame_corr_sig_gen = momntm_mu_sig_nn;
    momntm_mu_bc_frame_corr_sig_gen.Boost(X_bc_corr_sig_gen);
   
    float mu_en_bc_frame_corr_sig_gen = momntm_mu_bc_frame_corr_sig_gen.E();


    float bc_px_corr_sig_gen = momntm_bc_corr_sig_gen.Px();
    float bc_py_corr_sig_gen = momntm_bc_corr_sig_gen.Py();
    float bc_pz_corr_sig_gen = momntm_bc_corr_sig_gen.Pz();
    
    //resolution variables

    float bc_pt_sig_reso_nn = (bc_pt_sig_nn - bc_pt_sig_gen)/bc_pt_sig_gen;
    float bc_px_sig_reso_nn = (bc_px_sig_nn - bc_px_sig_gen)/bc_px_sig_gen;
    float bc_py_sig_reso_nn = (bc_py_sig_nn - bc_py_sig_gen)/bc_py_sig_gen;
    float bc_pz_sig_reso_nn = (bc_pz_sig_nn - bc_pz_sig_gen)/bc_pz_sig_gen;
    float m2_miss_sig_reso_nn = (m2_miss_sig_nn - m2_miss_sig_gen);
    float Q2_sig_reso_nn = (Q2_sig_nn - Q2_sig_gen)/Q2_sig_gen;
    float Pt_miss_sig_reso_nn = (Pt_miss_sig_nn - Pt_miss_sig_gen)/Pt_miss_sig_gen;

    float bc_pt_sig_reso_corr = (bc_pt_corr_sig_gen - bc_pt_sig_gen)/bc_pt_sig_gen;
    float bc_px_sig_reso_corr = (bc_px_corr_sig_gen - bc_px_sig_gen)/bc_px_sig_gen;
    float bc_py_sig_reso_corr = (bc_py_corr_sig_gen - bc_py_sig_gen)/bc_py_sig_gen;
    float bc_pz_sig_reso_corr = (bc_pz_corr_sig_gen - bc_pz_sig_gen)/bc_pz_sig_gen;
    float m2_miss_sig_reso_corr = (m2_miss_corr_sig_gen - m2_miss_sig_gen);
    float Q2_sig_reso_corr = (Q2_corr_sig_gen - Q2_sig_gen)/Q2_sig_gen;
    float Pt_miss_sig_reso_corr = (Pt_miss_corr_sig_gen - Pt_miss_sig_gen)/Pt_miss_sig_gen;

   
    

    //Filling the histograms

    //if ((abs(mu_eta_sig_nn) < 1.3 and mu_pt_sig_nn > 3.3) or (1.3 < abs(mu_eta_sig_nn) < 2.2 and mu_pt_sig_nn > 2.9) or (2.2 < abs(mu_eta_sig_nn) < 2.4 and mu_pt_sig_nn > 2.4)){ //acceptance cuts
    //if ((Trig_dimuon_sig_nn == true) and (Trig_JpsiTk_sig_nn == true)){ //trigger flags

        hmu_m2_miss_sig_gen->Fill(m2_miss_sig_gen);
	hmu_Q2_sig_gen->Fill(Q2_sig_gen);
	hmu_Pt_miss_sig_gen->Fill(Pt_miss_sig_gen);
	hmu_m2_miss_corr_sig_gen->Fill(m2_miss_corr_sig_gen);
	hmu_Q2_corr_sig_gen->Fill(Q2_corr_sig_gen);
	hmu_Pt_miss_corr_sig_gen->Fill(Pt_miss_corr_sig_gen);
	hmu_Pt_var_sig_gen->Fill(Pt_var_sig_gen);
	hmu_mu_pt_sig_gen->Fill(mu_pt_sig_gen);
	hmu_mu_eta_sig_gen->Fill(mu_eta_sig_gen);
	hmu_mu_phi_sig_gen->Fill(mu_phi_sig_gen);
	hmu_del_R_sig_gen->Fill(del_R_sig_gen);
	hmu_mu_en_bc_frame_sig_gen->Fill(mu_en_bc_frame_sig_gen);
	hmu_mu_en_jpsi_frame_sig_gen->Fill(mu_en_jpsi_frame_sig_gen);
	hmu_mu_en_bc_frame_corr_sig_gen->Fill(mu_en_bc_frame_corr_sig_gen);
	hmu_bc_pt_sig_gen->Fill(bc_pt_sig_gen);
	hmu_bc_pt_corr_sig_gen->Fill(bc_pt_corr_sig_gen);
	hmu_bc_pt_ratio_sig_gen->Fill(bc_pt_ratio_sig_gen);
	hprof_sig_gen->Fill(bc_pt_sig_gen,bc_pt_ratio_sig_gen,1);
	hmu_bc_pt_scatter_sig_gen->Fill(bc_pt_sig_gen,bc_pt_corr_sig_gen);
	hmu_m2_miss_jpsi_sig_gen->Fill(m2_miss_jpsi_sig_gen);
	
	
	hmu_m2_miss_sig_nn->Fill(m2_miss_sig_nn);
	hmu_Q2_sig_nn->Fill(Q2_sig_nn);
	hmu_Pt_miss_sig_nn->Fill(Pt_miss_sig_nn);
	hmu_Pt_var_sig_nn->Fill(Pt_var_sig_nn);
	hmu_mu_pt_sig_nn->Fill(mu_pt_sig_nn);
	hmu_mu_eta_sig_nn->Fill(mu_eta_sig_nn);
	hmu_mu_phi_sig_nn->Fill(mu_phi_sig_nn);
	hmu_del_R_sig_nn->Fill(del_R_sig_nn);
	hmu_mu_en_bc_frame_sig_nn->Fill(mu_en_bc_frame_sig_nn);
	hmu_mu_en_jpsi_frame_sig_nn->Fill(mu_en_jpsi_frame_sig_nn);
	hmu_bc_pt_sig_nn->Fill(bc_pt_sig_nn);
	hmu_bc_pt_sig_recco->Fill(bc_pt_sig_recco);
	hmu_bc_pt_ratio_sig_nn->Fill(bc_pt_ratio_sig_nn);
	hprof_sig_nn->Fill(bc_pt_sig_gen,bc_pt_ratio_sig_nn,1);
	hmu_bc_pt_scatter_sig_nn->Fill(bc_pt_sig_gen,bc_pt_sig_nn);
	hmu_m2_miss_jpsi_sig_nn->Fill(m2_miss_jpsi_sig_nn);
	hmu_m2_miss_sig_recco->Fill(m2_miss_sig_recco);
	hmu_Q2_sig_recco->Fill(Q2_sig_recco);
	hmu_Pt_miss_sig_recco->Fill(Pt_miss_sig_recco);


	hmu_bc_pt_sig_reso_nn->Fill(bc_pt_sig_reso_nn);
	hmu_bc_px_sig_reso_nn->Fill(bc_px_sig_reso_nn);
	hmu_bc_py_sig_reso_nn->Fill(bc_py_sig_reso_nn);
	hmu_bc_pz_sig_reso_nn->Fill(bc_pz_sig_reso_nn);
	hmu_m2_miss_sig_reso_nn->Fill(m2_miss_sig_reso_nn);
	hmu_Q2_sig_reso_nn->Fill(Q2_sig_reso_nn);
	hmu_Pt_miss_sig_reso_nn->Fill(Pt_miss_sig_reso_nn);
	hprof_bc_pt_sig_reso_nn->Fill(bc_pt_sig_gen,bc_pt_sig_reso_nn,1);
	hprof_bc_px_sig_reso_nn->Fill(bc_px_sig_gen,bc_px_sig_reso_nn,1);
	hprof_bc_py_sig_reso_nn->Fill(bc_py_sig_gen,bc_py_sig_reso_nn,1);
	hprof_bc_pz_sig_reso_nn->Fill(bc_pz_sig_gen,bc_pz_sig_reso_nn,1);
	hprof_m2_miss_sig_reso_nn->Fill(m2_miss_sig_gen,m2_miss_sig_reso_nn,1);
	hprof_Q2_sig_reso_nn->Fill(Q2_sig_gen,Q2_sig_reso_nn,1);
	hprof_Pt_miss_sig_reso_nn->Fill(Pt_miss_sig_gen,Pt_miss_sig_reso_nn,1);

	hmu_bc_pt_sig_reso_corr->Fill(bc_pt_sig_reso_corr);
	hmu_bc_px_sig_reso_corr->Fill(bc_px_sig_reso_corr);
	hmu_bc_py_sig_reso_corr->Fill(bc_py_sig_reso_corr);
	hmu_bc_pz_sig_reso_corr->Fill(bc_pz_sig_reso_corr);
	hmu_m2_miss_sig_reso_corr->Fill(m2_miss_sig_reso_corr);
	hmu_Q2_sig_reso_corr->Fill(Q2_sig_reso_corr);
	hmu_Pt_miss_sig_reso_corr->Fill(Pt_miss_sig_reso_corr);
	hprof_bc_pt_sig_reso_corr->Fill(bc_pt_sig_gen,bc_pt_sig_reso_corr,1);
	hprof_bc_px_sig_reso_corr->Fill(bc_px_sig_gen,bc_px_sig_reso_corr,1);
	hprof_bc_py_sig_reso_corr->Fill(bc_py_sig_gen,bc_py_sig_reso_corr,1);
	hprof_bc_pz_sig_reso_corr->Fill(bc_pz_sig_gen,bc_pz_sig_reso_corr,1);
	hprof_m2_miss_sig_reso_corr->Fill(m2_miss_sig_gen,m2_miss_sig_reso_corr,1);
	hprof_Q2_sig_reso_corr->Fill(Q2_sig_gen,Q2_sig_reso_corr,1);
	hprof_Pt_miss_sig_reso_corr->Fill(Pt_miss_sig_gen,Pt_miss_sig_reso_corr,1);

	hmu_m2_miss_scatter_sig_nn->Fill(m2_miss_sig_gen,m2_miss_sig_nn);
	hmu_Q2_scatter_sig_nn->Fill(Q2_sig_gen,Q2_sig_nn);
	hmu_Pt_miss_scatter_sig_nn->Fill(Pt_miss_sig_gen,Pt_miss_sig_nn);

	//}
	//}
  }  
  
  
   //Drawing Histograms for gen level

  /**
  //Missing mass Squared
  hmu_m2_miss_norm_gen->SetLineColor(kRed);
  hmu_m2_miss_norm_gen->GetXaxis()->SetTitle("m^{2}_{miss}[GeV^{2}]");
  hmu_m2_miss_norm_gen->GetYaxis()->SetTitle("a.u.");
  hmu_m2_miss_norm_gen->DrawNormalized();
  hmu_m2_miss_sig_gen->SetLineColor(kBlue);
  hmu_m2_miss_sig_gen->DrawNormalized("same");
  legend->AddEntry(hmu_m2_miss_norm_gen, "#mu channel", "l");
  legend->AddEntry(hmu_m2_miss_sig_gen, "#tau channel", "l");
  legend->Draw();
  c1->SaveAs("M2_miss_gen.png");

  
  //Squared four momentum transfered to lepton system
  hmu_Q2_norm_gen->SetLineColor(kRed);
  hmu_Q2_norm_gen->GetXaxis()->SetTitle("Q^{2}[GeV^{2}]");
  hmu_Q2_norm_gen->GetYaxis()->SetTitle("a.u.");
  hmu_Q2_norm_gen->DrawNormalized();
  hmu_Q2_sig_gen->SetLineColor(kBlue);
  hmu_Q2_sig_gen->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("Q2_gen.png");
  


  //Missing transverse momentum
  hmu_Pt_miss_norm_gen->SetLineColor(kRed);
  hmu_Pt_miss_norm_gen->GetXaxis()->SetTitle("P_{T}^{miss}[GeV]");
  hmu_Pt_miss_norm_gen->GetYaxis()->SetTitle("a.u.");
  hmu_Pt_miss_norm_gen->DrawNormalized();
  hmu_Pt_miss_sig_gen->SetLineColor(kBlue);
  hmu_Pt_miss_sig_gen->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("Pt_miss_gen.png");


  //Missing mass Squared corrected
  hmu_m2_miss_corr_norm_gen->SetLineColor(kRed);
  hmu_m2_miss_corr_norm_gen->GetXaxis()->SetTitle("m^{2}_{miss}[GeV^{2}]");
  hmu_m2_miss_corr_norm_gen->GetYaxis()->SetTitle("a.u.");
  hmu_m2_miss_corr_norm_gen->DrawNormalized();
  hmu_m2_miss_corr_sig_gen->SetLineColor(kBlue);
  hmu_m2_miss_corr_sig_gen->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("M2_miss_corr.png");

  
  //Squared four momentum transfered to lepton system corrected
  hmu_Q2_corr_norm_gen->SetLineColor(kRed);
  hmu_Q2_corr_norm_gen->GetXaxis()->SetTitle("Q^{2}[GeV^{2}]");
  hmu_Q2_corr_norm_gen->GetYaxis()->SetTitle("a.u.");
  hmu_Q2_corr_norm_gen->DrawNormalized();
  hmu_Q2_corr_sig_gen->SetLineColor(kBlue);
  hmu_Q2_corr_sig_gen->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("Q2_corr.png");
  


  //Missing transverse momentum corrected
  hmu_Pt_miss_corr_norm_gen->SetLineColor(kRed);
  hmu_Pt_miss_corr_norm_gen->GetXaxis()->SetTitle("P_{T}^{miss}[GeV]");
  hmu_Pt_miss_corr_norm_gen->GetYaxis()->SetTitle("a.u.");
  hmu_Pt_miss_corr_norm_gen->DrawNormalized();
  hmu_Pt_miss_corr_sig_gen->SetLineColor(kBlue);
  hmu_Pt_miss_corr_sig_gen->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("Pt_miss_corr.png");

 

  //Tranverse variable momentum
  hmu_Pt_var_norm_gen->SetLineColor(kRed);
  hmu_Pt_var_norm_gen->GetXaxis()->SetTitle("P_{T}^{var}[GeV]");
  hmu_Pt_var_norm_gen->GetYaxis()->SetTitle("a.u.");
  hmu_Pt_var_norm_gen->DrawNormalized();
  hmu_Pt_var_sig_gen->SetLineColor(kBlue);
  hmu_Pt_var_sig_gen->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("Pt_var_gen.png");
  
  //Delta R
  hmu_del_R_norm_gen->SetLineColor(kRed);
  hmu_del_R_norm_gen->GetXaxis()->SetTitle("#Delta R(#mu_{1},#mu_{2})");
  hmu_del_R_norm_gen->GetYaxis()->SetTitle("a.u.");
  hmu_del_R_norm_gen->DrawNormalized();
  hmu_del_R_sig_gen->SetLineColor(kBlue);
  hmu_del_R_sig_gen->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("del_R_gen.png");

  //Missing mass Squared for Jpsi
  hmu_m2_miss_jpsi_norm_gen->SetLineColor(kRed);
  hmu_m2_miss_jpsi_norm_gen->GetXaxis()->SetTitle("m^{2}_{miss}[GeV^{2}]");
  hmu_m2_miss_jpsi_norm_gen->GetYaxis()->SetTitle("a.u.");
  hmu_m2_miss_jpsi_norm_gen->DrawNormalized();
  hmu_m2_miss_jpsi_sig_gen->SetLineColor(kBlue);
  hmu_m2_miss_jpsi_sig_gen->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("M2_miss_jpsi_gen.png");



  //Unpaired muon tranverse momentum
  hmu_mu_pt_norm_gen->SetLineColor(kRed);
  hmu_mu_pt_norm_gen->GetXaxis()->SetTitle("P_{T}^{#mu}[GeV]");
  hmu_mu_pt_norm_gen->GetYaxis()->SetTitle("a.u.");
  hmu_mu_pt_norm_gen->DrawNormalized();
  hmu_mu_pt_sig_gen->SetLineColor(kBlue);
  hmu_mu_pt_sig_gen->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("mu_pt_gen.png");


  //Unpaired muon psuedorapidity
  hmu_mu_eta_norm_gen->SetLineColor(kRed);
  hmu_mu_eta_norm_gen->GetXaxis()->SetTitle("#eta^{#mu}");
  hmu_mu_eta_norm_gen->GetYaxis()->SetTitle("a.u.");
  hmu_mu_eta_norm_gen->DrawNormalized();
  hmu_mu_eta_sig_gen->SetLineColor(kBlue);
  hmu_mu_eta_sig_gen->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("mu_eta_gen.png");

  
  //Unpaired muon phi
  hmu_mu_phi_norm_gen->SetLineColor(kRed);
  hmu_mu_phi_norm_gen->GetXaxis()->SetTitle("#phi^{#mu}");
  hmu_mu_phi_norm_gen->GetYaxis()->SetTitle("a.u.");
  hmu_mu_phi_norm_gen->DrawNormalized();
  hmu_mu_phi_sig_gen->SetLineColor(kBlue);
  hmu_mu_phi_sig_gen->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("mu_phi_gen.png");

  
  //Energy of unpaired muon in Bc rest frame
  hmu_mu_en_bc_frame_norm_gen->SetLineColor(kRed);
  hmu_mu_en_bc_frame_norm_gen->GetXaxis()->SetTitle("E(#mu) in rest frame of B_{c} (GeV)");
  hmu_mu_en_bc_frame_norm_gen->GetYaxis()->SetTitle("a.u.");
  hmu_mu_en_bc_frame_norm_gen->DrawNormalized();
  hmu_mu_en_bc_frame_sig_gen->SetLineColor(kBlue);
  hmu_mu_en_bc_frame_sig_gen->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("E_bc_gen.png");

  //Energy of unpaired muon in Jpsi rest frame
  hmu_mu_en_jpsi_frame_norm_gen->SetLineColor(kRed);
  hmu_mu_en_jpsi_frame_norm_gen->GetXaxis()->SetTitle("E(#mu) in rest frame of J/(#Psi) (GeV)");
  hmu_mu_en_jpsi_frame_norm_gen->GetYaxis()->SetTitle("a.u.");
  hmu_mu_en_jpsi_frame_norm_gen->DrawNormalized();
  hmu_mu_en_jpsi_frame_sig_gen->SetLineColor(kBlue);
  hmu_mu_en_jpsi_frame_sig_gen->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("E_jpsi_gen.png");

  //Energy of unpaired muon in Bc rest frame corrected
  hmu_mu_en_bc_frame_corr_norm_gen->SetLineColor(kRed);
  hmu_mu_en_bc_frame_corr_norm_gen->GetXaxis()->SetTitle("E(#mu) in rest frame of B_{c} (GeV)");
  hmu_mu_en_bc_frame_corr_norm_gen->GetYaxis()->SetTitle("a.u.");
  hmu_mu_en_bc_frame_corr_norm_gen->DrawNormalized();
  hmu_mu_en_bc_frame_corr_sig_gen->SetLineColor(kBlue);
  hmu_mu_en_bc_frame_corr_sig_gen->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("E_bc_corr_gen.png");

  
  //Transverse momentum of Bc
  hmu_bc_pt_norm_gen->SetLineColor(kRed);
  hmu_bc_pt_norm_gen->GetXaxis()->SetTitle("P_{T}^{Bc} (GeV)");
  hmu_bc_pt_norm_gen->GetYaxis()->SetTitle("a.u.");
  hmu_bc_pt_norm_gen->DrawNormalized();
  hmu_bc_pt_sig_gen->SetLineColor(kBlue);
  hmu_bc_pt_sig_gen->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("Pt_Bc_gen.png");

  
  //Corrected transverse momentum of Bc
  hmu_bc_pt_corr_norm_gen->SetLineColor(kRed);
  hmu_bc_pt_corr_norm_gen->GetXaxis()->SetTitle("Corrected P_{T}^{Bc} (GeV)");
  hmu_bc_pt_corr_norm_gen->GetYaxis()->SetTitle("a.u.");
  hmu_bc_pt_corr_norm_gen->DrawNormalized();
  hmu_bc_pt_corr_sig_gen->SetLineColor(kBlue);
  hmu_bc_pt_corr_sig_gen->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("Corr_pt_bc.png");

  //Ratio of true and corrected transverse momentum of Bc
  hmu_bc_pt_ratio_norm_gen->SetLineColor(kRed);
  hmu_bc_pt_ratio_norm_gen->GetXaxis()->SetTitle("Ratio of corrected and true P_{T}^{Bc}");
  hmu_bc_pt_ratio_norm_gen->GetYaxis()->SetTitle("a.u.");
  hmu_bc_pt_ratio_norm_gen->DrawNormalized();
  hmu_bc_pt_ratio_sig_gen->SetLineColor(kBlue);
  hmu_bc_pt_ratio_sig_gen->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("Ratio_pt_bc_gen.png");
  

  //Scatter plot of Transverse momentum of Bc
  hmu_bc_pt_scatter_norm_gen->SetMarkerColor(kRed);
  hmu_bc_pt_scatter_norm_gen->GetXaxis()->SetTitle("True P_{T}^{Bc} (GeV)");
  hmu_bc_pt_scatter_norm_gen->GetYaxis()->SetTitle("Corrected P_{T}^{Bc} (GeV)");
  hmu_bc_pt_scatter_norm_gen->DrawNormalized();
  hmu_bc_pt_scatter_sig_gen->SetMarkerColor(kBlue);
  hmu_bc_pt_scatter_sig_gen->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("scatter_Pt_Bc_gen.png");


  //Profile histogram of transverse momentum of Bc
  hprof_norm_gen->SetLineColor(kRed);
  hprof_norm_gen->GetXaxis()->SetTitle("True P_{T}^{Bc} (GeV)");
  hprof_norm_gen->GetYaxis()->SetTitle("Ratio of corrected and true P_{T}^{Bc}");
  hprof_norm_gen->Draw();
  hprof_sig_gen->SetLineColor(kBlue);
  hprof_sig_gen->Draw("same");
  legend->Draw();
  c1->SaveAs("Profile_pt_bc_gen.png");


  //Drawing Histograms for nn output
  
  //Missing mass Squared
  hmu_m2_miss_norm_nn->SetLineColor(kRed);
  hmu_m2_miss_norm_nn->GetXaxis()->SetTitle("m^{2}_{miss}[GeV^{2}]");
  hmu_m2_miss_norm_nn->GetYaxis()->SetTitle("a.u.");
  hmu_m2_miss_norm_nn->DrawNormalized();
  hmu_m2_miss_sig_nn->SetLineColor(kBlue);
  hmu_m2_miss_sig_nn->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("M2_miss_nn.png");

  
  //Squared four momentum transfered to lepton system
  hmu_Q2_norm_nn->SetLineColor(kRed);
  hmu_Q2_norm_nn->GetXaxis()->SetTitle("Q^{2}[GeV^{2}]");
  hmu_Q2_norm_nn->GetYaxis()->SetTitle("a.u.");
  hmu_Q2_norm_nn->DrawNormalized();
  hmu_Q2_sig_nn->SetLineColor(kBlue);
  hmu_Q2_sig_nn->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("Q2_nn.png");
  

  
  //Missing transverse momentum
  hmu_Pt_miss_norm_nn->SetLineColor(kRed);
  hmu_Pt_miss_norm_nn->GetXaxis()->SetTitle("P_{T}^{miss}[GeV]");
  hmu_Pt_miss_norm_nn->GetYaxis()->SetTitle("a.u.");
  hmu_Pt_miss_norm_nn->DrawNormalized();
  hmu_Pt_miss_sig_nn->SetLineColor(kBlue);
  hmu_Pt_miss_sig_nn->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("Pt_miss_nn.png");
  
 
  //Transverse Momenta of B_{c} resolution
  hmu_bc_pt_norm_reso_nn->SetLineColor(kRed);
  hmu_bc_pt_norm_reso_nn->GetXaxis()->SetTitle("Resolution");
  hmu_bc_pt_norm_reso_nn->GetYaxis()->SetTitle("a.u.");
  hmu_bc_pt_norm_reso_nn->DrawNormalized();
  hmu_bc_pt_norm_reso_corr->SetLineColor(kBlue);
  hmu_bc_pt_norm_reso_corr->DrawNormalized("same");
  legend->AddEntry(hmu_bc_pt_norm_reso_nn, "nn_predicted", "l");
  legend->AddEntry(hmu_bc_pt_norm_reso_corr, "collinear approximation", "l");
  legend->Draw();
  c1->SaveAs("Bc_pt_norm_reso.png");

  //Transverse Momenta of B_{c} resolution
  hmu_bc_pt_sig_reso_nn->SetLineColor(kRed);
  hmu_bc_pt_sig_reso_nn->GetXaxis()->SetTitle("Resolution");
  hmu_bc_pt_sig_reso_nn->GetYaxis()->SetTitle("a.u.");
  hmu_bc_pt_sig_reso_nn->DrawNormalized();
  hmu_bc_pt_sig_reso_corr->SetLineColor(kBlue);
  hmu_bc_pt_sig_reso_corr->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("Bc_pt_sig_reso.png");

  //Transverse Momenta of B_{c} resolution
  hmu_bc_px_norm_reso_nn->SetLineColor(kRed);
  hmu_bc_px_norm_reso_nn->GetXaxis()->SetTitle("Resolution");
  hmu_bc_px_norm_reso_nn->GetYaxis()->SetTitle("a.u.");
  hmu_bc_px_norm_reso_nn->DrawNormalized();
  hmu_bc_px_norm_reso_corr->SetLineColor(kBlue);
  hmu_bc_px_norm_reso_corr->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("Bc_px_norm_reso.png");

  //Transverse Momenta of B_{c} resolution
  hmu_bc_px_sig_reso_nn->SetLineColor(kRed);
  hmu_bc_px_sig_reso_nn->GetXaxis()->SetTitle("Resolution");
  hmu_bc_px_sig_reso_nn->GetYaxis()->SetTitle("a.u.");
  hmu_bc_px_sig_reso_nn->DrawNormalized();
  hmu_bc_px_sig_reso_corr->SetLineColor(kBlue);
  hmu_bc_px_sig_reso_corr->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("Bc_px_sig_reso.png");

  //Transverse Momenta of B_{c} resolution
  hmu_bc_py_norm_reso_nn->SetLineColor(kRed);
  hmu_bc_py_norm_reso_nn->GetXaxis()->SetTitle("Resolution");
  hmu_bc_py_norm_reso_nn->GetYaxis()->SetTitle("a.u.");
  hmu_bc_py_norm_reso_nn->DrawNormalized();
  hmu_bc_py_norm_reso_corr->SetLineColor(kBlue);
  hmu_bc_py_norm_reso_corr->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("Bc_py_norm_reso.png");

  //Transverse Momenta of B_{c} resolution
  hmu_bc_py_sig_reso_nn->SetLineColor(kRed);
  hmu_bc_py_sig_reso_nn->GetXaxis()->SetTitle("Resolution");
  hmu_bc_py_sig_reso_nn->GetYaxis()->SetTitle("a.u.");
  hmu_bc_py_sig_reso_nn->DrawNormalized();
  hmu_bc_py_sig_reso_corr->SetLineColor(kBlue);
  hmu_bc_py_sig_reso_corr->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("Bc_py_sig_reso.png");


  //Transverse Momenta of B_{c} resolution
  hmu_bc_pz_norm_reso_nn->SetLineColor(kRed);
  hmu_bc_pz_norm_reso_nn->GetXaxis()->SetTitle("Resolution");
  hmu_bc_pz_norm_reso_nn->GetYaxis()->SetTitle("a.u.");
  hmu_bc_pz_norm_reso_nn->DrawNormalized();
  hmu_bc_pz_norm_reso_corr->SetLineColor(kBlue);
  hmu_bc_pz_norm_reso_corr->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("Bc_pz_norm_reso.png");

  //Transverse Momenta of B_{c} resolution
  hmu_bc_pz_sig_reso_nn->SetLineColor(kRed);
  hmu_bc_pz_sig_reso_nn->GetXaxis()->SetTitle("Resolution");
  hmu_bc_pz_sig_reso_nn->GetYaxis()->SetTitle("a.u.");
  hmu_bc_pz_sig_reso_nn->DrawNormalized();
  hmu_bc_pz_sig_reso_corr->SetLineColor(kBlue);
  hmu_bc_pz_sig_reso_corr->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("Bc_pz_sig_reso.png");

 
  //Missing mass Squared
  hmu_m2_miss_norm_reso_nn->SetLineColor(kRed);
  hmu_m2_miss_norm_reso_nn->GetXaxis()->SetTitle("Resolution");
  hmu_m2_miss_norm_reso_nn->GetYaxis()->SetTitle("a.u.");
  hmu_m2_miss_norm_reso_nn->DrawNormalized();
  hmu_m2_miss_norm_reso_corr->SetLineColor(kBlue);
  hmu_m2_miss_norm_reso_corr->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("M2_miss_norm_reso.png");

  //Missing mass Squared
  hmu_m2_miss_sig_reso_nn->SetLineColor(kRed);
  hmu_m2_miss_sig_reso_nn->GetXaxis()->SetTitle("Resolution");
  hmu_m2_miss_sig_reso_nn->GetYaxis()->SetTitle("a.u.");
  hmu_m2_miss_sig_reso_nn->DrawNormalized();
  hmu_m2_miss_sig_reso_corr->SetLineColor(kBlue);
  hmu_m2_miss_sig_reso_corr->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("M2_miss_sig_reso.png");

  //Missing mass Squared
  hmu_Q2_norm_reso_nn->SetLineColor(kRed);
  hmu_Q2_norm_reso_nn->GetXaxis()->SetTitle("Resolution");
  hmu_Q2_norm_reso_nn->GetYaxis()->SetTitle("a.u.");
  hmu_Q2_norm_reso_nn->DrawNormalized();
  hmu_Q2_norm_reso_corr->SetLineColor(kBlue);
  hmu_Q2_norm_reso_corr->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("Q2_norm_reso.png");

  //Missing mass Squared
  hmu_Q2_sig_reso_nn->SetLineColor(kRed);
  hmu_Q2_sig_reso_nn->GetXaxis()->SetTitle("Resolution");
  hmu_Q2_sig_reso_nn->GetYaxis()->SetTitle("a.u.");
  hmu_Q2_sig_reso_nn->DrawNormalized();
  hmu_Q2_sig_reso_corr->SetLineColor(kBlue);
  hmu_Q2_sig_reso_corr->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("Q2_sig_reso.png");

  //Missing mass Squared
  hmu_w_ratio_Q2_reso_nn->SetLineColor(kRed);
  hmu_w_ratio_Q2_reso_nn->GetXaxis()->SetTitle("w_ratio_Q2_resolution");
  hmu_w_ratio_Q2_reso_nn->GetYaxis()->SetTitle("a.u.");
  hmu_w_ratio_Q2_reso_nn->DrawNormalized();
  hmu_w_ratio_Q2_reso_corr->SetLineColor(kBlue);
  hmu_w_ratio_Q2_reso_corr->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("w_ratio_Q2_reso.png");

  //Missing mass Squared
  hmu_w_ratio_Q2_scatter_nn->SetMarkerColor(kRed);
  hmu_w_ratio_Q2_scatter_nn->GetXaxis()->SetTitle("w_ratio_true");
  hmu_w_ratio_Q2_scatter_nn->GetYaxis()->SetTitle("w_ratio_nn/corr");
  hmu_w_ratio_Q2_scatter_nn->Draw();
  hmu_w_ratio_Q2_scatter_corr->SetMarkerColor(kBlue);
  hmu_w_ratio_Q2_scatter_corr->Draw("same");
  legend->Draw();
  c1->SaveAs("w_ratio_Q2_scatter.png");  


  //Missing mass Squared
  hmu_Pt_miss_norm_reso_nn->SetLineColor(kRed);
  hmu_Pt_miss_norm_reso_nn->GetXaxis()->SetTitle("Resolution");
  hmu_Pt_miss_norm_reso_nn->GetYaxis()->SetTitle("a.u.");
  hmu_Pt_miss_norm_reso_nn->DrawNormalized();
  hmu_Pt_miss_norm_reso_corr->SetLineColor(kBlue);
  hmu_Pt_miss_norm_reso_corr->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("Pt_miss_norm_reso.png");

  //Missing mass Squared
  hmu_Pt_miss_sig_reso_nn->SetLineColor(kRed);
  hmu_Pt_miss_sig_reso_nn->GetXaxis()->SetTitle("Resolution");
  hmu_Pt_miss_sig_reso_nn->GetYaxis()->SetTitle("a.u.");
  hmu_Pt_miss_sig_reso_nn->DrawNormalized();
  hmu_Pt_miss_sig_reso_corr->SetLineColor(kBlue);
  hmu_Pt_miss_sig_reso_corr->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("Pt_miss_sig_reso.png");

  
  //Missing mass Squared
  hprof_bc_pt_norm_reso_nn->SetLineColor(kRed);
  hprof_bc_pt_norm_reso_nn->GetXaxis()->SetTitle("True bc_pt");
  hprof_bc_pt_norm_reso_nn->GetYaxis()->SetTitle("Resolution");
  hprof_bc_pt_norm_reso_nn->Draw();
  hprof_bc_pt_norm_reso_corr->SetLineColor(kBlue);
  hprof_bc_pt_norm_reso_corr->Draw("same");
  legend->Draw();
  c1->SaveAs("Profile_bc_pt_norm_reso.png");

  //Missing mass Squared
  hprof_bc_pt_sig_reso_nn->SetLineColor(kRed);
  hprof_bc_pt_sig_reso_nn->GetXaxis()->SetTitle("True bc_pt");
  hprof_bc_pt_sig_reso_nn->GetYaxis()->SetTitle("Resolution");
  hprof_bc_pt_sig_reso_nn->Draw();
  hprof_bc_pt_sig_reso_corr->SetLineColor(kBlue);
  hprof_bc_pt_sig_reso_corr->Draw("same");
  legend->Draw();
  c1->SaveAs("Profile_bc_pt_sig_reso.png");


  //Missing mass Squared
  hprof_bc_px_norm_reso_nn->SetLineColor(kRed);
  hprof_bc_px_norm_reso_nn->GetXaxis()->SetTitle("True bc_px");
  hprof_bc_px_norm_reso_nn->GetYaxis()->SetTitle("Resolution");
  hprof_bc_px_norm_reso_nn->Draw();
  hprof_bc_px_norm_reso_corr->SetLineColor(kBlue);
  hprof_bc_px_norm_reso_corr->Draw("same");
  legend->Draw();
  c1->SaveAs("Profile_bc_px_norm_reso.png");

  //Missing mass Squared
  hprof_bc_px_sig_reso_nn->SetLineColor(kRed);
  hprof_bc_px_sig_reso_nn->GetXaxis()->SetTitle("True bc_px");
  hprof_bc_px_sig_reso_nn->GetYaxis()->SetTitle("Resolution");
  hprof_bc_px_sig_reso_nn->Draw();
  hprof_bc_px_sig_reso_corr->SetLineColor(kBlue);
  hprof_bc_px_sig_reso_corr->Draw("same");
  legend->Draw();
  c1->SaveAs("Profile_bc_px_sig_reso.png");

  //Missing mass Squared
  hprof_bc_py_norm_reso_nn->SetLineColor(kRed);
  hprof_bc_py_norm_reso_nn->GetXaxis()->SetTitle("True bc_py");
  hprof_bc_py_norm_reso_nn->GetYaxis()->SetTitle("Resolution");
  hprof_bc_py_norm_reso_nn->Draw();
  hprof_bc_py_norm_reso_corr->SetLineColor(kBlue);
  hprof_bc_py_norm_reso_corr->Draw("same");
  legend->Draw();
  c1->SaveAs("Profile_bc_py_norm_reso.png");

  //Missing mass Squared
  hprof_bc_py_sig_reso_nn->SetLineColor(kRed);
  hprof_bc_py_sig_reso_nn->GetXaxis()->SetTitle("True bc_py");
  hprof_bc_py_sig_reso_nn->GetYaxis()->SetTitle("Resolution");
  hprof_bc_py_sig_reso_nn->Draw();
  hprof_bc_py_sig_reso_corr->SetLineColor(kBlue);
  hprof_bc_py_sig_reso_corr->Draw("same");
  legend->Draw();
  c1->SaveAs("Profile_bc_py_sig_reso.png");


  //Missing mass Squared
  hprof_bc_pz_norm_reso_nn->SetLineColor(kRed);
  hprof_bc_pz_norm_reso_nn->GetXaxis()->SetTitle("True bc_pz");
  hprof_bc_pz_norm_reso_nn->GetYaxis()->SetTitle("Resolution");
  hprof_bc_pz_norm_reso_nn->Draw();
  hprof_bc_pz_norm_reso_corr->SetLineColor(kBlue);
  hprof_bc_pz_norm_reso_corr->Draw("same");
  legend->Draw();
  c1->SaveAs("Profile_bc_pz_norm_reso.png");

  //Missing mass Squared
  hprof_bc_pz_sig_reso_nn->SetLineColor(kRed);
  hprof_bc_pz_sig_reso_nn->GetXaxis()->SetTitle("True bc_pz");
  hprof_bc_pz_sig_reso_nn->GetYaxis()->SetTitle("Resolution");
  hprof_bc_pz_sig_reso_nn->Draw();
  hprof_bc_pz_sig_reso_corr->SetLineColor(kBlue);
  hprof_bc_pz_sig_reso_corr->Draw("same");
  legend->Draw();
  c1->SaveAs("Profile_bc_pz_sig_reso.png");
  
 
  //Missing mass Squared
  hprof_m2_miss_norm_reso_nn->SetLineColor(kRed);
  hprof_m2_miss_norm_reso_nn->GetXaxis()->SetTitle("True m2_miss");
  hprof_m2_miss_norm_reso_nn->GetYaxis()->SetTitle("Resolution");
  hprof_m2_miss_norm_reso_nn->Draw();
  hprof_m2_miss_norm_reso_corr->SetLineColor(kBlue);
  hprof_m2_miss_norm_reso_corr->Draw("same");
  legend->Draw();
  c1->SaveAs("Profile_m2_miss_norm_reso.png");

  //Missing mass Squared
  hprof_m2_miss_sig_reso_nn->SetLineColor(kRed);
  hprof_m2_miss_sig_reso_nn->GetXaxis()->SetTitle("True m2_miss");
  hprof_m2_miss_sig_reso_nn->GetYaxis()->SetTitle("Resolution");
  hprof_m2_miss_sig_reso_nn->Draw();
  hprof_m2_miss_sig_reso_corr->SetLineColor(kBlue);
  hprof_m2_miss_sig_reso_corr->Draw("same");
  legend->Draw();
  c1->SaveAs("Profile_m2_miss_sig_reso.png");
  
  //Missing mass Squared
  hprof_Pt_miss_norm_reso_nn->SetLineColor(kRed);
  hprof_Pt_miss_norm_reso_nn->GetXaxis()->SetTitle("True Pt_miss");
  hprof_Pt_miss_norm_reso_nn->GetYaxis()->SetTitle("Resolution");
  hprof_Pt_miss_norm_reso_nn->Draw();
  hprof_Pt_miss_norm_reso_corr->SetLineColor(kBlue);
  hprof_Pt_miss_norm_reso_corr->Draw("same");
  legend->Draw();
  c1->SaveAs("Profile_Pt_miss_norm_reso.png");

  //Missing mass Squared
  hprof_Pt_miss_sig_reso_nn->SetLineColor(kRed);
  hprof_Pt_miss_sig_reso_nn->GetXaxis()->SetTitle("True Pt_miss");
  hprof_Pt_miss_sig_reso_nn->GetYaxis()->SetTitle("Resolution");
  hprof_Pt_miss_sig_reso_nn->Draw();
  hprof_Pt_miss_sig_reso_corr->SetLineColor(kBlue);
  hprof_Pt_miss_sig_reso_corr->Draw("same");
  legend->Draw();
  c1->SaveAs("Profile_Pt_miss_sig_reso.png");

  //Missing mass Squared
  hprof_Q2_norm_reso_nn->SetLineColor(kRed);
  hprof_Q2_norm_reso_nn->GetXaxis()->SetTitle("True Q2");
  hprof_Q2_norm_reso_nn->GetYaxis()->SetTitle("Resolution");
  hprof_Q2_norm_reso_nn->Draw();
  hprof_Q2_norm_reso_corr->SetLineColor(kBlue);
  hprof_Q2_norm_reso_corr->Draw("same");
  legend->Draw();
  //gPad->SetLogy();
  //c1->SaveAs("Profile_Q2_norm_reso_logy.png");
  c1->SaveAs("Profile_Q2_norm_reso.png");

  //Missing mass Squared
  hprof_Q2_sig_reso_nn->SetLineColor(kRed);
  hprof_Q2_sig_reso_nn->GetXaxis()->SetTitle("True Q2");
  hprof_Q2_sig_reso_nn->GetYaxis()->SetTitle("Resolution");
  hprof_Q2_sig_reso_nn->Draw();
  hprof_Q2_sig_reso_corr->SetLineColor(kBlue);
  hprof_Q2_sig_reso_corr->Draw("same");
  legend->Draw();
  //c1->SaveAs("Profile_Q2_sig_reso_logy.png");
  c1->SaveAs("Profile_Q2_sig_reso.png");

  
  //Scatter plot of Transverse momentum of Bc
  hmu_m2_miss_scatter_norm_nn->SetMarkerColor(kRed);
  hmu_m2_miss_scatter_norm_nn->GetXaxis()->SetTitle("True m_{miss}^{2} [GeV^{2}]");
  hmu_m2_miss_scatter_norm_nn->GetYaxis()->SetTitle("Predicted_nn m_{miss}^{2} [GeV^{2}]");
  hmu_m2_miss_scatter_norm_nn->Draw();
  hmu_m2_miss_scatter_sig_nn->SetMarkerColor(kBlue);
  hmu_m2_miss_scatter_sig_nn->Draw("same");
  legend->Draw();
  c1->SaveAs("scatter_m2_miss_nn.png");

  //Scatter plot of Transverse momentum of Bc
  hmu_Q2_scatter_norm_nn->SetMarkerColor(kRed);
  hmu_Q2_scatter_norm_nn->GetXaxis()->SetTitle("True Q^{2} [GeV^{2}]");
  hmu_Q2_scatter_norm_nn->GetYaxis()->SetTitle("Predicted_nn Q^{2} [GeV^{2}]");
  hmu_Q2_scatter_norm_nn->Draw();
  hmu_Q2_scatter_sig_nn->SetMarkerColor(kBlue);
  hmu_Q2_scatter_sig_nn->Draw("same");
  legend->Draw();
  c1->SaveAs("scatter_Q2_nn.png");

  //Scatter plot of Transverse momentum of Bc
  hmu_Pt_miss_scatter_norm_nn->SetMarkerColor(kRed);
  hmu_Pt_miss_scatter_norm_nn->GetXaxis()->SetTitle("True P_{T}^{miss} [GeV]");
  hmu_Pt_miss_scatter_norm_nn->GetYaxis()->SetTitle("Predicted_nn P_{T}^{miss} [GeV]");
  hmu_Pt_miss_scatter_norm_nn->Draw();
  hmu_Pt_miss_scatter_sig_nn->SetMarkerColor(kBlue);
  hmu_Pt_miss_scatter_sig_nn->Draw("same");
  legend->Draw();
  c1->SaveAs("scatter_Pt_miss_nn.png");
  **/
  
 
  /**
  //Missing mass Squared
  hmu_m2_miss_norm_recco->SetLineColor(kRed);
  hmu_m2_miss_norm_recco->GetXaxis()->SetTitle("m^{2}_{miss}[GeV^{2}]");
  hmu_m2_miss_norm_recco->GetYaxis()->SetTitle("a.u.");
  hmu_m2_miss_norm_recco->DrawNormalized();
  hmu_m2_miss_sig_recco->SetLineColor(kBlue);
  hmu_m2_miss_sig_recco->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("M2_miss_recco.png");

  
  //Squared four momentum transfered to lepton system
  hmu_Q2_norm_recco->SetLineColor(kRed);
  hmu_Q2_norm_recco->GetXaxis()->SetTitle("Q^{2}[GeV^{2}]");
  hmu_Q2_norm_recco->GetYaxis()->SetTitle("a.u.");
  hmu_Q2_norm_recco->DrawNormalized();
  hmu_Q2_sig_recco->SetLineColor(kBlue);
  hmu_Q2_sig_recco->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("Q2_recco.png");
  

  
  //Missing transverse momentum
  hmu_Pt_miss_norm_recco->SetLineColor(kRed);
  hmu_Pt_miss_norm_recco->GetXaxis()->SetTitle("P_{T}^{miss}[GeV]");
  hmu_Pt_miss_norm_recco->GetYaxis()->SetTitle("a.u.");
  hmu_Pt_miss_norm_recco->DrawNormalized();
  hmu_Pt_miss_sig_recco->SetLineColor(kBlue);
  hmu_Pt_miss_sig_recco->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("Pt_miss_recco.png");
  

  //Tranverse variable momentum
  hmu_Pt_var_norm_nn->SetLineColor(kRed);
  hmu_Pt_var_norm_nn->GetXaxis()->SetTitle("P_{T}^{var}[GeV]");
  hmu_Pt_var_norm_nn->GetYaxis()->SetTitle("a.u.");
  hmu_Pt_var_norm_nn->DrawNormalized();
  hmu_Pt_var_sig_nn->SetLineColor(kBlue);
  hmu_Pt_var_sig_nn->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("Pt_var_nn.png");
  
  //Delta R
  hmu_del_R_norm_nn->SetLineColor(kRed);
  hmu_del_R_norm_nn->GetXaxis()->SetTitle("#Delta R(#mu_{1},#mu_{2})");
  hmu_del_R_norm_nn->GetYaxis()->SetTitle("a.u.");
  hmu_del_R_norm_nn->DrawNormalized();
  hmu_del_R_sig_nn->SetLineColor(kBlue);
  hmu_del_R_sig_nn->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("del_R_nn.png");

  //Missing mass Squared for Jpsi
  hmu_m2_miss_jpsi_norm_nn->SetLineColor(kRed);
  hmu_m2_miss_jpsi_norm_nn->GetXaxis()->SetTitle("m^{2}_{miss}[GeV^{2}]");
  hmu_m2_miss_jpsi_norm_nn->GetYaxis()->SetTitle("a.u.");
  hmu_m2_miss_jpsi_norm_nn->DrawNormalized();
  hmu_m2_miss_jpsi_sig_nn->SetLineColor(kBlue);
  hmu_m2_miss_jpsi_sig_nn->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("M2_miss_jpsi_nn.png");


  //Unpaired muon tranverse momentum
  hmu_mu_pt_norm_nn->SetLineColor(kRed);
  hmu_mu_pt_norm_nn->GetXaxis()->SetTitle("P_{T}^{#mu}[GeV]");
  hmu_mu_pt_norm_nn->GetYaxis()->SetTitle("a.u.");
  hmu_mu_pt_norm_nn->DrawNormalized();
  hmu_mu_pt_sig_nn->SetLineColor(kBlue);
  hmu_mu_pt_sig_nn->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("mu_pt_nn.png");


  //Unpaired muon psuedorapidity
  hmu_mu_eta_norm_nn->SetLineColor(kRed);
  hmu_mu_eta_norm_nn->GetXaxis()->SetTitle("#eta^{#mu}");
  hmu_mu_eta_norm_nn->GetYaxis()->SetTitle("a.u.");
  hmu_mu_eta_norm_nn->DrawNormalized();
  hmu_mu_eta_sig_nn->SetLineColor(kBlue);
  hmu_mu_eta_sig_nn->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("mu_eta_nn.png");

  
  //Unpaired muon phi
  hmu_mu_phi_norm_nn->SetLineColor(kRed);
  hmu_mu_phi_norm_nn->GetXaxis()->SetTitle("#phi^{#mu}");
  hmu_mu_phi_norm_nn->GetYaxis()->SetTitle("a.u.");
  hmu_mu_phi_norm_nn->DrawNormalized();
  hmu_mu_phi_sig_nn->SetLineColor(kBlue);
  hmu_mu_phi_sig_nn->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("mu_phi_nn.png");

  
  //Energy of unpaired muon in Bc rest frame
  hmu_mu_en_bc_frame_norm_nn->SetLineColor(kRed);
  hmu_mu_en_bc_frame_norm_nn->GetXaxis()->SetTitle("E(#mu) in rest frame of B_{c} (GeV)");
  hmu_mu_en_bc_frame_norm_nn->GetYaxis()->SetTitle("a.u.");
  hmu_mu_en_bc_frame_norm_nn->DrawNormalized();
  hmu_mu_en_bc_frame_sig_nn->SetLineColor(kBlue);
  hmu_mu_en_bc_frame_sig_nn->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("E_bc_nn.png");

  //Energy of unpaired muon in Jpsi rest frame
  hmu_mu_en_jpsi_frame_norm_nn->SetLineColor(kRed);
  hmu_mu_en_jpsi_frame_norm_nn->GetXaxis()->SetTitle("E(#mu) in rest frame of J/(#Psi) (GeV)");
  hmu_mu_en_jpsi_frame_norm_nn->GetYaxis()->SetTitle("a.u.");
  hmu_mu_en_jpsi_frame_norm_nn->DrawNormalized();
  hmu_mu_en_jpsi_frame_sig_nn->SetLineColor(kBlue);
  hmu_mu_en_jpsi_frame_sig_nn->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("E_jpsi_nn.png");\

  
  //Transverse momentum of Bc
  hmu_bc_pt_norm_nn->SetLineColor(kRed);
  hmu_bc_pt_norm_nn->GetXaxis()->SetTitle("P_{T}^{Bc} (GeV)");
  hmu_bc_pt_norm_nn->GetYaxis()->SetTitle("a.u.");
  hmu_bc_pt_norm_nn->DrawNormalized();
  hmu_bc_pt_sig_nn->SetLineColor(kBlue);
  hmu_bc_pt_sig_nn->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("Pt_Bc_nn.png");

  
  //Ratio of true and predicted transverse momentum of Bc
  hmu_bc_pt_ratio_norm_nn->SetLineColor(kRed);
  hmu_bc_pt_ratio_norm_nn->GetXaxis()->SetTitle("Ratio of predicted_nn and true P_{T}^{Bc}");
  hmu_bc_pt_ratio_norm_nn->GetYaxis()->SetTitle("a.u.");
  hmu_bc_pt_ratio_norm_nn->DrawNormalized();
  hmu_bc_pt_ratio_sig_nn->SetLineColor(kBlue);
  hmu_bc_pt_ratio_sig_nn->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("Ratio_pt_bc_nn.png");
  

  //Profile histogram of transverse momentum of Bc
  hprof_norm_nn->SetLineColor(kRed);
  hprof_norm_nn->GetXaxis()->SetTitle("Predicted_nn P_{T}^{Bc} (GeV)");
  hprof_norm_nn->GetYaxis()->SetTitle("Ratio of predicted_nn and true P_{T}^{Bc}");
  hprof_norm_nn->Draw();
  hprof_sig_nn->SetLineColor(kBlue);
  hprof_sig_nn->Draw("same");
  legend->Draw();
  c1->SaveAs("Profile_pt_bc_nn.png");

  //Scatter plot of Transverse momentum of Bc
  hmu_bc_pt_scatter_norm_nn->SetMarkerColor(kRed);
  hmu_bc_pt_scatter_norm_nn->GetXaxis()->SetTitle("True P_{T}^{Bc} (GeV)");
  hmu_bc_pt_scatter_norm_nn->GetYaxis()->SetTitle("Predicted_nn P_{T}^{Bc} (GeV)");
  hmu_bc_pt_scatter_norm_nn->DrawNormalized();
  hmu_bc_pt_scatter_sig_nn->SetMarkerColor(kBlue);
  hmu_bc_pt_scatter_sig_nn->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("scatter_Pt_Bc_nn.png");


  **/
  //Overlap of 2 possibilities (true and Jona recipe)

  //Comparison of True Bc Pt, Jona recipe for normalization channel
  hmu_bc_pt_corr_norm_gen->SetLineColor(kGreen);
  hmu_bc_pt_corr_norm_gen->GetXaxis()->SetTitle("P_{T}^{Bc} (GeV)");
  hmu_bc_pt_corr_norm_gen->GetYaxis()->SetTitle("a.u.");
  hmu_bc_pt_corr_norm_gen->DrawNormalized();
  hmu_bc_pt_norm_gen->SetLineColor(kRed);
  hmu_bc_pt_norm_gen->DrawNormalized("same");
  legend->AddEntry(hmu_bc_pt_norm_gen, "True", "l");
  legend->AddEntry(hmu_bc_pt_corr_norm_gen, "collinear approximation", "l");
  legend->Draw();
  c1->SaveAs("Compare2_pt_bc_norm.png");

  //Comparison of True Bc Pt, Jona recipe for signal channel
  //Comparison of Bc Pt and Reconstructed Bc Pt
  hmu_bc_pt_corr_sig_gen->SetLineColor(kGreen);
  hmu_bc_pt_corr_sig_gen->GetXaxis()->SetTitle("P_{T}^{Bc} (GeV)");
  hmu_bc_pt_corr_sig_gen->GetYaxis()->SetTitle("a.u.");
  hmu_bc_pt_corr_sig_gen->DrawNormalized();
  hmu_bc_pt_sig_gen->SetLineColor(kRed);
  hmu_bc_pt_sig_gen->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("Compare2_pt_bc_sig.png");

  /**
  //Comparison of True m2 miss, Jona recipe for normalization channel
  hmu_m2_miss_corr_norm_gen->SetLineColor(kGreen);
  hmu_m2_miss_corr_norm_gen->GetXaxis()->SetTitle("m^{2}_{miss}[GeV^{2}]");
  hmu_m2_miss_corr_norm_gen->GetYaxis()->SetTitle("a.u.");
  hmu_m2_miss_corr_norm_gen->DrawNormalized();
  hmu_m2_miss_norm_gen->SetLineColor(kRed);
  hmu_m2_miss_norm_gen->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("Compare2_m2_miss_norm.png");

  //Comparison of True m2 miss, Jona recipe for signal channel 
  hmu_m2_miss_corr_sig_gen->SetLineColor(kGreen);
  hmu_m2_miss_corr_sig_gen->GetXaxis()->SetTitle("m^{2}_{miss}[GeV^{2}]");
  hmu_m2_miss_corr_sig_gen->GetYaxis()->SetTitle("a.u.");
  hmu_m2_miss_corr_sig_gen->DrawNormalized();
  hmu_m2_miss_sig_gen->SetLineColor(kRed);
  hmu_m2_miss_sig_gen->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("Compare2_m2_miss_sig.png");


  //Comparison of True Q2, Jona recipe for normalization channel
  hmu_Q2_corr_norm_gen->SetLineColor(kGreen);
  hmu_Q2_corr_norm_gen->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
  hmu_Q2_corr_norm_gen->GetYaxis()->SetTitle("a.u.");
  hmu_Q2_corr_norm_gen->DrawNormalized();
  hmu_Q2_norm_gen->SetLineColor(kRed);
  hmu_Q2_norm_gen->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("Compare2_Q2_norm.png");

  //Comparison of True Q2, Jona recipe for signal channel
  hmu_Q2_corr_sig_gen->SetLineColor(kGreen); 
  hmu_Q2_corr_sig_gen->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
  hmu_Q2_corr_sig_gen->GetYaxis()->SetTitle("a.u.");
  hmu_Q2_corr_sig_gen->DrawNormalized();
  hmu_Q2_sig_gen->SetLineColor(kRed);
  hmu_Q2_sig_gen->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("Compare2_Q2_sig.png");


  //Comparison of True pt miss, Jona recipe for normalization channel
  hmu_Pt_miss_corr_norm_gen->SetLineColor(kGreen);
  hmu_Pt_miss_corr_norm_gen->GetXaxis()->SetTitle("P_{T}^{miss} (GeV)");
  hmu_Pt_miss_corr_norm_gen->GetYaxis()->SetTitle("a.u.");
  hmu_Pt_miss_corr_norm_gen->DrawNormalized();
  hmu_Pt_miss_norm_gen->SetLineColor(kRed);
  hmu_Pt_miss_norm_gen->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("Compare2_Pt_miss_norm.png");

  //Comparison of True pt miss, Jona recipe for signal channel
  hmu_Pt_miss_corr_sig_gen->SetLineColor(kGreen);
  hmu_Pt_miss_corr_sig_gen->GetXaxis()->SetTitle("P_{T}^{miss} (GeV)");
  hmu_Pt_miss_corr_sig_gen->GetYaxis()->SetTitle("a.u.");
  hmu_Pt_miss_corr_sig_gen->DrawNormalized();
  hmu_Pt_miss_sig_gen->SetLineColor(kRed);
  hmu_Pt_miss_sig_gen->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("Compare2_Pt_miss_sig.png");

  
  
  //Overlap of all three possibilities
  
  //Comparison of True Bc Pt, Jona recipe and prediction from neural network for normalization channel
  hmu_bc_pt_norm_nn->SetLineColor(kBlack);
  hmu_bc_pt_norm_nn->GetXaxis()->SetTitle("P_{T}^{Bc} (GeV)");
  hmu_bc_pt_norm_nn->GetYaxis()->SetTitle("a.u.");
  hmu_bc_pt_norm_nn->DrawNormalized();
  hmu_bc_pt_corr_norm_gen->SetLineColor(kGreen);
  hmu_bc_pt_corr_norm_gen->DrawNormalized("same");
  hmu_bc_pt_norm_gen->SetLineColor(kRed);
  hmu_bc_pt_norm_gen->DrawNormalized("same");
  legend->AddEntry(hmu_bc_pt_norm_gen, "True", "l");
  legend->AddEntry(hmu_bc_pt_corr_norm_gen, "Collinear approximation", "l");
  legend->AddEntry(hmu_bc_pt_norm_nn, "nn_predicted", "l");
  legend->Draw();
  c1->SaveAs("Compare3_pt_bc_norm.pdf");

  //Comparison of True Bc Pt, Jona recipe and prediction from neural network for signal channel
  //Comparison of Bc Pt and Reconstructed Bc Pt
  hmu_bc_pt_sig_nn->SetLineColor(kBlack);
  hmu_bc_pt_sig_nn->GetXaxis()->SetTitle("P_{T}^{Bc} (GeV)");
  hmu_bc_pt_sig_nn->GetYaxis()->SetTitle("a.u.");
  hmu_bc_pt_sig_nn->DrawNormalized();
  hmu_bc_pt_corr_sig_gen->SetLineColor(kGreen);
  hmu_bc_pt_corr_sig_gen->DrawNormalized("same");
  hmu_bc_pt_sig_gen->SetLineColor(kRed);
  hmu_bc_pt_sig_gen->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("Compare3_pt_bc_sig.pdf");

  //Comparison of True m2 miss, Jona recipe and prediction from neural network for normalization channel
  hmu_m2_miss_norm_nn->SetLineColor(kBlack);
  hmu_m2_miss_norm_nn->GetXaxis()->SetTitle("m^{2}_{miss}[GeV^{2}]");
  hmu_m2_miss_norm_nn->GetYaxis()->SetTitle("a.u.");
  hmu_m2_miss_norm_nn->DrawNormalized();
  hmu_m2_miss_corr_norm_gen->SetLineColor(kGreen);
  hmu_m2_miss_corr_norm_gen->DrawNormalized("same");
  hmu_m2_miss_norm_gen->SetLineColor(kRed);
  hmu_m2_miss_norm_gen->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("Compare3_m2_miss_norm.pdf");

  //Comparison of True m2 miss, Jona recipe and prediction from neural network for signal channel
  hmu_m2_miss_sig_nn->SetLineColor(kBlack);
  hmu_m2_miss_sig_nn->GetXaxis()->SetTitle("m^{2}_{miss}[GeV^{2}]");
  hmu_m2_miss_sig_nn->GetYaxis()->SetTitle("a.u.");
  hmu_m2_miss_sig_nn->DrawNormalized();
  hmu_m2_miss_corr_sig_gen->SetLineColor(kGreen);
  hmu_m2_miss_corr_sig_gen->DrawNormalized("same");
  hmu_m2_miss_sig_gen->SetLineColor(kRed);
  hmu_m2_miss_sig_gen->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("Compare3_m2_miss_sig.pdf");


  //Comparison of True Q2, Jona recipe and prediction from neural network for normalization channel
  hmu_Q2_norm_nn->SetLineColor(kBlack);
  hmu_Q2_norm_nn->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
  hmu_Q2_norm_nn->GetYaxis()->SetTitle("a.u.");
  hmu_Q2_norm_nn->DrawNormalized();
  hmu_Q2_corr_norm_gen->SetLineColor(kGreen);
  hmu_Q2_corr_norm_gen->DrawNormalized("same");
  hmu_Q2_norm_gen->SetLineColor(kRed);
  hmu_Q2_norm_gen->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("Compare3_Q2_norm.pdf");

  //Comparison of True Q2, Jona recipe and prediction from neural network for signal channel
  hmu_Q2_sig_nn->SetLineColor(kBlack);
  hmu_Q2_sig_nn->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
  hmu_Q2_sig_nn->GetYaxis()->SetTitle("a.u.");
  hmu_Q2_sig_nn->DrawNormalized();
  hmu_Q2_corr_sig_gen->SetLineColor(kGreen);
  hmu_Q2_corr_sig_gen->DrawNormalized("same");
  hmu_Q2_sig_gen->SetLineColor(kRed);
  hmu_Q2_sig_gen->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("Compare3_Q2_sig.pdf");

  
  //Comparison of True pt miss, Jona recipe and prediction from neural network for normalization channel
  hmu_Pt_miss_norm_nn->SetLineColor(kBlack);
  hmu_Pt_miss_norm_nn->GetXaxis()->SetTitle("P_{T}^{miss} (GeV)");
  hmu_Pt_miss_norm_nn->GetYaxis()->SetTitle("a.u.");
  hmu_Pt_miss_norm_nn->DrawNormalized();
  hmu_Pt_miss_corr_norm_gen->SetLineColor(kGreen);
  hmu_Pt_miss_corr_norm_gen->DrawNormalized("same");
  hmu_Pt_miss_norm_gen->SetLineColor(kRed);
  hmu_Pt_miss_norm_gen->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("Compare3_Pt_miss_norm.pdf");

  //Comparison of True pt miss, Jona recipe and prediction from neural network for signal channel
  hmu_Pt_miss_sig_nn->SetLineColor(kBlack);
  hmu_Pt_miss_sig_nn->GetXaxis()->SetTitle("P_{T}^{miss} (GeV)");
  hmu_Pt_miss_sig_nn->GetYaxis()->SetTitle("a.u.");
  hmu_Pt_miss_sig_nn->DrawNormalized();
  hmu_Pt_miss_corr_sig_gen->SetLineColor(kGreen);
  hmu_Pt_miss_corr_sig_gen->DrawNormalized("same");
  hmu_Pt_miss_sig_gen->SetLineColor(kRed);
  hmu_Pt_miss_sig_gen->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("Compare3_Pt_miss_sig.pdf");
  


  
  hmu_m2_miss_norm_nn->Add(hmu_m2_miss_sig_nn);
  hmu_ratio_m2_miss_nn = hmu_m2_miss_sig_nn/hmu_m2_miss_norm_nn;
  hmu_ratio_m2_miss_nn->SetLineColor(kBlack);
  hmu_ratio_m2_miss_nn->Draw();
  c1->SaveAs("Ratio_m2_miss_nn");
  **/

  //Closing the files
  //f1->Close();
  //f2->Close();
  f3->Close();
  f4->Close();
   
}
