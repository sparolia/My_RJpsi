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

void Macro(){

  //declaring all the required variable for normalization channel
  float bc_pt_norm_gen=0,bc_eta_norm_gen=0,bc_phi_norm_gen=0,bc_mass_norm_gen=0;
  float mu_pt_norm_gen=0,mu_eta_norm_gen=0,mu_phi_norm_gen=0;
  float mu1_pt_norm_gen=0,mu1_eta_norm_gen=0,mu1_phi_norm_gen=0;
  float mu2_pt_norm_gen=0,mu2_eta_norm_gen=0,mu2_phi_norm_gen=0;
  float jpsi_pt_norm_gen=0,jpsi_eta_norm_gen=0,jpsi_phi_norm_gen=0,jpsi_mass_norm_gen=0;
  float Trig_dimuon_norm_gen=0,Trig_JpsiTk_norm_gen=0;
  float munu_pt_norm_gen=0,munu_eta_norm_gen=0,munu_phi_norm_gen=0;

  //declaring all the required variable for signal channel
  float bc_pt_sig_gen=0,bc_eta_sig_gen=0,bc_phi_sig_gen=0,bc_mass_sig_gen=0;
  float mu_pt_sig_gen=0,mu_eta_sig_gen=0,mu_phi_sig_gen=0;
  float mu1_pt_sig_gen=0,mu1_eta_sig_gen=0,mu1_phi_sig_gen=0;
  float mu2_pt_sig_gen=0,mu2_eta_sig_gen=0,mu2_phi_sig_gen=0;
  float jpsi_pt_sig_gen=0,jpsi_eta_sig_gen=0,jpsi_phi_sig_gen=0,jpsi_mass_sig_gen=0;
  float Trig_dimuon_sig_gen=0,Trig_JpsiTk_sig_gen=0;
  float taunu_pt_sig_gen=0,taunu_eta_sig_gen=0,taunu_phi_sig_gen=0;
  float taunu2_pt_sig_gen=0,taunu2_eta_sig_gen=0,taunu2_phi_sig_gen=0;
  float munu_pt_sig_gen=0,munu_eta_sig_gen=0,munu_phi_sig_gen=0;
  
  TFile *f1 = TFile::Open("bc_jpsi_mu_nu_gen_v2.root"); //opening the norm_genalization root file
  TTree *tree_norm_gen = (TTree*)f1->Get("tree"); //reading the tree from the ntuple
  TFile *f2 = TFile::Open("bc_jpsi_tau_nu_gen_v2.root"); //opening the sig_gennal root file
  TTree *tree_sig_gen = (TTree*)f2->Get("tree"); //reading the tree from the ntuple


  //reading variables from the normalization tree
  tree_norm_gen->SetBranchAddress("bc_pt",&bc_pt_norm_gen); 
  tree_norm_gen->SetBranchAddress("bc_eta",&bc_eta_norm_gen); 
  tree_norm_gen->SetBranchAddress("bc_phi",&bc_phi_norm_gen); 
  tree_norm_gen->SetBranchAddress("bc_mass",&bc_mass_norm_gen); 

  tree_norm_gen->SetBranchAddress("mu_pt",&mu_pt_norm_gen); 
  tree_norm_gen->SetBranchAddress("mu_eta",&mu_eta_norm_gen); 
  tree_norm_gen->SetBranchAddress("mu_phi",&mu_phi_norm_gen);

  tree_norm_gen->SetBranchAddress("mu1_pt",&mu1_pt_norm_gen); 
  tree_norm_gen->SetBranchAddress("mu1_eta",&mu1_eta_norm_gen); 
  tree_norm_gen->SetBranchAddress("mu1_phi",&mu1_phi_norm_gen);

  tree_norm_gen->SetBranchAddress("mu2_pt",&mu2_pt_norm_gen); 
  tree_norm_gen->SetBranchAddress("mu2_eta",&mu2_eta_norm_gen); 
  tree_norm_gen->SetBranchAddress("mu2_phi",&mu2_phi_norm_gen);

  tree_norm_gen->SetBranchAddress("jpsi_pt",&jpsi_pt_norm_gen); 
  tree_norm_gen->SetBranchAddress("jpsi_eta",&jpsi_eta_norm_gen); 
  tree_norm_gen->SetBranchAddress("jpsi_phi",&jpsi_phi_norm_gen);
  tree_norm_gen->SetBranchAddress("jpsi_mass",&jpsi_mass_norm_gen);

  tree_norm_gen->SetBranchAddress("trigDimuon0",&Trig_dimuon_norm_gen);
  tree_norm_gen->SetBranchAddress("trigJpsiTk",&Trig_JpsiTk_norm_gen);

  tree_norm_gen->SetBranchAddress("munu_pt",&munu_pt_norm_gen); 
  tree_norm_gen->SetBranchAddress("munu_eta",&munu_eta_norm_gen); 
  tree_norm_gen->SetBranchAddress("munu_phi",&munu_phi_norm_gen);

  //reading variables from the signal tree
  tree_sig_gen->SetBranchAddress("bc_pt",&bc_pt_sig_gen); 
  tree_sig_gen->SetBranchAddress("bc_eta",&bc_eta_sig_gen); 
  tree_sig_gen->SetBranchAddress("bc_phi",&bc_phi_sig_gen); 
  tree_sig_gen->SetBranchAddress("bc_mass",&bc_mass_sig_gen);

  tree_sig_gen->SetBranchAddress("mu_pt",&mu_pt_sig_gen); 
  tree_sig_gen->SetBranchAddress("mu_eta",&mu_eta_sig_gen); 
  tree_sig_gen->SetBranchAddress("mu_phi",&mu_phi_sig_gen);

  tree_sig_gen->SetBranchAddress("mu1_pt",&mu1_pt_sig_gen); 
  tree_sig_gen->SetBranchAddress("mu1_eta",&mu1_eta_sig_gen); 
  tree_sig_gen->SetBranchAddress("mu1_phi",&mu1_phi_sig_gen);

  tree_sig_gen->SetBranchAddress("mu2_pt",&mu2_pt_sig_gen); 
  tree_sig_gen->SetBranchAddress("mu2_eta",&mu2_eta_sig_gen); 
  tree_sig_gen->SetBranchAddress("mu2_phi",&mu2_phi_sig_gen);

  tree_sig_gen->SetBranchAddress("jpsi_pt",&jpsi_pt_sig_gen); 
  tree_sig_gen->SetBranchAddress("jpsi_eta",&jpsi_eta_sig_gen); 
  tree_sig_gen->SetBranchAddress("jpsi_phi",&jpsi_phi_sig_gen); 
  tree_sig_gen->SetBranchAddress("jpsi_mass",&jpsi_mass_sig_gen);

  tree_sig_gen->SetBranchAddress("trigDimuon0",&Trig_dimuon_sig_gen);
  tree_sig_gen->SetBranchAddress("trigJpsiTk",&Trig_JpsiTk_sig_gen);

  tree_sig_gen->SetBranchAddress("taunu_pt",&taunu_pt_sig_gen); 
  tree_sig_gen->SetBranchAddress("taunu_eta",&taunu_eta_sig_gen); 
  tree_sig_gen->SetBranchAddress("taunu_phi",&taunu_phi_sig_gen);

  tree_sig_gen->SetBranchAddress("taunu2_pt",&taunu2_pt_sig_gen); 
  tree_sig_gen->SetBranchAddress("taunu2_eta",&taunu2_eta_sig_gen); 
  tree_sig_gen->SetBranchAddress("taunu2_phi",&taunu2_phi_sig_gen);

  tree_sig_gen->SetBranchAddress("munu_pt",&munu_pt_sig_gen); 
  tree_sig_gen->SetBranchAddress("munu_eta",&munu_eta_sig_gen); 
  tree_sig_gen->SetBranchAddress("munu_phi",&munu_phi_sig_gen);
  

  //Defining Canvas, removing statistics and setting the legend
  auto c1 = new TCanvas("c", "c", 600,500);
  //gStyle->SetOptStat(0);
  auto* legend = new TLegend(0.7, 0.3, 0.9, 0.5);
   

  //booking the histograms
  TH1D* hmu_m2_miss_norm_gen = new TH1D("hmu_m2_miss_norm_gen","Missing mass squared",80,0,8); 
  TH1D* hmu_m2_miss_sig_gen = new TH1D("hmu_m2_miss_sig_gen","Missing mass squared",80,0,8);
  TH1D* hmu_Q2_norm_gen = new TH1D("hmu_Q2_norm_gen","Squared four momentum transfered to lepton system",100,0,10);
  TH1D* hmu_Q2_sig_gen = new TH1D("hmu_Q2_sig_gen","Squared four momentum transfered to lepton system",100,0,10);
  TH1D* hmu_Pt_miss_norm_gen = new TH1D("hmu_Pt_miss_norm_gen","Missing transverse momentum",200,0,20);
  TH1D* hmu_Pt_miss_sig_gen = new TH1D("hmu_Pt_miss_sig_gen","Missing transverse momentum",200,0,20);
  TH1D* hmu_m2_miss_corr_norm_gen = new TH1D("hmu_m2_miss_corr_norm_gen","Missing mass squared corrected",400,-20,20); 
  TH1D* hmu_m2_miss_corr_sig_gen = new TH1D("hmu_m2_miss_corr_sig_gen","Missing mass squared corrected",400,-20,20);
  TH1D* hmu_Q2_corr_norm_gen = new TH1D("hmu_Q2_corr_norm_gen","Squared four momentum transfered to lepton system corrected",100,0,10);
  TH1D* hmu_Q2_corr_sig_gen = new TH1D("hmu_Q2_corr_sig_gen","Squared four momentum transfered to lepton system corrected",100,0,10);
  TH1D* hmu_Pt_miss_corr_norm_gen = new TH1D("hmu_Pt_miss_corr_norm_gen","Missing transverse momentum corrected",500,0,50);
  TH1D* hmu_Pt_miss_corr_sig_gen = new TH1D("hmu_Pt_miss_corr_sig_gen","Missing transverse momentum corrected",500,0,50);
  TH1D* hmu_Pt_var_norm_gen = new TH1D("hmu_Pt_var_norm_gen","Tranverse variable momentum",600,-30,30);
  TH1D* hmu_Pt_var_sig_gen = new TH1D("hmu_Pt_var_sig_gen","Tranverse variable momentum",600,-30,30);
  TH1D* hmu_del_R_norm_gen = new TH1D("hmu_del_R_norm_gen","Delta R",200,0,2);
  TH1D* hmu_del_R_sig_gen = new TH1D("hmu_del_R_sig_gen","Delta R",200,0,2);
  TH1D* hmu_mu_pt_norm_gen = new TH1D("hmu_mu_pt_norm_gen","Unpaired muon tranverse momentum",200,0,20);
  TH1D* hmu_mu_pt_sig_gen = new TH1D("hmu_mu_pt_sig_gen","Unpaired muon tranverse momentum",200,0,20);
  TH1D* hmu_mu_px_norm_gen = new TH1D("hmu_mu_px_norm_gen","Unpaired muon momentum in x-direction",200,0,50);
  TH1D* hmu_mu_px_sig_gen = new TH1D("hmu_mu_px_sig_gen","Unpaired muon momentum in x-direction",200,0,50);
  TH1D* hmu_mu_py_norm_gen = new TH1D("hmu_mu_py_norm_gen","Unpaired muon momentum in y-direction",200,0,50);
  TH1D* hmu_mu_py_sig_gen = new TH1D("hmu_mu_py_sig_gen","Unpaired muon momentum in y-direction",200,0,50);
  TH1D* hmu_mu_pz_norm_gen = new TH1D("hmu_mu_pz_norm_gen","Unpaired muon momentum in z-direction",200,0,50);
  TH1D* hmu_mu_pz_sig_gen = new TH1D("hmu_mu_pz_sig_gen","Unpaired muon momentum in z-direction",200,0,50);
  TH1D* hmu_mu_eta_norm_gen = new TH1D("hmu_mu_eta_norm_gen","Unpaired muon psuedorapidity",60,-3,3);
  TH1D* hmu_mu_eta_sig_gen = new TH1D("hmu_mu_eta_sig_gen","Unpaired muon psuedorapidity",60,-3,3);
  TH1D* hmu_mu_phi_norm_gen = new TH1D("hmu_mu_phi_norm_gen","Unpaired muon phi",100,-5,5);
  TH1D* hmu_mu_phi_sig_gen = new TH1D("hmu_mu_phi_sig_gen","Unpaired muon phi",100,-5,5);
  TH1D* hmu_mu_en_bc_frame_norm_gen = new TH1D("hmu_mu_en_bc_frame_norm_gen","Energy of Unpaired muon in Bc rest frame",40,0,4);
  TH1D* hmu_mu_en_bc_frame_sig_gen = new TH1D("hmu_mu_en_bc_frame_sig_gen","Energy of Unpaired muon in Bc rest frame",40,0,4);
  TH1D* hmu_mu_en_jpsi_frame_norm_gen = new TH1D("hmu_mu_en_jpsi_frame_norm_gen","Energy of Unpaired muon in Jpsi rest frame",50,0,5);
  TH1D* hmu_mu_en_jpsi_frame_sig_gen = new TH1D("hmu_mu_en_jpsi_frame_sig_gen","Energy of Unpaired muon in Jpsi rest frame",50,0,5);
  TH1D* hmu_bc_pt_norm_gen = new TH1D("hmu_bc_pt_norm_gen","Transverse momentum of B_{c}",100,0,100);
  TH1D* hmu_bc_pt_sig_gen = new TH1D("hmu_bc_pt_sig_gen","Transverse momentum of B_{c}",100,0,100);
  TH1D* hmu_bc_pt_ratio_norm_gen = new TH1D("hmu_bc_pt_ratio_norm_gen","Ratio of true and corrected transverse momentum of B_{c}",100,0,2);
  TH1D* hmu_bc_pt_ratio_sig_gen = new TH1D("hmu_bc_pt_ratio_sig_gen","Ratio of true and corrected transverse momentum of B_{c}",100,0,2);
  TH1D* hmu_bc_pt_corr_norm_gen = new TH1D("hmu_bc_pt_corr_norm_gen","Corrected transverse momentum of B_{c}",100,0,100);
  TH1D* hmu_bc_pt_corr_sig_gen = new TH1D("hmu_bc_pt_corr_sig_gen","Corrected transverse momentum of B_{c}",100,0,100);
  auto hprof_norm_gen  = new TProfile("hprof_norm_gen","Profile of Transverse momentum of B_{c}",100,0,100,0,1000);
  auto hprof_sig_gen  = new TProfile("hprof_sig_gen","Profile of Transverse momentum of B_{c}",100,0,100,0,1000);
  TH2D* hmu_bc_pt_scatter_norm_gen = new TH2D("hmu_bc_pt_scatter_norm_gen","Scatter plot of transverse momentum of B_{c}",100,0,100,100,0,100);
  TH2D* hmu_bc_pt_scatter_sig_gen = new TH2D("hmu_bc_pt_scatter_sig_gen","Scatter plot of transverse momentum of B_{c}",100,0,100,100,0,100);
  TH1D* hmu_recco_bc_pt_norm_gen = new TH1D("hmu_recco_bc_pt_norm_gen","Reconstructed transverse momentum of B_{c}",100,0,100);
  TH1D* hmu_recco_bc_pt_sig_gen = new TH1D("hmu_recco_bc_pt_sig_gen","Reconstructed transverse momentum of B_{c}",100,0,100);
  
 


  //ParticleMass muonMass = 0.10565837;
  float mu_mass = 0.113428;

  //loop over events for norm_genalized channel
  for(int i=0; i<tree_norm_gen->GetEntries();i++){
    tree_norm_gen->GetEntry(i);
    
    //TLorentzVector class for calculating four momenta
    TLorentzVector momntm_bc_norm_gen;
    momntm_bc_norm_gen.SetPtEtaPhiM(bc_pt_norm_gen, bc_eta_norm_gen,bc_phi_norm_gen,bc_mass_norm_gen);
    TLorentzVector momntm_mu1_norm_gen;
    momntm_mu1_norm_gen.SetPtEtaPhiM(mu1_pt_norm_gen, mu1_eta_norm_gen,mu1_phi_norm_gen,mu_mass);
    TLorentzVector momntm_mu2_norm_gen;
    momntm_mu2_norm_gen.SetPtEtaPhiM(mu2_pt_norm_gen, mu2_eta_norm_gen,mu2_phi_norm_gen,mu_mass);
    TLorentzVector momntm_mu_norm_gen;
    momntm_mu_norm_gen.SetPtEtaPhiM(mu_pt_norm_gen, mu_eta_norm_gen,mu_phi_norm_gen,mu_mass);
    TLorentzVector momntm_jpsi_norm_gen;
    momntm_jpsi_norm_gen.SetPtEtaPhiM(jpsi_pt_norm_gen, jpsi_eta_norm_gen,jpsi_phi_norm_gen,jpsi_mass_norm_gen);

    //calculating the kinematic variables
    float m2_miss_norm_gen =(momntm_bc_norm_gen - momntm_mu1_norm_gen - momntm_mu2_norm_gen - momntm_mu_norm_gen).M2();
    float Q2_norm_gen = (momntm_bc_norm_gen - momntm_mu1_norm_gen - momntm_mu2_norm_gen).M2();
    float Pt_miss_norm_gen = (bc_pt_norm_gen - mu1_pt_norm_gen - mu2_pt_norm_gen - mu_pt_norm_gen);
    float Pt_var_norm_gen = (jpsi_pt_norm_gen - mu_pt_norm_gen );
    //float del_R_norm_gen = sqrt(pow((mu1_phi_norm_gen - mu2_phi_norm_gen),2)+pow((mu1_eta_norm_gen - mu2_eta_norm_gen),2));
    float del_R_norm_gen = momntm_mu1_norm_gen.DeltaR(momntm_mu2_norm_gen);
    TLorentzVector muon_system_norm_gen = (momntm_mu1_norm_gen+momntm_mu2_norm_gen+momntm_mu_norm_gen);
    //cout << muon_system_mass_norm_gen ;
    //float muon_system_pt_norm_gen = mu_pt_norm_gen + mu1_pt_norm_gen + mu2_pt_norm_gen;
    float bc_pt_corr_norm_gen = (bc_mass_norm_gen * muon_system_norm_gen.Pt())/ (muon_system_norm_gen.M()); //corrected Pt according to recipe in Jona thesis
    float bc_pt_ratio_norm_gen = bc_pt_corr_norm_gen/bc_pt_norm_gen;
    
    TLorentzVector momntm_bc_corr_norm_gen;
    momntm_bc_corr_norm_gen.SetPtEtaPhiM(bc_pt_corr_norm_gen, muon_system_norm_gen.Eta(), muon_system_norm_gen.Phi(),momntm_bc_norm_gen.M());
    float Q2_corr_norm_gen = (momntm_bc_corr_norm_gen - momntm_mu1_norm_gen - momntm_mu2_norm_gen).M2();
    float m2_miss_corr_norm_gen = (momntm_bc_corr_norm_gen - momntm_mu1_norm_gen - momntm_mu2_norm_gen - momntm_mu_norm_gen).M2();
    float Pt_miss_corr_norm_gen = (bc_pt_corr_norm_gen - mu1_pt_norm_gen - mu2_pt_norm_gen - mu_pt_norm_gen);

    TVector3 X_bc_norm_gen = -momntm_bc_corr_norm_gen.BoostVector();
    TVector3 X_jpsi_norm_gen = -momntm_jpsi_norm_gen.BoostVector();
    TLorentzVector momntm_mu_bc_frame_norm_gen = momntm_mu_norm_gen;
    momntm_mu_bc_frame_norm_gen.Boost(X_bc_norm_gen);
    TLorentzVector momntm_mu_jpsi_frame_norm_gen = momntm_mu_norm_gen;
    momntm_mu_jpsi_frame_norm_gen.Boost(X_jpsi_norm_gen);

    float mu_en_bc_frame_norm_gen = momntm_mu_bc_frame_norm_gen.E();
    float mu_en_jpsi_frame_norm_gen = momntm_mu_jpsi_frame_norm_gen.E();

    float recco_bc_pt_norm_gen = mu1_pt_norm_gen + mu2_pt_norm_gen + mu_pt_norm_gen + munu_pt_norm_gen;
    float mu_px_norm_gen = momntm_mu_norm_gen.Px();
    float mu_py_norm_gen = momntm_mu_norm_gen.Py();
    float mu_pz_norm_gen = momntm_mu_norm_gen.Pz();
        
    //Filling the histograms

    if ((abs(mu_eta_norm_gen) < 1.3 and mu_pt_norm_gen > 3.3) or (1.3 < abs(mu_eta_norm_gen) < 2.2 and mu_pt_norm_gen > 2.9) or (2.2 < abs(mu_eta_norm_gen) < 2.4 and mu_pt_norm_gen > 2.4)){ //acceptance cuts
      if ((Trig_dimuon_norm_gen == false) and (Trig_JpsiTk_norm_gen == true)){ //trigger flags
	hmu_m2_miss_norm_gen->Fill(m2_miss_norm_gen);
	hmu_Q2_norm_gen->Fill(Q2_norm_gen);
	hmu_Pt_miss_norm_gen->Fill(Pt_miss_norm_gen);
	hmu_m2_miss_corr_norm_gen->Fill(m2_miss_corr_norm_gen);
	hmu_Q2_corr_norm_gen->Fill(Q2_corr_norm_gen);
	hmu_Pt_miss_corr_norm_gen->Fill(Pt_miss_corr_norm_gen);
	hmu_Pt_var_norm_gen->Fill(Pt_var_norm_gen);
	hmu_mu_pt_norm_gen->Fill(mu_pt_norm_gen);
	hmu_mu_px_norm_gen->Fill(mu_px_norm_gen);
	hmu_mu_py_norm_gen->Fill(mu_py_norm_gen);
	hmu_mu_pz_norm_gen->Fill(mu_pz_norm_gen);
	hmu_mu_eta_norm_gen->Fill(mu_eta_norm_gen);
	hmu_mu_phi_norm_gen->Fill(mu_phi_norm_gen);
	hmu_del_R_norm_gen->Fill(del_R_norm_gen);
	hmu_mu_en_bc_frame_norm_gen->Fill(mu_en_bc_frame_norm_gen);
	hmu_mu_en_jpsi_frame_norm_gen->Fill(mu_en_jpsi_frame_norm_gen);
	hmu_bc_pt_norm_gen->Fill(bc_pt_norm_gen);
	hmu_bc_pt_corr_norm_gen->Fill(bc_pt_corr_norm_gen);
	hmu_bc_pt_ratio_norm_gen->Fill(bc_pt_ratio_norm_gen);
	hprof_norm_gen->Fill(bc_pt_norm_gen,bc_pt_ratio_norm_gen,1);
	hmu_bc_pt_scatter_norm_gen->Fill(bc_pt_norm_gen,bc_pt_corr_norm_gen);
	hmu_recco_bc_pt_norm_gen->Fill(recco_bc_pt_norm_gen);
      }
    }
  }

  //loop over events for sig_gennal channel
  for(int i=0; i<tree_sig_gen->GetEntries();i++){
    tree_sig_gen->GetEntry(i);

    //TLorentzVector class for calculating four momenta
    TLorentzVector momntm_bc_sig_gen;
    momntm_bc_sig_gen.SetPtEtaPhiM(bc_pt_sig_gen, bc_eta_sig_gen,bc_phi_sig_gen,bc_mass_sig_gen);
    TLorentzVector momntm_mu1_sig_gen;
    momntm_mu1_sig_gen.SetPtEtaPhiM(mu1_pt_sig_gen, mu1_eta_sig_gen,mu1_phi_sig_gen,mu_mass);
    TLorentzVector momntm_mu2_sig_gen;
    momntm_mu2_sig_gen.SetPtEtaPhiM(mu2_pt_sig_gen, mu2_eta_sig_gen,mu2_phi_sig_gen,mu_mass);
    TLorentzVector momntm_mu_sig_gen;
    momntm_mu_sig_gen.SetPtEtaPhiM(mu_pt_sig_gen, mu_eta_sig_gen,mu_phi_sig_gen,mu_mass);
    TLorentzVector momntm_jpsi_sig_gen;
    momntm_jpsi_sig_gen.SetPtEtaPhiM(jpsi_pt_sig_gen, jpsi_eta_sig_gen,jpsi_phi_sig_gen,jpsi_mass_sig_gen);

    //calculating the kinematic variables
    float m2_miss_sig_gen = (momntm_bc_sig_gen - momntm_mu1_sig_gen - momntm_mu2_sig_gen - momntm_mu_sig_gen).M2();
    float Q2_sig_gen = (momntm_bc_sig_gen - momntm_mu1_sig_gen - momntm_mu2_sig_gen).M2();
    float Pt_miss_sig_gen = (bc_pt_sig_gen-mu1_pt_sig_gen-mu2_pt_sig_gen-mu_pt_sig_gen);
    float Pt_var_sig_gen = (jpsi_pt_sig_gen - mu_pt_sig_gen );
    //float del_R_sig_gen = sqrt(pow((mu1_phi_sig_gen - mu2_phi_sig_gen),2)+pow((mu1_eta_sig_gen - mu2_eta_sig_gen),2));
    float del_R_sig_gen = momntm_mu1_sig_gen.DeltaR(momntm_mu2_sig_gen);
    TLorentzVector muon_system_sig_gen = (momntm_mu1_sig_gen+momntm_mu2_sig_gen+momntm_mu_sig_gen);
    //cout << muon_system_mass_sig_gen ;
    //float muon_system_pt_sig_gen = mu_pt_sig_gen + mu1_pt_sig_gen + mu2_pt_sig_gen;
    float bc_pt_corr_sig_gen = (bc_mass_sig_gen * muon_system_sig_gen.Pt())/ (muon_system_sig_gen.M());
    float bc_pt_ratio_sig_gen = bc_pt_corr_sig_gen/bc_pt_sig_gen;

    TLorentzVector momntm_bc_corr_sig_gen;
    momntm_bc_corr_sig_gen.SetPtEtaPhiM(bc_pt_corr_sig_gen, muon_system_sig_gen.Eta(), muon_system_sig_gen.Phi(), momntm_bc_sig_gen.M());
    float Q2_corr_sig_gen = (momntm_bc_corr_sig_gen - momntm_mu1_sig_gen - momntm_mu2_sig_gen).M2();
    float m2_miss_corr_sig_gen = (momntm_bc_corr_sig_gen - momntm_mu1_sig_gen - momntm_mu2_sig_gen - momntm_mu_sig_gen).M2();
    float Pt_miss_corr_sig_gen = (bc_pt_corr_sig_gen - mu1_pt_sig_gen - mu2_pt_sig_gen - mu_pt_sig_gen);
   
    TVector3 X_bc_sig_gen = -momntm_bc_corr_sig_gen.BoostVector();
    TVector3 X_jpsi_sig_gen = -momntm_jpsi_sig_gen.BoostVector();
    TLorentzVector momntm_mu_bc_frame_sig_gen = momntm_mu_sig_gen;
    momntm_mu_bc_frame_sig_gen.Boost(X_bc_sig_gen);
    TLorentzVector momntm_mu_jpsi_frame_sig_gen = momntm_mu_sig_gen;
    momntm_mu_jpsi_frame_sig_gen.Boost(X_jpsi_sig_gen);
    
    float mu_en_bc_frame_sig_gen = momntm_mu_bc_frame_sig_gen.E();
    float mu_en_jpsi_frame_sig_gen = momntm_mu_jpsi_frame_sig_gen.E();

    float recco_bc_pt_sig_gen = mu1_pt_sig_gen + mu2_pt_sig_gen + mu_pt_sig_gen + taunu_pt_sig_gen + taunu2_pt_sig_gen + munu_pt_sig_gen;
    float mu_px_sig_gen = momntm_mu_sig_gen.Px();
    float mu_py_sig_gen = momntm_mu_sig_gen.Py();
    float mu_pz_sig_gen = momntm_mu_sig_gen.Pz();
    
    //Filling the histograms
    if ((abs(mu_eta_sig_gen) < 1.3 and mu_pt_sig_gen > 3.3) or (1.3 < abs(mu_eta_sig_gen) < 2.2 and mu_pt_sig_gen > 2.9) or (2.2 < abs(mu_eta_sig_gen) < 2.4 and mu_pt_sig_gen > 2.4)){ //acceptance cuts
      if ((Trig_dimuon_sig_gen == false) and (Trig_JpsiTk_sig_gen == true)){ //trigger flags
	hmu_m2_miss_sig_gen->Fill(m2_miss_sig_gen);
	hmu_Q2_sig_gen->Fill(Q2_sig_gen);
	hmu_Pt_miss_sig_gen->Fill(Pt_miss_sig_gen);
	hmu_m2_miss_corr_sig_gen->Fill(m2_miss_corr_sig_gen);
	hmu_Q2_corr_sig_gen->Fill(Q2_corr_sig_gen);
	hmu_Pt_miss_corr_sig_gen->Fill(Pt_miss_corr_sig_gen);
	hmu_Pt_var_sig_gen->Fill(Pt_var_sig_gen);
	hmu_mu_pt_sig_gen->Fill(mu_pt_sig_gen);
	hmu_mu_px_sig_gen->Fill(mu_px_sig_gen);
	hmu_mu_py_sig_gen->Fill(mu_py_sig_gen);
	hmu_mu_pz_sig_gen->Fill(mu_pz_sig_gen);
	hmu_mu_eta_sig_gen->Fill(mu_eta_sig_gen);
	hmu_mu_phi_sig_gen->Fill(mu_phi_sig_gen);
	hmu_del_R_sig_gen->Fill(del_R_sig_gen);
	hmu_mu_en_bc_frame_sig_gen->Fill(mu_en_bc_frame_sig_gen);
	hmu_mu_en_jpsi_frame_sig_gen->Fill(mu_en_jpsi_frame_sig_gen);
	hmu_bc_pt_sig_gen->Fill(bc_pt_sig_gen);
	hmu_bc_pt_corr_sig_gen->Fill(bc_pt_corr_sig_gen);
	hmu_bc_pt_ratio_sig_gen->Fill(bc_pt_ratio_sig_gen);
	hprof_sig_gen->Fill(bc_pt_sig_gen,bc_pt_ratio_sig_gen,1);
	hmu_bc_pt_scatter_sig_gen->Fill(bc_pt_sig_gen,bc_pt_corr_sig_gen);
	hmu_recco_bc_pt_sig_gen->Fill(recco_bc_pt_sig_gen);
      }
    }
  }  
  
  //Drawing Histograms

  
  //Missing mass Squared
  //float norm_gen = hmu_m2_miss_norm_gen->GetEntries();
  hmu_m2_miss_norm_gen->SetLineColor(kRed);
  hmu_m2_miss_norm_gen->GetXaxis()->SetTitle("m^{2}_{miss}[GeV^{2}]");
  hmu_m2_miss_norm_gen->GetYaxis()->SetTitle("a.u.");
  hmu_m2_miss_norm_gen->DrawNormalized();
  hmu_m2_miss_sig_gen->SetLineColor(kBlue);
  hmu_m2_miss_sig_gen->DrawNormalized("same");
  legend->AddEntry(hmu_m2_miss_norm_gen, "#mu channel", "l");
  legend->AddEntry(hmu_m2_miss_sig_gen, "#tau channel", "l");
  legend->Draw();
  c1->SaveAs("M2_miss.png");

  
  //Squared four momentum transfered to lepton system
  hmu_Q2_norm_gen->SetLineColor(kRed);
  hmu_Q2_norm_gen->GetXaxis()->SetTitle("Q^{2}[GeV^{2}]");
  hmu_Q2_norm_gen->GetYaxis()->SetTitle("a.u.");
  hmu_Q2_norm_gen->DrawNormalized();
  hmu_Q2_sig_gen->SetLineColor(kBlue);
  hmu_Q2_sig_gen->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("Q2.png");
  


  //Missing transverse momentum
  hmu_Pt_miss_norm_gen->SetLineColor(kRed);
  hmu_Pt_miss_norm_gen->GetXaxis()->SetTitle("P_{T}^{miss}[GeV]");
  hmu_Pt_miss_norm_gen->GetYaxis()->SetTitle("a.u.");
  hmu_Pt_miss_norm_gen->DrawNormalized();
  hmu_Pt_miss_sig_gen->SetLineColor(kBlue);
  hmu_Pt_miss_sig_gen->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("Pt_miss.png");


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
  c1->SaveAs("Pt_var.png");
  
  //Delta R
  hmu_del_R_norm_gen->SetLineColor(kRed);
  hmu_del_R_norm_gen->GetXaxis()->SetTitle("#Delta R(#mu_{1},#mu_{2})");
  hmu_del_R_norm_gen->GetYaxis()->SetTitle("a.u.");
  hmu_del_R_norm_gen->DrawNormalized();
  hmu_del_R_sig_gen->SetLineColor(kBlue);
  hmu_del_R_sig_gen->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("del_R.png");


  //Unpaired muon tranverse momentum
  hmu_mu_pt_norm_gen->SetLineColor(kRed);
  hmu_mu_pt_norm_gen->GetXaxis()->SetTitle("P_{T}^{#mu}[GeV]");
  hmu_mu_pt_norm_gen->GetYaxis()->SetTitle("a.u.");
  hmu_mu_pt_norm_gen->DrawNormalized();
  hmu_mu_pt_sig_gen->SetLineColor(kBlue);
  hmu_mu_pt_sig_gen->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("mu_pt.png");


  //Unpaired muon momentum in x-direction
  hmu_mu_px_norm_gen->SetLineColor(kRed);
  hmu_mu_px_norm_gen->GetXaxis()->SetTitle("P_{x}^{#mu}[GeV]");
  hmu_mu_px_norm_gen->GetYaxis()->SetTitle("a.u.");
  hmu_mu_px_norm_gen->DrawNormalized();
  hmu_mu_px_sig_gen->SetLineColor(kBlue);
  hmu_mu_px_sig_gen->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("mu_px.png");


  //Unpaired muon momentum in y-direction
  hmu_mu_py_norm_gen->SetLineColor(kRed);
  hmu_mu_py_norm_gen->GetXaxis()->SetTitle("P_{y}^{#mu}[GeV]");
  hmu_mu_py_norm_gen->GetYaxis()->SetTitle("a.u.");
  hmu_mu_py_norm_gen->DrawNormalized();
  hmu_mu_py_sig_gen->SetLineColor(kBlue);
  hmu_mu_py_sig_gen->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("mu_py.png");


  //Unpaired muon momentum in z-direction
  hmu_mu_pz_norm_gen->SetLineColor(kRed);
  hmu_mu_pz_norm_gen->GetXaxis()->SetTitle("P_{z}^{#mu}[GeV]");
  hmu_mu_pz_norm_gen->GetYaxis()->SetTitle("a.u.");
  hmu_mu_pz_norm_gen->DrawNormalized();
  hmu_mu_pz_sig_gen->SetLineColor(kBlue);
  hmu_mu_pz_sig_gen->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("mu_pz.png");

  
  //Unpaired muon psuedorapidity
  hmu_mu_eta_norm_gen->SetLineColor(kRed);
  hmu_mu_eta_norm_gen->GetXaxis()->SetTitle("#eta^{#mu}");
  hmu_mu_eta_norm_gen->GetYaxis()->SetTitle("a.u.");
  hmu_mu_eta_norm_gen->DrawNormalized();
  hmu_mu_eta_sig_gen->SetLineColor(kBlue);
  hmu_mu_eta_sig_gen->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("mu_eta.png");

  
  //Unpaired muon phi
  hmu_mu_phi_norm_gen->SetLineColor(kRed);
  hmu_mu_phi_norm_gen->GetXaxis()->SetTitle("#phi^{#mu}");
  hmu_mu_phi_norm_gen->GetYaxis()->SetTitle("a.u.");
  hmu_mu_phi_norm_gen->DrawNormalized();
  hmu_mu_phi_sig_gen->SetLineColor(kBlue);
  hmu_mu_phi_sig_gen->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("mu_phi.png");

  
  //Energy of unpaired muon in Bc rest frame
  hmu_mu_en_bc_frame_norm_gen->SetLineColor(kRed);
  hmu_mu_en_bc_frame_norm_gen->GetXaxis()->SetTitle("E(#mu) in rest frame of B_{c} (GeV)");
  hmu_mu_en_bc_frame_norm_gen->GetYaxis()->SetTitle("a.u.");
  hmu_mu_en_bc_frame_norm_gen->DrawNormalized();
  hmu_mu_en_bc_frame_sig_gen->SetLineColor(kBlue);
  hmu_mu_en_bc_frame_sig_gen->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("E_bc.png");

  //Energy of unpaired muon in Jpsi rest frame
  hmu_mu_en_jpsi_frame_norm_gen->SetLineColor(kRed);
  hmu_mu_en_jpsi_frame_norm_gen->GetXaxis()->SetTitle("E(#mu) in rest frame of J/(#Psi) (GeV)");
  hmu_mu_en_jpsi_frame_norm_gen->GetYaxis()->SetTitle("a.u.");
  hmu_mu_en_jpsi_frame_norm_gen->DrawNormalized();
  hmu_mu_en_jpsi_frame_sig_gen->SetLineColor(kBlue);
  hmu_mu_en_jpsi_frame_sig_gen->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("E_jpsi.png");\

  
  //Transverse momentum of Bc
  hmu_bc_pt_norm_gen->SetLineColor(kRed);
  hmu_bc_pt_norm_gen->GetXaxis()->SetTitle("P_{T}^{Bc} (GeV)");
  hmu_bc_pt_norm_gen->GetYaxis()->SetTitle("a.u.");
  hmu_bc_pt_norm_gen->DrawNormalized();
  hmu_bc_pt_sig_gen->SetLineColor(kBlue);
  hmu_bc_pt_sig_gen->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("Pt_Bc.png");

  
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
  hmu_bc_pt_ratio_norm_gen->GetXaxis()->SetTitle("Ratio of true and corrected P_{T}^{Bc} (GeV)");
  hmu_bc_pt_ratio_norm_gen->GetYaxis()->SetTitle("a.u.");
  hmu_bc_pt_ratio_norm_gen->DrawNormalized();
  hmu_bc_pt_ratio_sig_gen->SetLineColor(kBlue);
  hmu_bc_pt_ratio_sig_gen->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("Ratio_pt_bc.png");
  

  //Scatter plot of Transverse momentum of Bc
  hmu_bc_pt_scatter_norm_gen->SetMarkerColor(kRed);
  hmu_bc_pt_scatter_norm_gen->GetXaxis()->SetTitle("P_{T}^{Bc} (GeV)");
  hmu_bc_pt_scatter_norm_gen->GetYaxis()->SetTitle("Corrected P_{T}^{Bc} (GeV)");
  hmu_bc_pt_scatter_norm_gen->DrawNormalized();
  hmu_bc_pt_scatter_sig_gen->SetMarkerColor(kBlue);
  hmu_bc_pt_scatter_sig_gen->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("scatter_Pt_Bc.png");


  //Profile histogram of transverse momentum of Bc
  hprof_norm_gen->SetLineColor(kRed);
  hprof_norm_gen->GetXaxis()->SetTitle("P_{T}^{Bc} (GeV)");
  hprof_norm_gen->GetYaxis()->SetTitle("Ratio of true and corrected P_{T}^{Bc}");
  hprof_norm_gen->Draw();
  hprof_sig_gen->SetLineColor(kBlue);
  hprof_sig_gen->Draw("same");
  legend->Draw();
  c1->SaveAs("Profile_pt_bc.png");
  
  
  //Comparison of Bc Pt and corrected Bc Pt
  hmu_bc_pt_norm_gen->SetLineColor(kRed);
  hmu_bc_pt_norm_gen->GetXaxis()->SetTitle("Corrected P_{T}^{Bc} (GeV)");
  hmu_bc_pt_norm_gen->GetYaxis()->SetTitle("a.u.");
  hmu_bc_pt_norm_gen->DrawNormalized();
  hmu_bc_pt_sig_gen->SetLineColor(kBlue);
  hmu_bc_pt_sig_gen->DrawNormalized("same");
  hmu_bc_pt_corr_norm_gen->SetLineColor(kGreen);
  hmu_bc_pt_corr_norm_gen->DrawNormalized("same");
  hmu_bc_pt_corr_sig_gen->SetLineColor(kYellow);
  hmu_bc_pt_corr_sig_gen->DrawNormalized("same");
  legend->AddEntry(hmu_bc_pt_corr_norm_gen, "#mu channel_recco", "l");
  legend->AddEntry(hmu_bc_pt_corr_sig_gen, "#tau channel_recco", "l");
  legend->Draw();
  c1->SaveAs("Compare_pt_bc.png");


  //Comparison of missing mass square and corrected missing mass square
  hmu_m2_miss_norm_gen->SetLineColor(kRed);
  hmu_m2_miss_norm_gen->GetXaxis()->SetTitle("Corrected m_{miss}^{2} (GeV^{2})");
  hmu_m2_miss_norm_gen->GetYaxis()->SetTitle("a.u.");
  hmu_m2_miss_norm_gen->DrawNormalized();
  hmu_m2_miss_sig_gen->SetLineColor(kBlue);
  hmu_m2_miss_sig_gen->DrawNormalized("same");
  hmu_m2_miss_corr_norm_gen->SetLineColor(kGreen);
  hmu_m2_miss_corr_norm_gen->DrawNormalized("same");
  hmu_m2_miss_corr_sig_gen->SetLineColor(kYellow);
  hmu_m2_miss_corr_sig_gen->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("Compare_m2_miss.png");


  //Comparison of Squared four momentum transfered to lepton system and corrected Squared four momentum transfered to lepton system
  hmu_Q2_norm_gen->SetLineColor(kRed);
  hmu_Q2_norm_gen->GetXaxis()->SetTitle("Corrected Q^{2} (GeV^{2})");
  hmu_Q2_norm_gen->GetYaxis()->SetTitle("a.u.");
  hmu_Q2_norm_gen->DrawNormalized();
  hmu_Q2_sig_gen->SetLineColor(kBlue);
  hmu_Q2_sig_gen->DrawNormalized("same");
  hmu_Q2_corr_norm_gen->SetLineColor(kGreen);
  hmu_Q2_corr_norm_gen->DrawNormalized("same");
  hmu_Q2_corr_sig_gen->SetLineColor(kYellow);
  hmu_Q2_corr_sig_gen->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("Compare_Q2.png");


  //Comparison of Missing transverse momentum and corrected Missing transverse momentum
  hmu_Pt_miss_norm_gen->SetLineColor(kRed);
  hmu_Pt_miss_norm_gen->GetXaxis()->SetTitle("Corrected P_{T}^{miss} (GeV)");
  hmu_Pt_miss_norm_gen->GetYaxis()->SetTitle("a.u.");
  hmu_Pt_miss_norm_gen->DrawNormalized();
  hmu_Pt_miss_sig_gen->SetLineColor(kBlue);
  hmu_Pt_miss_sig_gen->DrawNormalized("same");
  hmu_Pt_miss_corr_norm_gen->SetLineColor(kGreen);
  hmu_Pt_miss_corr_norm_gen->DrawNormalized("same");
  hmu_Pt_miss_corr_sig_gen->SetLineColor(kYellow);
  hmu_Pt_miss_corr_sig_gen->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("Compare_Pt_miss.png");


  //Comparison of Bc Pt and Reconstructed Bc Pt
  hmu_bc_pt_norm_gen->SetLineColor(kRed);
  hmu_bc_pt_norm_gen->GetXaxis()->SetTitle("P_{T}^{Bc} (GeV)");
  hmu_bc_pt_norm_gen->GetYaxis()->SetTitle("a.u.");
  hmu_bc_pt_norm_gen->DrawNormalized();
  hmu_bc_pt_sig_gen->SetLineColor(kBlue);
  hmu_bc_pt_sig_gen->DrawNormalized("same");
  hmu_recco_bc_pt_norm_gen->SetLineColor(kGreen);
  hmu_recco_bc_pt_norm_gen->DrawNormalized("same");
  hmu_recco_bc_pt_sig_gen->SetLineColor(kYellow);
  hmu_recco_bc_pt_sig_gen->DrawNormalized("same");
  legend->Draw();
  c1->SaveAs("Compare_recco_pt_bc.png");


  //Closing the files
  f1->Close();
  f2->Close();
  
}
