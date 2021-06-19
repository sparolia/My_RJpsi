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

void PV_reso(){

  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain
  
  //Opening File
  TFile *f1 = TFile::Open("RootupleBcTo3Mu_muonChannel.root"); //opening the normalization root file for nn output
  TTree *tree = (TTree*)f1->Get("rootuple/ntuple"); //reading the tree from the ntuple
 
  
  //Declaring leaf Variables

  vector<double> *allPrimaryVertexX;
  vector<double> *allPrimaryVertexY;
  vector<double> *allPrimaryVertexZ;
  vector<double> *primaryVertexX;
  vector<double> *primaryVertexY;
  vector<double> *primaryVertexZ;
  vector<short> *triggerMatchDimuon0;
  vector<short> *triggerMatchJpsiTk;
  TVector3 *gen_b_vtx;

  //List of Branches
  TBranch *b_allPrimaryVertexX;
  TBranch *b_allPrimaryVertexY;
  TBranch *b_allPrimaryVertexZ;
  TBranch *b_primaryVertexX;
  TBranch *b_primaryVertexY;
  TBranch *b_primaryVertexZ;
  TBranch *b_triggerMatchDimuon0;
  TBranch *b_triggerMatchJpsiTk;
  TBranch *b_gen_b_vtx;

 
  //Initializing values of the variable pointers
  allPrimaryVertexX = 0;
  allPrimaryVertexY = 0;
  allPrimaryVertexZ = 0;
  primaryVertexX = 0;
  primaryVertexY = 0;
  primaryVertexZ = 0;
  triggerMatchDimuon0 = 0;
  triggerMatchJpsiTk = 0;
  gen_b_vtx = 0;

  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  
  //reading variables tree
  fChain->SetBranchAddress("gen_b_vtx",&gen_b_vtx,&b_gen_b_vtx);
  fChain->SetBranchAddress("allPrimaryVertexX",&allPrimaryVertexX,&b_allPrimaryVertexX);
  fChain->SetBranchAddress("allPrimaryVertexY",&allPrimaryVertexY,&b_allPrimaryVertexY);
  fChain->SetBranchAddress("allPrimaryVertexZ",&allPrimaryVertexZ,&b_allPrimaryVertexZ);
  fChain->SetBranchAddress("primaryVertexX",&primaryVertexX,&b_primaryVertexX);
  fChain->SetBranchAddress("primaryVertexY",&primaryVertexY,&b_primaryVertexY);
  fChain->SetBranchAddress("primaryVertexZ",&primaryVertexZ,&b_primaryVertexZ);
  fChain->SetBranchAddress("triggerMatchJpsiTk",&triggerMatchJpsiTk,&b_triggerMatchJpsiTk);
  fChain->SetBranchAddress("triggerMatchDimuon0",&triggerMatchDimuon0,&b_triggerMatchDimuon0);

  //Defining Canvas and setting the legend
  auto c1 = new TCanvas("c", "c", 500,500);
  gStyle->SetOptStat(111111);
  auto* legend = new TLegend(0.0, 0.7, 0.2, 0.9);


  //Declaring Histograms
  
  TH1D* hmu_reso_PV_z_DiMu = new TH1D("hmu_reso_z_PV_DiMu","Resolution of Primary Vertex Z coordinate for Dimuon trigger only",100,-0.04,0.04);
  TH1D* hmu_reso_PV_z_JpsiTk = new TH1D("hmu_reso_PV_z_JpsiTk","Resolution of Primary Vertex Z coordinate for Jpsi Tk trigger only",100,-0.04,0.04);
  TH1D* hmu_reso_PV_z_DiMu_JpsiTk = new TH1D("hmu_reso_PV_z_DiMu_JpsiTk","Resolution of Primary Vertex Z coordinate for Dimuon and Jpsi trigger",100,-0.04,0.04);

  TH1D* hmu_z_recco_true = new TH1D("hmu_z_recco_true","Z distance between recco and true PV ",100,-0.1,0.1);
  TH1D* hmu_z_recco_closest = new TH1D("hmu_z_recco_closest","Z distance between recco and closest PV ",100,-0.1,0.1);
  TH1D* hmu_z_closest_true = new TH1D("hmu_z_closest_true","Z distance between closest and true PV ",100,-0.1,0.1);
  TH1D* hmu_z_selected_true = new TH1D("hmu_z_selected_true","Z distance between recco and true PV when it is closest",100,-0.04,0.04);

  TH1D* hmu_closest_dis = new TH1D("hmu_closest_dis","Distribution of closest PV",100,-15,15);
  TH1D* hmu_recco_dis = new TH1D("hmu_recco_dis","Distribution of Reconstructed PV",100,-15,15);
  TH1D* hmu_true_dis = new TH1D("hmu_true_dis","Distribution of True PV",100,-15,15);
  
  double reso_PV_z = 0, true_PV_z = 0, recco_PV_z = 0;
  int count = 0,count_pass = 0;

  //Starting the loop over all events
  for(int i=0; i<tree->GetEntries();i++){
    tree->GetEntry(i);

    if((triggerMatchJpsiTk->at(0)) or (triggerMatchDimuon0->at(0))){ //Checking the trigger condition and reading the variables of interest
      count += 1;
      true_PV_z = gen_b_vtx->z();
      recco_PV_z = primaryVertexZ->at(0);
      reso_PV_z = true_PV_z - recco_PV_z;
    
      if((!triggerMatchJpsiTk->at(0)) and (triggerMatchDimuon0->at(0))){ //Only DiMuon trigger
	hmu_reso_PV_z_DiMu->Fill(reso_PV_z);
      }

      if((triggerMatchJpsiTk->at(0)) and (!triggerMatchDimuon0->at(0))){ //Only JpsiTk trigger
	hmu_reso_PV_z_JpsiTk->Fill(reso_PV_z);
      }
      
      if((triggerMatchJpsiTk->at(0)) and (triggerMatchDimuon0->at(0))){ //Dimuon AND Jpsi trigger
	hmu_reso_PV_z_DiMu_JpsiTk->Fill(reso_PV_z);
      }
      
      double closest_PV_z = 0.0, eps = 100.0;
    
      for(int j=0; j< allPrimaryVertexZ->size(); j++){ //Loop over all primary vertices in a single event
	if(abs(allPrimaryVertexZ->at(j)- true_PV_z) < eps){
	  eps = abs(allPrimaryVertexZ->at(j)- true_PV_z);
	  closest_PV_z = allPrimaryVertexZ->at(j); //Selecting the closest PV to true PV
	}
      }
      if(closest_PV_z - recco_PV_z == 0){ //Checking if the recco PV is the closest one.
	count_pass +=1;
	hmu_z_selected_true->Fill(true_PV_z - recco_PV_z);
      }
      else
	{
	  hmu_z_recco_closest->Fill(closest_PV_z - recco_PV_z);
	  hmu_z_recco_true->Fill(true_PV_z - recco_PV_z);
	  hmu_z_closest_true->Fill(closest_PV_z - true_PV_z);
	}
      hmu_closest_dis->Fill(closest_PV_z);
      hmu_recco_dis->Fill(recco_PV_z);
      hmu_true_dis->Fill(true_PV_z);
    }
  }
  cout<<tree->GetEntries()<<endl;
  cout<<count<<endl;
  cout<<count_pass<<endl;
  
  //Resolution of Primary Vertex Z for Dimuon Trigger only
  hmu_reso_PV_z_DiMu->SetLineColor(kRed);
  hmu_reso_PV_z_DiMu->GetXaxis()->SetTitle("Resolution in Z (cm) of PV");
  hmu_reso_PV_z_DiMu->GetYaxis()->SetTitle("Events");
  hmu_reso_PV_z_DiMu->Draw();
  c1->SetLogy();
  c1->SaveAs("reso_PV_Dimu.png");

  //Resolution of Primary Vertex Z for JpsiTk Trigger only
  hmu_reso_PV_z_JpsiTk->SetLineColor(kRed);
  hmu_reso_PV_z_JpsiTk->GetXaxis()->SetTitle("Resolution in Z (cm) of PV");
  hmu_reso_PV_z_JpsiTk->GetYaxis()->SetTitle("Events");
  hmu_reso_PV_z_JpsiTk->Draw();
  c1->SaveAs("reso_PV_JpsiTk.png");

  //Resolution of Primary Vertex Z for Dimuon and JpsiTk Trigger
  hmu_reso_PV_z_DiMu_JpsiTk->SetLineColor(kRed);
  hmu_reso_PV_z_DiMu_JpsiTk->GetXaxis()->SetTitle("Resolution in Z (cm) of PV");
  hmu_reso_PV_z_DiMu_JpsiTk->GetYaxis()->SetTitle("Events");
  hmu_reso_PV_z_DiMu_JpsiTk->Draw();
  c1->SaveAs("reso_PV_Dimu_JpsiTk.png");

  //Resolution of Z for recco and true Primary Vertex
  hmu_z_recco_true->SetLineColor(kRed);
  hmu_z_recco_true->GetXaxis()->SetTitle("Distance in Z (cm)");
  hmu_z_recco_true->GetYaxis()->SetTitle("Events");
  hmu_z_recco_true->Draw();
  c1->SaveAs("z_recco_true.png");

  //Resolution of Z for recco and closest Primary Vertex
  hmu_z_recco_closest->SetLineColor(kRed);
  hmu_z_recco_closest->GetXaxis()->SetTitle("Distance in Z (cm)");
  hmu_z_recco_closest->GetYaxis()->SetTitle("Events");
  hmu_z_recco_closest->Draw();
  c1->SaveAs("z_recco_closest.png");

  //Resolution of Z for closest and true Primary Vertex
  hmu_z_closest_true->SetLineColor(kRed);
  hmu_z_closest_true->GetXaxis()->SetTitle("Distance in Z (cm)");
  hmu_z_closest_true->GetYaxis()->SetTitle("Events");
  hmu_z_closest_true->Draw();
  c1->SaveAs("z_closest_true.png");

  //Resolution of Z for selected and true Primary Vertex
  hmu_z_selected_true->SetLineColor(kRed);
  hmu_z_selected_true->GetXaxis()->SetTitle("Distance in Z (cm)");
  hmu_z_selected_true->GetYaxis()->SetTitle("Events");
  hmu_z_selected_true->Draw();
  c1->SaveAs("z_selected_true.png");

  //Distribution of Closest Primary Vertex
  hmu_closest_dis->SetLineColor(kRed);
  hmu_closest_dis->GetXaxis()->SetTitle("Distance in Z (cm)");
  hmu_closest_dis->GetYaxis()->SetTitle("Events");
  hmu_closest_dis->Draw();
  c1->SaveAs("closest_dis.png");

  //Distribution of Reconstructed Primary Vertex
  hmu_recco_dis->SetLineColor(kRed);
  hmu_recco_dis->GetXaxis()->SetTitle("Distance in Z (cm)");
  hmu_recco_dis->GetYaxis()->SetTitle("Events");
  hmu_recco_dis->Draw();
  c1->SaveAs("recco_dis.png");

  //Distribution of True Primary Vertex
  hmu_true_dis->SetLineColor(kRed);
  hmu_true_dis->GetXaxis()->SetTitle("Distance in Z (cm)");
  hmu_true_dis->GetYaxis()->SetTitle("Events");
  hmu_true_dis->Draw();
  c1->SaveAs("true_dis.png");
  
  f1->Close();
}
