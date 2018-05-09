#if !defined(__CINT__) || defined(__MAKECINT__)
#include <stdio.h>

#include <TCanvas.h>
#include <TTree.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLine.h>
#include <TList.h>
#include <TSystem.h>
#include <TGrid.h>
#include <TString.h>
#include <TStopwatch.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TStopwatch.h>
#endif

void create_TH3_from_tree(int RUN_NUMBER = 246984){
  //============================================================================
  //SET ADDRESS
  //============================================================================
  TStopwatch *clock1 = new TStopwatch();
  clock1 -> Start();
  //char *PATH_IN = "../GRID_FILES/";
  char *PATH_IN = "../../../../../PbPb_2015_TREE/"; //PATH FOR MAC
  char FILE_NAME_IN[400];
  sprintf(FILE_NAME_IN,"%s/Tree_%i.root",PATH_IN,RUN_NUMBER);
  //sprintf(FILE_NAME_IN,"Tree_%i.root",RUN_NUMBER);
  //============================================================================
  //CREATE HISTOGRAMS
  //============================================================================
  double PI = TMath::Pi();
  TH3D *hMassCostPhiHE_0pt2_2m = new TH3D("hMassCostPhiHE_0pt2_2m","hMassCostPhiHE_0pt2_2m",100,-1,1,100,-PI,PI,120,2,5);
  TH3D *hMassCostPhiCS_0pt2_2m = new TH3D("hMassCostPhiCS_0pt2_2m","hMassCostPhiCS_0pt2_2m",100,-1,1,100,-PI,PI,120,2,5);
  TH3D *hMassCostPhiHE_2pt6_2m = new TH3D("hMassCostPhiHE_2pt6_2m","hMassCostPhiHE_2pt6_2m",100,-1,1,100,-PI,PI,120,2,5);
  TH3D *hMassCostPhiCS_2pt6_2m = new TH3D("hMassCostPhiCS_2pt6_2m","hMassCostPhiCS_2pt6_2m",100,-1,1,100,-PI,PI,120,2,5);
  TH3D *hMassCostPhiHE_6pt12_2m = new TH3D("hMassCostPhiHE_6pt12_2m","hMassCostPhiHE_6pt12_2m",100,-1,1,100,-PI,PI,120,2,5);
  TH3D *hMassCostPhiCS_6pt12_2m = new TH3D("hMassCostPhiCS_6pt12_2m","hMassCostPhiCS_6pt12_2m",100,-1,1,100,-PI,PI,120,2,5);
  //============================================================================
  //OPENING THE FILE
  //============================================================================
  TChain *chain = new TChain("PbPbTree");
  Long_t *dummy1 = 0, *dummy2 = 0, *dummy3 = 0, *dummy4 = 0;
  if(gSystem -> GetPathInfo(FILE_NAME_IN,dummy1,dummy2,dummy3,dummy4) == 0){
  printf("Opening %s\n",FILE_NAME_IN);
  chain -> Add(FILE_NAME_IN);
  //============================================================================
  //SET TREE VARIABLES
  //============================================================================
  char TrigClass[200];
  Float_t PercV0M, PercCL0, PercCL1;
  Int_t NMuons, NTracklets, NContributors;
  Double_t Vertex[3];
  Double_t Pt[300],E[300], Px[300], Py[300], Pz[300], Y[300], Eta[300];
  Double_t TrackChi2[300], MatchTrigChi2[300], DCA[300], RAtAbsEnd[300];
  Int_t Charge[300], MatchTrig[300];
  Int_t NDimu;
  Double_t DimuPt[3000], DimuPx[3000], DimuPy[3000],DimuPz[3000], DimuY[3000];
  Double_t DimuMass[3000];
  Int_t DimuCharge[3000], DimuMatch[3000];
  Int_t DimuMu[3000][2];
  Double_t CostHE[3000], PhiHE[3000], CostCS[3000], PhiCS[3000];
  UInt_t inpmask;
  Bool_t IsPhysSelected;

  chain -> SetBranchAddress("FiredTriggerClasses",TrigClass);
  //chain -> SetBranchAddress("inpmask",&inpmask);
  chain -> SetBranchAddress("NMuons",&NMuons);
  //chain -> SetBranchAddress("NContributors",&NContributors);
  //chain -> SetBranchAddress("NTracklets",&NTracklets);
  chain -> SetBranchAddress("Vertex",Vertex);
  chain -> SetBranchAddress("PercentV0M",&PercV0M);
  //chain -> SetBranchAddress("PercentCL0",&PercCL0);
  //chain -> SetBranchAddress("PercentCL1",&PercCL1);
  chain -> SetBranchAddress("Pt",Pt);
  chain -> SetBranchAddress("E",E);
  chain -> SetBranchAddress("Px",Px);
  chain -> SetBranchAddress("Py",Py);
  chain -> SetBranchAddress("Pz",Pz);
  chain -> SetBranchAddress("Y",Y);
  chain -> SetBranchAddress("Eta",Eta);
  chain -> SetBranchAddress("MatchTrig",MatchTrig);
  //chain -> SetBranchAddress("TrackChi2",TrackChi2);
  chain -> SetBranchAddress("MatchTrigChi2",MatchTrigChi2);
  //chain -> SetBranchAddress("DCA",DCA);
  chain -> SetBranchAddress("Charge",Charge);
  chain -> SetBranchAddress("RAtAbsEnd",RAtAbsEnd);
  chain -> SetBranchAddress("NDimu",&NDimu);
  chain -> SetBranchAddress("DimuPt",DimuPt);
  chain -> SetBranchAddress("DimuPx",DimuPx);
  chain -> SetBranchAddress("DimuPy",DimuPy);
  chain -> SetBranchAddress("DimuPz",DimuPz);
  chain -> SetBranchAddress("DimuY",DimuY);
  chain -> SetBranchAddress("DimuMass",DimuMass);
  chain -> SetBranchAddress("DimuCharge",DimuCharge);
  chain -> SetBranchAddress("DimuMatch",DimuMatch);
  chain -> SetBranchAddress("DimuMu",DimuMu);
  chain -> SetBranchAddress("CostHE",CostHE);
  chain -> SetBranchAddress("PhiHE",PhiHE);
  chain -> SetBranchAddress("CostCS",CostCS);
  chain -> SetBranchAddress("PhiCS",PhiCS);
  chain -> SetBranchAddress("IsPhysSelected",&IsPhysSelected);

  //============================================================================
  //SET Cost & Phi RANGES
  //============================================================================
  /*const int N_ang_bins = 10;
  double Cost_max = 2, Cost_min = -2;
  double Cost_bin_width = (Cost_max - Cost_min)/N_ang_bins;
  double Cost_bin_max[N_ang_bins],Cost_bin_min[N_ang_bins],Cost_bin_central[N_ang_bins];
  for(int i = 0;i < N_ang_bins;i++){
    Cost_bin_max[i] = Cost_min + (i + 1)*Cost_bin_width;
    Cost_bin_min[i] = Cost_min + i*Cost_bin_width;
    Cost_bin_central[i] = Cost_bin_min[i] + (Cost_bin_width/2);
    //cout << Cost_bin_max[i] << " " << Cost_bin_min[i] << endl;
  }

  //double Phi_max = TMath::Pi(), Phi_min = -TMath::Pi();
  double Phi_max = 4, Phi_min = -4;
  double Phi_bin_width = (Phi_max - Phi_min)/N_ang_bins;
  double Phi_bin_max[N_ang_bins],Phi_bin_min[N_ang_bins],Phi_bin_central[N_ang_bins];
  for(int i = 0;i < N_ang_bins;i++){
    Phi_bin_max[i] = Phi_min + (i + 1)*Phi_bin_width;
    Phi_bin_min[i] = Phi_min + i*Phi_bin_width;
    Phi_bin_central[i] = Phi_bin_min[i] + (Phi_bin_width/2);
    //cout << Phi_bin_max[i] << " " << Phi_bin_min[i] << endl;
  }

  TH1D *hMassOS_2m_matrix[N_ang_bins][N_ang_bins];
  char hMassOS_2m_matrix_name[20];

  for(int i = 0;i < N_ang_bins;i++){
    for(int j = 0;j < N_ang_bins;j++){
      sprintf(hMassOS_2m_matrix_name,"hMassOS_2m_%i_%i",i,j);
      hMassOS_2m_matrix[i][j] = new TH1D(hMassOS_2m_matrix_name,hMassOS_2m_matrix_name,100,2.,5.);
    }
  }*/
  //============================================================================
  //FILLING HISTOS
  //============================================================================
  Int_t NEntries = chain -> GetEntries();
  for(int i = 0;i < NEntries;i++){
    printf("%i -> %i : %3.2f%\r",NEntries,i,(double) i/NEntries*100);
    chain -> GetEntry(i);
    for(int k = 0;k < NDimu;k++){

      //Double_t Eta_Mu0 = (Double_t) Eta[DimuMu[k][0]];
      //Double_t Eta_Mu1 = (Double_t) Eta[DimuMu[k][1]];
      //Double_t RAbs_Mu0 = (Double_t) RAtAbsEnd[DimuMu[k][0]];
      //Double_t RAbs_Mu1 = (Double_t) RAtAbsEnd[DimuMu[k][1]];

      if(IsPhysSelected){
        TString Trigger = TrigClass;
        Bool_t TriggerSelected = kFALSE;
        if(Trigger.Contains("CMUL7-B-NOPF-MUFAST")) TriggerSelected = kTRUE;
        if(DimuY[k] > -4. && DimuY[k] < -2.5){
          if(TriggerSelected){
            if(DimuMatch[k] == 2){
              if(DimuMass[k] > 2 && DimuMass[k] < 5){
                if(DimuPt[k] > 0 && DimuPt[k] <= 2){
                  hMassCostPhiHE_0pt2_2m -> Fill(CostHE[k],PhiHE[k],DimuMass[k]);
                  hMassCostPhiCS_0pt2_2m -> Fill(CostCS[k],PhiCS[k],DimuMass[k]);
                }
                if(DimuPt[k] > 2 && DimuPt[k] <= 6){
                  hMassCostPhiHE_2pt6_2m -> Fill(CostHE[k],PhiHE[k],DimuMass[k]);
                  hMassCostPhiCS_2pt6_2m -> Fill(CostCS[k],PhiCS[k],DimuMass[k]);
                }
                if(DimuPt[k] > 6 && DimuPt[k] <= 12){
                  hMassCostPhiHE_6pt12_2m -> Fill(CostHE[k],PhiHE[k],DimuMass[k]);
                  hMassCostPhiCS_6pt12_2m -> Fill(CostCS[k],PhiCS[k],DimuMass[k]);
                }
              }
            }
          }
        }
      }
    }
  }
  printf("%s\n");
  clock1 -> Stop();
  clock1 -> Print();

  //============================================================================
  //SAVE FILES
  //============================================================================
  TStopwatch *clock2 = new TStopwatch(); //<-
  clock2 -> Start();

  //char *PATH_OUT = "GRID_FILES/ANGULAR_DISTRIBUTIONS/DIMUON_MASS2_CORRECT";
  char *PATH_OUT = "HISTOS_FOR_SIGNAL_EXTRACTION";
  char FILE_NAME_OUT[400];
  sprintf(FILE_NAME_OUT,"%s/HistosFromTree_%i.root",PATH_OUT,RUN_NUMBER);
  //sprintf(FILE_NAME_OUT,"HistosFromTree_Luca_%i.root",RUN_NUMBER);
  TFile *file_out = new TFile(FILE_NAME_OUT,"RECREATE");
  file_out -> cd();
  hMassCostPhiHE_0pt2_2m -> Write();
  hMassCostPhiCS_0pt2_2m -> Write();
  hMassCostPhiHE_2pt6_2m -> Write();
  hMassCostPhiCS_2pt6_2m -> Write();
  hMassCostPhiHE_6pt12_2m -> Write();
  hMassCostPhiCS_6pt12_2m -> Write();
  clock2 -> Stop();
  clock2 -> Print(); //<-
  /*for(int i = 0;i < N_ang_bins;i++){
    for(int j = 0;j < N_ang_bins;j++){
        hMassOS_2m_matrix[i][j] -> Write();
    }
  }*/
  //============================================================================
  //DRAW ANGULAR DISTRIBUTIONS
  //============================================================================
  /*bool print_ang_distrib = kFALSE;
  if(print_ang_distrib){
    TH2D *h_CostHE = new TH2D("h_CostHE","h_CostHE",600,-2,2,1000,0,300);
    TH2D *h_PhiHE = new TH2D("h_PhiHE","h_PhiHE",600,-4,4,1000,0,300);
    TH2D *h_CostCS = new TH2D("h_CostCS","h_CostCS",600,-2,2,1000,0,300);
    TH2D *h_PhiCS = new TH2D("h_PhiCS","h_PhiCS",600,-4,4,1000,0,300);
    TCanvas *c_ang_distrib = new TCanvas("c_ang_distrib","c_ang_distrib",20,20,600,600);
    c_ang_distrib -> Divide(2,2);

    c_ang_distrib -> cd(1);
    h_CostHE -> Draw();
    hCostHE_2m -> Draw("same");

    c_ang_distrib -> cd(2);
    h_PhiHE -> Draw();
    hPhiHE_2m -> Draw("same");

    c_ang_distrib -> cd(3);
    h_CostCS -> Draw();
    hCostCS_2m -> Draw("same");

    c_ang_distrib -> cd(4);
    h_PhiCS -> Draw();
    hPhiCS_2m -> Draw("same");
  }*/
  }
}
