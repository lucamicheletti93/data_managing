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

void create_tree_from_tree(int RUN_NUMBER = 246763){
  //============================================================================
  //SET ADDRESS
  //============================================================================
  TStopwatch *clock1 = new TStopwatch();
  clock1 -> Start();
  //char *PATH_IN = "../GRID_FILES/";
  char *PATH_IN = "../../../PbPb_2015_TREE//"; //PATH FOR MAC
  char FILE_NAME_IN[400];
  sprintf(FILE_NAME_IN,"%s/Tree_%i.root",PATH_IN,RUN_NUMBER);
  //sprintf(FILE_NAME_IN,"Tree_%i.root",RUN_NUMBER);
  double PI = TMath::Pi();
  //============================================================================
  //OPENING THE FILE
  //============================================================================
  TChain *chain = new TChain("PbPbTree");
  Long_t *dummy1 = 0, *dummy2 = 0, *dummy3 = 0, *dummy4 = 0;
  if(gSystem -> GetPathInfo(FILE_NAME_IN,dummy1,dummy2,dummy3,dummy4) == 0){
  printf("Opening %s\n",FILE_NAME_IN);
  chain -> Add(FILE_NAME_IN);
  chain -> ls();
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
  // tree with essential info for the unbinned fit
  //============================================================================
  Double_t DimuPt_unb;
  Double_t DimuY_unb;
  Double_t DimuMass_unb;
  Double_t CostHE_unb, PhiHE_unb, CostCS_unb, PhiCS_unb;

  TTree *output_tree = new TTree("output_tree","output_tree");
  output_tree -> Branch("DimuPt_unb",&DimuPt_unb,"DimuPt_unb/D");
  output_tree -> Branch("DimuY_unb",&DimuY_unb,"DimuY_unb/D");
  output_tree -> Branch("DimuMass_unb",&DimuMass_unb,"DimuMass_unb/D");
  output_tree -> Branch("CostHE_unb",&CostHE_unb,"CostHE_unb/D");
  output_tree -> Branch("PhiHE_unb",&PhiHE_unb,"PhiHE_unb/D");
  output_tree -> Branch("CostCS_unb",&CostCS_unb,"CostCS_unb/D");
  output_tree -> Branch("PhiCS_unb",&PhiCS_unb,"PhiCS_unb/D");

  //============================================================================
  //FILLING HISTOS
  //============================================================================

  Int_t NEntries = chain -> GetEntries();
  for(int i = 0;i < NEntries;i++){
    printf("%i -> %i : %3.2f%\r",NEntries,i,(double) i/NEntries*100);
    chain -> GetEntry(i);
    for(int k = 0;k < NDimu;k++){

      if(IsPhysSelected){
        TString Trigger = TrigClass;
        Bool_t TriggerSelected = kFALSE;
        if(Trigger.Contains("CMUL7-B-NOPF-MUFAST")) TriggerSelected = kTRUE;
        if(DimuY[k] > -4. && DimuY[k] < -2.5){
          if(TriggerSelected){
            if(DimuMatch[k] == 2){
              if(DimuMass[k] > 2 && DimuMass[k] < 5){
                if(DimuPt[k] > 2 && DimuPt[k] <= 6){
                  DimuPt_unb = DimuPt[k];
                  DimuY_unb = DimuY[k];
                  DimuMass_unb = DimuMass[k];
                  CostHE_unb = CostHE[k];
                  PhiHE_unb = TMath::Abs(PhiHE[k]);
                  CostCS_unb = CostCS[k];
                  PhiCS_unb = TMath::Abs(PhiCS[k]);
                  output_tree -> Fill();
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
  TStopwatch *clock2 = new TStopwatch();
  clock2 -> Start();

  char *PATH_OUT = "UNBINNED_TREES";
  char FILE_NAME_OUT[400];
  sprintf(FILE_NAME_OUT,"%s/File_for_UnbinnedFit_%i.root",PATH_OUT,RUN_NUMBER);
  TFile *file_out = new TFile(FILE_NAME_OUT,"RECREATE");
  file_out -> cd();
  output_tree -> Write();
  clock2 -> Stop();
  clock2 -> Print();
}
}
