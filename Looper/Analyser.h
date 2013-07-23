#pragma once
#include <iostream>
#include "TTree.h"
#include <stdio.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
using namespace std;
// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class TreeLooper {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   UInt_t          SANSFrameIndex;
   Int_t           nC0;
   Int_t           nC1;
   Int_t           nC2;
   Int_t           nC3;
   Int_t           nC4;
   Int_t           nC5;
   Int_t           nC6;
   Int_t           nC7;
   Int_t           nC8;
   Int_t           nC9;
   Int_t           nC10;
   Int_t           nC11;
   Int_t           nC12;
   Int_t           nC13;
   Int_t           nC14;
   Int_t           nC15;
   Int_t           nNot6417;
   Int_t           nTotal;
   Float_t         C0[64];   //[nHitsChnl0]
   Float_t         C1[64];   //[nHitsChnl1]
   Float_t         C2[64];   //[nHitsChnl2]
   Float_t         C3[64];   //[nHitsChnl3]
   Float_t         C4[64];   //[nHitsChnl4]
   Float_t         C5[8];   //[nHitsChnl5]
   Float_t         C6[8];   //[nHitsChnl6]
   Float_t         C7[8];   //[nHitsChnl7]
   Float_t         C8[64];   //[nHitsChnl8]
   Float_t         C9[256];   //[nHitsChnl9]
   Float_t         C10[16];   //[nHitsChnl10]
   Float_t         C11[16];   //[nHitsChnl11]
   Float_t         C12[64];   //[nHitsChnl12]
   Float_t         C13[16];   //[nHitsChnl13]
   Float_t         C14[16];   //[nHitsChnl14]
   Float_t         C15[16];   //[nHitsChnl15]
   Float_t         CNot64[256];   //[nHitsChnl17]
   Int_t           CNot64Pxl1[256];   //[nHitsChnl17]
   Int_t           CNot64Pxl2[256];   //[nHitsChnl17]
   UInt_t          BNL105FrameIndex;
   Int_t           nBNL;
   Float_t         BNL[512];   //[nHitsChnl16]
   Float_t         BNLx[512];   //[nHitsChnl16]
   Float_t         BNLy[512];   //[nHitsChnl16]
   UInt_t          PulseShapeFrameIndex;
   UInt_t          PulseID;
   Float_t         Integral;
   Int_t           FWHM;
   Int_t           MaximumIndex;
   Int_t           RiseStart;
   Float_t         Maximum;
   Int_t           PulseVectorLength;
   Float_t         PulseVector[1024];   //[PulseVectorLength]
   Int_t           PulseVectorIndex[1024];   //[PulseVectorLength]
   Double_t        AproxTime;

   // List of branches
   TBranch        *b_SANSFrameIndex;   //!
   TBranch        *b_nHitsChnl0;   //!
   TBranch        *b_nHitsChnl1;   //!
   TBranch        *b_nHitsChnl2;   //!
   TBranch        *b_nHitsChnl3;   //!
   TBranch        *b_nHitsChnl4;   //!
   TBranch        *b_nHitsChnl5;   //!
   TBranch        *b_nHitsChnl6;   //!
   TBranch        *b_nHitsChnl7;   //!
   TBranch        *b_nHitsChnl8;   //!
   TBranch        *b_nHitsChnl9;   //!
   TBranch        *b_nHitsChnl10;   //!
   TBranch        *b_nHitsChnl11;   //!
   TBranch        *b_nHitsChnl12;   //!
   TBranch        *b_nHitsChnl13;   //!
   TBranch        *b_nHitsChnl14;   //!
   TBranch        *b_nHitsChnl15;   //!
   TBranch        *b_nHitsChnl17;   //!
   TBranch        *b_nHitsChnl18;   //!
   TBranch        *b_C0;   //!
   TBranch        *b_C1;   //!
   TBranch        *b_C2;   //!
   TBranch        *b_C3;   //!
   TBranch        *b_C4;   //!
   TBranch        *b_C5;   //!
   TBranch        *b_C6;   //!
   TBranch        *b_C7;   //!
   TBranch        *b_C8;   //!
   TBranch        *b_C9;   //!
   TBranch        *b_C10;   //!
   TBranch        *b_C11;   //!
   TBranch        *b_C12;   //!
   TBranch        *b_C13;   //!
   TBranch        *b_C14;   //!
   TBranch        *b_C15;   //!
   TBranch        *b_CNot64;   //!
   TBranch        *b_CNot64Pxl1;   //!
   TBranch        *b_CNot64Pxl2;   //!
   TBranch        *b_BNL105FrameIndex;   //!
   TBranch        *b_nHitsChnl16;   //!
   TBranch        *b_BNL;   //!
   TBranch        *b_BNLx;   //!
   TBranch        *b_BNLy;   //!
   TBranch        *b_PulseShapeFrameIndex;   //!
   TBranch        *b_PulseID;   //!
   TBranch        *b_Integral;   //!
   TBranch        *b_FWHM;   //!
   TBranch        *b_MaximumIndex;   //!
   TBranch        *b_RiseStart;   //!
   TBranch        *b_Maximum;   //!
   TBranch        *b_PulseVectorLength;   //!
   TBranch        *b_PulseVector;   //!
   TBranch        *b_PulseVectorIndex;   //!
   TBranch        *b_TimeStampPlusFrameTime;   //!

   // List of branches
   TBranch        *b_hits;   //!

   TreeLooper(TTree *tree=0);
   virtual ~TreeLooper();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop_RunStatistics(double rotation=0);
   virtual void     Loop_BNL_FlatField(double rotation,double rotxlowcut,double rotxhighcut,double rotylowcut,double rotyhighcut,int *Flucturating_Bin_GoesHigh,double LimitMultiplyer=1.,int NumberOfRunTimeBins=1000);
   virtual void     Loop_BNL_Analysis(double rotation, double xSignalLow,double xSignalHigh,double ySignalLow,double ySignalHigh,char *flatfieldFileName=0);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};


TreeLooper::TreeLooper(TTree *tree) : fChain(0){
	Init(tree);
}
TreeLooper::~TreeLooper()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}
Int_t TreeLooper::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t TreeLooper::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}
void TreeLooper::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);
   if(tree->GetListOfBranches()->FindObject("SANSFrameIndex")==0){
	   cout<<"****************************************************************************"<<endl;
	   cout<<"****************************************************************************"<<endl;
	   cout<<"******************************WARNING:: No SANS data!!!*********************"<<endl;
	   cout<<"****************************************************************************"<<endl;
	   cout<<"****************************************************************************"<<endl;
   } else {
   fChain->SetBranchAddress("SANSFrameIndex", &SANSFrameIndex, &b_SANSFrameIndex);
   fChain->SetBranchAddress("nC0", &nC0, &b_nHitsChnl0);
   fChain->SetBranchAddress("nC1", &nC1, &b_nHitsChnl1);
   fChain->SetBranchAddress("nC2", &nC2, &b_nHitsChnl2);
   fChain->SetBranchAddress("nC3", &nC3, &b_nHitsChnl3);
   fChain->SetBranchAddress("nC4", &nC4, &b_nHitsChnl4);
   fChain->SetBranchAddress("nC5", &nC5, &b_nHitsChnl5);
   fChain->SetBranchAddress("nC6", &nC6, &b_nHitsChnl6);
   fChain->SetBranchAddress("nC7", &nC7, &b_nHitsChnl7);
   fChain->SetBranchAddress("nC8", &nC8, &b_nHitsChnl8);
   fChain->SetBranchAddress("nC9", &nC9, &b_nHitsChnl9);
   fChain->SetBranchAddress("nC10", &nC10, &b_nHitsChnl10);
   fChain->SetBranchAddress("nC11", &nC11, &b_nHitsChnl11);
   fChain->SetBranchAddress("nC12", &nC12, &b_nHitsChnl12);
   fChain->SetBranchAddress("nC13", &nC13, &b_nHitsChnl13);
   fChain->SetBranchAddress("nC14", &nC14, &b_nHitsChnl14);
   fChain->SetBranchAddress("nC15", &nC15, &b_nHitsChnl15);
   fChain->SetBranchAddress("nNot6417", &nNot6417, &b_nHitsChnl17);
   fChain->SetBranchAddress("nTotal", &nTotal, &b_nHitsChnl18);
   fChain->SetBranchAddress("C0", C0, &b_C0);
   fChain->SetBranchAddress("C1", C1, &b_C1);
   fChain->SetBranchAddress("C2", C2, &b_C2);
   fChain->SetBranchAddress("C3", C3, &b_C3);
   fChain->SetBranchAddress("C4", C4, &b_C4);
   fChain->SetBranchAddress("C5", &C5, &b_C5);
   fChain->SetBranchAddress("C6", &C6, &b_C6);
   fChain->SetBranchAddress("C7", &C7, &b_C7);
   fChain->SetBranchAddress("C8", C8, &b_C8);
   fChain->SetBranchAddress("C9", C9, &b_C9);
   fChain->SetBranchAddress("C10", C10, &b_C10);
   fChain->SetBranchAddress("C11", C11, &b_C11);
   fChain->SetBranchAddress("C12", C12, &b_C12);
   fChain->SetBranchAddress("C13", C13, &b_C13);
   fChain->SetBranchAddress("C14", C14, &b_C14);
   fChain->SetBranchAddress("C15", C15, &b_C15);
   fChain->SetBranchAddress("CNot64", CNot64, &b_CNot64);
   fChain->SetBranchAddress("CNot64Pxl1", CNot64Pxl1, &b_CNot64Pxl1);
   fChain->SetBranchAddress("CNot64Pxl2", CNot64Pxl2, &b_CNot64Pxl2);
   }
   if(tree->GetListOfBranches()->FindObject("BNL105FrameIndex")==0){
		cout<<"****************************************************************************"<<endl;
		cout<<"****************************************************************************"<<endl;
		cout<<"****************************WARNING:: No BNL105 data!!!*********************"<<endl;
		cout<<"****************************************************************************"<<endl;
		cout<<"****************************************************************************"<<endl;
   } else {
   fChain->SetBranchAddress("BNL105FrameIndex", &BNL105FrameIndex, &b_BNL105FrameIndex);
   fChain->SetBranchAddress("nBNL", &nBNL, &b_nHitsChnl16);
   fChain->SetBranchAddress("BNL", BNL, &b_BNL);
   fChain->SetBranchAddress("BNLx", BNLx, &b_BNLx);
   fChain->SetBranchAddress("BNLy", BNLy, &b_BNLy);
   }
   if(tree->GetListOfBranches()->FindObject("PulseShapeFrameIndex")==0){
		cout<<"****************************************************************************"<<endl;
		cout<<"****************************************************************************"<<endl;
		cout<<"****************************WARNING:: No Pulse data!!!**********************"<<endl;
		cout<<"****************************************************************************"<<endl;
		cout<<"****************************************************************************"<<endl;
   } else {
   fChain->SetBranchAddress("PulseShapeFrameIndex", &PulseShapeFrameIndex, &b_PulseShapeFrameIndex);
   fChain->SetBranchAddress("PulseID", &PulseID, &b_PulseID);
   fChain->SetBranchAddress("Integral", &Integral, &b_Integral);
   fChain->SetBranchAddress("FWHM", &FWHM, &b_FWHM);
   fChain->SetBranchAddress("MaximumIndex", &MaximumIndex, &b_MaximumIndex);
   fChain->SetBranchAddress("RiseStart", &RiseStart, &b_RiseStart);
   fChain->SetBranchAddress("Maximum", &Maximum, &b_Maximum);
   fChain->SetBranchAddress("PulseVectorLength", &PulseVectorLength, &b_PulseVectorLength);
   fChain->SetBranchAddress("PulseVector", PulseVector, &b_PulseVector);
   fChain->SetBranchAddress("PulseVectorIndex", PulseVectorIndex, &b_PulseVectorIndex);
   fChain->SetBranchAddress("AproxTime", &AproxTime, &b_TimeStampPlusFrameTime);
   }

   Notify();
}
Bool_t TreeLooper::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}
void TreeLooper::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t TreeLooper::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
