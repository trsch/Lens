#include <iostream>
#include <math.h>
#include "TTree.h"
#include "TFile.h"
#include "TF1.h"
#include "TF2.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TLine.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TSCBasicInlineFunctions.h"
#include "TSCMathFunctions.h"
#include "Analyser.h"
#include "TRandom2.h"
#include <string>
using namespace std;

void Looper(char *filename="45mmWater.root"){
	TFile *f=new TFile(filename);
	TTree *T;
	f->GetObject("T",T);
	TreeLooper *O=new TreeLooper(T);
	O->Loop();
}

/*
int main(){
	Looper();
	return 0;
}
*/

	/*
	sprintf(chnlname[0],"chnl 0:  2D SANS detector");
	sprintf(chnlname[1],"chnl 1:  2D SANS detector");
	sprintf(chnlname[2],"chnl 2:  2D SANS detector");
	sprintf(chnlname[3],"chnl 3:  2D SANS detector");
	sprintf(chnlname[4],"chnl 4:  2D SANS detector");
	sprintf(chnlname[5],"Chnl 5: no connection");
	sprintf(chnlname[6],"Chnl 6: no connection");
	sprintf(chnlname[7],"chnl 7:  SANS BKGD pencil detector");
	sprintf(chnlname[8],"chnl 8: SANS beam monitor");
	sprintf(chnlname[9],"chnl 9: 40 Hz clock from timing system (elapsed number of frames)");
	sprintf(chnlname[10],"chnl 10: IUPTC BKGD pencil detector");
	sprintf(chnlname[11],"chnl 11: decimated Q-pulses (Q/10)");
	sprintf(chnlname[12],"chnl 12: 3He (~10pct efficiencty) spectrum detector");
	sprintf(chnlname[13],"chnl 13: GS-20 scintillation detector for emission time spectrometer");
	sprintf(chnlname[14],"chnl 14 RERP-2 BKGD pencil detector");
	sprintf(chnlname[15],"chnl 15 To from Q box (frames with beam present)");
	sprintf(chnlname[16],"The BNL105 Detector");
	sprintf(chnlname[17],"SANS Overflow chnls");
	sprintf(chnlname[18],"Total data counts");
	sprintf(chnlname[19],"this one should b empty!");
	*/

void TreeLooper::Loop(){
	//****************************************************************************
	//*************************Initializing***************************************
	//****************************************************************************

	TFile *SaveFile=new TFile("output.root","RECREATE");
	
	TH2D *BNL_2dPxl_uncut;
	BNL_2dPxl_uncut=new TH2D("BNL_2dPxl_uncut","",256,0,255,256,0,255);
	TSC::NameHist(BNL_2dPxl_uncut,1,"","Pxl x-index [#]","Pxl y-index [#]");
	BNL_2dPxl_uncut->Sumw2();

	TH3D *BNL_3dPxlTime_uncut;
	BNL_3dPxlTime_uncut=new TH3D("BNL_3dPxlTime_uncut","",256,0,255,256,0,255,250,0,25000);
	TSC::NameHist(BNL_3dPxlTime_uncut,1,"","Pxl x-index [#]","Pxl y-index [#]","Time [ns]");
	BNL_3dPxlTime_uncut->Sumw2();

	TH1D *BNL_Time_uncut;
	BNL_Time_uncut=new TH1D("BNL_time_uncut","",25000,0,50000);
	TSC::NameHist(BNL_Time_uncut,1,"","Pxl x-index [#]","Pxl y-index [#]");
	BNL_Time_uncut->Sumw2();

/*	
	"chnl 0:  2D SANS detector");
	"chnl 1:  2D SANS detector");
	"chnl 2:  2D SANS detector");
	"chnl 3:  2D SANS detector");
	"chnl 4:  2D SANS detector");
	"Chnl 5: no connection");
	"Chnl 6: no connection");
	"chnl 7:  SANS BKGD pencil detector");
	"chnl 8: SANS beam monitor");
	"chnl 9: 40 Hz clock from timing system (elapsed number of frames)");
	"chnl 10: IUPTC BKGD pencil detector");
	"chnl 11: decimated Q-pulses (Q/10)");
	"chnl 12: 3He (~10pct efficiencty) spectrum detector");
	"chnl 13: GS-20 scintillation detector for emission time spectrometer");
	"chnl 14 RERP-2 BKGD pencil detector");
	"chnl 15 To from Q box (frames with beam present)");
	"The BNL105 Detector");
	"SANS Overflow chnls");
	"Total data counts");
	"this one should b empty!");
	*/

	int itr; // iterator

	if (fChain == 0) return;
	cout<<"Starting loop"<<endl;
	Long64_t nentries = fChain->GetEntriesFast();
	cout<<"Looping over: "<<nentries<<" events"<<endl;
	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		if(jentry%10000==0)cout<<jentry<<"  "<<TSC::TheMatrix(60)<<endl;
		
		//**********************************************************************
		//*****************************BNL loop*********************************
		//**********************************************************************
		for(itr=0;itr<nBNL;itr++){
			BNL_2dPxl_uncut->Fill(BNLx[itr],BNLy[itr]);
			BNL_3dPxlTime_uncut->Fill(BNLx[itr],BNLy[itr],BNL[itr]);
			BNL_Time_uncut->Fill(BNL[itr]);
		}

	}
	cout<<"Loop done"<<endl;

	TCanvas *c1=new TCanvas(TSC::uname());
	BNL_2dPxl_uncut->Draw("colz");

	TCanvas *c2=new TCanvas(TSC::uname());
	BNL_Time_uncut->Draw();

	SaveFile->Write();
	//****************************************************************************
	//*************************Ending*********************************************
	//****************************************************************************
	cout<<"Drawing phase done"<<endl;
}
