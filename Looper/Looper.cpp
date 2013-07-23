#include <iostream>
#include <stdio.h>
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
#include "../TSC_Headers/TSCBasicInlineFunctions.h"
#include "../TSC_Headers/TSCMathFunctions.h"
#include "Analyser.h"
#include "TRandom2.h"
#include <string>
using namespace std;

	//Chnl Names global...
	char chnlname[20][150];
void DeclareChnlNames(){
	sprintf(chnlname[0],"0_2D_SANS_detector");
	sprintf(chnlname[1],"1_2D_SANS_detector");
	sprintf(chnlname[2],"2_2D_SANS_detector");
	sprintf(chnlname[3],"3_2D_SANS_detector");
	sprintf(chnlname[4],"4_2D_SANS_detector");
	sprintf(chnlname[5],"5_no_connection");
	sprintf(chnlname[6],"6_no_connection");
	sprintf(chnlname[7],"7_SANS_BKGD_pencil_detector");
	sprintf(chnlname[8],"8_SANS_beam_monitor");
	sprintf(chnlname[9],"9_40_Hz_clock_from_timing_system");
	sprintf(chnlname[10],"10_IUPTC_BKGD_pencil_detector");
	sprintf(chnlname[11],"11_decimated_Q-pulses_");
	sprintf(chnlname[12],"12_3He_10pct_efficiencty_spectrum_detector");
	sprintf(chnlname[13],"13_GS_20_scintillation_detector");
	sprintf(chnlname[14],"14_RERP_2_BKGD_pencil_detector");
	sprintf(chnlname[15],"15_To_from_Q_box");
	sprintf(chnlname[16],"16_The_BNL105_Detector");
	sprintf(chnlname[17],"17_SANS_Overflow_chnls");
	sprintf(chnlname[18],"18_Total_data_counts");
	sprintf(chnlname[19],"19_this_one_should_b_empty!");
}

void TreeLooper::Loop_RunStatistics(double rotation){
	//**********************************************************************
	//*****************************Initializing*****************************
	//**********************************************************************
	cout<<"Initializing..."<<endl;
	
	if (fChain == 0){cout<<"no tree..."<<endl;return;}

	//**********************************************************************
	//*****************************Initial loop - in needed*****************
	//**********************************************************************
	double par[7];
	par[0]=1e10;
	for(int i=1;i<6;i++)par[i]=0;

	cout<<"Initializing loop..."<<endl;
	bool HasPulse=false;
	bool HasBNL=false;
	bool HasSANS=false;
	if(fChain->GetListOfBranches()->FindObject("PulseID"))HasPulse=true;
	else cout<<"!!!!!!!!!!!!!!!!!!!!!!WARNING:: The data sample does not containg Pulse information: Run Time will be ajusted!!!"<<endl;
	if(fChain->GetListOfBranches()->FindObject("BNL105FrameIndex"))HasBNL=true;
	if(fChain->GetListOfBranches()->FindObject("SANSFrameIndex"))HasSANS=true;
	
	fChain->SetBranchStatus("*", 0);
	fChain->SetBranchStatus("Aprox*", 1);
	fChain->SetBranchStatus("*Frame*", 1);
	fChain->SetBranchStatus("PulseID", 1);
	{
	Long64_t nentries = fChain->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;
	cout<<"Preloop - Looping over: "<<nentries<<" events..."<<endl;
	for (Long64_t jentry=0; jentry<nentries;jentry+=5000){
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
			//fixing Aprox time - if needed:
		if(!HasPulse){
			if(HasBNL)AproxTime=BNL105FrameIndex/40.;
			else if(HasSANS)AproxTime=SANSFrameIndex/40.;
			else cout<<"ERRRORORORORORORO..."<<endl;
		}
		if(AproxTime!=0)par[0]=min(par[0],AproxTime);
		par[1]=max(par[1],AproxTime);
		if(HasPulse){
			par[3]=max(par[3],(double)PulseID);
		}
	}
	}
	fChain->SetBranchStatus("*", 1);


	// setting up stuf from the pre-loop
	double First_AproxTime=par[0];
	double Last_AproxTime=par[1];
	double First_PulseID=0;
	double Last_PulseID=par[3];
	double Low_InterestingPulseRange=50;
	double High_InterestingPulseRange=250;


	double BNL_rotation_angle=rotation;

	//**********************************************************************
	//*****************************Declaring output File********************
	//**********************************************************************

	TFile *SaveFile=new TFile("RunStatistics_output.root","RECREATE");


	//**********************************************************************
	//*****************************Declaring Histograms*********************
	//**********************************************************************
	int NumberOfRunTimeBins=1000;

	//**********************************************************************
	//************Histograms for BNL - req: BNL*****************************
	//**********************************************************************
	cout<<"Declaring Histograms..."<<endl;
	TH1D *BNL_xPxl_uncut;
	BNL_xPxl_uncut=new TH1D("BNL_xPxl_uncut","",256,0,256);
	TSC::NameHist(BNL_xPxl_uncut,1,"","Pxl x-index [#]","Entries");
	BNL_xPxl_uncut->Sumw2();

	TH1D *BNL_yPxl_uncut;
	BNL_yPxl_uncut=new TH1D("BNL_yPxl_uncut","",256,0,256);
	TSC::NameHist(BNL_yPxl_uncut,1,"","Pxl y-index [#]","Entries");
	BNL_yPxl_uncut->Sumw2();

	TH1D *BNL_xPxl_rot;
	BNL_xPxl_rot=new TH1D("BNL_xPxl_rot","",256,0,256);
	TSC::NameHist(BNL_xPxl_rot,1,"","Pxl x-index [#]","Entries");
	BNL_xPxl_rot->Sumw2();

	TH1D *BNL_yPxl_rot;
	BNL_yPxl_rot=new TH1D("BNL_yPxl_rot","",256,0,256);
	TSC::NameHist(BNL_yPxl_rot,1,"","Pxl y-index [#]","Entries");
	BNL_yPxl_rot->Sumw2();

	TH2D *BNL_2dPxl_uncut;
	BNL_2dPxl_uncut=new TH2D("BNL_2dPxl_uncut","",256,0,256,256,0,256);
	TSC::NameHist(BNL_2dPxl_uncut,1,"","Pxl x-index [#]","Pxl y-index [#]");
	BNL_2dPxl_uncut->Sumw2();
	TH2D *BNL_2dPxl_rot;
	BNL_2dPxl_rot=new TH2D("BNL_2dPxl_rot","",256,0,256,256,0,256);
	TSC::NameHist(BNL_2dPxl_rot,1,"","Pxl x-index [#]","Pxl y-index [#]");
	BNL_2dPxl_rot->Sumw2();

	TH3D *BNL_3dPxlTime_uncut;
	BNL_3dPxlTime_uncut=new TH3D("BNL_3dPxlTime_uncut","",256,0,256,256,0,256,250,0,25000);
	TSC::NameHist(BNL_3dPxlTime_uncut,1,"","Pxl x-index [#]","Pxl y-index [#]","Time [#mu s]");

	TH3D *BNL_3dPxlTime_rot;
	BNL_3dPxlTime_rot=new TH3D("BNL_3dPxlTime_rot","",256,0,256,256,0,256,250,0,25000);
	TSC::NameHist(BNL_3dPxlTime_rot,1,"","Pxl x-index [#]","Pxl y-index [#]","Time [#mu s]");

	TH1D *BNL_Time_uncut;
	BNL_Time_uncut=new TH1D("BNL_Time_uncut","",25000,0,50000);
	TSC::NameHist(BNL_Time_uncut,1,"","Time [#mu s]","Intensity [PDU]");
	BNL_Time_uncut->Sumw2();

	TH1D *BNL_InstantaniousCountrate_uncut;
	BNL_InstantaniousCountrate_uncut=new TH1D("BNL_InstantaniousCountrate_uncut","",1000,0,25000);
	TSC::NameHist(BNL_InstantaniousCountrate_uncut,1,"","Time [#mu s]","Instantanious Countrate [counts/s/frame]");
	BNL_InstantaniousCountrate_uncut->Sumw2();

	//**********************************************************************
	//************Histograms for P-Pulse - req: Pulse***********************
	//**********************************************************************

	TH2D *Pulse_AllPulse_FullPulse;
	Pulse_AllPulse_FullPulse=new TH2D("Pulse_AllPulse_FullPulse","",(Last_PulseID-First_PulseID+1),First_PulseID,Last_PulseID,1+High_InterestingPulseRange-Low_InterestingPulseRange,Low_InterestingPulseRange,High_InterestingPulseRange);
	TSC::NameHist(Pulse_AllPulse_FullPulse,1,"","Pulse ID [#]","Pulse entry [#]");


	//**********************************************************************
	//************Histograms for RunTimeStatistics - req: BNL,SANS,Pulse****
	//**********************************************************************

	TH2D *RunTimeStatistics_BNL_InstantaniousCountrate_overTime_uncut;
	RunTimeStatistics_BNL_InstantaniousCountrate_overTime_uncut=new TH2D("RunTimeStatistics_BNL_InstantaniousCountrate_overTime_uncut","",NumberOfRunTimeBins,First_AproxTime,Last_AproxTime,1000,0,25000);
	TSC::NameHist(RunTimeStatistics_BNL_InstantaniousCountrate_overTime_uncut,1,"","Run time [s]","Time [#mu s]","Instantanious Countrate [counts/s/frame]");
	RunTimeStatistics_BNL_InstantaniousCountrate_overTime_uncut->Sumw2();

	TH2D *RunTimeStatistics_BNL_Renormalizer_2d;
	RunTimeStatistics_BNL_Renormalizer_2d=new TH2D("RunTimeStatistics_BNL_Renormalizer_2d","",NumberOfRunTimeBins,First_AproxTime,Last_AproxTime,1000,0,25000);
	TSC::NameHist(RunTimeStatistics_BNL_Renormalizer_2d,1,"","Run time [s]","Time [#mu s]","Renomalization value [PDU]");

	TH1D *RunTimeStatistics_BNL_Renormalizer;
	RunTimeStatistics_BNL_Renormalizer=new TH1D("RunTimeStatistics_BNL_Renormalizer","",NumberOfRunTimeBins,First_AproxTime,Last_AproxTime);
	TSC::NameHist(RunTimeStatistics_BNL_Renormalizer,1,"","Run time [s]","Renomalization value [PDU]");

	TH1D *RunTimeStatistics_SANS_Renormalizer;
	RunTimeStatistics_SANS_Renormalizer=new TH1D("RunTimeStatistics_SANS_Renormalizer","",NumberOfRunTimeBins,First_AproxTime,Last_AproxTime);
	TSC::NameHist(RunTimeStatistics_SANS_Renormalizer,1,"","Run time [s]","Renomalization value [PDU]");

	TH1D *RunTimeStatistics_Pulse_Renormalizer;
	RunTimeStatistics_Pulse_Renormalizer=new TH1D("RunTimeStatistics_Pulse_Renormalizer","",NumberOfRunTimeBins,First_AproxTime,Last_AproxTime);
	TSC::NameHist(RunTimeStatistics_Pulse_Renormalizer,1,"","Run time [s]","Renomalization value [PDU]");

	TH1D *RunTimeStatistics_Pulse_Integral;
	RunTimeStatistics_Pulse_Integral=new TH1D("RunTimeStatistics_Pulse_Integral","",NumberOfRunTimeBins,First_AproxTime,Last_AproxTime);
	TSC::NameHist(RunTimeStatistics_Pulse_Integral,1,"","Run time [s]","Integral [PDU]");

	TH1D *RunTimeStatistics_Pulse_FWHM;
	RunTimeStatistics_Pulse_FWHM=new TH1D("RunTimeStatistics_Pulse_FWHM","",NumberOfRunTimeBins,First_AproxTime,Last_AproxTime);
	TSC::NameHist(RunTimeStatistics_Pulse_FWHM,1,"","Run time [s]","FWHM [# /Pulse]");

	TH1D *RunTimeStatistics_Pulse_RiseStart;
	RunTimeStatistics_Pulse_RiseStart=new TH1D("RunTimeStatistics_Pulse_RiseStart","",NumberOfRunTimeBins,First_AproxTime,Last_AproxTime);
	TSC::NameHist(RunTimeStatistics_Pulse_RiseStart,1,"","Run time [s]","Pulse Start [# /Pulse]");

	TH1D *RunTimeStatistics_Pulse_MaximumIndex;
	RunTimeStatistics_Pulse_MaximumIndex=new TH1D("RunTimeStatistics_Pulse_MaximumIndex","",NumberOfRunTimeBins,First_AproxTime,Last_AproxTime);
	TSC::NameHist(RunTimeStatistics_Pulse_MaximumIndex,1,"","Run time [s]","Maximum index [# /Pulse]");

	TH1D *RunTimeStatistics_Pulse_Maximum;
	RunTimeStatistics_Pulse_Maximum=new TH1D("RunTimeStatistics_Pulse_Maximum","",NumberOfRunTimeBins,First_AproxTime,Last_AproxTime);
	TSC::NameHist(RunTimeStatistics_Pulse_Maximum,1,"","Run time [s]","Maximum [# /Pulse]");

	TH1D *RunTimeStatistics_BNL_Activity;
	RunTimeStatistics_BNL_Activity=new TH1D("RunTimeStatistics_BNL_Activity","",NumberOfRunTimeBins,First_AproxTime,Last_AproxTime);
	TSC::NameHist(RunTimeStatistics_BNL_Activity,1,"","Run time [s]","Countrate [counts/s]");

	TH1D *RunTimeStatistics_chnlX_Activity[16];
	for(int i=0;i<16;i++){
		RunTimeStatistics_chnlX_Activity[i]=new TH1D(TSC::CHAR("RunTimeStatistics_chnlX_Activity",chnlname[i]),"",NumberOfRunTimeBins,First_AproxTime,Last_AproxTime);
		TSC::NameHist(RunTimeStatistics_chnlX_Activity[i],1,"","Run time [s]","Countrate [counts/s]");
	}
	//**********************************************************************
	//*****************************Declating counters and others************
	//**********************************************************************

	int Counter_NumberOfActiveBNLFrames=0;
	double BNLx_rot,BNLy_rot;
	fChain->SetBranchStatus("*",1);
	 // switching off some branches to save time...


	//**********************************************************************
	//*****************************Loop phase*******************************
	//**********************************************************************
	cout<<"Initializing loop..."<<endl;
	int itr; // iterator
	Long64_t nentries = 400000;
	//**********************do fast loop? - comment out next line:
	nentries = fChain->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;
	cout<<"Looping over: "<<nentries<<" events..."<<endl;
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		if(jentry%10000==0)cout<<jentry<<"  "<<TSC::TheMatrix(110)<<endl;  // Printing status


		//fixing Aprox time - if needed:
		if(!HasPulse){
			if(HasBNL)AproxTime=BNL105FrameIndex/40.;
			else if(HasSANS)AproxTime=SANSFrameIndex/40.;
		}

		//BNL Analysis:
		if(HasBNL){
			if(nBNL>0)Counter_NumberOfActiveBNLFrames++;
			RunTimeStatistics_BNL_Activity->Fill(AproxTime,nBNL);
			RunTimeStatistics_BNL_Renormalizer->Fill(AproxTime);
			//BNL loop:
			for(itr=0;itr<nBNL;itr++){
				//Rotating BNL:
				BNLx_rot=cos(BNL_rotation_angle*TSC::deg2rad)*BNLx[itr]-sin(BNL_rotation_angle*TSC::deg2rad)*BNLy[itr];
				BNLy_rot=cos(BNL_rotation_angle*TSC::deg2rad)*BNLy[itr]+sin(BNL_rotation_angle*TSC::deg2rad)*BNLx[itr];

				BNL_2dPxl_rot->Fill(BNLx_rot,BNLy_rot);
				BNL_3dPxlTime_rot->Fill(BNLx_rot,BNLy_rot,BNL[itr]);
				BNL_xPxl_rot->Fill(BNLx_rot);
				BNL_yPxl_rot->Fill(BNLy_rot);
				BNL_xPxl_uncut->Fill(BNLx[itr]);
				BNL_yPxl_uncut->Fill(BNLy[itr]);

				RunTimeStatistics_BNL_InstantaniousCountrate_overTime_uncut->Fill(AproxTime,BNL[itr]);
				BNL_2dPxl_uncut->Fill(BNLx[itr],BNLy[itr]);
				BNL_3dPxlTime_uncut->Fill(BNLx[itr],BNLy[itr],BNL[itr]);
				BNL_Time_uncut->Fill(BNL[itr]);
				BNL_InstantaniousCountrate_uncut->Fill(BNL[itr]);
			}
		}

		//Pulse Analysis
		if(HasPulse){
			if(PulseVectorLength>0){
				RunTimeStatistics_Pulse_Renormalizer->Fill(AproxTime);
				RunTimeStatistics_Pulse_Integral->Fill(AproxTime,Integral);
				RunTimeStatistics_Pulse_FWHM->Fill(AproxTime,FWHM);
				RunTimeStatistics_Pulse_RiseStart->Fill(AproxTime,RiseStart);
				RunTimeStatistics_Pulse_MaximumIndex->Fill(AproxTime,MaximumIndex);
				RunTimeStatistics_Pulse_Maximum->Fill(AproxTime,Maximum);
			}
			for(itr=0;itr<PulseVectorLength;itr++){
				Pulse_AllPulse_FullPulse->Fill(PulseID,itr,PulseVector[itr]);
			}
		}
		//SANS Analysis:
		if(HasSANS){
			RunTimeStatistics_SANS_Renormalizer->Fill(AproxTime);
			RunTimeStatistics_chnlX_Activity[0]->Fill(AproxTime,nC0);
			RunTimeStatistics_chnlX_Activity[1]->Fill(AproxTime,nC1);
			RunTimeStatistics_chnlX_Activity[2]->Fill(AproxTime,nC2);
			RunTimeStatistics_chnlX_Activity[3]->Fill(AproxTime,nC3);
			RunTimeStatistics_chnlX_Activity[4]->Fill(AproxTime,nC4);
			RunTimeStatistics_chnlX_Activity[5]->Fill(AproxTime,nC5);
			RunTimeStatistics_chnlX_Activity[6]->Fill(AproxTime,nC6);
			RunTimeStatistics_chnlX_Activity[7]->Fill(AproxTime,nC7);
			RunTimeStatistics_chnlX_Activity[8]->Fill(AproxTime,nC8);
			RunTimeStatistics_chnlX_Activity[9]->Fill(AproxTime,nC9);
			RunTimeStatistics_chnlX_Activity[10]->Fill(AproxTime,nC10);
			RunTimeStatistics_chnlX_Activity[11]->Fill(AproxTime,nC11);
			RunTimeStatistics_chnlX_Activity[12]->Fill(AproxTime,nC12);
			RunTimeStatistics_chnlX_Activity[13]->Fill(AproxTime,nC13);
			RunTimeStatistics_chnlX_Activity[14]->Fill(AproxTime,nC14);
			RunTimeStatistics_chnlX_Activity[15]->Fill(AproxTime,nC15);
		}
	}
	cout<<"Loop done"<<endl;
	//**********************************************************************
	//*****************************Analysis phase***************************
	//**********************************************************************
	cout<<"Analysing..."<<endl;
	if(Counter_NumberOfActiveBNLFrames!=0){
		BNL_InstantaniousCountrate_uncut->Scale(1/(Counter_NumberOfActiveBNLFrames*BNL_InstantaniousCountrate_uncut->GetBinWidth(1)*1e-6));
	}
	for(int i=0;i<RunTimeStatistics_BNL_Renormalizer_2d->GetNbinsX()+1;i++)
		for(int j=0;j<RunTimeStatistics_BNL_Renormalizer_2d->GetNbinsY()+1;j++)
			RunTimeStatistics_BNL_Renormalizer_2d->SetBinContent(i,j,RunTimeStatistics_BNL_Renormalizer->GetBinContent(i));

	RunTimeStatistics_BNL_InstantaniousCountrate_overTime_uncut->Divide(RunTimeStatistics_BNL_Renormalizer_2d);
	RunTimeStatistics_BNL_InstantaniousCountrate_overTime_uncut->Scale(1/(BNL_InstantaniousCountrate_uncut->GetBinWidth(1)*1e-6));

	RunTimeStatistics_Pulse_Integral->Divide(RunTimeStatistics_Pulse_Renormalizer);
	RunTimeStatistics_Pulse_FWHM->Divide(RunTimeStatistics_Pulse_Renormalizer);
	RunTimeStatistics_Pulse_RiseStart->Divide(RunTimeStatistics_Pulse_Renormalizer);
	RunTimeStatistics_Pulse_MaximumIndex->Divide(RunTimeStatistics_Pulse_Renormalizer);
	RunTimeStatistics_Pulse_Maximum->Divide(RunTimeStatistics_Pulse_Renormalizer);

	RunTimeStatistics_BNL_Activity->Scale(1/RunTimeStatistics_BNL_Activity->GetBinWidth(1));//Divide(RunTimeStatistics_BNL_Renormalizer);
	
	for(int i=0;i<16;i++)
		RunTimeStatistics_chnlX_Activity[i]->Scale(1/RunTimeStatistics_chnlX_Activity[i]->GetBinWidth(1));//Divide(RunTimeStatistics_SANS_Renormalizer);
	cout<<"Saving..."<<endl;
	SaveFile->Write();
	cout<<"Looper done"<<endl;
}
void Draw_RunStatistics(char *filename="RunStatistics_output.root"){
	//**********************************************************************
	//*****************************Initializing*****************************
	//**********************************************************************
	cout<<"Loading file and histograms..."<<endl;
	TFile *f=new TFile(filename);
	TH2D *BNL_2dPxl_uncut=(TH2D*)f->Get("BNL_2dPxl_uncut");
	//TH3D *BNL_3dPxlTime_uncut=(TH3D*)f->Get("BNL_3dPxlTime_uncut");

	TH1D *BNL_Time_uncut=(TH1D*)f->Get("BNL_Time_uncut");
	TH1D *BNL_InstantaniousCountrate_uncut=(TH1D*)f->Get("BNL_InstantaniousCountrate_uncut");
	TH2D *Pulse_AllPulse_FullPulse=(TH2D*)f->Get("Pulse_AllPulse_FullPulse");
	TH2D *RunTimeStatistics_BNL_InstantaniousCountrate_overTime_uncut=(TH2D*)f->Get("RunTimeStatistics_BNL_InstantaniousCountrate_overTime_uncut");
	TH2D *RunTimeStatistics_BNL_Renormalizer_2d=(TH2D*)f->Get("RunTimeStatistics_BNL_Renormalizer_2d");
	TH1D *RunTimeStatistics_BNL_Renormalizer=(TH1D*)f->Get("RunTimeStatistics_BNL_Renormalizer");
	TH1D *RunTimeStatistics_SANS_Renormalizer=(TH1D*)f->Get("RunTimeStatistics_SANS_Renormalizer");
	TH1D *RunTimeStatistics_Pulse_Renormalizer=(TH1D*)f->Get("RunTimeStatistics_Pulse_Renormalizer");
	TH1D *RunTimeStatistics_Pulse_Integral=(TH1D*)f->Get("RunTimeStatistics_Pulse_Integral");
	TH1D *RunTimeStatistics_Pulse_FWHM=(TH1D*)f->Get("RunTimeStatistics_Pulse_FWHM");
	TH1D *RunTimeStatistics_Pulse_RiseStart=(TH1D*)f->Get("RunTimeStatistics_Pulse_RiseStart");
	TH1D *RunTimeStatistics_Pulse_MaximumIndex=(TH1D*)f->Get("RunTimeStatistics_Pulse_MaximumIndex");
	TH1D *RunTimeStatistics_Pulse_Maximum=(TH1D*)f->Get("RunTimeStatistics_Pulse_Maximum");
	TH1D *RunTimeStatistics_BNL_Activity=(TH1D*)f->Get("RunTimeStatistics_BNL_Activity");
	TH1D *RunTimeStatistics_chnlX_Activity[16];
	TH2D *BNL_2dPxl_rot=(TH2D*)f->Get("BNL_2dPxl_rot");
	for(int i=0;i<16;i++){
		RunTimeStatistics_chnlX_Activity[i]=(TH1D*)f->Get(TSC::CHAR("RunTimeStatistics_chnlX_Activity",chnlname[i]));
	}


	//**********************************************************************
	//*****************************Drawing phase****************************
	//**********************************************************************
	TCanvas *c1=new TCanvas(TSC::uname());
	TSC::DrawName(BNL_2dPxl_uncut);
	BNL_2dPxl_uncut->Draw("colz");

	TCanvas *c2=new TCanvas(TSC::uname());
	TSC::DrawName(BNL_2dPxl_rot);
	BNL_2dPxl_rot->Draw("colz");

	TCanvas *c3=new TCanvas(TSC::uname());
	TSC::DrawName(	BNL_InstantaniousCountrate_uncut);
	BNL_InstantaniousCountrate_uncut->Draw("colz");

	const int I_RTStats_Max=3;
	int I_RTStats=0;
	TLegend *L_RTStats[I_RTStats_Max];
	TCanvas *C_RTStats=new TCanvas(TSC::uname(),"Run Time Statistics",1400,I_RTStats_Max*300);
	C_RTStats->Divide(1,I_RTStats_Max);
	
	I_RTStats++;
	C_RTStats->cd(I_RTStats);
	RunTimeStatistics_chnlX_Activity[0]->SetTitle("2d SANS Detector and BG detectors");
	RunTimeStatistics_chnlX_Activity[0]->GetYaxis()->SetRangeUser(0,TSC::FindMaximum(RunTimeStatistics_chnlX_Activity[0],RunTimeStatistics_chnlX_Activity[1],RunTimeStatistics_chnlX_Activity[2],RunTimeStatistics_chnlX_Activity[3],RunTimeStatistics_chnlX_Activity[4],RunTimeStatistics_chnlX_Activity[7],/*RunTimeStatistics_chnlX_Activity[10],*/RunTimeStatistics_chnlX_Activity[14]));
	RunTimeStatistics_chnlX_Activity[0]->SetLineColor(2);
	RunTimeStatistics_chnlX_Activity[1]->SetLineColor(4);
	RunTimeStatistics_chnlX_Activity[2]->SetLineColor(6);
	RunTimeStatistics_chnlX_Activity[3]->SetLineColor(8);
	RunTimeStatistics_chnlX_Activity[4]->SetLineColor(1);
	RunTimeStatistics_chnlX_Activity[7]->SetLineColor(7);
	RunTimeStatistics_chnlX_Activity[10]->SetLineColor(5);
	RunTimeStatistics_chnlX_Activity[14]->SetLineColor(3);
	L_RTStats[I_RTStats]=new TLegend(0.7,0.2,0.95,0.60);
	L_RTStats[I_RTStats]->AddEntry(RunTimeStatistics_chnlX_Activity[0],chnlname[0],"l");
	L_RTStats[I_RTStats]->AddEntry(RunTimeStatistics_chnlX_Activity[1],chnlname[1],"l");
	L_RTStats[I_RTStats]->AddEntry(RunTimeStatistics_chnlX_Activity[2],chnlname[2],"l");
	L_RTStats[I_RTStats]->AddEntry(RunTimeStatistics_chnlX_Activity[3],chnlname[3],"l");
	L_RTStats[I_RTStats]->AddEntry(RunTimeStatistics_chnlX_Activity[4],chnlname[4],"l");
	L_RTStats[I_RTStats]->AddEntry(RunTimeStatistics_chnlX_Activity[7],chnlname[7],"l");
	L_RTStats[I_RTStats]->AddEntry(RunTimeStatistics_chnlX_Activity[10],chnlname[10],"l");
	L_RTStats[I_RTStats]->AddEntry(RunTimeStatistics_chnlX_Activity[14],chnlname[14],"l");
	RunTimeStatistics_chnlX_Activity[0]->Draw();
	RunTimeStatistics_chnlX_Activity[1]->Draw("same");
	RunTimeStatistics_chnlX_Activity[2]->Draw("same");
	RunTimeStatistics_chnlX_Activity[3]->Draw("same");
	RunTimeStatistics_chnlX_Activity[4]->Draw("same");
	RunTimeStatistics_chnlX_Activity[7]->Draw("same");
	RunTimeStatistics_chnlX_Activity[10]->Draw("same");
	RunTimeStatistics_chnlX_Activity[14]->Draw("same");
	L_RTStats[I_RTStats]->Draw();
	
	I_RTStats++;
	C_RTStats->cd(I_RTStats);
	RunTimeStatistics_chnlX_Activity[8]->SetTitle("Proton and neutron beam monitors");
	RunTimeStatistics_chnlX_Activity[8]->GetYaxis()->SetRangeUser(0,TSC::FindMaximum(RunTimeStatistics_chnlX_Activity[8],RunTimeStatistics_chnlX_Activity[11],RunTimeStatistics_chnlX_Activity[12]));
	RunTimeStatistics_chnlX_Activity[8]->SetLineColor(2);
	RunTimeStatistics_chnlX_Activity[11]->SetLineColor(4);
	RunTimeStatistics_chnlX_Activity[12]->SetLineColor(6);
	L_RTStats[I_RTStats]=new TLegend(0.7,0.2,0.95,0.60);
	L_RTStats[I_RTStats]->AddEntry(RunTimeStatistics_chnlX_Activity[8],chnlname[8],"l");
	L_RTStats[I_RTStats]->AddEntry(RunTimeStatistics_chnlX_Activity[11],chnlname[11],"l");
	L_RTStats[I_RTStats]->AddEntry(RunTimeStatistics_chnlX_Activity[12],chnlname[12],"l");
	RunTimeStatistics_chnlX_Activity[8]->Draw();
	RunTimeStatistics_chnlX_Activity[11]->Draw("same");
	RunTimeStatistics_chnlX_Activity[12]->Draw("same");
	L_RTStats[I_RTStats]->Draw();
	

	I_RTStats++;
	C_RTStats->cd(I_RTStats);
	double scale=1./RunTimeStatistics_BNL_Activity->GetMaximum()*10;
	RunTimeStatistics_BNL_Activity->Scale(scale);
	RunTimeStatistics_chnlX_Activity[13]->SetTitle("SANS and BNL105 signal");
	RunTimeStatistics_chnlX_Activity[13]->GetYaxis()->SetRangeUser(0,TSC::FindMaximum(RunTimeStatistics_chnlX_Activity[13],RunTimeStatistics_BNL_Activity,RunTimeStatistics_Pulse_Integral));
	RunTimeStatistics_chnlX_Activity[13]->SetLineColor(2);
	RunTimeStatistics_Pulse_Integral->SetLineColor(6);
	RunTimeStatistics_BNL_Activity->SetLineColor(4);
	L_RTStats[I_RTStats]=new TLegend(0.7,0.2,0.95,0.60);
	L_RTStats[I_RTStats]->AddEntry(RunTimeStatistics_chnlX_Activity[13],"Time Focused Spectrometer","l");
	L_RTStats[I_RTStats]->AddEntry(RunTimeStatistics_BNL_Activity,"BNL105 - scaled down","l");
	L_RTStats[I_RTStats]->AddEntry(RunTimeStatistics_Pulse_Integral,"Pulse Integral sum [PDU]","l");
	RunTimeStatistics_chnlX_Activity[13]->Draw();
	RunTimeStatistics_BNL_Activity->DrawClone("same");
	RunTimeStatistics_Pulse_Integral->Draw("same");
	L_RTStats[I_RTStats]->Draw();
	RunTimeStatistics_BNL_Activity->Scale(1./scale);

	//****************************************************************************
	//*************************Ending*********************************************
	//****************************************************************************
	cout<<"Drawing phase done"<<endl;
}
void TreeLooper::Loop_BNL_FlatField(double rotation,double rotxlowcut,double rotxhighcut,double rotylowcut,double rotyhighcut,int *Flucturating_Bin_GoesHigh,double LimitMultiplyer,int NumberOfRunTimeBins){


	//**********************************************************************
	//*****************************Initializing*****************************
	//**********************************************************************
	cout<<"Initializing..."<<endl;
	
	if (fChain == 0){cout<<"no tree..."<<endl;return;}
	// setting initial conditions
	double BNL_rotation_angle=rotation;
	double RunTime_Low=1e10;
	double RunTime_High=0;


	//**********************************************************************
	//*****************************Initial loop - in needed*****************
	//**********************************************************************

	cout<<"Initializing loop..."<<endl;
	bool HasPulse=false;
	bool HasBNL=false;
	bool HasSANS=false;
	if(fChain->GetListOfBranches()->FindObject("PulseID"))HasPulse=true;
	else cout<<"!!!!!!!!!!!!!!!!!!!!!!WARNING:: The data sample does not containg Pulse information: Run Time will be ajusted!!!"<<endl;
	if(fChain->GetListOfBranches()->FindObject("BNL105FrameIndex"))HasBNL=true;
	if(fChain->GetListOfBranches()->FindObject("SANSFrameIndex"))HasSANS=true;






	fChain->SetBranchStatus("*", 0);
	fChain->SetBranchStatus("Aprox*", 1);
	fChain->SetBranchStatus("*Frame*", 1);
	{
	Long64_t nentries = fChain->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;
	cout<<"Preloop (1/2)- Looping over: "<<nentries<<" events..."<<endl;
	for (Long64_t jentry=0; jentry<nentries;jentry+=5000){
		Long64_t ientry = LoadTree(jentry);
		if(jentry%250000==0)cout<<jentry<<endl;  // Printing status
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
			//fixing Aprox time - if needed:
		if(!HasPulse){
			if(HasBNL)AproxTime=BNL105FrameIndex/40.;
			else if(HasSANS)AproxTime=SANSFrameIndex/40.;
			else cout<<"ERRRORORORORORORO..."<<endl;
		}
		if(AproxTime!=0)RunTime_Low=min(RunTime_Low,AproxTime);
		RunTime_High=max(RunTime_High,AproxTime);
	}
	}
	fChain->SetBranchStatus("*", 1);

	RunTime_High+=(RunTime_High-RunTime_Low)*0.025;


	//**********************************************************************
	//*****************************Declaring output File********************
	//**********************************************************************

	TFile *SaveFile=new TFile("FlatField_output.root","RECREATE");



	//**********************************************************************
	//*****************************Declaring Histograms*********************
	//**********************************************************************
	
	cout<<"Declaring Histograms..."<<endl;
	
	TH1D *BNL_HighBin_overTime;
	BNL_HighBin_overTime=new TH1D("BNL_HighBin_overTime","",NumberOfRunTimeBins,RunTime_Low,RunTime_High);
	TSC::NameHist(BNL_HighBin_overTime,1,"","RunTime [s]","Entries");
	BNL_HighBin_overTime->Sumw2();

	TH1D *BNL_LowBin_overTime;
	BNL_LowBin_overTime=new TH1D("BNL_LowBin_overTime","",NumberOfRunTimeBins,RunTime_Low,RunTime_High);
	TSC::NameHist(BNL_LowBin_overTime,1,"","RunTime [s]","Entries");
	BNL_LowBin_overTime->Sumw2();

	TH2D *Flatfield_UsedSurface;
	Flatfield_UsedSurface=new TH2D("Flatfield_UsedSurface","",256,0,256,256,0,256);
	TSC::NameHist(Flatfield_UsedSurface,1,"","Pxl x-index [#]","Pxl y-index [#]");

	TH2D *FlatfieldGenerator_X;
	FlatfieldGenerator_X=new TH2D("FlatfieldGenerator_X","",256,0,256,NumberOfRunTimeBins,RunTime_Low,RunTime_High);
	TSC::NameHist(FlatfieldGenerator_X,1,"","Pxl x-index [#]","RunTime [s]","Weight");
	
	TH2D *FlatfieldGenerator_Y;
	FlatfieldGenerator_Y=new TH2D("FlatfieldGenerator_Y","",256,0,256,NumberOfRunTimeBins,RunTime_Low,RunTime_High);
	TSC::NameHist(FlatfieldGenerator_Y,1,"","Pxl y-index  [#]","RunTime [s]","Weight");

	// Histograms for analysis phase
	TH1D *WhereHigh;
	WhereHigh=new TH1D("WhereHigh","",NumberOfRunTimeBins,RunTime_Low,RunTime_High);
	TSC::NameHist(WhereHigh,1,"","RunTime [s]","High or Low");

	TH1D *FlatfieldHigh_X;
	FlatfieldHigh_X=new TH1D("FlatfieldHigh_X","",256,0,256);
	TSC::NameHist(FlatfieldHigh_X,1,"","Pxl x-index [#]","1/Weight");
	
	TH1D *FlatfieldHigh_Y;
	FlatfieldHigh_Y=new TH1D("FlatfieldHigh_Y","",256,0,256);
	TSC::NameHist(FlatfieldHigh_Y,1,"","Pxl y-index  [#]","1/Weight");

	TH1D *FlatfieldLow_X;
	FlatfieldLow_X=new TH1D("FlatfieldLow_X","",256,0,256);
	TSC::NameHist(FlatfieldLow_X,1,"","Pxl x-index [#]","1/Weight");
	
	TH1D *FlatfieldLow_Y;
	FlatfieldLow_Y=new TH1D("FlatfieldLow_Y","",256,0,256);
	TSC::NameHist(FlatfieldLow_Y,1,"","Pxl y-index  [#]","1/Weight");

	//**********************************************************************
	//*****************************Declating counters and others************
	//**********************************************************************

	int Counter_NumberOfActiveBNLFrames=0;
	double BNLx_rot,BNLy_rot;
	 // switching off some branches to save time...
	fChain->SetBranchStatus("*",0);
	fChain->SetBranchStatus("nBNL",1);
	fChain->SetBranchStatus("BNL",1);
	fChain->SetBranchStatus("BNLx",1);
	fChain->SetBranchStatus("BNLy",1);
	fChain->SetBranchStatus("AproxTime",1);
	fChain->SetBranchStatus("BNL105FrameIndex*",1);
	
	//**********************************************************************
	//***************Loop phase - Fetching RunTime information**************
	//**********************************************************************
	cout<<"Initializing loop 2/2..."<<endl;
	int itr,jtr; // iterator
	int LengthOfHigh=0;
	if(Flucturating_Bin_GoesHigh!=0)LengthOfHigh=sizeof(Flucturating_Bin_GoesHigh)/sizeof(Flucturating_Bin_GoesHigh[0]);
	Long64_t nentries = fChain->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;
	cout<<"Looping over: "<<nentries<<" events..."<<endl;
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		if(jentry%10000==0)cout<<jentry<<"  "<<TSC::TheMatrix(110)<<endl;  // Printing status

		//fixing Aprox time - if needed:
		if(!HasPulse){
			if(HasBNL){AproxTime=BNL105FrameIndex/40.;}
			else {cout<<"ERROR:: no pulse and no BNL!"<<endl;break;}
		}
		for(itr=0;itr<nBNL;itr++){
			//if(LowToFLimit>BNL[itr])continue;
			BNLx_rot=cos(BNL_rotation_angle*TSC::deg2rad)*BNLx[itr]-sin(BNL_rotation_angle*TSC::deg2rad)*BNLy[itr];
			BNLy_rot=cos(BNL_rotation_angle*TSC::deg2rad)*BNLy[itr]+sin(BNL_rotation_angle*TSC::deg2rad)*BNLx[itr];
			if((BNLx_rot>rotxhighcut||BNLx_rot<rotxlowcut)||(BNLy_rot>rotyhighcut||BNLy_rot<rotylowcut))
				Flatfield_UsedSurface->Fill(BNLx[itr],BNLy[itr]);
			if(BNLy_rot>rotyhighcut||BNLy_rot<rotylowcut)
				FlatfieldGenerator_X->Fill(BNLx[itr],AproxTime);
			if(BNLx_rot>rotxhighcut||BNLx_rot<rotxlowcut)
				FlatfieldGenerator_Y->Fill(BNLy[itr],AproxTime);
			if(LengthOfHigh!=0){
				for(jtr=0;jtr<LengthOfHigh;jtr++){
					if(BNLx[itr]>=Flucturating_Bin_GoesHigh[jtr]&&BNLx[itr]<Flucturating_Bin_GoesHigh[jtr]+1){
						BNL_HighBin_overTime->Fill(AproxTime);
					}
				}
			}
			BNL_LowBin_overTime->Fill(AproxTime);
			/*for(jtr=0;jtr<LengthOfLow;jtr++){
				if(BNLx[itr]>=Flucturating_Bin_GoesLow[jtr]&&BNLx[itr]<Flucturating_Bin_GoesLow[jtr]+1){
					BNL_LowBin_overTime->Fill(AproxTime);
				}
			}*/
		}
	}
	cout<<"Loop done"<<endl;
	//**********************************************************************
	//*****************************Analysis phase***************************
	//**********************************************************************
	cout<<"Analysing..."<<endl;
	if(LengthOfHigh!=0){	
		TCanvas *c1=new TCanvas(TSC::uname());
		BNL_HighBin_overTime->SetTitle("Fluctutations over time (High bin)");
		BNL_HighBin_overTime->Draw();
	}
	
	TCanvas *c2=new TCanvas(TSC::uname());
	BNL_LowBin_overTime->SetTitle("Total fluctutations over time");
	BNL_LowBin_overTime->Draw();

	TCanvas *c3=new TCanvas(TSC::uname());
	TH1D *Flucturations=(TH1D*)(BNL_HighBin_overTime->Clone("Flucturations"));
	Flucturations->SetTitle("Fluctutations over time (High/Total)");
	Flucturations->Divide(BNL_LowBin_overTime);
	Flucturations->Draw("HIST");
	
	double HighBin_integral=0;
	double HighBin_errorsum=0;
	double HighBin_activbins=0;
	if(LengthOfHigh!=0){
		for(int i=1;i<BNL_HighBin_overTime->GetNbinsX()+1;i++){
			if(Flucturations->GetBinContent(i)!=0)
				HighBin_activbins++;
				HighBin_errorsum+=Flucturations->GetBinError(i);
				HighBin_integral+=Flucturations->GetBinContent(i);
		}
	}
	double highLowLimit=0;
	if(LengthOfHigh!=0)highLowLimit=HighBin_integral/HighBin_activbins*LimitMultiplyer;
	
	TLine *L=new TLine(0,highLowLimit,1e16,highLowLimit);
	L->SetLineColor(2);
	L->Draw("same");

	TCanvas *c31=new TCanvas(TSC::uname());
	Flatfield_UsedSurface->SetTitle("Hits generating the flatt field");
	Flatfield_UsedSurface->Draw("colz");

	double dummydouble;
	for(int i=1;i<BNL_HighBin_overTime->GetNbinsX()+1;i++){
		if(LengthOfHigh==0||highLowLimit>Flucturations->GetBinContent(i)){ // low <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
			for(int j=1;j<Flucturations->GetNbinsX()+1;j++){
				dummydouble=FlatfieldLow_X->GetBinContent(j)+FlatfieldGenerator_X->GetBinContent(j,i);
				FlatfieldLow_X->SetBinContent(j,dummydouble);
				dummydouble=FlatfieldLow_Y->GetBinContent(j)+FlatfieldGenerator_Y->GetBinContent(j,i);
				FlatfieldLow_Y->SetBinContent(j,dummydouble);
			}
		} else { // high <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
			for(int j=1;j<Flucturations->GetNbinsX()+1;j++){
				dummydouble=FlatfieldHigh_X->GetBinContent(j)+FlatfieldGenerator_X->GetBinContent(j,i);
				FlatfieldHigh_X->SetBinContent(j,dummydouble);
				dummydouble=FlatfieldHigh_Y->GetBinContent(j)+FlatfieldGenerator_Y->GetBinContent(j,i);
				FlatfieldHigh_Y->SetBinContent(j,dummydouble);
			}
			WhereHigh->SetBinContent(i,1);
		}
	}
	
	TCanvas *c4=new TCanvas(TSC::uname());
	WhereHigh->SetTitle("High flat field will be used in these ranges");
	WhereHigh->Draw();
	if(LengthOfHigh!=0)FlatfieldHigh_X->Scale(1./FlatfieldHigh_X->Integral()*256);
	if(LengthOfHigh!=0)FlatfieldHigh_Y->Scale(1./FlatfieldHigh_Y->Integral()*256);
	FlatfieldLow_X->Scale(1./FlatfieldLow_X->Integral()*256);
	FlatfieldLow_Y->Scale(1./FlatfieldLow_Y->Integral()*256);

	FlatfieldHigh_X->SetLineColor(2);
	FlatfieldLow_X->SetLineColor(4);
	FlatfieldHigh_Y->SetLineColor(2);
	FlatfieldLow_Y->SetLineColor(4);

	TCanvas *c41=new TCanvas(TSC::uname());
	FlatfieldLow_X->Draw();
	FlatfieldHigh_X->Draw("same");

	TCanvas *c42=new TCanvas(TSC::uname());
	FlatfieldLow_Y->Draw();
	FlatfieldHigh_Y->Draw("same");


	cout<<"Saving..."<<endl;

	SaveFile->Write();
	cout<<">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"<<endl;
	cout<<">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"<<endl;
	cout<<"Supply the BNL looper with the saved file to use the generated Flatfield: FlatField_output.root"<<endl;
	cout<<">>>>>>>>>>>Note: Be sure to use the same number of time bins!!!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"<<endl;
	cout<<">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"<<endl;
	cout<<">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"<<endl;
	cout<<"Looper done"<<endl;
	
}
void TreeLooper::Loop_BNL_Analysis(double rotation, double xSignalLow,double xSignalHigh,double ySignalLow,double ySignalHigh,char *flatfieldFileName){
	
	//**********************************************************************
	//*****************************Initializing*****************************
	//**********************************************************************
	if(!flatfieldFileName)cout<<"No flat field file found - try running TreeLooper::Loop_BNL_FlatField"<<endl;
	
	double BNL_rotation_angle=rotation;

	cout<<"Fetching flat field from: "<<flatfieldFileName<<"..."<<endl;
	TFile *f=new TFile(flatfieldFileName);
	TH1D *FFHy=(TH1D*)f->Get("FlatfieldHigh_Y");
	TH1D *FFLy=(TH1D*)f->Get("FlatfieldLow_Y");
	TH1D *FFHx=(TH1D*)f->Get("FlatfieldHigh_X");
	TH1D *FFLx=(TH1D*)f->Get("FlatfieldLow_X");
	TH1D *WhereHigh=(TH1D*)f->Get("WhereHigh");
	cout<<"Generating a single y flat field..."<<endl;
	TH1D *FFy=(TH1D*)(FFHy->Clone("Flatfield_Y"));
	FFy->Add(FFLy);
	FFy->Scale(.5);



	
	//**********************************************************************
	//*****************************Initial loop - in needed*****************
	//**********************************************************************

	cout<<"Initializing loop..."<<endl;
	bool HasPulse=false;
	bool HasBNL=false;
	bool HasSANS=false;
	if(fChain->GetListOfBranches()->FindObject("PulseID"))HasPulse=true;
	else cout<<"!!!!!!!!!!!!!!!!!!!!!!WARNING:: The data sample does not containg Pulse information: Run Time will be ajusted!!!"<<endl;
	if(fChain->GetListOfBranches()->FindObject("BNL105FrameIndex"))HasBNL=true;
	if(fChain->GetListOfBranches()->FindObject("SANSFrameIndex"))HasSANS=true;

	double RunTime_Low=1e10,RunTime_High=0;

	fChain->SetBranchStatus("*", 0);
	fChain->SetBranchStatus("Aprox*", 1);
	fChain->SetBranchStatus("*Frame*", 1);
	{
	Long64_t nentries = fChain->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;
	cout<<"Preloop (1/2)- Looping over: "<<nentries<<" events..."<<endl;
	for (Long64_t jentry=0; jentry<nentries;jentry+=5000){
		Long64_t ientry = LoadTree(jentry);
		if(jentry%250000==0)cout<<jentry<<endl;  // Printing status
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
			//fixing Aprox time - if needed:
		if(!HasPulse){
			if(HasBNL)AproxTime=BNL105FrameIndex/40.;
			else if(HasSANS)AproxTime=SANSFrameIndex/40.;
			else cout<<"ERRRORORORORORORO..."<<endl;
		}
		if(AproxTime!=0)RunTime_Low=min(RunTime_Low,AproxTime);
		RunTime_High=max(RunTime_High,AproxTime);
	}
	}
	fChain->SetBranchStatus("*", 1);
	RunTime_High+=(RunTime_High-RunTime_Low)*0.025;


	//**********************************************************************
	//*****************************Declaring output File********************
	//**********************************************************************

	cout<<"Making new file: BNL_Analysis_output.root and saving flat fields..."<<endl;
	TFile *SaveFile=new TFile("BNL_Analysis_output.root","RECREATE");
	FFy->Write();
	FFHx->Write();
	FFLx->Write();
	WhereHigh->Write();

	//**********************************************************************
	//*****************************Declaring Histograms*********************
	//**********************************************************************

	cout<<"Declaring Histograms..."<<endl;

	TH2D *BNL_xy[4][3];

	TH2D *BNL_xt[4][3];

	TH1D *BNL_x[4][3];
	
	TH1D *BNL_y[4][3];

	TH1D *BNL_ToF[4][3];



	char Endingname[100];

	for(int i=0;i<3;i++){
		for(int j=0;j<4;j++){

		if(i==0)sprintf(Endingname,"_HighAndLow");
		if(i==1)sprintf(Endingname,"_HighOnly");
		if(i==2)sprintf(Endingname,"_LowOnly");

		if(j==0)sprintf(Endingname,"%s_nocut",Endingname);
		if(j==1)sprintf(Endingname,"%s_rotated",Endingname);
		if(j==2)sprintf(Endingname,"%s_rotated_flatfield",Endingname);
		if(j==3)sprintf(Endingname,"%s_flatfield",Endingname);

		BNL_xy[j][i]=new TH2D(TSC::CHAR("BNL_xy",Endingname),"",350,-50,300,350,-50,300);
		TSC::NameHist(BNL_xy[j][i],1,"","x-index [#]","y-index [#]","Intensity [PDU]");
		BNL_xy[j][i]->Sumw2();

		BNL_xt[j][i]=new TH2D(TSC::CHAR("BNL_xt",Endingname),"",350,-50,300,1000,0,25000);
		TSC::NameHist(BNL_xt[j][i],1,"","x-index [#]","Time [s]","Intensity [PDU]");
		BNL_xt[j][i]->Sumw2();

		BNL_x[j][i]=new TH1D(TSC::CHAR("BNL_x",Endingname),"",350,-50,300);
		TSC::NameHist(BNL_x[j][i],1,"","x-index [#]","Intensity [PDU]");
		BNL_x[j][i]->Sumw2();

		BNL_y[j][i]=new TH1D(TSC::CHAR("BNL_y",Endingname),"",350,-50,300);
		TSC::NameHist(BNL_y[j][i],1,"","y-index [#]","Intensity [PDU]");
		BNL_y[j][i]->Sumw2();
		
		}
	}

	//**********************************************************************
	//*****************************Declating counters and others************
	//**********************************************************************

	int Counter_NumberOfActiveBNLFrames=0;
	double BNLx_rot,BNLy_rot;
	// switching off some branches to save time...
	fChain->SetBranchStatus("*",0);
	fChain->SetBranchStatus("nBNL",1);
	fChain->SetBranchStatus("BNL",1);
	fChain->SetBranchStatus("BNLx",1);
	fChain->SetBranchStatus("BNLy",1);
	fChain->SetBranchStatus("AproxTime",1);
	fChain->SetBranchStatus("BNL105FrameIndex*",1);

	int AproxtimeWhereHighIndex=1;
	bool NowIsHigh=(WhereHigh->GetBinContent(1));
	double flatfieldW=0;
	//**********************************************************************
	//***************Loop phase - Fetching RunTime information**************
	//**********************************************************************
	cout<<"Initializing loop 2/2..."<<endl;
	int itr,jtr; // iterator
	Long64_t nentries = fChain->GetEntriesFast();
	Long64_t nbytes = 0, nb = 0;
	cout<<"Looping over: "<<nentries<<" events..."<<endl;
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		if(jentry%5000==0)cout<<jentry<<"  "<<TSC::TheMatrix(110)<<endl;  // Printing status
		//fixing Aprox time - if needed:
		if(!HasPulse){
			if(HasBNL){AproxTime=BNL105FrameIndex/40.;}
			else {cout<<"ERROR:: no pulse and no BNL!"<<endl;break;}
		}
		// checking time vs "WhereHigh" - Figuring out if the broken wire was on or off...
		if(AproxTime>=WhereHigh->GetBinLowEdge(AproxtimeWhereHighIndex+1)){
			AproxtimeWhereHighIndex++;
			//if(NowIsHigh!=(WhereHigh->GetBinContent(AproxtimeWhereHighIndex)))cout<<"Swapping to: "<<WhereHigh->GetBinContent(AproxtimeWhereHighIndex)<<endl;
			NowIsHigh=(WhereHigh->GetBinContent(AproxtimeWhereHighIndex));
		}
		
		// Fetching flatfield weight at this point:


		for(itr=0;itr<nBNL;itr++){
			//cutting out strange hits
			//if(BNLx[itr]>254||BNLy[itr]>254||BNLx[itr]<0||BNLy[itr]<0)break;
			// calculating rotation: (fast)
			//BNLx_rot=BNL_rotation_angle*TSC::deg2rad*BNLy[itr];
			//BNLy_rot=BNL_rotation_angle*TSC::deg2rad*BNLx[itr];
			BNLx_rot=cos(BNL_rotation_angle*TSC::deg2rad)*BNLx[itr]-sin(BNL_rotation_angle*TSC::deg2rad)*BNLy[itr];
			BNLy_rot=cos(BNL_rotation_angle*TSC::deg2rad)*BNLy[itr]+sin(BNL_rotation_angle*TSC::deg2rad)*BNLx[itr];
			if(NowIsHigh){
				flatfieldW=1./(FFHx->GetBinContent(BNLx[itr]+1)*FFy->GetBinContent(BNLy[itr]+1));
			} else {
				flatfieldW=1./(FFLx->GetBinContent(BNLx[itr]+1)*FFy->GetBinContent(BNLy[itr]+1));
			}
			if(flatfieldW>300)flatfieldW=0; // this line is extreamly inportant for speed (due to some insane weights when only very few entries exists)
			// filling histograms:
			for(jtr=0;jtr<3;jtr++){
				if(jtr==1&&!NowIsHigh)continue; // only high if 1
				if(jtr==2&&NowIsHigh)continue;  // only low if 2
				
				BNL_xy[0][jtr]->Fill(BNLx[itr],BNLy[itr]);
				BNL_xy[1][jtr]->Fill(BNLx_rot,BNLy_rot);
				BNL_xy[2][jtr]->Fill(BNLx_rot,BNLy_rot,flatfieldW);
				BNL_xy[3][jtr]->Fill(BNLx[itr],BNLy[itr],flatfieldW);

				if(BNLx_rot>xSignalLow&&BNLx_rot<xSignalHigh){
					BNL_y[0][jtr]->Fill(BNLy[itr]);
					BNL_y[1][jtr]->Fill(BNLy_rot);
					BNL_y[2][jtr]->Fill(BNLy_rot,flatfieldW);
					BNL_y[3][jtr]->Fill(BNLy[itr],flatfieldW);
				}
				
				if(BNLy_rot>ySignalLow&&BNLy_rot<ySignalHigh){
					BNL_x[0][jtr]->Fill(BNLx[itr]);
					BNL_x[1][jtr]->Fill(BNLx_rot);
					BNL_x[2][jtr]->Fill(BNLx_rot,flatfieldW);
					BNL_x[3][jtr]->Fill(BNLx[itr],flatfieldW);

					BNL_xt[0][jtr]->Fill(BNLx[itr],BNL[itr]);
					BNL_xt[1][jtr]->Fill(BNLx_rot,BNL[itr]);
					BNL_xt[2][jtr]->Fill(BNLx_rot,BNL[itr],flatfieldW);
					BNL_xt[3][jtr]->Fill(BNLx[itr],BNL[itr],flatfieldW);
				}

				//if(flatfieldW>1000||flatfieldW<.00001)cout<<flatfieldW<<"  --  "<<FFHx->GetBinContent(BNLx[itr]+1)<<"<<>>"<<FFy->GetBinContent(BNLy[itr]+1)<<" :: "<<BNLx[itr]<<"<<>>"<<BNLy[itr]<<endl;
				//if(jentry%5000==0)cout<<nBNL<<","<<BNLx[itr]<<","<<BNLy[itr]<<","<<BNLx_rot<<","<<BNLy_rot<<","<<BNL<<","<<endl;

			}
			// checking time vs "WhereHigh" - Figuring out if the broken wire was on or off...
		}
	}
	cout<<"Loop done"<<endl;
	//**********************************************************************
	//*****************************Analysis phase***************************
	//**********************************************************************
	cout<<"Analysing..."<<endl;



	cout<<"Ploting..."<<endl;

	cout<<"Saving file..."<<endl;
	SaveFile->Write();


}
void Draw_BNL_Analysis(double xSignalLow,double xSignalHigh,double ySignalLow,double ySignalHigh,double xBGLow,double xBGHigh,double yBGLow,double yBGHigh,double VallyLow, double VallyHigh, double DistVallyShoulder,char *filename="BNL_Analysis_output.root"){
	//**********************************************************************
	//*****************************Initializing*****************************
	//**********************************************************************
	cout<<"Loading file and histograms..."<<endl;
	TFile *f=new TFile(filename);

	TH2D *BNL_xy[4][3];
	TH2D *BNL_xt[4][3];
	TH1D *BNL_x[4][3];
	TH1D *BNL_y[4][3];

	TH1D *FFy=(TH1D*)f->Get("Flatfield_Y");
	TH1D *FFHx=(TH1D*)f->Get("FlatfieldHigh_X");
	TH1D *FFLx=(TH1D*)f->Get("FlatfieldLow_X");
	char Endingname[100];
	for(int i=0;i<3;i++){
		for(int j=0;j<4;j++){

			if(i==0)sprintf(Endingname,"_HighAndLow");
			if(i==1)sprintf(Endingname,"_HighOnly");
			if(i==2)sprintf(Endingname,"_LowOnly");

			if(j==0)sprintf(Endingname,"%s_nocut",Endingname);
			if(j==1)sprintf(Endingname,"%s_rotated",Endingname);
			if(j==2)sprintf(Endingname,"%s_rotated_flatfield",Endingname);
			if(j==3)sprintf(Endingname,"%s_flatfield",Endingname);

			BNL_xy[j][i]=(TH2D*)f->Get(TSC::CHAR("BNL_xy",Endingname));
			BNL_xt[j][i]=(TH2D*)f->Get(TSC::CHAR("BNL_xt",Endingname));
			BNL_x[j][i]=(TH1D*)f->Get(TSC::CHAR("BNL_x",Endingname));
			BNL_y[j][i]=(TH1D*)f->Get(TSC::CHAR("BNL_y",Endingname));
		}
	}
	
	//**********************************************************************
	//*****************************Drawing phase****************************
	//**********************************************************************
	TCanvas *c0=new TCanvas(TSC::uname());
	TSC::DrawName((TH1D*)BNL_xy[2][0],"2.7 deg rotation");
	BNL_xy[2][0]->Draw("colz");
	TLine *lin1=new TLine(xSignalLow,-1000,xSignalLow,1000);
	TLine *lin2=new TLine(xSignalHigh,-1000,xSignalHigh,1000);
	TLine *lin3=new TLine(-10000,ySignalLow,10000,ySignalLow);
	TLine *lin4=new TLine(-10000,ySignalHigh,10000,ySignalHigh);
	lin1->Draw("same");
	lin2->Draw("same");
	lin3->Draw("same");
	lin4->Draw("same");

	

	TLine *lin51=new TLine(VallyLow,-1000,VallyLow,1000);
	TLine *lin52=new TLine(VallyHigh,-1000,VallyHigh,1000);
	TLine *lin53=new TLine(VallyLow-DistVallyShoulder,-1000,VallyLow-DistVallyShoulder,1000);
	TLine *lin54=new TLine(VallyHigh+DistVallyShoulder,-1000,VallyHigh+DistVallyShoulder,1000);
	lin51->SetLineColor(6);
	lin52->SetLineColor(6);
	lin53->SetLineColor(8);
	lin54->SetLineColor(8);
	lin51->Draw("same");
	lin52->Draw("same");
	lin53->Draw("same");
	lin54->Draw("same");

	TLine *lin11=new TLine(xBGLow,-1000,xBGLow,1000);
	TLine *lin21=new TLine(xBGHigh,-1000,xBGHigh,1000);
	TLine *lin31=new TLine(-10000,yBGLow,10000,yBGLow);
	TLine *lin41=new TLine(-10000,yBGHigh,10000,yBGHigh);
	lin11->SetLineColor(2);
	lin21->SetLineColor(2);
	lin31->SetLineColor(2);
	lin41->SetLineColor(2);
	lin11->Draw("same");
	lin21->Draw("same");
	lin31->Draw("same");
	lin41->Draw("same");

	TCanvas *c1=new TCanvas(TSC::uname());
	TLegend *l1=new TLegend(0.7,0.2,0.95,0.60);
	TSC::DrawName(BNL_x[0][0],"Flatfielding");
	BNL_x[0][0]->Scale(1/BNL_x[0][0]->Integral()*BNL_x[0][0]->GetNbinsX());
	BNL_x[0][0]->SetLineColor(1);
	l1->AddEntry(BNL_x[0][0],"no flatfield applied","l");

	BNL_x[3][0]->Scale(1/BNL_x[3][0]->Integral()*BNL_x[3][0]->GetNbinsX());
	BNL_x[3][0]->SetLineColor(2);
	l1->AddEntry(BNL_x[3][0],"Flatfield applied (good and bad)","l");	
	BNL_x[3][1]->Scale(1/BNL_x[3][1]->Integral()*BNL_x[3][1]->GetNbinsX());
	BNL_x[3][1]->SetLineColor(4);
	l1->AddEntry(BNL_x[3][1],"Flatfield applied (only bad data)","l");	
	BNL_x[3][2]->Scale(1/BNL_x[3][2]->Integral()*BNL_x[3][2]->GetNbinsX());
	BNL_x[3][2]->SetLineColor(6);
	l1->AddEntry(BNL_x[3][2],"Flatfield applied (only good data)","l");
	BNL_x[0][2]->Scale(1/BNL_x[0][2]->Integral()*BNL_x[0][2]->GetNbinsX());
	BNL_x[0][2]->SetLineColor(8);
	l1->AddEntry(BNL_x[0][2],"no flatfield applied (only good data)","l");

	FFHx->SetLineColor(7);
	l1->AddEntry(FFHx,"Flatfield - High","l");
	FFLx->SetLineColor(5);
	l1->AddEntry(FFLx,"Flatfield - Low","l");


	BNL_x[0][0]->Draw("hist");
	BNL_x[3][0]->Draw("same:hist");
	BNL_x[3][1]->Draw("same:hist");
	BNL_x[3][2]->Draw("same:hist");
	BNL_x[0][2]->Draw("same:hist");
	FFHx->Draw("same");
	FFLx->Draw("same");
	l1->Draw();

// At different times
	TCanvas *c2=new TCanvas(TSC::uname());
	TLegend *l2=new TLegend(0.7,0.2,0.95,0.60);
	TSC::DrawName(BNL_x[2][0],"Time dependent x-structure");
	BNL_x[2][0]->Scale(1./BNL_x[2][0]->Integral()*BNL_x[2][0]->GetNbinsX());
	BNL_x[2][0]->SetLineColor(1);

	BNL_xt[2][0]->GetYaxis()->SetRangeUser(0,500);
	TH1D *H0_500=(TH1D*)(BNL_xt[2][0]->ProjectionX("H0_500"));
	BNL_xt[2][0]->GetYaxis()->SetRangeUser(500,1500);
	TH1D *H500_1500=(TH1D*)(BNL_xt[2][0]->ProjectionX("H500_1500"));
	BNL_xt[2][0]->GetYaxis()->SetRangeUser(1500,5000);
	TH1D *H1500_5000=(TH1D*)(BNL_xt[2][0]->ProjectionX("H1500_5000"));
	BNL_xt[2][0]->GetYaxis()->SetRangeUser(5000,25000);
	TH1D *H5000_25000=(TH1D*)(BNL_xt[2][0]->ProjectionX("H5000_25000"));
	
	cout<<H1500_5000->Integral()<<endl;
	H0_500->Scale(1./H0_500->Integral()*H0_500->GetNbinsX());
	H500_1500->Scale(1./H500_1500->Integral()*H500_1500->GetNbinsX());
	H1500_5000->Scale(1./H1500_5000->Integral()*H1500_5000->GetNbinsX());
	H5000_25000->Scale(1./H5000_25000->Integral()*H5000_25000->GetNbinsX());



	H0_500->SetLineColor(2);
	H500_1500->SetLineColor(4);
	H1500_5000->SetLineColor(6);
	H5000_25000->SetLineColor(8);

	l2->AddEntry(BNL_x[2][0],"Any ToF","l");
	l2->AddEntry(H0_500,"ToF #epsilon [0;500]","l");
	l2->AddEntry(H500_1500,"ToF #epsilon [500;1500]","l");
	l2->AddEntry(H1500_5000,"ToF #epsilon [1500;5000]","l");
	l2->AddEntry(H5000_25000,"ToF #epsilon [5000;25000]","l");

	BNL_x[2][0]->GetYaxis()->SetRangeUser(0,TSC::FindMaximum(BNL_x[2][0],H0_500,H500_1500,H1500_5000,H5000_25000));

	/*
	BNL_x_rotated_flatfield[0]->Draw("hist");
	H0_500->Draw("hist:same");
	H500_1500->Draw("hist:same");
	H1500_5000->Draw("hist:same");
	H5000_25000->Draw("hist:same");
	*/
	BNL_x[2][0]->Draw("e");
	H0_500->Draw("e:same");
	H500_1500->Draw("e:same");
	H1500_5000->Draw("e:same");
	H5000_25000->Draw("e:same");
	
	l2->Draw();

	TCanvas *c3=new TCanvas(TSC::uname());
	TLegend *l3=new TLegend(0.7,0.2,0.95,0.60);
	
	BNL_x[2][0]->Scale(1./BNL_x[2][0]->Integral()*BNL_x[2][0]->GetNbinsX());
	BNL_x[2][0]->SetLineColor(1);
	BNL_xt[2][0]->GetYaxis()->SetRangeUser(10,25000);

	BNL_xt[2][0]->GetXaxis()->SetRangeUser(VallyLow,VallyHigh);
	TH1D *HVally=(TH1D*)(BNL_xt[2][0]->ProjectionY("HVally"));
	BNL_xt[2][0]->GetXaxis()->SetRangeUser(xSignalLow,VallyLow-DistVallyShoulder);
	TH1D *HLowShoulder=(TH1D*)(BNL_xt[2][0]->ProjectionY("HLowShoulder"));
	BNL_xt[2][0]->GetXaxis()->SetRangeUser(VallyHigh+DistVallyShoulder,xSignalHigh);
	TH1D *HHighShoulder=(TH1D*)(BNL_xt[2][0]->ProjectionY("HHighShoulder"));
	TH1D *HShoulders=(TH1D*)(BNL_xt[2][0]->ProjectionY("HShoulders"));

	HShoulders->Add(HLowShoulder);
	HVally->Scale(1./HVally->Integral()*HVally->GetNbinsX());
	HHighShoulder->Scale(1./HHighShoulder->Integral()*HHighShoulder->GetNbinsX());
	HLowShoulder->Scale(1./HLowShoulder->Integral()*HLowShoulder->GetNbinsX());
	HShoulders->Scale(1./HShoulders->Integral()*HShoulders->GetNbinsX());

	HVally->SetLineColor(2);
	HShoulders->SetLineColor(4);
	HLowShoulder->SetLineColor(6);
	HHighShoulder->SetLineColor(8);

	l2->AddEntry(HVally,"Vally","l");
	l2->AddEntry(HShoulders,"Shoulders","l");
	l2->AddEntry(HLowShoulder,"Left side shoulder","l");
	l2->AddEntry(HHighShoulder,"Right side shoulder","l");

	TSC::DrawName(HLowShoulder,"ToF Vally vs Shoulder");
	HLowShoulder->Draw("Hist");
	HHighShoulder->Draw("Hist:same");
	HShoulders->Draw("Hist:same");
	HVally->Draw("Hist:same");
	//****************************************************************************
	//*************************Ending*********************************************
	//****************************************************************************
	cout<<"Drawing phase done"<<endl;
}
void Looper(){
	gStyle->SetOptStat(0);
	DeclareChnlNames();
	/*
	TFile *f=new TFile("c:/WorkDir/LENS/runner/debug/45mmSlabs.root");
	TTree *T;
	f->GetObject("T",T);
	TreeLooper *O=new TreeLooper(T);

	int HighList[]={63,64,96,97};
	
	O->Loop_RunStatistics(2.7);
	O->Loop_BNL_FlatField(2.7,500,10,190,40,230,HighList,.88);
	
	//return;
	O->Loop_BNL_Analysis(2.7,35,170,70,200,"FlatField_output.root");
	Draw_BNL_Analysis(35,170,70,200,10,190,40,230,95,120,15);
	// do draw - Can be run over a outputfile from the loop...
	Draw_RunStatistics();
	

	*/
	// do loop - can be done fully compiled (90% of the root compiled time)
	TFile *f=new TFile("c:/WorkDir/LENS/runner/debug/polysi80_1.root");
	TTree *T;
	f->GetObject("T",T);
	TreeLooper *O=new TreeLooper(T);

	int *HighList=0;  // using a null pointer, as there are no flucturating wire in this run...
	
	//O->Loop_RunStatistics(2.6);
	//O->Loop_BNL_FlatField(2.6,8,188,40,220,HighList,1e16); // 1e16 is to insure that only one (Low) flatfield is used.
	
	O->Loop_BNL_Analysis(2.6,25,170,70,195,"FlatField_output.root");
	Draw_BNL_Analysis(25,170,70,195,8,188,40,220,50,95,15);
	// do draw - Can be run over a outputfile from the loop...
	//Draw_RunStatistics();
	
}
int main(){
	Looper();
	return 0;
}

