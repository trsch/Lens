#pragma once
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include <string>
#include "DataReader.h"
#include <time.h>
using namespace std;

 class BatchRunner{
	clock_t tim[100];
	bool SANS_NEOF;
	bool BNL105_NEOF;
	bool PulseShape_NEOF;
	void SansWorker();
	//void TreeWorker();
	//void PrintWorker();
	void PulseWorker();
	void BNLWorker();
	int skipnextBNL;
	int skipnextSANS;

	char AliasAproxTime0Input[100];
	int dumint;
	int dumjnt;
	bool DoBNLSkipper();
	bool DoSANSSkipper();
	string OldTimeStampChar[30];
	bool isFirstBatch;
	bool isSANS;
	bool isBNL105;
	bool isPulseShape;
	bool LookSANS;
	bool LookBNL105;
	bool LookPulseShape;
	NeutronEventReader *SANS; // This will be the SANS file Reader
	NeutronEventReader *BNL105; // This will be the BNL105 file Reader
	PulseShapeReader *PulseShape; // This will be the Pulse shape file Reader
	void SetTreeSANS();
	void SetTreeBNL105();
	void SetTreePulseShape();
	bool TestRunMode;
	int BNL105QueIndex;
	int BNL105NFilesInQue;
	string BNL105QueFileName[50];
	int SANSQueIndex;
	int SANSNFilesInQue;
	string SANSQueFileName[50];
	int PulseShapeQueIndex;
	int PulseShapeNFilesInQue;
	string PulseShapeQueFileName[50];
	bool LoadNextBNL105FromQue();
	bool LoadNextPulseShapeFromQue();
	bool LoadNextSANSFromQue();
	bool FlattFloatBNL;
	bool CalculateRotationMatrix;
	bool HasStartet;
	int fastdummyindex;
	bool DoSkipIfNoPulse;
	NeutronEventReader *FlatField;
	TH2 *FlatFieldHistFull;
	TH2 *FlatFieldHistExpanded;
	int FlatFieldArray[512];
	bool isFlatField;
	float FramesPerSec;
	int FrameCounter;
	int PauseForNFrames;
	TFile *File;
	TTree *T;
	TH2 *BNL2d;
	TH2 *BNL2d2500;
	void PrintStats();

  
	int EntrySkipperBNL_entry[4096];
	int EntrySkipperBNL_SkipN[4096];
	int EntrySkipperBNL_max;
	int EntrySkipperBNL_index;
	int PeriodicSkipperBNL_period;
	int EntrySkipperSANS_entry[4096];
	int EntrySkipperSANS_SkipN[4096];
	int EntrySkipperSANS_max;
	int EntrySkipperSANS_index;
	int PeriodicSkipperSANS_period;


public :

	BatchRunner(string OutPutFileName=0);
	~BatchRunner();



	void ImportFlatField(string FlatFiledFileName);

	void SetAcceleratorFrequency(float setFramesPerSec){
		FramesPerSec=setFramesPerSec;
	}

	void SetSkipBNLIfNoPulse(bool doDoSkipIfNoPulse=true){
		DoSkipIfNoPulse=true;
	}
	void SetEntrySkipperBNL(int setEntrySkipperBNL_entry,int setEntrySkipperBNL_SkipN){
		EntrySkipperBNL_entry[EntrySkipperBNL_max]=setEntrySkipperBNL_entry;
		EntrySkipperBNL_SkipN[EntrySkipperBNL_max]=setEntrySkipperBNL_SkipN;
		EntrySkipperBNL_max++;
	}
	
	void SetPeriodicSkipperBNL(int setPeriodicSkipperBNL_period){
		PeriodicSkipperBNL_period=setPeriodicSkipperBNL_period;
	}
	void SetEntrySkipperSANS(int setEntrySkipperSANS_entry,int setEntrySkipperSANS_SkipN){
		EntrySkipperSANS_entry[EntrySkipperSANS_max]=setEntrySkipperSANS_entry;
		EntrySkipperSANS_SkipN[EntrySkipperSANS_max]=setEntrySkipperSANS_SkipN;
		EntrySkipperSANS_max++;
	}

	void SetPeriodicSkipperSANS(int setPeriodicSkipperSANS_period){
		PeriodicSkipperSANS_period=setPeriodicSkipperSANS_period;
	}


	void setFloatFlatteningOfBNLPxl(bool doFlattFloatBNL=true){
		if(!HasStartet){
			FlattFloatBNL=doFlattFloatBNL;
		} else {cout<<"ERROR: BatchRunner::setFloatFlatteningOfBNLPxl called too late"<<endl;}
	}
	void setCalculateRotationMatrix(bool doCalculateRotationMatrix=true){
		if(!HasStartet){
			CalculateRotationMatrix=doCalculateRotationMatrix;
		} else {cout<<"ERROR: BatchRunner::setCalculateRotationMatrix called too late"<<endl;}
	}

	void SetTestRunMode(bool Mode=true){
		if(!HasStartet){
			TestRunMode=Mode;
		} else {cout<<"ERROR: BatchRunner::SetTestRunMode called too late"<<endl;}
	}





	void AddBNL105FileToQue(string QueThisFileName);
	void AddSANSFileToQue(string QueThisFileName);
	void AddPulseShapeFileToQue(string QueThisFileName);



	void RunBatch(int StartFrameSANS=0,int StartFrameBNL105=0,int StartFramePulseShape=0,string StartTime=0);
	void LoadBatch(string filenameSANS=0, string filenameBNL105=0, string filenamePulseShape=0);
};


BatchRunner::~BatchRunner(){

  if(CalculateRotationMatrix){
    cout<<"Calculating rotation matrix and center..."<<endl;
    T->SetAlias("BNLtheta","3.5*1.");
    T->SetAlias("BNLCenterY","110*1.");
    T->SetAlias("BNLCenterX","130*1.");
    T->SetAlias("BNLyrot","(BNLx-BNLCenterX)*sin(BNLtheta)+(BNLy-BNLCenterY)*cos(BNLtheta)+BNLCenterY");
    T->SetAlias("BNLxrot","(BNLx-BNLCenterX)*cos(BNLtheta)-(BNLy-BNLCenterY)*sin(BNLtheta)+BNLCenterX");

  }
  cout<<"Optimizing tree for speed and size"<<endl;
  T->OptimizeBaskets();
  cout<<"Closing file"<<endl;
  File->Write();
  File->Close();
  if(BNL105)delete BNL105;
  if(SANS)delete SANS;
  if(PulseShape)delete PulseShape;
	cout<<"BatchRunner finished!"<<endl;

}
BatchRunner::BatchRunner(string OutPutFileName){
	BNL105=0;
	SANS=0;
	PulseShape=0;
	isFlatField=false;
	skipnextBNL=0;
	skipnextSANS=0;
	AliasAproxTime0Input[0]=0;
	cout<<"Creating file"<<endl;
	HasStartet=false;
	FlattFloatBNL=false;
	CalculateRotationMatrix=true;
	EntrySkipperBNL_max=0;
	EntrySkipperBNL_index=0;
	PeriodicSkipperBNL_period=0;

	if(OutPutFileName.size()!=0){
	cout<<"here"<<endl;
		char *DummyFileName=new char[OutPutFileName.size()+1];
		DummyFileName[OutPutFileName.size()]=0;
		memcpy(DummyFileName,OutPutFileName.c_str(),OutPutFileName.size());
		File = new TFile(DummyFileName,"RECREATE");
	} else {
		
	cout<<"here"<<endl;

		File = new TFile("Output.root","RECREATE");
	}
	File->cd();
	T=new TTree("T","DataTree");

	FrameCounter=0;
	FramesPerSec=40;
	isFirstBatch=true;
	isSANS=false;
	isBNL105=false;
	isPulseShape=false;
	TestRunMode=false;
	cout<<"initilization done"<<endl;

}
bool BatchRunner::DoBNLSkipper(){
  if(LookPulseShape)if(PulseShape->Integral==-8){BNL105->NoRead0Entry();return true;}
  if(EntrySkipperBNL_max>EntrySkipperBNL_index){
    if(BNL105->FrameIndex>=EntrySkipperBNL_entry[EntrySkipperBNL_index]){
      if(BNL105->FrameIndex==EntrySkipperBNL_entry[EntrySkipperBNL_index]){
	BNL105->NoRead0Entry();
	EntrySkipperBNL_SkipN[EntrySkipperBNL_index]--;
	if(EntrySkipperBNL_SkipN[EntrySkipperBNL_index]>0){
	  EntrySkipperBNL_entry[EntrySkipperBNL_index]++;
	} else {
	  EntrySkipperBNL_index++;
	}
	return true;
      } else {cout<<"ERROR: Skipping algoritm failed in BatchRunner::DoBNLSkipper"<<endl;}
    }
    }
  if(PeriodicSkipperBNL_period!=0){
    if((BNL105->FrameIndex%PeriodicSkipperBNL_period)!=0)return false;
    BNL105->NoRead0Entry();
    return true;
  }
  return false;
}
bool BatchRunner::DoSANSSkipper(){
  if(skipnextSANS>=3){
    skipnextSANS--;
    SANS->NoRead0Entry();
    return true;
  }
  if(SANS->FrameEntryCounter[9]>=2){  
    skipnextSANS=SANS->FrameEntryCounter[9];
    SANS->NoRead0Entry();
    return true;
  }
  if(EntrySkipperSANS_max>EntrySkipperSANS_index){
    if(SANS->FrameIndex>=EntrySkipperSANS_entry[EntrySkipperSANS_index]){
      if(SANS->FrameIndex==EntrySkipperSANS_entry[EntrySkipperSANS_index]){
	SANS->NoRead0Entry();
	EntrySkipperSANS_SkipN[EntrySkipperSANS_index]--;
	if(EntrySkipperSANS_SkipN[EntrySkipperSANS_index]>0){
	  EntrySkipperSANS_entry[EntrySkipperSANS_index]++;
	} else {
	  EntrySkipperSANS_index++;
	}
	return true;
      } else {cout<<"ERROR: Skipping algoritm failed in BatchRunner::DoSANSSkipper"<<endl;}
    }
  }
  if(PeriodicSkipperSANS_period!=0){
     if((SANS->FrameIndex%PeriodicSkipperSANS_period)!=0)return false;
    SANS->NoRead0Entry();
    return true;
    }
  return false;
}
void BatchRunner::LoadBatch(string filenameSANS, string filenameBNL105, string filenamePulseShape){
  if(filenameSANS!="0"){
    SANSQueIndex=0;
    SANSNFilesInQue=0;
    if(!isSANS)SANS=new NeutronEventReader(filenameSANS,42); // 42 is the code for SANS instrument
    if(isSANS)SANS->LoadNewDataFile(filenameSANS,42,isSANS);
    if(!isSANS)SetTreeSANS();
    isSANS=true;
  } else {isSANS=false;}
  if(filenameBNL105!="0"){
    BNL105QueIndex=0;
    BNL105NFilesInQue=0;
    if(!isBNL105){
      BNL105=new NeutronEventReader(filenameBNL105,60); // 60 is the code for BNL105 - and it is 42 written in base 7
      SetTreeBNL105();
      BNL105->SetupDoRandomBinFlatening(FlattFloatBNL); // Control of this is important - very timeconsuming!!!
      if(isFlatField)BNL105->SetFlatFieldArray(FlatFieldArray);
      if(CalculateRotationMatrix){  // a simple algoritm for setting up some initial aliases for rotations (to be done in the destructor - not yet coded...).
	BNL2d=new TH2F("BNLPxlHits","BNLPxlHits",256,.5,256.5,256,.5,256.5);
	BNL2d2500=new TH2F("BNLColdPxlHits","BNLColdPxlHits",256,.5,256.5,256,.5,256.5);
      }
    }
    if(isBNL105)BNL105->LoadNewDataFile(filenameBNL105,60,isBNL105);
    isBNL105=true;
  } else {isBNL105=false;}
  
  if(filenamePulseShape!="0"){
    PulseShapeQueIndex=0;
    PulseShapeNFilesInQue=0;
    if(!isPulseShape)PulseShape=new PulseShapeReader(filenamePulseShape);
    if(isPulseShape)PulseShape->LoadNewDataFile(filenamePulseShape,isPulseShape);
    if(!isPulseShape)SetTreePulseShape();
    isPulseShape=true;
  } else {isPulseShape=false;}
}
void BatchRunner::ImportFlatField(string FlatFiledFileName){
  FlatField=new NeutronEventReader(FlatFiledFileName,60);
  FlatFieldHistFull=new TH2F("FlatFieldFull","FlatFieldFull",256,0,256,256,0,256);
  FlatFieldHistExpanded=new TH2F("FlatFieldExpanded","FlatFieldExpanded",256,0,256,256,0,256);
 
  for(int i=0;i<256;i++){
    FlatFieldArray[i]=0;
    FlatFieldArray[i+256]=0;
  }
  T->Branch("nBNLFlatField",&FlatField->FrameEntryCounter[16],"FlatField_nHitsChnl16/I");
  T->Branch("BNLFlatField",FlatField->DataPointVector[16],"FlatField_HitsChnl16[FlatField_nHitsChnl16]/F");
  if(!FlattFloatBNL){
    T->Branch("BNLxFlatField",FlatField->DataPoint16PxlId1,"FlatField_HitsChnl16Pxl1[FlatField_nHitsChnl16]/I");
    T->Branch("BNLyFlatField",FlatField->DataPoint16PxlId2,"FlatField_HitsChnl16Pxl2[FlatField_nHitsChnl16]/I");
  } else {
    T->Branch("BNLxFlatField",FlatField->DataPoint16PxlId1flattfloat,"FlatField_HitsChnl16Pxl1flattfloat[FlatField_nHitsChnl16]/F");
    T->Branch("BNLyFlatField",FlatField->DataPoint16PxlId2flattfloat,"FlatField_HitsChnl16Pxl2flattfloat[FlatField_nHitsChnl16]/F");
  }
  while(FlatField->GetNextFrameData()){ 
    for(fastdummyindex=0;fastdummyindex<FlatField->FrameEntryCounter[16];fastdummyindex++){
      if(FlatField->DataPoint16PxlId2[fastdummyindex]<256&&FlatField->DataPoint16PxlId1[fastdummyindex]<256){
	FlatFieldHistFull->Fill(FlatField->DataPoint16PxlId2[fastdummyindex],FlatField->DataPoint16PxlId1[fastdummyindex]);
	FlatFieldArray[FlatField->DataPoint16PxlId2[fastdummyindex]+256]++;
	FlatFieldArray[FlatField->DataPoint16PxlId1[fastdummyindex]]++;
      }
    }
    T->Fill();
  }
  FlatField->FrameEntryCounter[16]=0;
  
  for(int i=0;i<256;i++){
    for(int j=0;j<256;j++)FlatFieldHistExpanded->Fill(j,i,FlatFieldArray[i]*FlatFieldArray[j+256]);
  }
  //for(int i=0;i<256;i++){if(FlatFieldArray[i+256]!=0)FlatFieldArray[i+256]=1./FlatFieldArray[i+256];if(FlatFieldArray[i]!=0)FlatFieldArray[i]=1./FlatFieldArray[i];}    
  isFlatField=true;
  

}
void BatchRunner::SansWorker(){
	if(LookSANS)
      if(!DoSANSSkipper())
	if(SANS->GetNextFrameData()){}else{
	  LoadNextSANSFromQue();
	  if(LookSANS){if(!SANS->GetNextFrameData())SANS_NEOF=false;
	  }else {SANS_NEOF=false;}
	}
}
void BatchRunner::PulseWorker(){
   //tim[1]=clock();
		if(LookPulseShape)if(PulseShape->GetNextFrameData()){}else{LoadNextPulseShapeFromQue();}
   //tim[2]+=tim[1]-clock();
}
void BatchRunner::BNLWorker(){
	if(LookBNL105)
      if(!DoBNLSkipper())
	if(BNL105->GetNextFrameData()){
	if(CalculateRotationMatrix)
	    for(fastdummyindex=0;fastdummyindex<BNL105->FrameEntryCounter[16];fastdummyindex++){
	      BNL2d->Fill(BNL105->DataPoint16PxlId2[fastdummyindex],BNL105->DataPoint16PxlId1[fastdummyindex]);
	      if(BNL105->DataPointVector[16][fastdummyindex]>2500)
		BNL2d2500->Fill(BNL105->DataPoint16PxlId2[fastdummyindex],BNL105->DataPoint16PxlId1[fastdummyindex]);
	    }
	}else{
	  LoadNextBNL105FromQue();
	  if(LookBNL105){if(!BNL105->GetNextFrameData())BNL105_NEOF=false;
	  }else {BNL105_NEOF=false;}
	}
}
void BatchRunner::RunBatch(int StartFrameSANS,int StartFrameBNL105,int StartFramePulseShape,string StartTime){
  HasStartet=true;

  if(isPulseShape){
    PulseShape->SetFramesPerSec(FramesPerSec);
    if(StartTime!="0"){
      cout<<"Fast fowarding to start time... : "<<StartTime<<endl;
      PulseShape->SpoolToRightStartingTime(StartTime);
      //cout<<StartTime<<endl;
    } else {
      PulseShape->GetTimeStamp();
    }
    cout<<"Logging startet at time="<<PulseShape->TimeStampChar<<" : "<<PulseShape->TimeStamp<<" secunds after new year"<<endl;
    if(AliasAproxTime0Input[0]==0){
      //sprintf(AliasAproxTime0Input,"AproxTime-%ld",PulseShape->TimeStamp);
      //cout<<"-----------"<<AliasAproxTime0Input<<"-----------"<<endl;
      T->SetAlias("AproxTime0",AliasAproxTime0Input);
    }
  }


  if(isSANS&&StartFrameSANS!=0)SANS->FrameIndex=StartFrameSANS;
  if(isPulseShape&&StartFramePulseShape!=0)PulseShape->FrameIndex=StartFramePulseShape;
  if(isBNL105&&StartFrameBNL105!=0)BNL105->FrameIndex=StartFrameBNL105;
  SANS_NEOF=isSANS;
  BNL105_NEOF=isBNL105;
  PulseShape_NEOF=true;
  //double Integral=0;
 
  LookSANS=isSANS;
  LookPulseShape=isPulseShape;
  LookBNL105=isBNL105;
  int PrinterCounter=0;
  int Dummy__counter=0;
   while(SANS_NEOF||BNL105_NEOF){
	   Dummy__counter++;
	   if(Dummy__counter%25000==1)cout<<Dummy__counter-1<<endl;
	PulseWorker();
	BNLWorker();	
	SansWorker();
	T->Fill();
    PrinterCounter++;
	//if(PrinterCounter%10000==0)cout<<tim[2]<<" - "<<tim[3]<<" | "<<tim[4]<<endl;
    if(PrinterCounter%100000==0){PrintStats();if(TestRunMode)break;}
  }
  cout<<"The run ended at ";
  if(isPulseShape)cout<<"time: "<<PulseShape->TimeStampChar<<endl;
  if(isSANS)cout<<" - Frames(skipped) SANS   : "<<SANS->FrameIndex<<"("<<SANS->SkippedCounter<<") containing "<<SANS->CDataPointsTotal[18]<<" entries"<<endl;
  if(isBNL105)cout<<" - Frames(entries) BNL105 : "<<BNL105->FrameIndex<<"("<<BNL105->SkippedCounter<<") containing "<<BNL105->CDataPointsTotal[18]<<" entries"<<endl;
  if(isPulseShape)cout<<" - in "<<PulseShape->PulseID<<" pulses"<<endl;
  
}
void BatchRunner::AddSANSFileToQue(string QueThisFileName){
  if(SANSNFilesInQue>=50)cout<<"ERROR in BatchRunner::AddSANSFileToQue"<<endl;
  SANSQueFileName[SANSNFilesInQue]=QueThisFileName;
  SANSNFilesInQue++;
}
void BatchRunner::AddPulseShapeFileToQue(string QueThisFileName){
  if(PulseShapeNFilesInQue>=50)cout<<"ERROR in BatchRunner::AddPulseShapeFileToQue"<<endl;
  PulseShapeQueFileName[PulseShapeNFilesInQue]=QueThisFileName;
  PulseShapeNFilesInQue++;

}
void BatchRunner::AddBNL105FileToQue(string QueThisFileName){
  if(BNL105NFilesInQue>=50)cout<<"ERROR in BatchRunner::AddBNL105FileToQue"<<endl;
  BNL105QueFileName[BNL105NFilesInQue]=QueThisFileName;
  BNL105NFilesInQue++;
}
bool BatchRunner::LoadNextBNL105FromQue(){
  if(BNL105QueIndex>=BNL105NFilesInQue){cout<<"No new files in que"<<endl;LookBNL105=false;return false;}
  BNL105->LoadNewDataFile(BNL105QueFileName[BNL105QueIndex],60,isBNL105);
  BNL105QueIndex++;
  isBNL105=true;
  cout<<"File: "<<BNL105QueFileName[BNL105QueIndex]<<" was loaded"<<endl;
  return true;
}
bool BatchRunner::LoadNextPulseShapeFromQue(){
  if(PulseShapeQueIndex>=PulseShapeNFilesInQue){cout<<"No new files in que"<<endl;LookPulseShape=false;return false;}
  PulseShape->LoadNewDataFile(PulseShapeQueFileName[PulseShapeQueIndex],isPulseShape);
  PulseShapeQueIndex++;
  isPulseShape=true;
  cout<<"File: "<<PulseShapeQueFileName[BNL105QueIndex]<<" was loaded"<<endl;
  return true;
}
bool BatchRunner::LoadNextSANSFromQue(){
  if(SANSQueIndex>=SANSNFilesInQue){cout<<"No new files in que"<<endl;LookSANS=false;return false;}
  SANS->LoadNewDataFile(SANSQueFileName[SANSQueIndex],42,isSANS);
  SANSQueIndex++;
  isSANS=true;
  cout<<"File: "<<SANSQueFileName[BNL105QueIndex]<<" was loaded"<<endl;
  return true;
}
void BatchRunner::SetTreeSANS(){
  T->Branch("SANSFrameIndex",&SANS->FrameIndex,"SANSFrameIndex/i");
  T->Branch("nC0",&SANS->FrameEntryCounter[0],"nHitsChnl0/I");
  T->Branch("nC1",&SANS->FrameEntryCounter[1],"nHitsChnl1/I");
  T->Branch("nC2",&SANS->FrameEntryCounter[2],"nHitsChnl2/I");
  T->Branch("nC3",&SANS->FrameEntryCounter[3],"nHitsChnl3/I");
  T->Branch("nC4",&SANS->FrameEntryCounter[4],"nHitsChnl4/I");
  T->Branch("nC5",&SANS->FrameEntryCounter[5],"nHitsChnl5/I");
  T->Branch("nC6",&SANS->FrameEntryCounter[6],"nHitsChnl6/I");
  T->Branch("nC7",&SANS->FrameEntryCounter[7],"nHitsChnl7/I");
  T->Branch("nC8",&SANS->FrameEntryCounter[8],"nHitsChnl8/I");
  T->Branch("nC9",&SANS->FrameEntryCounter[9],"nHitsChnl9/I");
  T->Branch("nC10",&SANS->FrameEntryCounter[10],"nHitsChnl10/I");
  T->Branch("nC11",&SANS->FrameEntryCounter[11],"nHitsChnl11/I");
  T->Branch("nC12",&SANS->FrameEntryCounter[12],"nHitsChnl12/I");
  T->Branch("nC13",&SANS->FrameEntryCounter[13],"nHitsChnl13/I");
  T->Branch("nC14",&SANS->FrameEntryCounter[14],"nHitsChnl14/I");
  T->Branch("nC15",&SANS->FrameEntryCounter[15],"nHitsChnl15/I");
  //T->Branch("nHitsChnl16",&SANS->FrameEntryCounter[16],"nHitsChnl16/I");
  T->Branch("nNot6417",&SANS->FrameEntryCounter[17],"nHitsChnl17/I");
  T->Branch("nTotal",&SANS->FrameEntryCounter[18],"nHitsChnl18/I");
  
  T->Branch("C0",SANS->DataPointVector[0],"HitsChnl0[nHitsChnl0]/F");
  T->Branch("C1",SANS->DataPointVector[1],"HitsChnl1[nHitsChnl1]/F");
  T->Branch("C2",SANS->DataPointVector[2],"HitsChnl2[nHitsChnl2]/F");
  T->Branch("C3",SANS->DataPointVector[3],"HitsChnl3[nHitsChnl3]/F");
  T->Branch("C4",SANS->DataPointVector[4],"HitsChnl4[nHitsChnl4]/F");
  T->Branch("C5",SANS->DataPointVector[5],"HitsChnl5[nHitsChnl5]/F");
  T->Branch("C6",SANS->DataPointVector[6],"HitsChnl6[nHitsChnl6]/F");
  T->Branch("C7",SANS->DataPointVector[7],"HitsChnl7[nHitsChnl7]/F");
  T->Branch("C8",SANS->DataPointVector[8],"HitsChnl8[nHitsChnl8]/F");
  T->Branch("C9",SANS->DataPointVector[9],"HitsChnl9[nHitsChnl9]/F");
  T->Branch("C10",SANS->DataPointVector[10],"HitsChnl10[nHitsChnl10]/F");
  T->Branch("C11",SANS->DataPointVector[11],"HitsChnl11[nHitsChnl11]/F");
  T->Branch("C12",SANS->DataPointVector[12],"HitsChnl12[nHitsChnl12]/F");
  T->Branch("C13",SANS->DataPointVector[13],"HitsChnl13[nHitsChnl13]/F");
  T->Branch("C14",SANS->DataPointVector[14],"HitsChnl14[nHitsChnl14]/F");
  T->Branch("C15",SANS->DataPointVector[15],"HitsChnl15[nHitsChnl15]/F");
  //T->Branch("HitsChnl16",SANS->DataPointVector[16],"HitsChnl16[nHitsChnl16]/F");
  T->Branch("CNot64",SANS->DataPointVector[17],"HitsChnl17[nHitsChnl17]/F");
  T->Branch("CNot64Pxl1",SANS->DataPoint17PxlId1,"HitsChnl17Pxl1[nHitsChnl17]/I");
  T->Branch("CNot64Pxl2",SANS->DataPoint17PxlId2,"HitsChnl17Pxl2[nHitsChnl17]/I");
}
void BatchRunner::SetTreeBNL105(){
  T->Branch("BNL105FrameIndex",&BNL105->FrameIndex,"BNL105FrameIndex/i");
  T->Branch("nBNL",&BNL105->FrameEntryCounter[16],"nHitsChnl16/I");
  T->Branch("BNL",BNL105->DataPointVector[16],"HitsChnl16[nHitsChnl16]/F");
  if(!FlattFloatBNL){
    T->Branch("BNLx",BNL105->DataPoint16PxlId1,"HitsChnl16Pxl1[nHitsChnl16]/I");
    T->Branch("BNLy",BNL105->DataPoint16PxlId2,"HitsChnl16Pxl2[nHitsChnl16]/I");
  } else {
    T->Branch("BNLx",BNL105->DataPoint16PxlId1flattfloat,"HitsChnl16Pxl1flattfloat[nHitsChnl16]/F");
    T->Branch("BNLy",BNL105->DataPoint16PxlId2flattfloat,"HitsChnl16Pxl2flattfloat[nHitsChnl16]/F");
  }
  if(isFlatField){
    T->Branch("BNLxWeight",BNL105->FlatFieldx,"FlatFieldx[nHitsChnl16]/F");
    T->Branch("BNLyWeight",BNL105->FlatFieldy,"FlatFieldy[nHitsChnl16]/F");
  }
}
void BatchRunner::SetTreePulseShape(){
  T->Branch("PulseShapeFrameIndex",&PulseShape->FrameIndex,"PulseShapeFrameIndex/i");
  T->Branch("PulseID",&PulseShape->PulseID,"PulseID/i");
  T->Branch("Integral",&PulseShape->Integral,"Integral/F");
  T->Branch("FWHM",&PulseShape->FWHM,"FWHM/I"); 
  T->Branch("MaximumIndex",&PulseShape->MaximumIndex,"MaximumIndex/I"); 
  T->Branch("RiseStart",&PulseShape->RiseStart,"RiseStart/I"); 
  T->Branch("Maximum",&PulseShape->Maximum,"Maximum/F"); 
  T->Branch("PulseVectorLength",&PulseShape->PulseVectorLength,"PulseVectorLength/I");
  T->Branch("PulseVector",PulseShape->PulseVector,"PulseVector[PulseVectorLength]/F");
  T->Branch("PulseVectorIndex",PulseShape->PulseVectorIndex,"PulseVectorIndex[PulseVectorLength]/I");
  T->Branch("AproxTime",&PulseShape->TimeStampPlusFrameTime,"TimeStampPlusFrameTime/D");

}
void BatchRunner::PrintStats(){


  cout<<endl<<"Frame: ";
  
  if(isSANS&&!isBNL105){
    cout<<SANS->FrameIndex;
  } else if(!isSANS&&isBNL105){ 
    cout<<BNL105->FrameIndex;
  } else if (isSANS&&isBNL105){
    if(SANS->FrameIndex>=BNL105->FrameIndex){
      cout<<SANS->FrameIndex;
    }else {
      cout<<BNL105->FrameIndex;
    }
  } else {
    cout<<"NA";
  }
  
  if(isSANS)cout<<"\t SANS Entries: "<<SANS->CDataPointsTotal[18];
  if(isBNL105)cout<<"\t BNL105 Entries: "<<BNL105->CDataPointsTotal[18];
  if(isPulseShape)cout<<"\t PulseID: "<<PulseShape->PulseID<<" at: "<<PulseShape->TimeStampChar;
  cout<<endl;
  if(isPulseShape)cout<<"  PulseShape::  -Integral="<<PulseShape->Integral<<" \t-FWHM="<<PulseShape->FWHM<<"\t -StartIndex="<<PulseShape->RiseStart<<"\t -Maximum="<<PulseShape->Maximum<<" at="<<PulseShape->MaximumIndex<<endl;
  
}

