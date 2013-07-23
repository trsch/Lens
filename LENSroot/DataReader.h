#pragma once
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string>
#include "TRandom2.h"
//#define _FILE_OFFSET_BITS 64
using namespace std;


class NeutronEventReader{
  char dummy4charpointerobjectthingy[4];
  //bool MakeNextFrameEmpty;
  //BNL rotator
  //int nBNLtotal[256][256];
  //int nBNLtotalX[256];
  //int nBNLtotalY[256];
  //Float_t StatForRotationMinimum;
  //bool StatForRotation; // true by default
  

  // BNL internal Bin Flatening;
  bool DoFlattening;  // false by defaul - uses quite some ekstre time!
  float LastDataValue;
  int dumint;
  int dumjnt;


  Float_t randarray[1024];
  int randarrayIndex;

  // general stuff
  unsigned char DataChunkString[4097];

  int ReadEntry_ReadIndex;
  bool LoadNewDataChunk();
  bool isEOF;
  int dummyindex;
  string FileName;
  //string chnlname[20];
  FILE *DataFile;
  
  //unsigned char XXX[8];
  
  unsigned int FrameCounterThisFile;
  
  bool isAddToChain;
  
  bool ReadEntry();
  bool ReadNextFrame();
  
  bool isSANS;
  bool isBNL105;
  
  float CurrentDataValue;
  int CurrentDataChnl;
  int CurrentDataCtrlBin;
  int CurrentData2ndChnl;
  int CurrentData3rdChnl;
  int CurrentData32bitChnl;
  float FlatFieldArray[512];
  bool isFlatFieldWeight;

 public :
  
  NeutronEventReader(string InputFile, int Instrument,bool addtochain=false);
  ~NeutronEventReader();

  int SkippedCounter;
  //void SetupGetRotation(bool doStatForRotation=true){
  //  StatForRotation=doStatForRotation;
  //}
  //Float_t GetAvgRotation();
  
  // BNL internal Bin Flatening;
  void SetupDoRandomBinFlatening(bool doDoFlattening=true){
    DoFlattening=doDoFlattening;
  }
  bool NoRead0Entry();


  float FlatFieldx[4096];
  float FlatFieldy[4096];

  void SetFlatFieldArray(int inFlatFieldArray[512]){
    for(int i=0;i<512;i++){
      if(inFlatFieldArray[i]>0){FlatFieldArray[i]=1./inFlatFieldArray[i];}
      else {FlatFieldArray[i]=1;}
      //cout<<FlatFieldArray[i]<<endl;
    }
    isFlatFieldWeight=true;
  }

  float DataPoint16PxlId1flattfloat[4096];
  float DataPoint16PxlId2flattfloat[4096];

  // general stuff:
  
  unsigned int CDataPointsThisFile[20];
  unsigned int CDataPointsTotal[20];

  void LoadNewDataFile(string InputFile, int Instrument,bool addtochain=false);

  int FrameEntryCounter[20];
  float DataPointVector[20][4096]; // maximum entries in a single chnl per frame "short";
  int DataPoint17PxlId1[4096]; // SANS 2d plate
  int DataPoint17PxlId2[4096]; 
  int DataPoint16PxlId1[4096]; // BNL105 
  int DataPoint16PxlId2[4096]; 
  unsigned int FrameIndex;

  bool GetNextFrameData();
  void PrintStatus(bool Printfull=false);
 
};



class PulseShapeReader{
  int dummycounter;
  char ReadEntry_chararray[12];
  int ReadEntry_index;
  Float_t dumflo;
  int dumint;
  int dumjnt;
  Float_t FramesPerSec;
  bool isEOF;
  FILE *DataFile;
  string FileName;

  unsigned int FrameCounterThisFile;
  bool isAddToChain;

  bool ReadEntry();

  float FramesPerPulseShapeReading;

  double CurrentDataValue;
    
  float PulseClock;

  double TimeUntilNextPulse;
  
  double DummyTimeStamp;
  string DummyTimeStampChar;
  
  double NextTimeStamp;
  string NextTimeStampChar;



 public :
  void SetFramesPerSec(Float_t freq){
    FramesPerSec=freq;
    FramesPerPulseShapeReading=1./FramesPerSec;
  }
  PulseShapeReader(string InputFile);
  ~PulseShapeReader();

  bool uselastpulsestatistics;
  double ReadTimeStampToDouble(string InputString);
  void SpoolToRightStartingTime(string starttime);
  bool GetTimeStamp();
  

  string StartTimeChar;
  double StartTimeDouble;
 
  double TimeStamp;
  string TimeStampChar;



  double TimeStampPlusFrameTime;

  int PulseVectorLength;
  float PulseVector[10000];
  int PulseVectorIndex[10000];





  unsigned int FrameIndex;
  unsigned int PulseID;
  unsigned int PulseIDthisfile;
  float Integral;
  int FWHM;
  //Float_t LinaryTopSloape;
  int RiseStart;
  float Maximum;
  int MaximumIndex;

  bool GetNextFrameData();
  //  void PrintStatus(bool Printfull=false);
  void SetStartTime(string starttime);
  //  bool GetTimeStamp();
  void LoadNewDataFile(string InputFile,bool addtochain=false);
  
};


NeutronEventReader::~NeutronEventReader(){
	fclose(DataFile);//delete DataFile;
}
NeutronEventReader::NeutronEventReader(string InputFile,int Instrument,bool addtochain){
	DataChunkString[4096]=0;
	DoFlattening=false;
	SkippedCounter=0;
	isFlatFieldWeight=false;
	LoadNewDataFile(InputFile,Instrument,addtochain);
}
void NeutronEventReader::LoadNewDataFile(string InputFile,int Instrument,bool addtochain){
	//Resetting
	isEOF=true;

	if(Instrument==42){FileName=InputFile;}//sprintf(FileName,"%s : SANSFile",InputFile);}
	else if(Instrument==60){FileName=InputFile;}//sprintf(FileName,"%s : BNL105File",InputFile);}
	else {cout<<"Unknown instrument"<<endl;return;}
	isAddToChain=addtochain;
	for(dumint=0;dumint<20;dumint++)CDataPointsThisFile[dumint]=0;
	FrameCounterThisFile=0;
	if(addtochain){
		fclose(DataFile);
	} else {
		FrameIndex=0;
		for(dumint=0;dumint<20;dumint++)CDataPointsTotal[dumint]=0;
	}

	//Setting up
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

	isSANS=false;
	isBNL105=false;
	//=(char*)&InputFile;
	char *DummyFileName=new char[InputFile.size()+1];
	DummyFileName[InputFile.size()]=0;
	memcpy(DummyFileName,InputFile.c_str(),InputFile.size());
	DataFile=fopen(DummyFileName,"rb");
	if(!DataFile)cout<<">>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<"<<endl<<"Warning: "<<InputFile<<" Was not loaded!!!"<<endl<<">>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<"<<endl;
	//delete[InputFile.size()+1] DummyFileName;
	//cout<<Instrument<<" |"<<(int)Instrument[0]<<"|"<<endl;
	if(Instrument==42){isSANS=true;isEOF=false;}
	else if(Instrument==60){isBNL105=true;isEOF=false;}
	else {cout<<"ERROR: Unknown instrument code in NeutronEventReader::LoadNewDataFile"<<endl;}
	ReadEntry_ReadIndex=4096;
}
bool NeutronEventReader::LoadNewDataChunk(){
	if(isEOF)return false;
	 ReadEntry_ReadIndex=-1;
	if(DoFlattening)gRandom->RndmArray(1024,randarray);
	if(!fread(DataChunkString, 4096, 1, DataFile)){
		cout<<"!!!!!!!!!!read fail!!!!!!!!!!!!"<<endl;
		isEOF=true;
		ReadEntry_ReadIndex=4096;
		return false;
	}
	if(!feof(DataFile))return true;
	isEOF=true;
	ReadEntry_ReadIndex=4096;
	return false;
}
bool NeutronEventReader::ReadEntry(){
	randarrayIndex++;
	if(ReadEntry_ReadIndex>4088){
		randarrayIndex=0;
		if(!LoadNewDataChunk())
		{cout<<FileName<<" is EOF! after "<<CDataPointsThisFile[18]<<" read entries, at frame: "<<FrameIndex<<". Last Frame was dumped!"<<endl;isEOF=true;for(dumint=0;dumint<20;dumint++)FrameEntryCounter[dumint]=0;return false;}
	}
	//CurrentDataValue=0.1*(float)*(unsigned int*)&DataChunkString[ReadEntry_ReadIndex+1];
	//CurrentDataValue=.1*(DataChunkString[ReadEntry_ReadIndex+1]+256*DataChunkString[ReadEntry_ReadIndex+2]+65536*DataChunkString[ReadEntry_ReadIndex+3]+16777216*DataChunkString[ReadEntry_ReadIndex+4]);
	CurrentDataValue=0.1*(DataChunkString[ReadEntry_ReadIndex+1]+256*DataChunkString[ReadEntry_ReadIndex+2]+65536*DataChunkString[ReadEntry_ReadIndex+3]);
	//+(16777216*DataChunkString[ReadEntry_ReadIndex+4])*); // this one is dumped as a float only contains 3*8 bits for the number itself and 8 bit for sign and order of magnitude

	CurrentDataChnl=DataChunkString[ReadEntry_ReadIndex+5];
	CurrentData2ndChnl=DataChunkString[ReadEntry_ReadIndex+6];
	CurrentData3rdChnl=DataChunkString[ReadEntry_ReadIndex+7];
	CurrentDataCtrlBin=DataChunkString[ReadEntry_ReadIndex+8];
	ReadEntry_ReadIndex+=8;
	CDataPointsThisFile[18]++;
	CDataPointsTotal[18]++;
	return true;
}
bool NeutronEventReader::NoRead0Entry(){
	//cout<<"Skipped: "<<FrameIndex<<" "<<isBNL105<<endl;
	//for(dumint=0;dumint<howmany;dumint++){
	SkippedCounter++;
	FrameIndex++;
	FrameCounterThisFile++;
	for(dumjnt=0;dumjnt<20;dumjnt++)FrameEntryCounter[dumjnt]=0;
	//}
	return true;
}
bool NeutronEventReader::ReadNextFrame(){

	if(isEOF)return false;
	//if(MakeNextFrameEmpty){if(FrameIndex%100==0)cout<<"skipped"<<endl;MakeNextFrameEmpty=false;return NoRead0Entry();}

	if(CurrentDataValue<0||CDataPointsThisFile[18]==0){
		if(!ReadEntry()){return false;      
		for(dumint=0;dumint<20;dumint++)FrameEntryCounter[dumint]=0;
		return false;
		}
	}

	LastDataValue=0;
	for(dumint=0;dumint<20;dumint++)FrameEntryCounter[dumint]=0;
	FrameIndex++;
	FrameCounterThisFile++;
	while(CurrentDataValue>=(LastDataValue-100)){
		//if(FrameEntryCounter[9]>=2){MakeNextFrameEmpty=true;break;}
		FrameEntryCounter[18]++;
		LastDataValue=CurrentDataValue;
		if(isBNL105){
			if(CurrentDataCtrlBin==0&&CurrentData3rdChnl<2){ // pxl dim 256*304
				CDataPointsThisFile[16]++;
				CDataPointsTotal[16]++;
				DataPointVector[16][FrameEntryCounter[16]]=CurrentDataValue; 
				DataPoint16PxlId1[FrameEntryCounter[16]]=CurrentDataChnl;
				DataPoint16PxlId2[FrameEntryCounter[16]]=CurrentData2ndChnl+256*CurrentData3rdChnl;
				if(isFlatFieldWeight){
					FlatFieldx[FrameEntryCounter[16]]=FlatFieldArray[CurrentDataChnl];
					FlatFieldy[FrameEntryCounter[16]]=FlatFieldArray[CurrentData2ndChnl+256];
					//if(DataPoint16PxlId1[FrameEntryCounter[16]]==42)cout<<FlatFieldxArray[FrameEntryCounter[16]]<<" | "<<FlatFieldyArray[FrameEntryCounter[16]]<<endl;
				}
	
				if(DoFlattening){ 
					DataPoint16PxlId1flattfloat[FrameEntryCounter[16]]=CurrentDataChnl+randarray[randarrayIndex];
					DataPoint16PxlId2flattfloat[FrameEntryCounter[16]]=CurrentData2ndChnl+randarray[randarrayIndex+512];
				}
				FrameEntryCounter[16]++;
			} else {
				CDataPointsThisFile[19]++;
				CDataPointsTotal[19]++;
				FrameEntryCounter[19]++;
			}

		} else if(CurrentDataCtrlBin==64&&CurrentDataChnl<16){
			DataPointVector[CurrentDataChnl][FrameEntryCounter[CurrentDataChnl]]=CurrentDataValue;
			FrameEntryCounter[CurrentDataChnl]++;
			CDataPointsThisFile[CurrentDataChnl]++;
			CDataPointsTotal[CurrentDataChnl]++;
		} else {
			DataPointVector[17][FrameEntryCounter[17]]=CurrentDataValue; 
			DataPoint17PxlId1[FrameEntryCounter[17]]=CurrentDataChnl;
			DataPoint17PxlId2[FrameEntryCounter[17]]=CurrentData2ndChnl+256*CurrentData3rdChnl;
			FrameEntryCounter[17]++;
		}

		if(!ReadEntry()){
			for(dumint=0;dumint<20;dumint++)FrameEntryCounter[dumint]=0;
			return false;
		}
	}
	return true;
}
void NeutronEventReader::PrintStatus(bool Printfull){
	if(Printfull)cout<<endl<<"-----------------STATUS--------------------"<<endl;
	cout<<CDataPointsThisFile[18]<<" entries have been read in "<<FileName<<", which is "<<FrameCounterThisFile<<" frames."<<endl;
	if(isAddToChain)cout<<CDataPointsTotal[18]<<" entries have been read, which is "<<FrameIndex<<" frames."<<endl;
	if(Printfull){

		if(isAddToChain){
		}
	}
}
bool NeutronEventReader::GetNextFrameData(){
	return ReadNextFrame();
}
//***************************************************************
//***************************************************************
//***************************************************************
//***************************************************************
//***************************************************************
//***************************************************************
//*********************Pulse Shape reader************************
//***************************************************************
//***************************************************************
//***************************************************************
//***************************************************************
//***************************************************************
//***************************************************************


PulseShapeReader::~PulseShapeReader(){
	fclose(DataFile);//delete DataFile;
}


double PulseShapeReader::ReadTimeStampToDouble(string InputString){

	/*
	cout<<"<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>"<<endl;
	cout<<"<<<<<<<<<<<<<<<<<<<<<<WARNING!!!>>>>>>>>>>>>>>>>>>>>>>>"<<endl;
	cout<<"<<<<<<<<<<ReadTimeStampToDouble was called!>>>>>>>>>>>>"<<endl;
	cout<<"<<<Dosnt work over daychanges 12/31->1/1 or 2/29->3/1>>"<<endl;
	cout<<"<<<<Also dayligth hour changes will pose a problem>>>>>"<<endl;
	//cout<<"<<<<<<<<<Time was hacked by Troels Schonfeldt>>>>>>>>>>"<<endl;
	cout<<"<<<<<<<<<<<<<<<<<<<<<<WARNING!!!>>>>>>>>>>>>>>>>>>>>>>>"<<endl;
	cout<<"<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>"<<endl;
	*/


  int monthnumber=0;
  
  if(InputString[0]=='J'||InputString[0]=='j'){
    if(InputString[1]=='a'){monthnumber=1;
    } else if(InputString[2]=='l'){monthnumber=7;
    } else if(InputString[2]=='n'){monthnumber=6;
    } else {cout<<"Wrong Time Format: type 1 -|"<<InputString<<"|"<<endl;return false;
    }
  } else if(InputString[0]=='F'||InputString[0]=='f'){monthnumber=2; 
  } else if(InputString[0]=='M'||InputString[0]=='m'){
    if(InputString[1]=='a'){monthnumber=3;
    } else if(InputString[1]=='e'){monthnumber=5;
    } else {cout<<"Wrong Time Format: type 2 -|"<<InputString<<"|"<<endl;return false;
    }
  } else if(InputString[0]=='A'||InputString[0]=='a'){
    if(InputString[1]=='p'){monthnumber=4;
    } else if(InputString[1]=='u'){monthnumber=8;
    } else {cout<<"Wrong Time Format: type 3 -|"<<InputString<<"|"<<endl;return false;
    }
  } else if(InputString[0]=='S'||InputString[0]=='s'){monthnumber=9;
  } else if(InputString[0]=='O'||InputString[0]=='o'){monthnumber=10;
  } else if(InputString[0]=='N'||InputString[0]=='n'){monthnumber=11;
  } else if(InputString[0]=='D'||InputString[0]=='d'){monthnumber=12;
  } else {cout<<"Wrong Time Format: type 4 -|"<<InputString<<"|"<<endl;return false;
  }
  char yearpart[3]={InputString[9],InputString[10],0};
  char daypart[3]={InputString[4],InputString[5],0};
  Float_t daynumber=strtod(daypart,NULL);
  Float_t yearnumber=strtod(yearpart,NULL);
  
  if(monthnumber>1)daynumber+=31;
  if(monthnumber>2){
    daynumber+=28;
    if((((int)yearnumber)%4)==0)daynumber++;
  }
  if(monthnumber>3)daynumber+=31;
  if(monthnumber>4)daynumber+=30;
  if(monthnumber>5)daynumber+=31;
  if(monthnumber>6)daynumber+=30;
  if(monthnumber>7)daynumber+=31;
  if(monthnumber>8)daynumber+=31;
  if(monthnumber>9)daynumber+=30;
  if(monthnumber>10)daynumber+=31;
  if(monthnumber>11)daynumber+=30;
  
  char hourpart[3]={InputString[12],InputString[13],0};  
  char minutespart[3]={InputString[15],InputString[16],0};  
  char secundpart[7]={InputString[18],InputString[19],InputString[20],InputString[21],InputString[22],InputString[23],0};
  //cout<<InputString<<"  |||  "<<((daynumber*24+strtod(hourpart,NULL))*60+strtod(minutespart,NULL))*60+strtod(secundpart,NULL)<<endl;
  return ((daynumber*24+strtod(hourpart,NULL))*60+strtod(minutespart,NULL))*60+strtod(secundpart,NULL);
}

PulseShapeReader::PulseShapeReader(string InputFile){
  FramesPerSec=40;                            // We are running on 40 Hz
  FramesPerPulseShapeReading=1./FramesPerSec;  
  LoadNewDataFile(InputFile,false);
  for(dumint=0;dumint<10000;dumint++)PulseVectorIndex[dumint]=dumint;
}

void PulseShapeReader::LoadNewDataFile(string InputFile,bool addtochain){
 
  isEOF=false;
  FileName=InputFile;
  isAddToChain=addtochain;
  PulseIDthisfile=0;
  if(addtochain){
    fclose(DataFile);//delete DataFile;
  } else {
    PulseID=0;
    FrameIndex=0;
  }
  char *DummyFileName=new char[InputFile.size()+1];
  DummyFileName[InputFile.size()]=0;
  memcpy(DummyFileName,InputFile.c_str(),InputFile.size());
  DataFile=fopen(DummyFileName,"r");
  if(!DataFile){cout<<"Warning: "<<InputFile<<" or "<<DummyFileName<<" Was not loaded!!"<<endl;isEOF=true;
  } else {isEOF=false;}
  //delete[InputFile.size()+1] DummyFileName;
}


void PulseShapeReader::SetStartTime(string starttime){
  dummycounter=0;
  StartTimeChar=starttime;
  StartTimeDouble=ReadTimeStampToDouble(StartTimeChar);
}


void PulseShapeReader::SpoolToRightStartingTime(string starttime){

if(starttime!=0){
    SetStartTime(starttime);
  }
  ReadEntry_index=0;

  GetTimeStamp();
   while(TimeStamp<StartTimeDouble){
     while(fgetc(DataFile)!=10);
     GetTimeStamp();
     if(isEOF){cout<<"The starting point of the file was never found by PulseShapeReader::SpoolToRightStartingTime. Therefore no pulse statistics from this file will be used!"<<endl;return;}
  }
  PulseClock=0;
  NextTimeStamp=TimeStamp;
  NextTimeStampChar=TimeStampChar;
}

bool PulseShapeReader::ReadEntry(){
  for(ReadEntry_index=0;ReadEntry_index<12;ReadEntry_index++){
    ReadEntry_chararray[ReadEntry_index]=fgetc(DataFile);
    if(feof(DataFile)){cout<<FileName<<" is EOF (in ::ReadEntry)!"<<endl;
      isEOF=true;
      return false;
    }else if(ReadEntry_chararray[ReadEntry_index]==44){
      break;
    }else if(ReadEntry_chararray[ReadEntry_index]==10){
      return false;
    }
  }
  CurrentDataValue=strtod(ReadEntry_chararray,NULL);
  return true;
}

bool PulseShapeReader::GetTimeStamp(){
  TimeStampChar="";
  //DummyTimeStampChar=TimeStampChar;
  for(ReadEntry_index=0;ReadEntry_index<24;ReadEntry_index++){
    TimeStampChar+=fgetc(DataFile);
  }
  if(feof(DataFile)){cout<<FileName<<" is EOF (in ::ReadEntry)! after pulse nr: "<<PulseID<<" at time="<<DummyTimeStampChar<<endl;
    if(uselastpulsestatistics){cout<<" If no new Pulse file is loaded - the last pulse statictics will be used for the rest of this run! (From frame: "<<FrameIndex<<")"<<endl;
    } else {
      cout<<" If no new Pulse file is loaded - an empty pulses will be entered for the rest of the run!"<<endl;
      Integral=0;
      FWHM=0; 
      MaximumIndex=0; 
      RiseStart=0;
      Maximum=0;
      PulseVectorLength=0;
    }
    isEOF=true;
    return false;
  }
  TimeStamp=ReadTimeStampToDouble(TimeStampChar);
 	return true;
}


bool PulseShapeReader::GetNextFrameData(){
  PulseVectorLength=0;
  if(isEOF){
    return false;
  }

  PulseClock+=FramesPerPulseShapeReading;
  TimeStampPlusFrameTime+=FramesPerPulseShapeReading;
  FrameIndex++;
  dummycounter++;
  //cout<<NextTimeStamp-TimeStamp<<"  "<<PulseClock<<endl;
  if(dummycounter>5000){
    cout<<PulseClock<<" - "<<FrameIndex<<" | "<<FramesPerPulseShapeReading<<" || "<<NextTimeStamp<<" , "<<TimeStamp<<" = "<<NextTimeStamp-TimeStamp<<endl;
    cout<<NextTimeStampChar<<" || "<<TimeStampChar<<"  "<<PulseID<<endl;
  }
  TimeUntilNextPulse=NextTimeStamp-TimeStamp;
  if(PulseClock>=TimeUntilNextPulse){
    TimeStampPlusFrameTime=TimeStamp+PulseClock;
    dummycounter=0;
    //cout<<FrameIndex<<" "<<PulseClock<<" | "<<(NextTimeStamp-TimeStamp)<<endl;
    PulseClock-=((NextTimeStamp-TimeStamp));
    //if(-1>PulseClock||1<PulseClock)cout<<FrameIndex<<" "<<PulseClock<<" | "<<(NextTimeStamp-TimeStamp)<<endl;
    TimeStamp=NextTimeStamp;
    TimeStampChar=NextTimeStampChar;
    Integral=0;
    Maximum=0;
    MaximumIndex=0;
    PulseID++;
    PulseIDthisfile++;
    while(ReadEntry()){
      PulseVector[PulseVectorLength]=CurrentDataValue;
      if(CurrentDataValue>Maximum){
	Maximum=CurrentDataValue;
	MaximumIndex=PulseVectorLength;
      }
    PulseVectorLength++;
    }
    dumflo=Maximum*0.5;
    for(RiseStart=0;RiseStart<PulseVectorLength;RiseStart++){

      if(PulseVector[RiseStart]>dumflo){
	for(FWHM=0;FWHM<PulseVectorLength-RiseStart;FWHM++){
	  if(PulseVector[FWHM+RiseStart]<dumflo)break;
	}
	if(RiseStart>5&&FWHM>5&&RiseStart+5<PulseVectorLength){
	  for(dumint=RiseStart-5;dumint<RiseStart+5+FWHM;dumint++){
	    Integral+=PulseVector[dumint];
	  }
	} else {Integral=-8;}
	break;
      }
    }
    
    //Junk code for time...
    DummyTimeStamp=TimeStamp;
    DummyTimeStampChar=TimeStampChar;
    GetTimeStamp();    
    NextTimeStampChar=TimeStampChar;
    NextTimeStamp=TimeStamp;
    TimeStamp=DummyTimeStamp;
    TimeStampChar=DummyTimeStampChar;
    if(!isEOF)return true;
  }  
  return true;
}


//**************************************************************
//**************************************************************
//**************************************************************
//**************************************************************
//**************************************************************
//**************************************************************
//**************************************************************
//**************************************************************
//**************************************************************
//**************************************************************


