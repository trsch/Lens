
//***********************************
// READ runner45mmWater() for an example of use
//***********************************

//***********************************
// Setting up to run c++:
//***********************************
#include <iostream>
#include "TTree.h"
#include "TROOT.h" 
#include "TFile.h"
#include <stdio.h>

using namespace std;

//**********************************************************************
// Includes needed by this program:
//**********************************************************************
#include "BatchRunner.h"  // used for compiling using a libery, in root or root indepentently
//#include "BatchRunner.cxx" // used for "root -l runner.cpp" with no libery compiling 
//#include "DataReader.cxx" // used for "root -l runner.cpp" with no libery compiling 
char __DataPathName[500];
char __GlobalCHARCharArray[500];
char *CHAR(const char *c1,const char *c2, const char *c3=0, const char *c4=0, const char *c5=0, const char *c6=0, const char *c7=0){
sprintf(__GlobalCHARCharArray,"%s%s\0",c1,c2);
if(c3)sprintf(__GlobalCHARCharArray,"%s%s\0",__GlobalCHARCharArray,c3);
if(c4)sprintf(__GlobalCHARCharArray,"%s%s\0",__GlobalCHARCharArray,c4);
if(c5)sprintf(__GlobalCHARCharArray,"%s%s\0",__GlobalCHARCharArray,c5);
if(c6)sprintf(__GlobalCHARCharArray,"%s%s\0",__GlobalCHARCharArray,c6);
if(c7)sprintf(__GlobalCHARCharArray,"%s%s\0",__GlobalCHARCharArray,c7);
return __GlobalCHARCharArray;
}
void runner45mmWater(bool inTestmode=false);
void runner45mmSlabs(bool inTestmode=false);
void runnerBGRun(bool inTestmode=false);

void runnerOld_polysi300(bool inTestmode=false);
void runnerOld_polysi80_1(bool inTestmode=false);
void runnerOld_polysi80_2(bool inTestmode=false);
void runnerOld_polymono300(bool inTestmode=false);
void runnerOld_polymono80(bool inTestmode=false);
void runnerOld_polyvac300_1(bool inTestmode=false);
void runnerOld_polyvac300_2(bool inTestmode=false);

void runner(){
       sprintf(__DataPathName,getenv("LENSDATA"));

	
	runnerOld_polysi300();
	runnerOld_polysi80_1();
	runnerOld_polysi80_2();
	runnerOld_polymono300();
	runnerOld_polymono80();
	runnerOld_polyvac300_1();
	runnerOld_polyvac300_2();
	
	runnerBGRun();
	runner45mmWater();
	runner45mmSlabs();
  
	/*  // example of how to run a batch:
	BatchRunner *BR;
	BR=new BatchRunner();
	BR->setCalculateRotationMatrix(true);
	BR->setFloatFlatteningOfBNLPxl(true);
	BR->SetTestRunMode(false);
	BR->LoadBatch(0,"../BNLTEST_15_662.dat",0);
	BR->RunBatch();
	cout<<"ending program"<<endl;
	delete BR;
	cout<<"Program done"<<endl;
	return;
	*/
}


void runner45mmWater(bool inTestmode){

	cout<<"<<<Runner: Initializing 45mmWater>>>"<<endl;

	BatchRunner *BR;
	//**********************************************************************
	// RunBatch(SANS_File_name , BNL105_File_name , PulseShape_File_name , StartTime)
	// StartTime must be provided in the format: "MMM-dd-yyyy hh:mm:ss.mmm" 
	//   where the last mmm is millisec and the first MMM is 3 letters month tag fx: "Jul" or "Dec"
	//**********************************************************************

	//**********************************************************************
	// Setup phase.
	//**********************************************************************
	// initializing the BatchRunner
	BR=new BatchRunner("45mmWater.root");

	// Setting up for advanced smearing of pxl hits in BNL detector
	// if not called or set to false the BNL data is clean/raw
	// else (if set true) each pxl bin is smeared with + random float [0;1[
	//   this avoids binning-bios when rotating the 2d pxl detector view.
	//   but at the cost of:
	//           1: more storrage space used by TFile
	//           2: much more caLculation time used by BatchRunner/DataReader
	//           3: a tiny bit (~10%) more reading time when reading BNL pxl hits
	//               from the final TFile/TTree when analysing.

	BR->setFloatFlatteningOfBNLPxl(true); // false by default

	// setting up to calculate a roration matrix, increases calculation time in
	//   DataReader/BatchRunner, but eases life while analysing
	BR->setCalculateRotationMatrix(true); // true by default


	// Set to run a small test run, if you just want a small overview TFile :)
	//    i (TSC) used this alot while debugging
	BR->SetTestRunMode(inTestmode); // false by default


	/* // you can set up lists of frames that will be entered into the datafile as 0 entries not reading data
	// This is a way to syncronize the BNL and SANS DataFile
	// the example below skips every 3143 SANS frame periodically, along with 800 BNL frames from frame nr 1000;
	BR->SetPeriodicSkipperSANS(3143);
	BR->SetEntrySkipperBNL(1000,800);
	in this case i just run the autosyncronizer BR->SetSkipBNLIfNoPulse(true); (this does not have a SANS version...) 
	*/


	BR->SetSkipBNLIfNoPulse(true);  // skips BNL if theres no pulseshape while the PulseShape file is not EOF;
	BR->SetAcceleratorFrequency(40.107);  // set to run at 40.107 Hz (is 40 by default too...) must be an int (use skippers to tune...)

	// Reading BNL datafile and defining weights for flat fielding (this takes up much disk space - but is fairly fast)
	//    This is not strictly needed as it can fairly easily be done using a BG TTree and the Data TTree inside the analysis
	//    But it does speed up and ease alot in this way - which belive to be better than HDD space...
	//    Also: It is rotation independent - which reduces the abount of human brain activity needed when analyzing :)
	cout<<"<<<Runner:Background statistics>>>"<<endl;
	//BR->ImportFlatField(CHAR(__DataPathName,"2012-12-NBL2/BNL105/BNLTEST_15_659.dat");


	//**********************************************************************
	// Loading and running phase.
	//**********************************************************************

	// loading the first batch of files:
	cout<<"<<<Runner: Processing first run>>>"<<endl;
	BR->LoadBatch(CHAR(__DataPathName,"2012-12-NBL2/SANS/SANS_5332_neutron_event.dat"),
		CHAR(__DataPathName,"2012-12-NBL2/BNL105/BNLTEST_15_668.dat"),
		CHAR(__DataPathName,"2012-12-NBL2/ProtonPulses/tds_out_Dec-07-2012-11.14.47.479.csv"));

	// adding another pulse file to que: (note that this will start right away - therefor no syncronizing)
	BR->AddPulseShapeFileToQue(CHAR(__DataPathName,"2012-12-NBL2/ProtonPulses/tds_out_Dec-08-2012-09.25.07.347.csv"));
	//BR->AddSANSFileToQue(char* QueThisFileName)
	//BR->AddBNL105FileToQue(char* QueThisFileName)


	// running the first batch, from starting frame (SANS,BNL,PulseShape) 0,0,0 and from starting time "x":
	BR->RunBatch(0,0,0,"Dec-07-2012 11:35:00.000");

	cout<<"<<<Runner: Processing secund run>>>"<<endl;
	// loading the next batch: (note that the next run will start simultaniously, therefor this is a syncronizing point)
	BR->LoadBatch(CHAR(__DataPathName,"2012-12-NBL2/SANS/SANS_5333_neutron_event.dat"),
		CHAR(__DataPathName,"2012-12-NBL2/BNL105/BNLTEST_15_669.dat"),
		CHAR(__DataPathName,"2012-12-NBL2/ProtonPulses/tds_out_Dec-08-2012-09.25.07.347.csv"));
	BR->AddPulseShapeFileToQue(CHAR(__DataPathName,"2012-12-NBL2/ProtonPulses/tds_out_Dec-08-2012-14.42.17.638.csv"));
	// running the next batch, from starting frame (SANS,BNL,PulseShape) 0,0,0 and from starting time "x":
	BR->RunBatch(5000000,5000000,5000000,"Dec-08-2012 13:09:00.000");


	// and now the pinhole part:
	cout<<"<<<Runner: Processing third run (pinhole part)>>>"<<endl;
	BR->LoadBatch(CHAR(__DataPathName,"2012-12-NBL2/SANS/SANS_5334_neutron_event.dat"),
		CHAR(__DataPathName,"2012-12-NBL2/BNL105/BNLTEST_15_670.dat"),
		CHAR(__DataPathName,"2012-12-NBL2/ProtonPulses/tds_out_Dec-08-2012-14.42.17.638.csv"));
	BR->AddPulseShapeFileToQue(CHAR(__DataPathName,"2012-12-NBL2/ProtonPulses/tds_out_Dec-09-2012-14.42.17.638.csv"));
	BR->RunBatch(10000000,10000000,10000000,"Dec-09-2012 13:47:00.000");

	cout<<"<<<Runner: Ending run>>>"<<endl;
	delete BR;
	cout<<"<<<Runner: Done!"<<endl;
}
void runner45mmSlabs(bool inTestmode){
	BatchRunner *BR;
	cout<<"<<<Runner: Initializing 45mm Slabs>>>"<<endl;
	BR=new BatchRunner("45mmSlabs.root");

	BR->setFloatFlatteningOfBNLPxl(true); 
	BR->setCalculateRotationMatrix(true); 
	BR->SetTestRunMode(inTestmode);
	BR->SetSkipBNLIfNoPulse(true); 
	BR->SetAcceleratorFrequency(40.107);


	cout<<"<<<Runner:Background statistics>>>"<<endl;
	//BR->ImportFlatField("CHAR(__DataPathName,"2012-12-NBL2/BNL105/BNLTEST_15_659.dat");
	cout<<"<<<Runner: Processing first run (pinhole part 1)>>>"<<endl;
	BR->LoadBatch(CHAR(__DataPathName,"2012-12-NBL2/SANS/SANS_5335_neutron_event.dat"),
		CHAR(__DataPathName,"2012-12-NBL2/BNL105/BNLTEST_15_671.dat"),
		CHAR(__DataPathName,"2012-12-NBL2/ProtonPulses/tds_out_Dec-10-2012-11.19.31.435.csv"));
	BR->RunBatch(0,0,0,"Dec-10-2012 11:54:00.000");

	cout<<"<<<Runner: Processing secund run (pinhole part 2)>>>"<<endl;
	BR->LoadBatch(CHAR(__DataPathName,"2012-12-NBL2/SANS/SANS_5336_neutron_event.dat"),
		CHAR(__DataPathName,"2012-12-NBL2/BNL105/BNLTEST_15_672.dat"),
		CHAR(__DataPathName,"2012-12-NBL2/ProtonPulses/tds_out_Dec-10-2012-11.19.31.435.csv"));
	BR->RunBatch(1500000,1500000,1500000,"Dec-10-2012 13:39:00.000");
  
	cout<<"<<<Runner: Processing third run >>>"<<endl;
	BR->LoadBatch(CHAR(__DataPathName,"2012-12-NBL2/SANS/SANS_5337_neutron_event.dat"),
		CHAR(__DataPathName,"2012-12-NBL2/BNL105/BNLTEST_15_673.dat"),
		CHAR(__DataPathName,"2012-12-NBL2/ProtonPulses/tds_out_Dec-10-2012-11.19.31.435.csv"));
	BR->RunBatch(5000000,5000000,5000000,"Dec-10-2012 15:06:00.000");
  
	cout<<"<<<Runner: Processing fourth run >>>"<<endl;
	BR->LoadBatch(CHAR(__DataPathName,"2012-12-NBL2/SANS/SANS_5338_neutron_event.dat"),
		CHAR(__DataPathName,"2012-12-NBL2/BNL105/BNLTEST_15_674.dat"),
		CHAR(__DataPathName,"2012-12-NBL2/ProtonPulses/tds_out_Dec-10-2012-11.19.31.435.csv"));
	BR->RunBatch(15000000,15000000,15000000,"Dec-11-2012 13:06:00.000");


	cout<<"<<<Runner: Processing fifth run >>>"<<endl;
	BR->LoadBatch(CHAR(__DataPathName,"2012-12-NBL2/SANS/SANS_5340_neutron_event.dat"),
		CHAR(__DataPathName,"2012-12-NBL2/BNL105/BNLTEST_15_676.dat"),
		CHAR(__DataPathName,"2012-12-NBL2/ProtonPulses/tds_out_Dec-10-2012-11.19.31.435.csv"));
	BR->RunBatch(25000000,25000000,25000000,"Dec-12-2012 09:39:00.000");

	cout<<"<<<Runner: Ending run>>>"<<endl;
	delete BR;
	cout<<"<<<Runner: Done!"<<endl;

}
void runnerBGRun(bool inTestmode){
	BatchRunner *BR;
	cout<<"<<<Runner: Initializing BG run>>>"<<endl;
	BR=new BatchRunner("BGRun.root");
	BR->setFloatFlatteningOfBNLPxl(true); 
	BR->setCalculateRotationMatrix(true); 
	BR->SetTestRunMode(inTestmode);
	BR->SetSkipBNLIfNoPulse(false); 
	BR->SetAcceleratorFrequency(40.107);

	cout<<"<<<Runner: Processing Background run>>>"<<endl;
	BR->LoadBatch("0",
		CHAR(__DataPathName,"2012-12-NBL2/BNL105/BNLTEST_15_659.dat"),
		"0");

	BR->RunBatch(0,0,0,"0");

	cout<<"<<<Runner: Ending run>>>"<<endl;
	delete BR;
	cout<<"<<<Runner: Done!"<<endl;

}
void runnerOld_polysi300(bool inTestmode){
	BatchRunner *BR;
	cout<<"<<<Runner: Initializing Poly Si 300>>>"<<endl;
	BR=new BatchRunner("polysi300.root");
	BR->setFloatFlatteningOfBNLPxl(true); 
	BR->setCalculateRotationMatrix(true); 
	BR->SetTestRunMode(inTestmode);
	BR->SetSkipBNLIfNoPulse(true); 
	BR->SetAcceleratorFrequency(40.107);

	cout<<"<<<Runner: processing run>>>"<<endl;
	BR->LoadBatch("0",
		CHAR(__DataPathName,"2012-06-NBL2/Datafiles/BNLTEST_15_645.dat"),
		"0");
	BR->AddBNL105FileToQue(CHAR(__DataPathName,"2012-06-NBL2/Datafiles/BNLTEST_15_646.dat"));
	BR->AddBNL105FileToQue(CHAR(__DataPathName,"2012-06-NBL2/Datafiles/BNLTEST_15_647.dat"));
	BR->AddBNL105FileToQue(CHAR(__DataPathName,"2012-06-NBL2/Datafiles/BNLTEST_15_648.dat"));
	BR->RunBatch(0,0,0,"0");

	cout<<"<<<Runner: Ending run>>>"<<endl;
	delete BR;
	cout<<"<<<Runner: Done!"<<endl;

}
void runnerOld_polysi80_1(bool inTestmode){
	BatchRunner *BR;
	cout<<"<<<Runner: Initializing Poly Si 80>>>"<<endl;
	BR=new BatchRunner("polysi80_1.root");
	BR->setFloatFlatteningOfBNLPxl(true); 
	BR->setCalculateRotationMatrix(true); 
	BR->SetTestRunMode(inTestmode);
	BR->SetSkipBNLIfNoPulse(true); 
	BR->SetAcceleratorFrequency(40.107);

	cout<<"<<<Runner: processing run>>>"<<endl;
	BR->LoadBatch("0",
		CHAR(__DataPathName,"2012-06-NBL2/Datafiles/BNLTEST_15_660.dat"),
		"0");
	BR->AddBNL105FileToQue(CHAR(__DataPathName,"2012-06-NBL2/Datafiles/BNLTEST_15_661.dat"));
	BR->AddBNL105FileToQue(CHAR(__DataPathName,"2012-06-NBL2/Datafiles/BNLTEST_15_662.dat"));
	BR->AddBNL105FileToQue(CHAR(__DataPathName,"2012-06-NBL2/Datafiles/BNLTEST_15_663.dat"));
	BR->AddBNL105FileToQue(CHAR(__DataPathName,"2012-06-NBL2/Datafiles/BNLTEST_15_664.dat"));	
	BR->RunBatch(0,0,0,"0");

	cout<<"<<<Runner: Ending run>>>"<<endl;
	delete BR;
	cout<<"<<<Runner: Done!"<<endl;
}
void runnerOld_polysi80_2(bool inTestmode){
	BatchRunner *BR;
	cout<<"<<<Runner: Initializing Poly Si 80 - 2>>>"<<endl;
	BR=new BatchRunner("polysi80_2.root");
	BR->setFloatFlatteningOfBNLPxl(true); 
	BR->setCalculateRotationMatrix(true); 
	BR->SetTestRunMode(inTestmode);
	BR->SetSkipBNLIfNoPulse(true); 
	BR->SetAcceleratorFrequency(40.107);

	cout<<"<<<Runner: processing run>>>"<<endl;
	BR->LoadBatch("0",
		CHAR(__DataPathName,"2012-06-NBL2/Datafiles/BNLTEST_15_660.dat"),
		"0");
	BR->RunBatch(0,0,0,"0");

	cout<<"<<<Runner: Ending run>>>"<<endl;
	delete BR;
	cout<<"<<<Runner: Done!"<<endl;
}
void runnerOld_polymono300(bool inTestmode){
	BatchRunner *BR;
	cout<<"<<<Runner: Initializing Poly Mono 300>>>"<<endl;
	BR=new BatchRunner("polymono300.root");
	BR->setFloatFlatteningOfBNLPxl(true); 
	BR->setCalculateRotationMatrix(true); 
	BR->SetTestRunMode(inTestmode);
	BR->SetSkipBNLIfNoPulse(true); 
	BR->SetAcceleratorFrequency(40.107);

	cout<<"<<<Runner: processing run>>>"<<endl;
	BR->LoadBatch("0",
		CHAR(__DataPathName,"2012-06-NBL2/Datafiles/BNLTEST_15_671.dat"),
		"0");
	BR->AddBNL105FileToQue(CHAR(__DataPathName,"2012-06-NBL2/Datafiles/BNLTEST_15_672.dat"));
	BR->RunBatch(0,0,0,"0");

	cout<<"<<<Runner: Ending run>>>"<<endl;
	delete BR;
	cout<<"<<<Runner: Done!"<<endl;
}
void runnerOld_polymono80(bool inTestmode){
	BatchRunner *BR;
	cout<<"<<<Runner: Initializing Poly mono 80>>>"<<endl;
	BR=new BatchRunner("polymone80.root");
	BR->setFloatFlatteningOfBNLPxl(true); 
	BR->setCalculateRotationMatrix(true); 
	BR->SetTestRunMode(inTestmode);
	BR->SetSkipBNLIfNoPulse(true); 
	BR->SetAcceleratorFrequency(40.107);

	cout<<"<<<Runner: processing run>>>"<<endl;
	BR->LoadBatch("0",
		CHAR(__DataPathName,"2012-06-NBL2/Datafiles/BNLTEST_15_676.dat"),
		"0");
	BR->AddBNL105FileToQue(CHAR(__DataPathName,"2012-06-NBL2/Datafiles/BNLTEST_15_676.dat"));
	BR->AddBNL105FileToQue(CHAR(__DataPathName,"2012-06-NBL2/Datafiles/BNLTEST_15_677.dat"));
	BR->AddBNL105FileToQue(CHAR(__DataPathName,"2012-06-NBL2/Datafiles/BNLTEST_15_678.dat"));
	BR->RunBatch(0,0,0,"0");

	cout<<"<<<Runner: Ending run>>>"<<endl;
	delete BR;
	cout<<"<<<Runner: Done!"<<endl;
};
void runnerOld_polyvac300_1(bool inTestmode){
	BatchRunner *BR;
	cout<<"<<<Runner: Initializing poly vac 300>>>"<<endl;
	BR=new BatchRunner("polyvac300_1.root");
	BR->setFloatFlatteningOfBNLPxl(true); 
	BR->setCalculateRotationMatrix(true); 
	BR->SetTestRunMode(inTestmode);
	BR->SetSkipBNLIfNoPulse(true); 
	BR->SetAcceleratorFrequency(40.107);

	cout<<"<<<Runner: processing run>>>"<<endl;
	BR->LoadBatch("0",
		CHAR(__DataPathName,"2012-06-NBL2/Datafiles/BNLTEST_15_681.dat"),
		"0");
	BR->AddBNL105FileToQue(CHAR(__DataPathName,"2012-06-NBL2/Datafiles/BNLTEST_15_682.dat"));

	BR->RunBatch(0,0,0,"0");

	cout<<"<<<Runner: Ending run>>>"<<endl;
	delete BR;
	cout<<"<<<Runner: Done!"<<endl;
}
void runnerOld_polyvac300_2(bool inTestmode){  
	BatchRunner *BR;
	cout<<"<<<Runner: Initializing poly vac 300 - 2>>>"<<endl;
	BR=new BatchRunner("polyvac300_2.root");
	BR->setFloatFlatteningOfBNLPxl(true); 
	BR->setCalculateRotationMatrix(true); 
	BR->SetTestRunMode(inTestmode);
	BR->SetSkipBNLIfNoPulse(true); 
	BR->SetAcceleratorFrequency(40.107);

	cout<<"<<<Runner: processing run>>>"<<endl;
	BR->LoadBatch("0",
		CHAR(__DataPathName,"2012-06-NBL2/Datafiles/BNLTEST_15_688.dat"),
		"0");
	BR->RunBatch(0,0,0,"0");

	cout<<"<<<Runner: Ending run>>>"<<endl;
	delete BR;
	cout<<"<<<Runner: Done!"<<endl;
}

// a trick for c++ compiling:
int main(){
	runner();
	return 0;
}
