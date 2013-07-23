#pragma once
#include <iostream>
//#include <stdio.h>
//#include <stdlib.h>
//#include "string.h"
#include "TTree.h"
#include "TFile.h"
#include <math.h>
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TRandom2.h"
using namespace std;

double pi=3.1415926535897932384626433832795;




class TSC {
public :
static unsigned int __GlobalUnameInteger;
static char __GlobalunameCharArray[1000];
static char __GlobalCHARCharArray[1000];
static char __GlobalTheMatrixCharArray[1000];
static bool __GlobalTheMatrixBoolArray[1000];
static double TSCTitlesize;
static double TSCLablesize;
static double TSCTitleOffset;
static double TSCyTitlesize;
static double TSCyLablesize;
static double TSCyTitleOffset;
static double TSCNamesize;
static double TSCNameOffset;
static double TSCLinewidth;
static double pi;
static double hbar; // eV*s
static double h; // eV*s
static double nmass; // eV
static double hovernmass;
static double eVperJ;
static double K_a; //1/mol
static double K_b;  //eV/K
static double q;
static double c;
static double deg2rad;  //rad/deg
static double rad2deg;  //deg/rad


static void DrawName(TH1 *h1,char *name=0){
   if(name==0){
     h1->SetTitle(h1->GetName());
   }
   else h1->SetTitle(name);
}
static int round(double x){
   return (int)(x+0.5);
}
static double FindMaximum(TH1 *a1,TH1 *a2=0,TH1 *a3=0,TH1 *a4=0,TH1 *a5=0,TH1 *a6=0,TH1 *a7=0,TH1 *a8=0,TH1 *a9=0){
	double maximum=a1->GetMaximum();
	if(a2!=0)maximum=max(maximum,a2->GetMaximum());
	if(a3!=0)maximum=max(maximum,a3->GetMaximum());
	if(a4!=0)maximum=max(maximum,a4->GetMaximum());
	if(a5!=0)maximum=max(maximum,a5->GetMaximum());
	if(a6!=0)maximum=max(maximum,a6->GetMaximum());
	if(a7!=0)maximum=max(maximum,a7->GetMaximum());
	if(a8!=0)maximum=max(maximum,a8->GetMaximum());
	if(a9!=0)maximum=max(maximum,a9->GetMaximum());
	return maximum;
}
static double FindMinimum(TH1 *a1,TH1 *a2=0,TH1 *a3=0,TH1 *a4=0,TH1 *a5=0,TH1 *a6=0,TH1 *a7=0,TH1 *a8=0,TH1 *a9=0){
	double minimum=a1->GetMinimum(1e-38);
	if(a2!=0)minimum=min(minimum,a2->GetMinimum(1e-38));
	if(a3!=0)minimum=min(minimum,a3->GetMinimum(1e-38));
	if(a4!=0)minimum=min(minimum,a4->GetMinimum(1e-38));
	if(a5!=0)minimum=min(minimum,a5->GetMinimum(1e-38));
	if(a6!=0)minimum=min(minimum,a6->GetMinimum(1e-38));
	if(a7!=0)minimum=min(minimum,a7->GetMinimum(1e-38));
	if(a8!=0)minimum=min(minimum,a8->GetMinimum(1e-38));
	if(a9!=0)minimum=min(minimum,a9->GetMinimum(1e-38));
	return minimum;
}
static char *CHAR(const char *c1,const char *c2=0, const char *c3=0, const char *c4=0, const char *c5=0, const char *c6=0, const char *c7=0){
sprintf(__GlobalCHARCharArray,"%s\0",c1);
if(c2)sprintf(__GlobalCHARCharArray,"%s%s\0",c1,c2);
if(c3)sprintf(__GlobalCHARCharArray,"%s%s\0",__GlobalCHARCharArray,c3);
if(c4)sprintf(__GlobalCHARCharArray,"%s%s\0",__GlobalCHARCharArray,c4);
if(c5)sprintf(__GlobalCHARCharArray,"%s%s\0",__GlobalCHARCharArray,c5);
if(c6)sprintf(__GlobalCHARCharArray,"%s%s\0",__GlobalCHARCharArray,c6);
if(c7)sprintf(__GlobalCHARCharArray,"%s%s\0",__GlobalCHARCharArray,c7);
sprintf(__GlobalCHARCharArray,"%s\0",__GlobalCHARCharArray);
return __GlobalCHARCharArray;
}
static char *uname(){
	__GlobalUnameInteger++;
	sprintf(__GlobalunameCharArray,"uname%d\0",__GlobalUnameInteger);
	return __GlobalunameCharArray;
}
static char *uname(int dummynumber){
	sprintf(__GlobalunameCharArray,"uname_dummy%d\0",dummynumber);
	return __GlobalunameCharArray;
}
static TLine *XLine(int col,double x,double ymin,double ymax){
	TLine *newline=new TLine(x,ymin,x,ymax);
	newline->SetLineColor(col);
	newline->SetLineWidth(3);
	return newline;
}
static TH1D *FastDraw1d(TTree *T,char* var,char*binning=0,char *cut=0,char *opt=0,int entries=100000000000,int start=0,TCanvas *dumcanvas=0){
	bool setDelete=true;
	if(dumcanvas){setDelete=false;
	} else {dumcanvas=new TCanvas();
	}
	char dummyhistname[100];
	sprintf(dummyhistname,"%s",uname());
	char setopt[100]={"\0"};
	if(opt)sprintf(setopt,"%s\0",opt);
	cout<<CHAR(var,">>",dummyhistname,"(",binning,")")<<endl;
	if(binning){T->Draw(CHAR(var,">>",dummyhistname,"(",binning,")"),cut,setopt,entries,start);
	} else {T->Draw(CHAR(var,">>",dummyhistname),cut,setopt,entries,start);
	}
	TH1D *hist=(TH1D*)(gDirectory->Get(dummyhistname)->Clone());
	((TH1D*)gDirectory->Get(dummyhistname))->Delete();
	if(setDelete)delete dumcanvas;
	return hist;
}
static TH2D *FastDraw2d(TTree *T,char* var,char*binning=0,char *cut=0,char *opt=0,int entries=100000000000,int start=0,TCanvas *dumcanvas=0){
	bool setDelete=true;
	if(dumcanvas){setDelete=false;
	} else {dumcanvas=new TCanvas();
	}
	char setopt[100];
	sprintf(setopt,"colz\0");
	if(opt)sprintf(setopt,"colz %s\0",opt);
	char dummyhistname[100];
	sprintf(dummyhistname,"%s\0",uname());
	
	cout<<CHAR(var,">>",dummyhistname,"(",binning,")")<<endl;
	if(binning){T->Draw(CHAR(var,">>",dummyhistname,"(",binning,")"),cut,setopt,entries,start);
	} else {T->Draw(CHAR(var,">>",dummyhistname),cut,setopt,entries,start);
	}
	cout<<"fetched"<<endl;
	TH2D *hist=(TH2D*)(gDirectory->Get(dummyhistname)->Clone());
	((TH2D*)gDirectory->Get(dummyhistname))->Delete();
	if(setDelete)delete dumcanvas;
	return hist;
}
static TTree *TreeLoader(TFile *f,const char *file,const char *path=""){
	cout<<"Loading: "<<CHAR(path,file,".root")<<" ..."<<endl;
	f= new TFile(CHAR(path,file,".root"));
	if(!f){cout<<"ERROR! <<<<<<<<<<<<<<<<<<<<LOADING FAILED<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"<<endl;return 0;}
	return (TTree*)f->Get("T")->Clone(uname());
}
static void NameHist(TH1D *h,int color,char *n=0, char *a1=0,char *a2=0,TLegend *l=0,char *hn=0,char *legopt=0){
	h->SetTitle(n);
	h->SetTitleSize(TSCNamesize);
	h->SetTitleSize(TSCNamesize);
	if(a1)h->GetXaxis()->SetTitle(a1);
	if(a2)h->GetYaxis()->SetTitle(a2);
	h->SetLineColor(color);
	if(l){
		if(legopt)l->AddEntry(h, hn, legopt);
		if(!legopt)l->AddEntry(h, hn, "l");
	}
	h->GetXaxis()->SetLabelSize(TSCLablesize);
	h->GetXaxis()->SetTitleOffset(TSCTitleOffset);
	h->GetXaxis()->SetTitleSize(TSCTitlesize);
	h->GetYaxis()->SetLabelSize(TSCyLablesize);
	h->GetYaxis()->SetTitleOffset(TSCyTitleOffset);
	h->GetYaxis()->SetTitleSize(TSCyTitlesize);
	h->SetLineWidth(TSCLinewidth);
	h->SetLineStyle(1);
}
static void NameHist(TH1D *h,int color,TLegend *l=0,char *hn=0,char *legopt=0){
	NameHist(h,color,0, 0,0,l,hn,legopt);
}
static void NameHist(TH2D *h,int color,char *n=0, char *a1=0,char *a2=0,char *a3=0,TLegend *l=0,char *hn=0,char *legopt=0){
	h->SetTitle(n);
	h->SetTitleSize(TSCNamesize);
	h->SetTitleSize(TSCNamesize);
	if(a1)h->GetXaxis()->SetTitle(a1);
	if(a2)h->GetYaxis()->SetTitle(a2);
	if(a3)h->GetZaxis()->SetTitle(a3);
	h->SetLineColor(color);
	if(l){
		if(legopt)l->AddEntry(h, hn, legopt);
		if(!legopt)l->AddEntry(h, hn, "l");		
	}
	h->GetXaxis()->SetLabelSize(TSCLablesize);
	h->GetXaxis()->SetTitleOffset(TSCTitleOffset);
	h->GetXaxis()->SetTitleSize(TSCTitlesize);
	h->GetYaxis()->SetLabelSize(TSCyLablesize);
	h->GetYaxis()->SetTitleOffset(TSCyTitleOffset);
	h->GetYaxis()->SetTitleSize(TSCyTitlesize);
}
static void NameHist(TH2D *h,int color,TLegend *l=0,char *hn=0,char *legopt=0){
	NameHist(h,color,0, 0,0,0,l,hn,legopt);
}
static void NameHist(TH3D *h,int color,char *n=0, char *a1=0,char *a2=0,char *a3=0,TLegend *l=0,char *hn=0,char *legopt=0){
	h->SetTitle(n);
	h->SetTitleSize(TSCNamesize);
	h->SetTitleSize(TSCNamesize);
	if(a1)h->GetXaxis()->SetTitle(a1);
	if(a2)h->GetYaxis()->SetTitle(a2);
	if(a3)h->GetZaxis()->SetTitle(a3);
	h->SetLineColor(color);
	if(l){
		if(legopt)l->AddEntry(h, hn, legopt);
		if(!legopt)l->AddEntry(h, hn, "l");		
	}
	h->GetXaxis()->SetLabelSize(TSCLablesize);
	h->GetXaxis()->SetTitleOffset(TSCTitleOffset);
	h->GetXaxis()->SetTitleSize(TSCTitlesize);
	h->GetYaxis()->SetLabelSize(TSCLablesize);
	h->GetYaxis()->SetTitleOffset(TSCTitleOffset);
	h->GetYaxis()->SetTitleSize(TSCTitlesize);
	h->GetZaxis()->SetLabelSize(TSCLablesize);
	h->GetZaxis()->SetTitleOffset(TSCTitleOffset);
	h->GetZaxis()->SetTitleSize(TSCTitlesize);


}
static void NameHist(TH3D *h,int color,TLegend *l=0,char *hn=0,char *legopt=0){
	NameHist(h,color,0, 0,0,0,l,hn,legopt);
}
static char *TheMatrix(int length){
	if(length>=1000)return 0;
	
	for(int i=0;i<length;i++){
		//cout<<gRandom->Rndm()<<endl;//<0.1)__GlobalTheMatrixCharArray[0]=56;
		if(__GlobalTheMatrixBoolArray[i])if(gRandom->Rndm()<=0.125)__GlobalTheMatrixBoolArray[i]=false;
		if(!__GlobalTheMatrixBoolArray[i])if(gRandom->Rndm()<=0.025)__GlobalTheMatrixBoolArray[i]=true;
		
		if(__GlobalTheMatrixBoolArray[i]){
			__GlobalTheMatrixCharArray[i]=round(gRandom->Rndm()*95+32);
		} else __GlobalTheMatrixCharArray[i]=32;
		//for(int i=1;i<30;i++){__GlobalTheMatrixCharArray[i]=i;cout<<i<<":"<<__GlobalTheMatrixCharArray[i]<<"|"<<endl;}
	}
	return __GlobalTheMatrixCharArray;
}
};
//initializers
double TSC::pi=3.1415926535897932384626433832795;
double TSC::c=299792458.; // m/s
double TSC::hbar=6.58211928e-16; // eV*s
double TSC::h=4.135667516e-15; // eV*s
double TSC::nmass=939565378.; // eV
double TSC::hovernmass=TSC::h/TSC::nmass*TSC::c*TSC::c; // s
double TSC::eVperJ=6.24150934e18;
double TSC::K_a=6.02214129; //1/mol
double TSC::K_b=8.6173324e-5;  //eV/K
double TSC::q=1.60217657e-19;  //eV/K
double TSC::deg2rad=TSC::pi/180.;  //rad/deg
double TSC::rad2deg=180./TSC::pi;  //deg/rad

char TSC::__GlobalCHARCharArray[1000];
char TSC::__GlobalunameCharArray[1000];
char TSC::__GlobalTheMatrixCharArray[1000];
bool TSC::__GlobalTheMatrixBoolArray[1000];
unsigned int TSC::__GlobalUnameInteger=0;
double TSC::TSCTitlesize=0.05;
double TSC::TSCLablesize=0.035;
double TSC::TSCTitleOffset=.85;
double TSC::TSCyTitlesize=0.05;
double TSC::TSCyLablesize=0.035;
double TSC::TSCyTitleOffset=.85;
double TSC::TSCNamesize=0.1;
double TSC::TSCNameOffset=1;
double TSC::TSCLinewidth=2;


class Neutron{
	public:
	bool SwichOff_DoYouKnowThatThisMightNotBeANeutron;
	bool UseAllParticles;
	double history,id,x,y,z,wx,wy,k,energy,time,weight;
	double i;
	double JGP;
	double JC;
	double IPT;
	double surface;
		
	double wz;
	double lambda;
	double l;
	double e;
	double loge;
	double logl;
	double ModeratorHight;
	double rsqr;
	double radtodeg;
	double r;
	double xsign;
	double ysign;
	double theta;
	double ycenter;
	double BigRsqr;
	double BigR;
	double phi;
	double wr;
	double wxsign;
	double wtheta0;
	double wthetaDiff;
	double wtheta;
	double wphi;
	int extractionindex;
	double wthetaDiff_center;
	double wtheta_center;
	void Reset(){
		extractionindex=-1;
		history=0;
		id=0;
		x=0;
		y=0;
		z=0;
		wx=0;
		wy=0;
		k=0;
		energy=0;
		time=0;
		weight=0;
		wthetaDiff_center=0;
		i=0;
		JGP=0;
		JC=0;
		IPT=0;
		surface=0;
		wz=0;
		lambda=0;
		l=0;
		e=0;
		loge=0;
		logl=0;
		ModeratorHight=0;
		rsqr=0;
		radtodeg=0;
		r=0;
		xsign=0;
		ysign=0;
		theta=0;
		ycenter=0;
		BigRsqr=0;
		BigR=0;
		phi=0;
		wr=0;
		wxsign=0;
		wtheta0=0;
		wthetaDiff=0;
		wtheta=0;
		wphi=0;
		wtheta_center=0;
	}
	void LoadNeutron(double inhistory,double inid,double inx,double iny,double inz,double inwx,double inwy,double ink,double inenergy,double intime,double inweight){

		id=inid;
		i=TSC::round(fabs(id/1E+6));
		JGP=TSC::round(i/200.0);
		JC=TSC::round(i/100.0) + 2*JGP;
		IPT=i-100*JC+200*JGP;

		if(IPT!=1){
			if(!UseAllParticles){
				if(SwichOff_DoYouKnowThatThisMightNotBeANeutron){
					return;
				} else {
					cout<<"WARNING: A non-neutron was used as a neutron, turn off this Warning:"<<endl;
					cout<<" - Neutron::SwichOff_DoYouKnowThatThisMightNotBeANeutron=true;"<<endl;
					cout<<" or set:"<<endl;cout<<" - Neutron::UseAllParticles=true"<<endl;
				}
			}
		}

		history=inhistory;
		x=inx;
		y=iny;
		z=inz;
		wx=inwx;
		wy=inwy;
		k=ink;
		energy=inenergy;
		time=intime;
		weight=inweight;


		surface=fabs(id)-1000000;
		//if(1-wx*wx-wy*wy<0)cout<<"ERROR TSC::SSW  ( 1-wx*wx-wy*wy<0 )"<<endl;
		wz=sqrt(1-wx*wx-wy*wy) * id/fabs(id);

		lambda=0.000286299/sqrt(energy);
		l=lambda;
		e=energy;
		loge=log10(energy);
		logl=log10(lambda);
		ModeratorHight=13.;
		rsqr=x*x+z*z;
		radtodeg=57.29577;
		r=sqrt(rsqr);
		xsign=(-2.*(x<0)+1.);
		ysign=(-2.*(y<0)+1.);
		theta=(x<0)*360.+xsign*(acos(z/r)*radtodeg);
		ycenter=ysign*ModeratorHight-y;
		BigRsqr=x*x+ycenter*ycenter+z*z;
		BigR=sqrt(x*x+ycenter*ycenter+z*z);
		phi=asin(ycenter/BigR)*radtodeg;
		wr=sqrt(wx*wx+wz*wz);
		wxsign=(-2.*(wx<0)+1.);
		wtheta0=(wx<0)*360.+wxsign*(acos(wz/wr)*radtodeg);
		wthetaDiff=(wtheta0-theta);
		wtheta=(wthetaDiff<-180)*360-(wthetaDiff>180)*360+wthetaDiff;
		wphi=asin(wy)*radtodeg;
		extractionindex=(theta-30)/5-(theta>90)-(theta>155)*10-(theta>265);
		wthetaDiff_center=(wtheta0-(extractionindex*5.+32.5+(extractionindex>11)*5.+(extractionindex>23)*50.+(extractionindex>35)*5.));
		//cout<<theta<<"  KK  "<<wtheta0<<" :: "<<wthetaDiff_center<<"  --  "<<wtheta<<endl;
		wtheta_center=(wthetaDiff_center<-180)*360-(wthetaDiff_center>180)*360+wthetaDiff_center;

		if(extractionindex<0||extractionindex>47){
			if(theta<=330.00001&&theta>=329){extractionindex=47;
			} else if(theta>=-0.00001&&theta<=1){extractionindex=0;
			}else {
			cout<<"strange extraction index : "<<extractionindex<<history<<"   -   "<<theta<<endl;
			}
		}
	}
	double GetExtractionIndexLimit(int inindex){
		return inindex*5.+30.+(inindex>11)*5.+(inindex>23)*50.+(inindex>35)*5.;
	}
	Neutron(){
		SwichOff_DoYouKnowThatThisMightNotBeANeutron=false;
		UseAllParticles=false;
		Reset();
	}
	
	void UseAllIPT(bool inputbool=true){UseAllParticles=inputbool;}// switch if you dont want the warning!
	Neutron(double inhistory,double inid,double inx,double iny,double inz,double inwx,double inwy,double ink,double inenergy,double intime,double inweight){
		SwichOff_DoYouKnowThatThisMightNotBeANeutron=false;
		UseAllParticles=false;
		LoadNeutron(inhistory,inid,inx,iny,inz,inwx,inwy,ink,inenergy,intime,inweight);
	}
	~Neutron(){}
};
void GetThoseNumbersGeometryLimits(double x,double z,double phi){
double XX=sin(phi*pi/180.)*2-x;
double ZZ=cos(phi*pi/180.)*2-z;
double R=sqrt(ZZ*ZZ+XX*XX);
cout<<"X: "<<.12*XX/R+x<<endl;
cout<<"Z: "<<.12*ZZ/R+z<<endl;
}
double ColdThermalThetaLimit(double *theta,double *par){
	double limx,limy,sign,cosang,rsqr,r;
	double R=2.;
	double Rsqr=4.;
	if(theta[0]<92.5){
			limx=.0891;
			limy=-.0128;
			sign=-1;
			rsqr=limx*limx+limy*limy;r=sqrt(rsqr);
			cosang=cos(theta[0]/180*pi-acos(limy/r));
		} else if(theta[0]<180){
			limx=.026;
			limy=-.0862;
			sign=-1;
			rsqr=limx*limx+limy*limy;r=sqrt(rsqr);
			cosang=cos(theta[0]/180*pi-acos(limy/r));
		} else if(theta[0]>267.5){
			limx=-.0891;
			limy=-.0128;
			sign=1;
			rsqr=limx*limx+limy*limy;r=sqrt(rsqr);
			cosang=cos(theta[0]/180*pi-(2*pi-acos(limy/r)));
		} else {
			limx=-.026;
			limy=-.0862;
			sign=1;
			rsqr=limx*limx+limy*limy;r=sqrt(rsqr);
			cosang=cos(theta[0]/180*pi-(2*pi-acos(limy/r)));
		}
		double over=Rsqr-r*R*cosang;
		double under=sqrt(rsqr+Rsqr-2*r*R*cosang)*R;
		return sign*acos(over/under)*180/pi;
	}
double ReflectorColdThetaLimit(double *theta,double *par){
	double limx,limy,sign,cosang,rsqr,r;
	double R=2.;
	double Rsqr=4.;
	if(theta[0]<92.5){
		limx=.0334;
		limy=.0836;
		sign=1;
		rsqr=limx*limx+limy*limy;r=sqrt(rsqr);
		cosang=cos(theta[0]/180*pi-acos(limy/r));
	} else if(theta[0]<180){
		limx=.0899;
		limy=.005;
		sign=1;
		rsqr=limx*limx+limy*limy;r=sqrt(rsqr);
		cosang=cos(theta[0]/180*pi-acos(limy/r));
	} else if(theta[0]>267.5){
		limx=-.0334;
		limy=.0836;
		sign=-1;
		rsqr=limx*limx+limy*limy;r=sqrt(rsqr);
		cosang=cos(theta[0]/180*pi-(2*pi-acos(limy/r)));
	} else {
		limx=-.0899;
		limy=.005;
		sign=-1;
		rsqr=limx*limx+limy*limy;r=sqrt(rsqr);
		cosang=cos(theta[0]/180*pi-(2*pi-acos(limy/r)));
	}
	double over=Rsqr-r*R*cosang;
	double under=sqrt(rsqr+Rsqr-2*r*R*cosang)*R;
	return sign*acos(over/under)*180/pi;
}
double ThermalReflectorThetaLimit(double *theta,double *par){
	double limx,limy,sign,cosang,rsqr,r;
	double R=2.;
	double Rsqr=4.;
	if(theta[0]<92.5){
		
		limx= .1603;
		limy= -.1055;
		sign=-1;
		rsqr=limx*limx+limy*limy;r=sqrt(rsqr);
		cosang=cos(theta[0]/180.*pi-acos(limy/r));
		//cosang=cos(x[0]/180*pi-asin(limx/r));
	} else if(theta[0]<180){
		limx= .0279;
		limy= -.1899;
		sign=-1;
		rsqr=limx*limx+limy*limy;r=sqrt(rsqr);
		cosang=cos(theta[0]/180.*pi-(2.*pi-acos(limy/r)));
	} else if(theta[0]>267.5){
		limx= -.1603;
		limy= -.1055;
		sign=1;
		rsqr=limx*limx+limy*limy;r=sqrt(rsqr);
		cosang=cos(theta[0]/180.*pi-(2.*pi-acos(limy/r)));
	} else {
		limx= -.0279;
		limy= -.1899;
		sign=1;
		rsqr=limx*limx+limy*limy;r=sqrt(rsqr);
		cosang=cos(theta[0]/180.*pi-acos(limy/r));
	}
	double over=Rsqr-r*R*cosang;
	double under=sqrt(rsqr+Rsqr-2*r*R*cosang)*R;
	return sign*acos(over/under)*180/pi;
}
double ThermalReflectorThetaUpperLimit(double *theta,double *par){
	double limx,limy,sign,cosang,rsqr,r;
	double R=2.;
	double Rsqr=4.;
	if(theta[0]<92.5){
		limx=.2801;
		limy= -.0986;
		sign=-1;
		rsqr=limx*limx+limy*limy;r=sqrt(rsqr);
		cosang=cos(theta[0]/180.*pi-acos(limy/r));
		//cosang=cos(x[0]/180*pi-asin(limx/r));
	} else if(theta[0]<180){
		limx= .02896;
		limy= -.29707;
		sign=-1;
		rsqr=limx*limx+limy*limy;r=sqrt(rsqr);
		cosang=cos(theta[0]/180.*pi-(2.*pi-acos(limy/r)));
	} else if(theta[0]>267.5){
		limx=-.2801;
		limy= -.0986;
		sign=1;
		rsqr=limx*limx+limy*limy;r=sqrt(rsqr);
		cosang=cos(theta[0]/180.*pi-(2.*pi-acos(limy/r)));
	} else {
		limx= -.02896;
		limy= -.29707;
		sign=1;
		rsqr=limx*limx+limy*limy;r=sqrt(rsqr);
		cosang=cos(theta[0]/180.*pi-acos(limy/r));
	}
	double over=Rsqr-r*R*cosang;
	double under=sqrt(rsqr+Rsqr-2*r*R*cosang)*R;
	return sign*acos(over/under)*180/pi;
}
double ReflectorColdThetaUpperLimit(double *theta,double *par){
	double limx,limy,sign,cosang,rsqr,r;
	double R=2.;
	double Rsqr=4.;
	if(theta[0]<92.5){
		limx=.0941;
		limy=.1871;
		sign=1;
		rsqr=limx*limx+limy*limy;r=sqrt(rsqr);
		cosang=cos(theta[0]/180*pi-acos(limy/r));
	} else if(theta[0]<180){
		limx=.2094;
		limy=.0063;
		sign=1;
		rsqr=limx*limx+limy*limy;r=sqrt(rsqr);
		cosang=cos(theta[0]/180*pi-acos(limy/r));
	} else if(theta[0]>267.5){
		limx=-.0941;
		limy=.1871;
		sign=-1;
		rsqr=limx*limx+limy*limy;r=sqrt(rsqr);
		cosang=cos(theta[0]/180*pi-(2*pi-acos(limy/r)));
	} else {
		limx=-.2094;
		limy=.0063;
		sign=-1;
		rsqr=limx*limx+limy*limy;r=sqrt(rsqr);
		cosang=cos(theta[0]/180*pi-(2*pi-acos(limy/r)));
	}
	double over=Rsqr-r*R*cosang;
	double under=sqrt(rsqr+Rsqr-2*r*R*cosang)*R;
	return sign*acos(over/under)*180/pi;
}




double AlphaColdTherm(double *x,double *par){
	double limx,limy,sign,cosang,rsqr,r;
	double R=2.;
	double Rsqr=4.;
		if(x[0]<92.5){
		limx=.0891;
		limy=-.0128;
		sign=-1;
		rsqr=limx*limx+limy*limy;r=sqrt(rsqr);
		cosang=cos(x[0]/180*pi-acos(limy/r));
	} else if(x[0]<180){
		limx=.026;
		limy=-.0862;
		sign=-1;
		rsqr=limx*limx+limy*limy;r=sqrt(rsqr);
		cosang=cos(x[0]/180*pi-acos(limy/r));
	} else if(x[0]>267.5){
		limx=-.0891;
		limy=-.0128;
		sign=1;
		rsqr=limx*limx+limy*limy;r=sqrt(rsqr);
		cosang=cos(x[0]/180*pi-(2*pi-acos(limy/r)));
	} else {
		limx=-.026;
		limy=-.0862;
		sign=1;
		rsqr=limx*limx+limy*limy;r=sqrt(rsqr);
		cosang=cos(x[0]/180*pi-(2*pi-acos(limy/r)));
	}
	double over=Rsqr-r*R*cosang;
	double under=sqrt(rsqr+Rsqr-2*r*R*cosang)*R;
	return sign*acos(over/under)*180/pi;

}
double AlphaRefCold(double *x,double *par){
	double limx,limy,sign,cosang,rsqr,r;
	double R=2.;
	double Rsqr=4.;
		if(x[0]<92.5){
		limx=.0334;
		limy=.0836;
		sign=1;
		rsqr=limx*limx+limy*limy;r=sqrt(rsqr);
		cosang=cos(x[0]/180*pi-acos(limy/r));
	} else if(x[0]<180){
		limx=.0899;
		limy=.005;
		sign=1;
		rsqr=limx*limx+limy*limy;r=sqrt(rsqr);
		cosang=cos(x[0]/180*pi-acos(limy/r));
	} else if(x[0]>267.5){
		limx=-.0334;
		limy=.0836;
		sign=-1;
		rsqr=limx*limx+limy*limy;r=sqrt(rsqr);
		cosang=cos(x[0]/180*pi-(2*pi-acos(limy/r)));
	} else {
		limx=-.0899;
		limy=.005;
		sign=-1;
		rsqr=limx*limx+limy*limy;r=sqrt(rsqr);
		cosang=cos(x[0]/180*pi-(2*pi-acos(limy/r)));
	}
	double over=Rsqr-r*R*cosang;
	double under=sqrt(rsqr+Rsqr-2*r*R*cosang)*R;
	return sign*acos(over/under)*180/pi;

}
double AlphaRefTherm(double *x,double *par){
	double limx,limy,sign,cosang,rsqr,r;
	double R=2.;
	double Rsqr=4.;

	if(x[0]<92.5){
		limx= .1603;
		limy= -.1055;
		sign=-1;
		rsqr=limx*limx+limy*limy;r=sqrt(rsqr);
		cosang=cos(x[0]/180.*pi-acos(limy/r));
		//cosang=cos(x[0]/180*pi-asin(limx/r));
	} else if(x[0]<180){
		limx= -.0279;
		limy= -.1899;
		sign=-1;
		rsqr=limx*limx+limy*limy;r=sqrt(rsqr);
		cosang=cos(x[0]/180.*pi-acos(limy/r));
	} else if(x[0]>267.5){
		limx= -.1603;
		limy= -.1055;
		sign=1;
		rsqr=limx*limx+limy*limy;r=sqrt(rsqr);
		cosang=cos(x[0]/180.*pi-(2.*pi-acos(limy/r)));
	} else {
		limx= .0279;
		limy= -.1899;
		sign=1;
		rsqr=limx*limx+limy*limy;r=sqrt(rsqr);
		cosang=cos(x[0]/180.*pi-(2.*pi-acos(limy/r)));
	}
	//cout<<acos(cosang)-x[0]/180*pi<<endl;
	double over=Rsqr-r*R*cosang;
	double under=sqrt(rsqr+Rsqr-2*r*R*cosang)*R;
	//cout<<cosang<<" :-: "<<over<<" :-: "<<under<<endl;
	return sign*acos(over/under)*180/pi;
	
	/*
	double xcosysin=limx*cos(x[0]/180.*3.pi)+limy*sin(x[0]/180.*pi);
	double R=2.;
	double Rsqr=4.;
	return sign*180/pi*acos((R-xcosysin)/sqrt(limx*limx+limy*limy-2*R*xcosysin+Rsqr));
*/
	//cout<<(2.-limy)<<"  ::  "<<sqrt(limx*limx+limy*limy)<<endl;
	//return x[0]-180/pi*acos((2.-limx)/sqrt(limx*limx+limy*limy));
	
}

class SSW{
	public:
	static void SetSSWAliases(TTree *T){
		T->SetAlias("i","TMath::Nint(TMath::Abs(id/1E+6))");
		T->SetAlias("JGP","-TMath::Nint(i/200.0)");
		T->SetAlias("JC","TMath::Nint(i/100.0) + 2*JGP");
		T->SetAlias("IPT","i-100*JC+200*JGP");
		T->SetAlias("surface","TMath::Abs(id)-1000000");
		T->SetAlias("wz","TMath::Sqrt(TMath::Max(0, 1-wx*wx-wy*wy)) * id/TMath::Abs(id)");

		T->SetAlias("lambda","0.000286299/sqrt(energy*1.)");
		T->SetAlias("l","0.000286299/sqrt(energy*1.)");
		T->SetAlias("e","energy*1.");
		T->SetAlias("loge","log10(energy*1.)");
		T->SetAlias("logl","log10(lambda*1.)");
		T->SetAlias("ModeratorHight","13.*1");
		T->SetAlias("rsqr","x*x+z*z");
		T->SetAlias("radtodeg","1.*57.29577");
		T->SetAlias("r","sqrt(rsqr*1.)");
		T->SetAlias("xsign","(-2.*(x<0)+1.)");
		T->SetAlias("ysign","(-2.*(y<0)+1.)");
		T->SetAlias("theta","(x<0)*360.+xsign*(acos(z/r)*radtodeg)");
		T->SetAlias("ycenter","ysign*ModeratorHight-y");
		T->SetAlias("BigRsqr","x*x+ycenter*ycenter+z*z");
		T->SetAlias("BigR","sqrt(x*x+ycenter*ycenter+z*z)");
		T->SetAlias("phi","asin(ycenter/BigR)*radtodeg");
		T->SetAlias("wr","sqrt(wx*wx+wz*wz)");
		T->SetAlias("wxsign","(-2.*(wx<0)+1.)");
		T->SetAlias("wtheta0","(wx<0)*360.+wxsign*(acos(wz/wr)*radtodeg)");
		T->SetAlias("wthetaDiff","(wtheta0-theta)");
		T->SetAlias("wtheta","(wthetaDiff<-180)*360-(wthetaDiff>180)*360+wthetaDiff");
		T->SetAlias("wphi","asin(wy)*radtodeg");
	}
	static double history,id,x,y,z,wx,wy,k,energy,time,weight;

	static double i;
	static double JGP;
	static double JC;
	static double IPT;
	static double surface;
		
	static double wz;
	static double lambda;
	static double l;
	static double e;
	static double loge;
	static double logl;
	static double ModeratorHight;
	static double rsqr;
	static double radtodeg;
	static double r;
	static double xsign;
	static double ysign;
	static double theta;
	static double ycenter;
	static double BigRsqr;
	static double BigR;
	static double phi;
	static double wr;
	static double wxsign;
	static double wtheta0;
	static double wthetaDiff;
	static double wtheta;
	static double wphi;
	static void LoadNeutron(double inhistory,double inid,double inx,double iny,double inz,double inwx,double inwy,double ink,double inenergy,double intime,double inweight){
		history=inhistory;
		id=inid;
		x=inx;
		y=iny;
		z=inz;
		wx=inwx;
		wy=inwy;
		k=ink;
		energy=inenergy;
		time=intime;
		weight=inweight;

		i=TSC::round(fabs(id/1E+6));
		JGP=TSC::round(i/200.0);
		JC=TSC::round(i/100.0) + 2*JGP;
		IPT=i-100*JC+200*JGP;
		surface=fabs(id)-1000000;
		//if(1-wx*wx-wy*wy<0)cout<<"ERROR TSC::SSW  ( 1-wx*wx-wy*wy<0 )"<<endl;
		wz=sqrt(1-wx*wx-wy*wy) * id/fabs(id);

		lambda=0.000286299/sqrt(energy*1.);
		l=0.000286299/sqrt(energy*1.);
		e=energy*1.;
		loge=log10(energy*1.);
		logl=log10(lambda*1.);
		ModeratorHight=13.*1;
		rsqr=x*x+z*z;
		radtodeg=1.*57.29577;
		r=sqrt(rsqr);
		xsign=(-2.*(x<0)+1.);
		ysign=(-2.*(y<0)+1.);
		theta=(x<0)*360.+xsign*(acos(z/r)*radtodeg);
		ycenter=ysign*ModeratorHight-y;
		BigRsqr=x*x+ycenter*ycenter+z*z;
		BigR=sqrt(x*x+ycenter*ycenter+z*z);
		phi=asin(ycenter/BigR)*radtodeg;
		wr=sqrt(wx*wx+wz*wz);
		wxsign=(-2.*(wx<0)+1.);
		wtheta0=(wx<0)*360.+wxsign*(acos(wz/wr)*radtodeg);
		wthetaDiff=(wtheta0-theta);
		wtheta=(wthetaDiff<-180)*360-(wthetaDiff>180)*360+wthetaDiff;
		wphi=asin(wy)*radtodeg;
	}
};


double SSW::history=0;
double SSW::id=0;
double SSW::x=0;
double SSW::y=0;
double SSW::z=0;
double SSW::wx=0;
double SSW::wy=0;
double SSW::k=0;
double SSW::energy=0;
double SSW::time=0;
double SSW::weight=0;

double SSW::i=0;
double SSW::JGP=0;
double SSW::JC=0;
double SSW::IPT=0;
double SSW::surface=0;
double SSW::wz=0;
double SSW::lambda=0;
double SSW::l=0;
double SSW::e=0;
double SSW::loge=0;
double SSW::logl=0;
double SSW::ModeratorHight=0;
double SSW::rsqr=0;
double SSW::radtodeg=0;
double SSW::r=0;
double SSW::xsign=0;
double SSW::ysign=0;
double SSW::theta=0;
double SSW::ycenter=0;
double SSW::BigRsqr=0;
double SSW::BigR=0;
double SSW::phi=0;
double SSW::wr=0;
double SSW::wxsign=0;
double SSW::wtheta0=0;
double SSW::wthetaDiff=0;
double SSW::wtheta=0;
double SSW::wphi=0;


