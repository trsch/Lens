#pragma once
#include <iostream>
#include "TF1.h"
#include "TF2.h"
#include "TF3.h"
#include <cmath>
#include "TSCBasicInlineFunctions.h"
using namespace std;


double V2L(double speed){  //ToF in micro secunds
	return TSC::hovernmass/speed*1e10; // Å (e10*e-6)
}
double E2T(double E){ //MeV to Kelvin
	return 11604.5e6*E;
}
double T2E(double T){ //MeV to Kelvin
	return T/11604.5e6;
}
double E2L(double E){ // MeV to Angstrom
	return 0.000286299/sqrt(E);
}
double dldE(double E){
	return 0.0001431495/(sqrt(E)*E);
}
double L2E(double l){ // Angstrom to MeV
	return 0.000286299*0.000286299/(l*l);
}
double dEdl(double l){
	return 0.000572598/(l*l*l);
}
double T2L(double T){
	return E2L(T2E(T));
}
double L2T(double T){
	return E2T(L2E(T));
}
double ToF2L(double *ToF,double *par=0){  //ToF in micro secunds
// par[0] = dist moderator to detector
	if(!par){
		par=new double [1];
		par[0]=4.;
	}
	return V2L(par[0]/ToF[0]*1e-6); // Å
}

/*
double erf(double x)
{
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;

    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x);

    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

    return sign*y;
}
*/
double Maxwellian(double *l,double *T){
	if(T[0]!=0&&l[0]>0){
		double aOlsqr=949./(T[0]*l[0]*l[0]);
		return 2.*aOlsqr*aOlsqr/l[0]*exp(-aOlsqr);
	}
	return 0;
}
double ScaledMaxwellian(double *l,double *par){
	if(par[1]!=0&&l[0]>0){
		double aOlsqr=949./(par[1]*l[0]*l[0]);
		return par[0]*2.*aOlsqr*aOlsqr/l[0]*exp(-aOlsqr);
	}
	return 0;
}
double Joining(double *x,double *par){ // speed,joinvalue
	return 1/(1+exp(par[1]*(x[0]-par[0])));
}
double Slowingdown(double *l,double *par){
	if(1+exp(par[0]*l[0]-par[1])!=0&&l[0]>0){
	return 1/(1+exp(par[0]*l[0]-par[1]))/(l[0]*l[0]*l[0]);
	}
	return 0;
}
double ScaledSlowingdown(double *l,double *par){
	if(1+exp(par[1]*l[0]-par[2])!=0&&l[0]>0){
		return par[0]/(1+exp(par[1]*l[0]-par[2]))/l[0]*l[0]*l[0];
	}
	return 0;
}
/*
double BispectralFolding(double *l,double *par){ //erf std, cutoff center,decay
	return 1-(+.5+.5*erf(par[0]*(l[0]-par[1])));
	
}
double BispectralFolding_thermal(double *l,double *par){ //erf std, cutoff center,decay
	return 0.9*BispectralFolding(l,par)+.1*exp(-par[2]*l[0]);
//	)+.1*exp(-par[2]*l[0]);
}

double BispectralFolding_cold(double *l,double *par){ //erf std, cutoff center,decay
	//return 0.96*(1-BispectralFolding(l,par));
	return 0.97*(1-BispectralFolding(l,par))-.1*(1-BispectralFolding(l,par))*exp(-par[2]*l[0]);
}
*/
double Joinedexponential(double *l, double *par){
	return (1-1/
		pow((1+exp(par[1]*(l[0]-par[0]))),2.5)
		)*exp(-par[2]*(l[0]-par[0]));
}

double ParaSwitch(double *e,double *par){ // cut, beta, delta
	if(e[0]>=par[0]){
		double x=par[1]*(e[0]-par[0]);
		return 1-exp(-x)*(1+x+.5*x*x);
	}
	return 0;
}
double ParaCut(double *e,double *par){
	if(e[0]>=par[0]){
		double x=par[1]*(e[0]-par[0]);
		return 1+exp(-x)*(1+x+.5*x*x);
	}
	return 0;
}

double ShortPulse(double time, double *par){	//2 par: par[0]=tau, par[1]=eta (pulseshape-parameter)
	if(par[0]!=0&&(1-par[1])!=0){
		double tauinvers=1/par[0];
		return (exp(-time*tauinvers)-exp(-par[1]*tauinvers*time))*par[1]/(par[1]-1)*tauinvers;
	}
	if(par[0]==0)return 0;
	cout<<"undefined-shortpulse"<<endl;
	return 0;
}
double LongPulseIExp(double t,double *par){ //[0]=tau, [1]=d;
	if(t<=par[1])return par[0]*(1-exp(-t/par[0]));
	return par[0]*(exp(par[1]/par[0]-t/par[0])-exp(-t/par[0]));
}
double LongPulse(double *t,double *par){ //[0]=tau,[1]=eta,[2]=d;
	if(par[0]!=0&&par[1]!=0&&par[1]!=1){
		double passpar1[2]={par[0],par[2]};
		double passpar2[2]={par[0]/par[1],par[2]};
		return (LongPulseIExp(t[0],passpar1)-LongPulseIExp(t[0],passpar2))*par[1]/(par[1]-1)/par[0]/par[2];
	}
	//if(par[1]==0)cout<<"strange"<<endl;
	//cout<<"undefined-Longpulse"<<endl;
	return 0;
}
double DoubleLongPulse(double *t,double *par){ //[0]=tau_pulse,tau_Reflector,[1]=eta,[2]=d;
	double otherpar[3]={par[0],par[2],par[3]};
	return (LongPulse(t,otherpar)+LongPulse(t,&par[1]))*.5; // Semi normalized!!!
}


// for long pulses:
double Decoupled(double *lt, double *par){ //IntencityM ( temperture : tau, eta, d) : IntencityS (cutoffspeed, cutoffPosition : tau, eta, d) 
	return par[0]*Maxwellian(&lt[0],&par[1])*LongPulse(&lt[1],&par[2])+par[5]*Joining(&lt[0],&par[6])*LongPulse(&lt[1],&par[7]);
}
double Coupled(double *lt, double *par){ //IntencityM ( temperture : tau_pulse,tau_reflector, eta, d) : IntencityS (cutoffspeed, cutoffPosition : tau, eta, d) 
	
	return par[0]*Maxwellian(&lt[0],&par[1])*(DoubleLongPulse(&lt[1],&par[2]))+par[6]*Joining(&lt[0],&par[7])*LongPulse(&lt[1],&par[8]);
}
double Coupled_l(double *l, double *par){ //IntencityM ( temperture : tau_pulse,tau_reflector, eta, d) : IntencityS (cutoffspeed, cutoffPosition : tau, eta, d) 
	double lt[2]={l[0],par[5]};
	return Coupled(lt,par);
}
double Decoupled_l(double *l, double *par){ //IntencityM ( temperture : tau, eta, d) : IntencityS (cutoffspeed, cutoffPosition : tau, eta, d) 
	double lt[2]={l[0],par[4]};
	return Decoupled(lt,par);
}


double MaxMaxSlow(double *l, double *par){
	return par[0]*Maxwellian(l,&par[1])+par[2]*Maxwellian(l,&par[3])+par[4]/((1+exp(par[5]*l[0]-par[6]))*l[0])+par[7]/((1+exp(par[8]*l[0]-par[9]))*l[0]);
}



double LeakingParaH2(double *l, double *par){
	return (1-1/
		pow((1+exp(par[1]*(l[0]-par[0]))),par[2])
		)*exp(-par[3]*l[0]);
}


double BispectralParaH2_ver2(double *l, double *par){
	return par[0]*Maxwellian(l,&par[1])*pow(l[0],par[2])  // Maxwellian
		+par[3]*Joinedexponential(l,&par[4])          // ParaHydrogen model
		+par[7]/((1+exp(par[8]*(l[0]-par[9])))*l[0])  // slowing Down
		+par[10]/((1+exp(par[11]*(l[0]-par[12])))*l[0]); // Slowing Down

}


double BispectralParaH2(double *l, double *par){
	double partherm[2]={par[7],par[9]};
	return par[0]*Maxwellian(l,&par[1])*pow(l[0],par[2])+par[3]*Joinedexponential(l,&par[4])
		+par[7]/((1+exp(par[8]*(l[0]-par[9])))*l[0])+par[10]/((1+exp(par[11]*(l[0]-par[12])))*l[0]);

}

double TSC_BispectralParaH2Coupled(double *l, double *par){
	double partherm[2]={par[7],par[9]};
	return par[0]*Maxwellian(l,&par[1])*pow(l[0],par[10])+par[2]*Joinedexponential(l,&par[3])+par[6]/(l[0])*(.55*Joining(l,partherm)+.45*Joining(l,&par[8]));
}
double TSC_BispectralParaH2Coupled_Slowingdown(double *l, double *par){
	double partherm[2]={par[7],par[9]};
	return par[6]/(l[0])*(.55*Joining(l,partherm)+.45*Joining(l,&par[8]));
}
double TSC_BispectralParaH2Coupled_Maxwellian(double *l, double *par){
	double partherm[2]={par[7],par[9]};
	return par[0]*Maxwellian(l,&par[1])*pow(l[0],par[10])+par[6]/(l[0])*.55*Joining(l,partherm);
}
double TSC_BispectralParaH2Coupled_Joinedexponential(double *l, double *par){
	return par[2]*Joinedexponential(l,&par[3])+par[6]/(l[0])*.45*Joining(l,&par[8]);
}

double BispectralParaH2Coupled(double *l, double *par){
	return par[0]*Maxwellian(l,&par[1])+par[2]*Joinedexponential(l,&par[3])+par[6]*Slowingdown(l,&par[7]);
}
double BispectralParaH2Coupled_Slowingdown(double *l, double *par){
	return par[6]*Slowingdown(l,&par[7]);
}
double BispectralParaH2Coupled_Maxwellian(double *l, double *par){
	return par[0]*Maxwellian(l,&par[1])+.5*par[6]*Slowingdown(l,&par[7]);
}
double BispectralParaH2Coupled_Joinedexponential(double *l, double *par){
	return par[2]*Joinedexponential(l,&par[3])+.5*par[6]*Slowingdown(l,&par[7]);
}
double ParaH2Coupled(double *l,double *par){ //Im : Temperature :: Is : cutvalue,cutspeed : absorption : cutvalue,cutspeed,bumpsize : epithermalScale,EpitermalPower
	//          Maxwellian
	double e[1]={L2E(l[0])*1e6};
	if(l[0]!=0&&e[0]!=0)
	return par[0]*Maxwellian(l,&par[1])
		//      Switching to energy:
		+ dEdl(l[0])
		//      ParaH2 Leakage cutoff
		* par[2]*ParaCut(e,&par[3])*exp(par[5]/sqrt(e[0])) 
		//           Epithermal part
		* dEdl(l[0])*ParaSwitch(e,&par[6])*pow(par[9]/e[0],1-par[10])
		// * par[2]  // this must be an error...
		;
	return 0;
}




double BGSpectrum(double *l,double *par){
// [0]normalizer;[1]time;[2]tau;[3]n;[4]d;[5]bg;
	double passpar[3]={par[2],par[3],par[4]};
	double passvar[1]={par[1]};
		
		return par[0]*LongPulse(passvar,passpar)*Slowingdown(l,&par[5]);
}
double BGSpectrum_Time(double *t,double *par){
	double l[1]={par[0]};
	par[0]=t[0];
	return BGSpectrum(t,par);
}
double SingleSpectrum(double *l,double *par){
	//0normalizer;1time;2Tempreture;3tau;4n;5d
	double passpar[3]={par[3],par[4],par[5]};
	double maxpasspar[1]={par[2]};
	double passvar[1]={par[1]};
	return par[0]*Maxwellian(l,maxpasspar)*LongPulse(passvar,passpar);
}
double SingleSpectrum_Time(double *t,double *par){
	double l[1]={par[0]};
	par[0]=t[0];
	return SingleSpectrum(t,par);
}


double Maxwellian_trunk(double *l,double *T){
	if(T[0]!=0&&l[0]>0){
		double a=949./T[0];
		double lsqr=l[0]*l[0];
		return 2.*a*a/lsqr/lsqr/l[0]*exp(-a/lsqr/lsqr);
	}
	return 0;
}
double SpectrumESS(double *l,double *par){
	if(l[0]!=0)return par[0]*Slowingdown(l,&par[1])+par[2]*Maxwellian(l,&par[3])
		+par[4]*Maxwellian(l,&par[5]);
		//+par[4]*Maxwellian_trunc(l,&par[5]);		//*(1+(par[6]*(l[0]-par[7])*(l[0]-par[7])*exp(-(l[0]-par[7])*(l[0]-par[7])/par[8]/par[8])));
	return 0;
}
double SpectrumESS_trunk(double *l,double *par){

	if(l[0]!=0)return par[0]*Slowingdown(l,&par[1])+par[2]*Maxwellian(l,&par[3])/l[0]//*pow(l[0],par[6])
		//+par[4]*Maxwellian(l,&par[5]);
		+par[4]*Maxwellian_trunk(l,&par[5])*pow(l[0],par[7]);				//*(1+(par[6]*(l[0]-par[7])*(l[0]-par[7])*exp(-(l[0]-par[7])*(l[0]-par[7])/par[8]/par[8])));
	return 0;
}
