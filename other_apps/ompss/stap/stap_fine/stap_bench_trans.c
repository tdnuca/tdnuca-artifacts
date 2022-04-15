/******************************************************************************/
// Source File: stap_bench_trans.c
// (Copyright © THALES 2010 All rights reserved) 
// THE PROGRAM IS PROVIDED ON AN "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS
// OF ANY KIND, EITHER EXPRESS OR IMPLIED INCLUDING, WITHOUT LIMITATION, 
// ANY WARRANTIES ON ITS, NON-INFRINGEMENT, MERCHANTABILITY, SECURED, 
// INNOVATIVE OR RELEVANT NATURE, FITNESS FOR A PARTICULAR PURPOSE OR 
// COMPATIBILITY WITH ANY EQUIPMENT OR SOFTWARE.
//
// Authors: 	Eric Lenormand (eric.lenormand@thalesgroup.com),
//		Remi Barrere (remi.barrere@thalesgroup.com),
//	   	Sami Yehia (sami.yehia@thalesgroup.com)
//	   	Teodora Petrisor (claudia-teodora.petrisor@thalesgroup.com)
//
// Date: Jan,8,2010
//
/******************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include <sys/time.h>
#include <time.h>

/**********************************
 * System parameters
 ********************************* */

#define Pi 3.141592653589793238
#define ATL_Lambda 0.03
#define ATL_Tr 0.0006  // 5.0E-4

#define ATL_Nrg 95			// number of range gates
#define ATL_tf 16			// filter length (16-tap)
#define ATL_Nrgf   ATL_Nrg-ATL_tf+1	
#define ATL_ntt 5			// number of pulses in sub-window
#define ATL_nsa 5			// number of antennas (=Nant)	
#define ATL_nrec 32			// number of pulses in a burst
#define ATL_tfac_L_V 5
#define ATL_tfac_L_RG 7
#define ATL_tfac_K   10.0
#define ATL_tfac_AVPWK  10000.0

#define ATL_Nth 20			// number of "Theta" points assumptions
#define ATL_Nv 20			// number of "velocity" points assumptions

#define ATL_Nrg_enlarged  ATL_Nrgf +  ATL_tfac_L_RG -1   // 86
#define ATL_Nv_enlarged ATL_Nv +  ATL_tfac_L_V -1   // 24

#define ATL_CarrierSpeed 200
#define ATL_BeamAngle 90
#define ATL_BeamWidth 5
#define ATL_SubarraySpacing 0.2
#define ATL_TargetDistance 12000
#define ATL_Nplots 100
#define ATL_rg_size  15




/************************************************************************/
// Define complex number
typedef struct {
     float re;
     float im;
} Cplfloat;
/************************************************************************/



Cplfloat filtre[ATL_tf];
Cplfloat Steervect[ATL_Nth][ATL_Nv][ATL_ntt*ATL_nsa];
Cplfloat tab_coh[ATL_Nv][ATL_nrec-ATL_ntt+1];
Cplfloat out_pulse[ATL_nsa][ATL_nrec][ATL_Nrgf];
Cplfloat in_apply[ATL_Nrgf][ATL_nrec][ATL_nsa];
Cplfloat in_cov[ATL_Nrgf][ATL_nsa][ATL_nrec];
Cplfloat out_avcov[ATL_Nrgf][ATL_ntt*ATL_nsa][ATL_ntt*ATL_nsa];
Cplfloat mat_inv[ATL_Nrgf][ATL_ntt*ATL_nsa][ATL_ntt*ATL_nsa];
Cplfloat Fil[ATL_Nrgf][ATL_Nth][ATL_Nv][ATL_ntt*ATL_nsa];
Cplfloat out_apply[ATL_Nrgf][ATL_Nth][ATL_Nv][ATL_nrec-ATL_ntt+1];
Cplfloat out_int_dop[ATL_Nrgf][ATL_Nth][ATL_Nv];
Cplfloat in_avpow[ATL_Nth][ATL_Nrgf][ATL_Nv];
Cplfloat out_avpow[ATL_Nth];
Cplfloat out_facedges[ATL_Nth][ATL_Nrg_enlarged][ATL_Nv_enlarged];
Cplfloat out_face1[ATL_Nth][ATL_Nrgf][ATL_Nv];
Cplfloat in_maxpow[ATL_Nrgf][ATL_Nv][ATL_Nth];
Cplfloat out_maxpow[ATL_Nrgf][ATL_Nv][ATL_Nth];
Cplfloat in_correctionV[ATL_Nrgf][ATL_Nth][ATL_Nv];
Cplfloat out_correctionV[ATL_Nrgf][ATL_Nth][ATL_Nv];



/************************************************************************/

double Reprodconj(double XR, double XI,double YR, double YI){
     return( XR*YR +XI*YI);
}

double ModSquare (Cplfloat Z){
     return(Reprodconj((double)Z.re, (double)Z.im,(double) Z.re,(double)Z.im));
}



/************************************************************************/
// type conversion functions


int Vfloat2vint(float V,int Nv, float lambda, float Tr){
	
     float Vmax = lambda/(4.*Tr);	
     float Vgate = 2.*Vmax/(float)Nv;
     int v;
     if(V>=0.){
	  v = (int)(1.E-5+ V/Vgate);
     }else{
	  v = (int)(1.E-5+ (V+2*Vmax)/Vgate);
     }
     if(v<0) v=0;
     if(v>Nv) v=Nv;
     return(v);
}


/************************************************************************/

float vint2Vfloat (int v,int Nv, float lambda, float Tr)
{
     float Vmax = lambda/(4.*Tr);
     float V;
     V= Vmax * (2.*v+1.)/(1.*Nv);
     if(v>=Nv/2){
	  V=V-2.*Vmax;
     }
     return(V);
}

/************************************************************************/


float vint2Vfloat_1 (int v,int Nv, float Vmax)
{
     float V;
     V= Vmax * (2.*v+1.)/(1.*Nv);
     if(v>=Nv/2){
	  V=V-2.*Vmax;
     }
     return(V);
}

/************************************************************************/

float thint2THfloat(int th, int Nth, float Beamwidth){
     float Beamwidth_rad = Beamwidth*Pi/180.;
     return(Beamwidth_rad*(-.5 + (2.*th+1.)/(2.*Nth)));
}

/************************************************************************/

int THfloat2thint(float TH , int Nth, float Beamwidth){
     float Beamwidth_rad = Beamwidth*Pi/180.;
     float Thgate = Beamwidth_rad/(float)Nth;
     int th; 
     th = (int) ( .5 +(TH +Beamwidth_rad/2.)/ Thgate);
     if(th<0) th=0;
     if (th>Nth) th=Nth;
     return(th);
}	


/************************************************************************/



void RecordPlots( char* name, int NbPlots, Cplfloat Plots[NbPlots], float rg_size, 
		  int Nv,float Vmax,int Nth, float Beamwidth){

     FILE *fp;
     int pl, rg, v, th;
     float Distance, Velocity, Angle;
     fp= fopen(name, "w");
     if(fp){
	  for (pl=0; pl<NbPlots; pl++) {
	       Cplfloat PLOT = Plots[pl];
	       if (PLOT.im!=0.) {
	            rg = (int)(PLOT.re/10000. +.1);
	            Distance = rg_size * rg;
	            v = (int)( (PLOT.re - 10000.*rg)/100. +.1);
	            Velocity = vint2Vfloat_1 (v, Nv,Vmax); 
	            th = (int)(PLOT.re - 10000.*rg -100.*v +.1);
	            Angle = thint2THfloat(th, Nth, Beamwidth);
	            Angle = 90. - Angle*(180./Pi);
	            fprintf(fp, " %2d:  Distance  %e  (%2d)     Velocity %e  ( %2d)     Angle %e ( %2d) _  Peak  %e\n", pl, Distance,rg,Velocity, v,Angle, th, PLOT.im);
	       }
	  }
	  fclose(fp);
     }  	   
}	



/************************************************************************/
// Read inputs from file InputStimuli.txt 

int HasNextBurst(FILE *fp, float *Tr ){
     int status =0;	
     char line[256];
     int chara;
     int i;
	
     //First significant line should contain pulse repetition period, "Tr=",
     // and the float value of Tr
	
     int HasEqual=0; 
     int EndOfFile =0;
     while( (!HasEqual)&&(!EndOfFile)){ 
	  HasEqual= 0;
	  i=0;
	  while(i<255) {
	       line[i]=fgetc(fp);
	       chara= (int)(line[i]);
	       if((chara<0) || (chara>=255)){
		    EndOfFile=1;
		    break;
	       }
	       if(line[i]=='=')
		    HasEqual=1;
	       if(line[i]=='\n') break;
	       if(line[i]!=' ') i++;
	  }
	  line[i+1]='\0';
     }
     if(HasEqual){
	  status=1;
	  strtok(line,"=");
	  sscanf(strtok(NULL,"\n"), "%f", Tr);
     }
     return status;
}


int GetNextBurst(FILE *fp, int Nrg, int Nsa, int Npul, Cplfloat Tab[Npul][Nsa][Nrg]) {


     float X;
     int pul, sa,rg;


     for(pul=0; pul<Npul;pul++) {
	  for(sa=0; sa<Nsa;sa++) {
	       for(rg=0; rg<Nrg;rg++) {
		    fscanf(fp, "%f", &X);
		    Tab[pul][sa][rg].re = X;
		    fscanf(fp, "%f", &X);
		    Tab[pul][sa][rg].im = X;
	       }
	  }
     }

     return 0;
}

/************************************************************************/

#pragma omp task out(tab_coh[0;dim])
void init_tab_coh_1(int iv, int nv, int dim, float Tr, 
		    Cplfloat tab_coh[dim]) 
{
     int idop;
     float V, time;

     float Vmax= ATL_Lambda/(4.* Tr);
     V= Vmax * (2.*iv+1.)/(1.*nv);
     
     if (iv>=nv/2) {
	  V=V-2*Vmax;
     }
     for (idop=0; idop<dim; idop++) {
	  time = 2 * Pi * idop * 2 * V * Tr * (1./ATL_Lambda);
	  tab_coh[idop].re= cos(time);
	  tab_coh[idop].im= -sin(time);
     }
}

void init_tab_coh(int nv, int dim, float Tr, Cplfloat tab_coh[nv][dim]) 
{
     int iv;
     
     for (iv=0; iv<nv; iv++) {
	  init_tab_coh_1(iv, nv, dim, Tr, tab_coh[iv]);
     }
}

/************************************************************************/

void init_chirp(int tf, Cplfloat Chirp[tf]) {
     int i=0;
     Chirp[i].re=1.;
     Chirp[i].im=0.;
     i++;
     Chirp[i].re=0.93461;
     Chirp[i].im=0.;
     i++;
     Chirp[i].re=0.832174452;
     Chirp[i].im=0.;
     i++;
     Chirp[i].re=0.617829847;
     Chirp[i].im=0.;
     i++;
     Chirp[i].re=0.309016994;
     Chirp[i].im=0.;
     i++;
     Chirp[i].re=-0.382683432;
     Chirp[i].im=0.;
     i++;
     Chirp[i].re=-0.891615484;
     Chirp[i].im=0.;
     i++;
     Chirp[i].re=-0.66535522;
     Chirp[i].im=0.;
     i++;
     Chirp[i].re=0.251076308;
     Chirp[i].im=0.;
     i++;
     Chirp[i].re=0.747688;
     Chirp[i].im=0.;
     i++;
     Chirp[i].re=-1.75436E-11;
     Chirp[i].im=0.;
     i++;
     Chirp[i].re=-0.623073334;
     Chirp[i].im=0.;
     i++;
     Chirp[i].re=0.173822059;
     Chirp[i].im=0.;
     i++;
     Chirp[i].re=0.380202983;
     Chirp[i].im=0.;
     i++;
     Chirp[i].re=-0.416087226;
     Chirp[i].im=0.;
     i++;
     Chirp[i].re=0.143506287;
     Chirp[i].im=0.;
     i++;
    
}

//  -----------------  stap_functions.c -----------------------------------------

/************************************************************************/
//#pragma omp task input(th, v, nth, nv, ntt, nsa, SubarraySpacing, Lambda, BeamAngle, CarrierSpeed, BeamWidth, TargetDistance, Tr) output(Steervect)
void Calc_steervect_1(int th, int v, int nth, int nv, int ntt, int nsa,
		    float SubarraySpacing,
		    float Lambda,
		    float BeamAngle,
		    float CarrierSpeed,
		    float BeamWidth,
		    float TargetDistance,
		    float Tr,
		    Cplfloat Steervect[ntt*nsa]) 
{
     int ant, rec;
     float V;
     double AngleCarrier_rad= BeamAngle*Pi/180.;

     double theta;
     double Kdeph = 2.*Pi/Lambda;

     double angldir, delta_dist, deph, deph_a, deph_v;

     theta = thint2THfloat(th, nth, BeamWidth);	  

     V = vint2Vfloat(v, nv, Lambda, Tr);

     for (rec=0; rec<ntt; rec++) {
	  for (ant=0; ant<nsa; ant++) {
	       angldir= AngleCarrier_rad + theta+ (SubarraySpacing*ant
						   - CarrierSpeed*Tr*rec)
		    *sin(AngleCarrier_rad)/TargetDistance;
	       delta_dist = -V*Tr*rec;
	       deph_a = Kdeph* SubarraySpacing*ant*cos(angldir);
	       deph_v = Kdeph*2* delta_dist;
	       deph=deph_a + deph_v;
	       
	       Steervect[ant + nsa*rec].re = cos(deph);
	       Steervect[ant + nsa*rec].im = sin(deph);
	       
	  }
     }

}
#pragma omp task out(Steervect[0;nv][0;ntt*nsa])
void Calc_steervect_2(int th, int nth, int nv, int ntt, int nsa,
		    float SubarraySpacing,
		    float Lambda,
		    float BeamAngle,
		    float CarrierSpeed,
		    float BeamWidth,
		    float TargetDistance,
		    float Tr,
		    Cplfloat Steervect[nv][ntt*nsa])
{
          int v;
	  for (v=0; v<nv; v++) {
	       Calc_steervect_1(th, v, nth, nv, ntt, nsa,
				SubarraySpacing, Lambda, BeamAngle,
				CarrierSpeed, BeamWidth, TargetDistance,
				Tr, Steervect[v]);
	  }

}

void Calc_steervect(int nth, int nv, int ntt, int nsa,
		    float SubarraySpacing,
		    float Lambda,
		    float BeamAngle,
		    float CarrierSpeed,
		    float BeamWidth,
		    float TargetDistance,
		    float Tr,
		    Cplfloat Steervect[nth][nv][ntt*nsa]) 
{
     int th;

     for (th=0; th<nth; th++) {
          Calc_steervect_2(th, nth, nv, ntt, nsa,
				SubarraySpacing, Lambda, BeamAngle,
				CarrierSpeed, BeamWidth, TargetDistance,
				Tr, Steervect[th]);
/*	  for (v=0; v<nv; v++) {
	       Calc_steervect_1(th, v, nth, nv, ntt, nsa,
				SubarraySpacing, Lambda, BeamAngle,
				CarrierSpeed, BeamWidth, TargetDistance,
				Tr, Steervect[th][v]);
	  }*/
     }
}

/************************************************************************/
//#pragma omp task input(nrg, ptrin, tf, ptrfiltre) output(ptrout)
void Pulse_Comp_1(int nrg, Cplfloat ptrin[nrg], int tf, 
		Cplfloat ptrfiltre[tf],
		Cplfloat ptrout[nrg-tf+1]) {

     int i, j;
     float R, I;
     for (i=0; i<nrg-tf+1; i++) {

	  R = 0.0;
	  I = 0.0;
	  for (j=0; j<tf; j++) {
	       R += ptrin[i+j].re * ptrfiltre[j].re - ptrin[i+j].im
                    * ptrfiltre[j].im;
	       I += ptrin[i+j].im * ptrfiltre[j].re + ptrin[i+j].re
                    * ptrfiltre[j].im;
	  }
	  ptrout[i].re = R;
	  ptrout[i].im = I;
     }
}


#pragma omp task in(ptrin[0;nrec][0;nrg], ptrfiltre[0;tf]) out(ptrout[0;nrec][0;nrg-tf+1])
void Pulse_Comp(int nrec, int nrg, Cplfloat ptrin[nrec][nrg], int tf, 
		Cplfloat ptrfiltre[tf],
		Cplfloat ptrout[nrec][nrg-tf+1]) {
          int k;

	  for (k=0; k<nrec; k++) {
	       Pulse_Comp_1(nrg, ptrin[k], tf, ptrfiltre, ptrout[k]);
	  }
}

/************************************************************************/
#pragma omp task in(In[0;nsa][0;nrec]) out(Out[0;ntt*nsa][0;ntt*nsa])
void CovAvCov(int nrec, int nsa, int ntt, 
	      Cplfloat In[nsa][nrec],
	      Cplfloat Out[ntt*nsa][ntt*nsa]) {

     int avlength = nrec - ntt +1;

     int r0, r1, a0, a1, l;
     Cplfloat S;
     for (a0=0; a0<nsa; a0++) {
	  for (a1=0; a1<nsa; a1++) {
	       for (r0=0; r0<ntt; r0++) {
		    for (r1=0; r1<ntt; r1++) {
			 S.re=0.;
			 S.im=0.;
			 for (l=0; l<avlength; l++) {
			      S.re += In[a0][r0+l].re * In[a1][r1+l].re
				   + In[a0][r0+l].im * In[a1][r1+l].im;
			      S.im += In[a0][r0+l].im * In[a1][r1+l].re
				   - In[a0][r0+l].re * In[a1][r1+l].im;
			 }
			 Out[a0+nsa*r0][a1+nsa*r1].re = S.re;
			 Out[a0+nsa*r0][a1+nsa*r1].im = S.im;
		    }
	       }
	  }
     }
}


/************************************************************************/
#pragma omp task in(In[0;ntt*nsa][0;ntt*nsa]) out(Out[0;ntt*nsa][0;ntt*nsa])
void Mat_Invert(int ntt, int nsa, Cplfloat In[ntt*nsa][ntt*nsa],
		Cplfloat Out[ntt*nsa][ntt*nsa]) {

     double inv[nsa*ntt][2*nsa*ntt][2];

     double pivot[2], coef[2];

     float re, im;

     int i=0, j=0, k=0, l=0;

     for (i=0; i<ntt*nsa; i++) {
	  for (j=0; j<ntt*nsa; j++) {
	       inv[i][j][0] = (double) In[i][j].re;
	       inv[i][j][1] = (double) In[i][j].im;
	       if (i == j) {
		    inv[i][j+nsa*ntt][0] =1.0;
		    inv[i][j+nsa*ntt][1] =0.0;
	       } else {
		    inv[i][j+nsa*ntt][0]=0.0;
		    inv[i][j+nsa*ntt][1] =0.0;
	       }
	  }
     }

     for (i=0; i<nsa*ntt; i++) {
	  pivot[0]=inv[i][i][0];
	  pivot[1]=inv[i][i][1];

	  if (pivot[0] == 0.) {
	       printf("\n Pivot nul re = %f , im = %f\n", pivot[0], pivot[1]);
	       exit(0);
	  }

	  for (j=i; j<2*nsa*ntt; j++) {
	       re = inv[i][j][0];
	       im = inv[i][j][1];
	       inv[i][j][0] = (re * pivot[0] + im * pivot[1])/(pivot[0] * pivot[0]
							       + pivot[1] * pivot[1]);
	       inv[i][j][1] = (im * pivot[0] - re * pivot[1])/(pivot[0] * pivot[0]
							       + pivot[1] * pivot[1]);
	  }

	  for (k=0; k<nsa*ntt; k++) {
	       if (i!=k) {
		    coef[0] = inv[k][i][0];
		    coef[1] = inv[k][i][1];
		    for (l=i; l<2*nsa*ntt; l++) {
			 inv[k][l][0] -= (coef[0] * inv[i][l][0] - coef[1]
					  * inv[i][l][1]);
			 inv[k][l][1] -= (coef[0] * inv[i][l][1] + coef[1]
					  * inv[i][l][0]);
		    }
	       }
	  }
     }

     for (i=0; i<nsa*ntt; i++) {
	  for (j=0; j<nsa*ntt; j++) {
	       Out[i][j].re = (float) inv[i][j+nsa*ntt][0];
	       Out[i][j].im = (float) inv[i][j+nsa*ntt][1];
	  }
     }
}




/************************************************************************/
//#pragma omp task input(ntt, nsa, Nv, CovInv, SteerVect) output(Fil)
void Calc_Filter_1(int ntt, int nsa, int Nv,
		   Cplfloat CovInv[ntt*nsa][ntt*nsa], 
		   Cplfloat SteerVect[ntt*nsa],
		   Cplfloat Fil[ntt*nsa]) 
{
     int i, j, length_filt;
     Cplfloat X, Y, Z;

     Cplfloat W[50];
     length_filt =ntt*nsa;

     for (i=0; i<length_filt; i++) {
	  Z.re=0.;
	  Z.im=0.;
	  for (j=0; j<length_filt; j++) {
	       X.re = SteerVect[j].re;
	       X.im = SteerVect[j].im;
	       Y.re= CovInv[j][i].re;
	       Y.im= CovInv[j][i].im;
	       Z.re+= X.re*Y.re + X.im*Y.im;
	       Z.im += -X.im*Y.re + X.re*Y.im;
	  }
	  W[i].re=Z.re;
	  W[i].im=Z.im;
     }
     Z.re=0.;
     Z.im=0.;
     for (i=0; i<length_filt; i++) {
	  X.re = SteerVect[i].re;
	  X.im = SteerVect[i].im;
	  Z.re += W[i].re*X.re -W[i].im*X.im;
     }
     for (i=0; i<length_filt; i++) {
	  Fil[i].re= W[i].re/Z.re;
	  Fil[i].im= W[i].im/Z.re;;
     }
}

#pragma omp task in(CovInv[0;ntt*nsa][0;ntt*nsa], SteerVect[0;Nv][0;ntt*nsa]) out(Fil[0;Nv][0;ntt*nsa])
void Calc_Filter(int ntt, int nsa, int Nv,
		 Cplfloat CovInv[ntt*nsa][ntt*nsa], 
		 Cplfloat SteerVect[Nv][ntt*nsa],
		 Cplfloat Fil[Nv][ntt*nsa]) 
{

     int v;
 
     for (v=0; v<Nv; v++) {
	  Calc_Filter_1(ntt, nsa, Nv, CovInv, SteerVect[v], Fil[v]);
     }
}

/************************************************************************/
void Apply_Filter(int Nrec, int Ntt, int Na, Cplfloat Sig[Nrec][Na],
		  Cplfloat Fil[Ntt*Na], Cplfloat Out[Nrec-Ntt+1]) {

     int i, a, r;
     float R, I;

     for (i=0; i<Nrec-Ntt+1; i++) {
	  R=0.;
	  I=0.;
	  for (r=0; r<Ntt; r++) {
	       for (a=0; a<Na; a++) {
		    R+= Sig[i+r][a].re * Fil[a+Na*r].re + Sig[i+r][a].im
			 * Fil[a+Na*r].im;
		    I+= -Sig[i+r][a].re * Fil[a+Na*r].im + Sig[i+r][a].im
			 * Fil[a+Na*r].re;
	       }
	  }
	  Out[i].re = R;
	  Out[i].im = I;
     }
}

#pragma omp task in(Sig[0;Nrec][0;Na], Fil[0;Nv][0;Ntt*Na]) out(Out[0;Nv][0;Nrec-Ntt+1])
void Apply_Filter_task(int Nv, int Nrec, int Ntt, int Na, 
		       Cplfloat Sig[Nrec][Na],
		       Cplfloat Fil[Nv][Ntt*Na],
		       Cplfloat Out[Nv][Nrec-Ntt+1])
{
     int l;

     for (l=0; l<Nv; l++) {
	  Apply_Filter(Nrec, Ntt, Na, Sig, Fil[l], Out[l]);
     }
}

/************************************************************************/
#pragma omp task in(ptrin1[0;nv][0;ndop], ptrin2[0;nv][0;ndop]) out(ptrout[0;nv])
void Int_Dop(int nv, int ndop, Cplfloat ptrin1[nv][ndop],
	     Cplfloat ptrin2[nv][ndop],        Cplfloat ptrout[nv]) {
     int i, j;
     float R, I, X, Y;

     for (i=0; i<nv; i++) {
	  R = 0.0;
	  I = 0.0;
	  for (j=0; j<ndop; j++) {
	       X= ptrin1[i][j].re * ptrin2[i][j].re + ptrin1[i][j].im
                    * ptrin2[i][j].im;
	       Y= -ptrin1[i][j].re * ptrin2[i][j].im + ptrin1[i][j].im
                    * ptrin2[i][j].re;
	       R+= X;
	       I+=Y;
	  }
	  ptrout[i].re = R;
	  ptrout[i].im = I;
     }
}


/************************************************************************/

#pragma omp task in(ptrin[0;Nrg][0;Nv]) out(Pow[0;1])
void average_power_1(int Nth, int Nrg, int Nv, 
		     Cplfloat ptrin[Nrg][Nv],
		     Cplfloat Pow[1])
{
     int v, rg;
     double PP=0.;
     for (rg=0; rg<Nrg; rg++) {
	  for (v=0; v<Nv; v++) {
	       PP += ptrin[rg][v].re *ptrin[rg][v].re
		    +ptrin[rg][v].im *ptrin[rg][v].im;
	  }
     }
     Pow->re= (float)(PP/((float)(Nv*Nrg)));
     Pow->im= 0.;
}

void average_power(int Nth, int Nrg, int Nv, 
		   Cplfloat ptrin[Nth][Nrg][Nv],
		   Cplfloat Pow[Nth]) {

     int th;

     for (th=0; th<Nth; th++) {
	  average_power_1(Nth, Nrg, Nv, ptrin[th], &Pow[th]);
     }
}

/************************************************************************/

#pragma omp task in(In[0;Nrg][0;Nv]) out(TOut[0;Nv-tfac_L_V+1])
void tfac_1(int rg, int Nth, int Nrg, int Nv, 
	    Cplfloat In[Nrg][Nv], 
	    int tfac_L_V, int tfac_L_RG, float tfac_K, float tfac_AVPWK,
	    Cplfloat TOut[Nv-tfac_L_V+1]) 
{
     int midRg, midV, nb_neighb;
     double AvPower;
     float val;
     int v, xv, xrg;
     double S, Sc, T;
     double sqrAvPow;
     int IsMax;

     midRg= (tfac_L_RG - 1)/2;
     midV= (tfac_L_V - 1)/2;
     nb_neighb = tfac_L_RG*tfac_L_V-1;

     AvPower= 0.;

     AvPower = (double)(In[0][0].re);
     sqrAvPow = sqrt(AvPower);

     for (v=0; v<=Nv -tfac_L_V; v++) {
	  
	  Sc=ModSquare(In[rg+midRg][v+midV]);
	  IsMax=1;
	  S=0.;
	  for (xrg=0; xrg<2*midRg+1; xrg++) {
	       for (xv=0; xv< 2*midV+1; xv++) {
		    if ( (xrg!=midRg)||(xv!=midV)) {
			 T=ModSquare(In[rg+xrg][v+xv]);
			 if (Sc<T) {
			      IsMax=0;
			 }
			 S+= T;
		    }
	       }
	  }
	  
	  if ((IsMax)&&(Sc>tfac_AVPWK* AvPower) &&(Sc> tfac_K*S
						   /(double)nb_neighb)) {
	       val = 1.*IsMax;
	       
	  } else {
	       val=0.;
	  }
	  TOut[v].re=val*In[rg+midRg][v+midV].re/sqrAvPow;
	  TOut[v].im=val*In[rg+midRg][v+midV].im/sqrAvPow;     
     }
}

void tfac(int Nth, int Nrg, int Nv, 
	  Cplfloat In[Nth][Nrg][Nv], 
	  int tfac_L_V, int tfac_L_RG, float tfac_K, float tfac_AVPWK,
	  Cplfloat TOut[Nth][Nrg-tfac_L_RG+1][Nv-tfac_L_V+1]) 
{
     int th, rg;

     for (th=0; th<Nth; th++) {
	  for (rg=0; rg<=Nrg - tfac_L_RG; rg++) {
	       tfac_1(rg, Nth, Nrg, Nv, In[th], 
		      tfac_L_V, tfac_L_RG, tfac_K, tfac_AVPWK,
		      TOut[th][rg]);
	  }
     }
}

/************************************************************************/

#pragma omp task in(In[0;Nv], Pow[0;1]) out(TOut[0;NvLarge])
void tfac_add_edges_1(int rg, int Nrg, int Nv, Cplfloat In[Nv],
		      Cplfloat Pow[1], int NrgLarge, int NvLarge,
		      Cplfloat TOut[NvLarge]) 
{
     int v;
     int sidesv, sidesRg;
     float sqpow;

     sqpow = (float) sqrt(Pow[0].re);

     sidesv = (NvLarge-Nv)/2;
     sidesRg = (NrgLarge-Nrg)/2;
     
     for (v=0; v<NvLarge; v++) {
	  if ( (rg<sidesRg) || (rg>=Nrg+sidesRg) ||
	       (v<sidesv) || (v>=Nv+sidesv)) {

	       TOut[v].re = sqpow;
	       TOut[v].im = 0;
	  } else {
	       TOut[v].re = In[v-sidesv].re;
	       TOut[v].im = In[v-sidesv].im;
	  }
     }
}

void tfac_add_edges(int Nth, int Nrg, int Nv, Cplfloat In[Nth][Nrg][Nv],
		    Cplfloat Pow[Nth], int NrgLarge, int NvLarge,
		    Cplfloat TOut[Nth][NrgLarge][NvLarge]) {
     
     int th, rg;
     int sidesRg = (NrgLarge-Nrg)/2;
     
     for (th=0; th<Nth; th++) {
	  for (rg=0; rg<NrgLarge; rg++) {
	       if ( (rg<sidesRg) ||(rg>=Nrg+sidesRg) ) {
		    tfac_add_edges_1(rg, Nrg, Nv, In[th][rg], &Pow[th],
				     NrgLarge, NvLarge, TOut[th][rg]);
	       } 
	       else {
		    tfac_add_edges_1(rg, Nrg, Nv, In[th][rg-sidesRg], &Pow[th],
				     NrgLarge, NvLarge, TOut[th][rg]);
		    
	       }
	  }
     }
}


/************************************************************************/
void MaxPower(int Nth, Cplfloat In[Nth], Cplfloat Out[Nth]) {

     double Pow, PowMax;
     int thmax;
     int th;
     thmax=0;

     PowMax=0.;
     for (th=0; th<Nth; th++) {
	  Pow = In[th].re*In[th].re+In[th].im*In[th].im;
	  if (Pow>PowMax) {
	       PowMax=Pow;
	       thmax=th;
	  }
     }
     for (th=0; th<Nth; th++) {
	  if (th!=thmax) {
	       Out[th].re=0.;
	       Out[th].im=0.;
	  } else {
	       Out[th].re=In[th].re;
	       Out[th].im=In[th].im;
	  }
     }
}

#pragma omp task in(In[0;nv][0;Nth]) out(Out[0;nv][0;Nth])
void MaxPower_task(int Nth, int nv, 
		   Cplfloat In[nv][Nth], Cplfloat Out[nv][Nth]) 
{
     int k;
     
     for (k=0; k<nv; k++) {
	  MaxPower(Nth, In[k], Out[k]);
     }
}

/************************************************************************/

#pragma omp task in(In[0;Nth][0;Nv]) out(Out[0;Nth][0;Nv])
void CorrecV(int Nth, int Nv, Cplfloat In[Nth][Nv], Cplfloat Out[Nth][Nv],
	     float Beamwidth, float AngleCarrier, float Tr,
	     float lambda, float Vcarrier) {

     int v, th, vcor;
     float AngleCarrier_rad, V, TH, dVcar;
     AngleCarrier_rad=AngleCarrier*Pi/180.;

     memset(Out, 0, Nth*Nv*sizeof(Cplfloat));

     for (th=0; th<Nth; th++) {
	  TH = thint2THfloat(th, Nth, Beamwidth);
	  TH += AngleCarrier_rad;
	  dVcar= Vcarrier*cos(TH);
	  for (v=0; v<Nv; v++) {
	       V= vint2Vfloat(v, Nv, lambda, Tr);
	       V -= dVcar;

	       vcor = Vfloat2vint(V, Nv, lambda, Tr);

	       Out[th][vcor].re = In[th][v].re;
	       Out[th][vcor].im = In[th][v].im;
	  }
     }
}

/************************************************************************/


void Mat2Plot(int Nth, int Nv, int Nrg, int nbplots,
	      Cplfloat In[Nrg][Nth][Nv], Cplfloat Out[nbplots]) {

     int th, v, rg;
     int pl=0;
     int i;
     double pow;

     for (rg=0; rg<Nrg; rg++) {
	  for (v=0; v<Nv; v++) {
	       for (th=0; th<Nth; th++) {
		    pow= ModSquare(In[rg][th][v]);
		    if ( (pow >1.E-5)&&(pl<nbplots)) {
			 Out[pl].re = (float)( 10000*rg+100*v+ th);
			 Out[pl].im = (float) pow;
			 pl++;
		    }
	       }
	  }
     }

     for (i=pl; i<nbplots; i++) {
	  Out[i].re = 0.;
	  Out[i].im = 0.;
     }
} 


/************************************************************************/
void RemAmbiguity(int nbplots, int Nraf, Cplfloat In[Nraf][nbplots], float pTr[Nraf],
		  int Nrg, int Nv, int Nth, float rg_size, float lambda, Cplfloat Out[4*nbplots]) {

    
     int pl, raf,i;
     int th, v, newv0, rg;

     memset(Out, 0, 4*nbplots*sizeof(Cplfloat));

     float TAB_RGVTH[Nrg][2*Nv][Nth];
     int TAB_NBHIT[Nrg][2*Nv][Nth];

     for (rg=0; rg<Nrg; rg++) {
	  for (v=0; v<2*Nv; v++) {
	       for (th=0; th<Nth; th++) {
		    TAB_RGVTH[rg][v][th]=0.;
		    TAB_NBHIT[rg][v][th]=0;
	       }
	  }
     }


     float V;
     float Trmin= pTr[0];
     for(raf=0; raf<Nraf;raf++){
	  if(pTr[raf]<Trmin){
	       Trmin=pTr[raf];
	  }
     }
     float Vmax0 = lambda/(4.*Trmin);
    	
    
     Cplfloat PLOT, NEW_PLOT; 
     for (raf=0; raf<Nraf; raf++) {
	  float Trec = pTr[raf];
	  for (pl=0; pl<nbplots; pl++) {
	       PLOT = In[raf][pl];
	       if (PLOT.im!=0.) {
		    rg = (int)(PLOT.re/10000. +.1);
		    v = (int)( (PLOT.re - 10000.*rg)/100. +.1);
		    th = (int)(PLOT.re - 10000.*rg -100.*v +.1);               
		    V = vint2Vfloat(v, Nv, lambda, Trec);
		    float Vmax= lambda/(4*Trec);
		    for(i=-2;i<3;i++){
			 if( (V + 2*i*Vmax>-2*Vmax0) && (V + 2*i*Vmax<2*Vmax0)){
			      newv0 = Vfloat2vint(V + 2*i*Vmax, 2*Nv, lambda, Trmin/2.);
			      TAB_RGVTH[rg][newv0][th] += PLOT.im;
			      TAB_NBHIT[rg][newv0][th] += 1;
			 }
		    }
	       }
	  }
     }
    
     int RG_WINDOW=1;
     int V_WINDOW=1;
     int TH_WINDOW=1;
     int a, b, c;
     pl=0;
 
     for (rg=RG_WINDOW; rg<Nrg-RG_WINDOW; rg++) {
	  for (v=V_WINDOW; v<2*Nv-V_WINDOW; v++) {
	       for (th=TH_WINDOW; th<Nth-TH_WINDOW; th++) {
		    if(TAB_NBHIT[rg][v][th]>0){
			 int isMax=1;
			 int nb_hits=TAB_NBHIT[rg][v][th];
			 float pow = TAB_RGVTH[rg][v][th];
			 float add_pow=0.;
			 for (a=-RG_WINDOW; a<=RG_WINDOW; a++) {
			      for (b=-V_WINDOW; b<=V_WINDOW; b++) {
				   for (c=-TH_WINDOW; c<=TH_WINDOW; c++) {
	                        	float pow_a= TAB_RGVTH[rg+a][v+b][th+c];
	                        	int h = TAB_NBHIT[rg+a][v+b][th+c];
					if (((a!=0)||(b!=0)||(c!=0))&&(h!=0)){
					     nb_hits+=h;
					     add_pow+= pow_a;
					     if (pow_a >pow) {
						  isMax=0;
					     }
					}
				   }
			      }
			 }
			 if ( (isMax>0) &&(nb_hits>= .8*Nraf)) {
			      NEW_PLOT.re = 10000.*rg + 100.*v + 1.*th;
			      NEW_PLOT.im = pow + add_pow;
			      Out[pl] = NEW_PLOT;
			      pl++;
			 }
		    }
	       }
	  }
     }
 
     while (pl<2*nbplots) {
	  NEW_PLOT.re = 0.;
	  NEW_PLOT.im = 0.;
	  Out[pl] = NEW_PLOT;
	  pl++;
     }
    

}





/************************************************************************/
void X_2(int nsa, int nrec, int dim3, Cplfloat a[nsa][nrec][dim3],
	 Cplfloat b[dim3][nsa][nrec]) {
     int j, k, l;

     for (j=0; j<nsa; j++) {
	  for (k=0; k<nrec; k++) {
	       for (l=0; l<dim3; l++) {
		    b[l][j][k].re = a[j][k][l].re;
		    b[l][j][k].im = a[j][k][l].im;
	       }
	  }
     }
}

/************************************************************************/
#pragma omp task in(a[0;dim2][0;dim3]) out(b[0;dim3][0;dim2])
void X_1_1(int dim2, int dim3, 
	   Cplfloat a[dim2][dim3],
	   Cplfloat b[dim3][dim2]) 
{
     int k, l;

     for (k=0; k<dim2; k++) {
	  for (l=0; l<dim3; l++) {
	       b[l][k].re = a[k][l].re;
	       b[l][k].im = a[k][l].im;
	  }
     }
}

// parallelize transpose of internal 2D arrays
void X_1(int dim1, int dim2, int dim3, 
	 Cplfloat a[dim1][dim2][dim3],
	 Cplfloat b[dim1][dim3][dim2]) 
{
     int j;

     for (j=0; j<dim1; j++) {
	  X_1_1(dim2, dim3, a[j], b[j]);
     }
}


/************************************************************************/
void turn5(int dim1, int ntt, int nsa, Cplfloat a[dim1][ntt][nsa][ntt][nsa],
	   Cplfloat b[nsa][ntt][nsa][ntt][dim1]) {
     int i, j, k, l, t;

     for (t=0; t<dim1; t++) {
	  for (i=0; i<ntt; i++) {
	       for (j=0; j<nsa; j++) {
		    for (k=0; k<ntt; k++) {
			 for (l=0; l<nsa; l++) {
			      b[i][j][k][l][t].re = a[t][i][j][k][l].re;
			      b[i][j][k][l][t].im = a[t][i][j][k][l].im;
			 }
		    }
	       }
	  }
     }
}

/************************************************************************/

#pragma omp task in(a[0;dim3]) out(b[0;dim3])
void X_3_1(int dim3, 
	   Cplfloat a[dim3],
	   Cplfloat b[dim3]) 
{
     int k;
     for (k=0; k<dim3; k++) {
	  b[k].re = a[k].re;
	  b[k].im = a[k].im;
     }
}

void X_3(int dim1, int dim2, int dim3, 
	 Cplfloat a[dim1][dim2][dim3],
	 Cplfloat b[dim2][dim1][dim3]) 
{
     int i, j;

     for (i=0; i<dim1; i++) {
	  for (j=0; j<dim2; j++) {
	       X_3_1(dim3, a[i][j], b[j][i]);
	  }
     }
}

/************************************************************************/
void X_4(int dim1, int dim2, int dim3, Cplfloat a[dim1][dim2][dim3],
	 Cplfloat b[dim2][dim3][dim1]) {
     int i, j, k;

     for (i=0; i<dim1; i++) {
	  for (j=0; j<dim2; j++) {
	       for (k=0; k<dim3; k++) {
		    b[j][k][i].re = a[i][j][k].re;
		    b[j][k][i].im = a[i][j][k].im;
	       }
	  }
     }
}

/************************************************************************/

#pragma omp task in(a[0;dim1][0;dim2]) out(b[0;dim2][0;dim1])
void X_5(int dim1, int dim2, Cplfloat a[dim1][dim2], Cplfloat b[dim2][dim1]) {
     int i, j;
     for (i=0; i<dim1; i++) {
	  for (j=0; j<dim2; j++) {
	       b[i][j].re = a[j][i].re;
	       b[i][j].im = a[j][i].im;
	  }
     }
}




/************************************************************************/
/***    Bursts processing  ******************************************************************/
/************************************************************************/


void trt_burst(int nrg, int tf, int nv, int nth, int nsa, int ntt,
	       int nrec, int nplots, 
	       int Tfac_SizeV, int Tfac_sizeRG, int Tfac_K, int Tfac_WK,
	       float SubarraySpacing, float Lambda, 
	       float BeamAngle, float CarrierSpeed, float BeamWidth, 
	       float TargetDistance,float Tr,
	       Cplfloat in_pulse[nsa][nrec][nrg],
	       Cplfloat out_plot[nplots]) {

     int nrgf = nrg-tf+1;
     int NrgLarge=  nrgf + Tfac_sizeRG -1; 
     int NvLarge= nv + Tfac_SizeV -1;
     int j, k, l;


//     init_steervect(nth, nv, ntt, nsa, Steervect);
     Calc_steervect(nth, nv, ntt, nsa,SubarraySpacing, Lambda, 
		    BeamAngle,CarrierSpeed,BeamWidth, TargetDistance,Tr, Steervect); 	
     
     // yoav: moved outside the function
//     init_chirp(tf, filtre);
     init_tab_coh(nv, nrec-ntt+1, Tr,tab_coh);

     for (j=0; j<nsa; j++) {
          Pulse_Comp(nrec,nrg, in_pulse[j], tf, filtre, out_pulse[j]);
/*	  for (k=0; k<nrec; k++) {
	       Pulse_Comp_1(nrg, in_pulse[j][k], tf, filtre, out_pulse[j][k]);
	  }*/
     }

// need to finish calculating out_pulse before the 3D transpose in X_2
// need to finish calculating tab_coh before using it later in Int_Dop
#pragma omp taskwait

     // SEQUENTIAL FUNCTION
     X_2(nsa, nrec, nrgf, out_pulse, in_cov); 

     X_1(nrgf, nsa, nrec, in_cov, in_apply);

     for (j=0; j<nrgf; j++) {
	  CovAvCov(nrec, nsa, ntt, in_cov[j], out_avcov[j]);
     }

     for (j=0; j<nrgf; j++) {
	  Mat_Invert(ntt, nsa, out_avcov[j], mat_inv[j]);

     }

     for (j=0; j<nrgf; j++) {
	  for (k=0; k<nth; k++) {   	
	       Calc_Filter(ntt, nsa, nv, mat_inv[j], Steervect[k], Fil[j][k]);
	  }
     }

     for (j=0; j<nrgf; j++) {
	  for (k=0; k<nth; k++) {
	       
	       Apply_Filter_task(nv, nrec, ntt, nsa, in_apply[j], Fil[j][k],
				 out_apply[j][k]);
	  }
     }

     // input dependencies in 2D, output dependency in 1D
     for (j=0; j<nrgf; j++) {
	  for (k=0; k<nth; k++) {
	       Int_Dop(nv, nrec-ntt+1, out_apply[j][k], tab_coh,
		       out_int_dop[j][k]);
	  }
     }
 
     // task input dependencies in 1D
     X_3(nrgf, nth, nv, out_int_dop, in_avpow);

     // task inputs are 2D (out of 3D). outputs at 0D point (out of 1D)
     average_power(nth, nrgf, nv, in_avpow, out_avpow);

     // task input deps at 0D (point in out_avpow). 
     // outputs at 1D (of 3D in out_facedges)
     tfac_add_edges(nth, nrgf, nv, in_avpow, out_avpow, 
		    NrgLarge, NvLarge, out_facedges);

// need to wait for 1D to finish so we can operate on 2D
#pragma omp taskwait  
 
     // operate on input/output at 2D granularity (out of 3D)
     tfac(nth, NrgLarge, NvLarge, out_facedges, 
	  Tfac_SizeV, Tfac_sizeRG, Tfac_K, Tfac_WK, out_face1);

// merge out_face1 before sequential X_4
#pragma omp taskwait  
 
     // REMAINS SEQUENTIAL
     X_4(nth, nrgf, nv, out_face1, in_maxpow);

     for (j=0; j<nrgf; j++) {
	  // external task operates on 2D arrays (of 3D in in_maxpow,
	  // out_maxpow)
	  MaxPower_task(nth, nv, in_maxpow[j], out_maxpow[j]);
     }

     for (j=0; j<nrgf; j++) {
	  // input/output are 2D arrays
	  X_5(nth, nv, out_maxpow[j], in_correctionV[j]);
     }

     for (j=0; j<nrgf; j++) {
	  // task inputs/outputs are 2D arrays
	  CorrecV(nth, nv, in_correctionV[j], out_correctionV[j], BeamWidth,
		  BeamAngle, Tr, Lambda, CarrierSpeed);
     }

#pragma omp taskwait  

     Mat2Plot(nth, nv, nrgf, nplots, out_correctionV, out_plot);
}


/************************************************************************/
/*** MAIN ***************************************************************/
/************************************************************************/
int main() {


     int Nraf_max=10;
     int Nraf=0;
     int i=0, j=0 /*remove j - arico*/;
     
     struct timeval start, stop;

     Cplfloat out_plot[Nraf_max][ATL_Nplots];
     Cplfloat amb_out[ATL_Nplots*4];
     Cplfloat in_pulse[Nraf_max][ATL_nsa][ATL_nrec][ATL_Nrg];
    
     float Tr[Nraf_max];
     FILE *fp;

     char *filename="InputStimuli.txt";

     //Open input file
     if( (fp=fopen(filename,"rb"))==NULL) {
	  fprintf(stderr,"couldn't open \"%s\"!\n",filename);
	  return 1;
     }

     float Trn=0.;
     while( HasNextBurst(fp, &Trn)){		
	  GetNextBurst(fp, ATL_Nrg, ATL_nsa, ATL_nrec, in_pulse[Nraf]) ;
	  Tr[Nraf]=Trn;
	  Nraf++;
     }

     // yoav: moved here from inside trt_burst
     init_chirp(ATL_tf, filtre);

     gettimeofday(&start, NULL);

     for (j=0; j<20; j++) { // arico: repeat 20 times to make it long enough for benchmarking 
     for (i=0; i<Nraf; i++) {
          trt_burst(ATL_Nrg, ATL_tf, ATL_Nv, ATL_Nth, ATL_nsa,ATL_ntt, 
		    ATL_nrec, ATL_Nplots, 
		    ATL_tfac_L_V, ATL_tfac_L_RG, ATL_tfac_K , ATL_tfac_AVPWK, 
		    ATL_SubarraySpacing, ATL_Lambda, 
		    ATL_BeamAngle,ATL_CarrierSpeed,ATL_BeamWidth, 
		    ATL_TargetDistance,/*ATL_Tr*/ Tr[i],
		    in_pulse[i], out_plot[i]);
     }
     }
#pragma omp taskwait
     gettimeofday(&stop, NULL);
     
     {
	     long secs, usecs;
             unsigned long elapsed;
	     secs = (long)(stop.tv_sec - start.tv_sec);
	     if(stop.tv_usec >= start.tv_usec) {
		     usecs = (long)(stop.tv_usec - start.tv_usec);
	     } else {
		     --secs;
		     usecs = 1000000 + (long)stop.tv_usec - (long)start.tv_usec;
	     }
	     
	     printf("Parallel section: %ld.%2ld seconds\n", secs, usecs);

             elapsed = 1000000 * (stop.tv_sec - start.tv_sec);
             elapsed += stop.tv_usec - start.tv_usec;
             printf("par_sec_time_us:%lu\n",elapsed);
     }

     RemAmbiguity(ATL_Nplots, Nraf, out_plot, Tr, ATL_Nrg-ATL_tf+1, ATL_Nv,ATL_Nth, ATL_rg_size, ATL_Lambda, amb_out);

     float Trmin= Tr[0];
     for(i=0; i<Nraf; i++){
	  if(Tr[i]<Trmin){
	       Trmin=Tr[i];
	  }
     }

     RecordPlots( "PlotsOut.txt", ATL_Nplots*4,amb_out, ATL_rg_size, 2*ATL_Nv,ATL_Lambda/(2.*Trmin),ATL_Nth,ATL_BeamWidth  );
    
     printf("STAP - end of program.\n# file 'PlotsOut.txt' generated.\n");
    
     return(0);
    
} 
