#include "Rate.hpp"
#include <algorithm>
#include <sstream>
#include "gsl/gsl_sf_bessel.h"
#include <boost/math/special_functions/bessel.hpp>

#define epsrel_CB 1.e-4
#define epsabs_CB 1.e-12
// This File stores all amplitude and rate calculations

//====================
//Calculation of Rate/Collision term
//====================

//For decay

/*
void Rate::GetProdRate1to2(long double &rate)
{
    long double tw = 0.0l;
    long double PSdm = 0.0l;
    long double PSsm = 0.0l, PSe = 0.0l, PSmu = 0.0l, PStau = 0.0l;
    if (m_Zp > 2.*m_h2)
    {
        // Zp -> X X
        PSdm = 1.;
    }
    rate = m_Zp/(12.l*M_PI)*PSdm*Get_NZp(temp);
}

long double Rate::Width_h2()
{
    long double tw = 0.0l;
    long double PSdm = 0.0l;
    long double PSsm = 0.0l, PSe = 0.0l;
    if (m_Zp > 2.*m_h2)
    {
        // Zp -> X X
        PSdm = 1.;
    }
    if (m_Zp > 2.*me)
    {
        // Zp -> e+ e-
        PSe = 1.;
    }
    PSsm = PSe;
    tw = 12.0l;
    return tw;
}

*/

//For 2to2 process
// id = 1: theta theta -> h h
// id = 2: theta theta -> b B
// id = 3: theta theta -> t T
// id = 4: theta theta -> s S
// id = 5: theta theta -> c C
// id = 6: theta theta -> d D
// id = 7: theta theta -> u U
// id = 8: theta theta -> l L
// id = 9: theta theta -> m M
// id = 10: theta theta -> e E
// id = 11: theta theta -> W+ W-
// id = 12: theta theta -> Z Z

/*
#include <iostream>
#include <fstream>
#include <string>
 
using namespace std;
string[] amplitudes(){
   ifstream file;
   
   //count the lines in file to get size of array
   string unused;
   int num_lines=0;
   
   file.open("/home/parracho/Desktop/Project/Cuba-4.2.2/Ultralight/src/var.txt"); //open a file to perform read operation using file object
   while( getline(file, unused)){
   	num_lines++;	
   }
   file.close();
   string amplitude[num_lines];
   file.open("/home/parracho/Desktop/Project/Cuba-4.2.2/Ultralight/src/var.txt"); //open a file to perform read operation using file object
   int counter = 0;
   while(getline(file, amplitude[counter])){ //read data from file object and put it into string.
      	counter++;
      }
   file.close(); //close the file object;
   return amplitude;
}
*/


long double Rate::SqAmp(long double s, long double cc, int id)
{
    long double amp = 0.l;
    long double num=0.l, den=0.l;

    
	if (id==1){
	num=pow(2*hhh*hHH*mW*SW*pow(mh2,2) - 3*g*HHH*SW*pow(mH,2)*(pow(mH,2) - 2*pow(mh2,2))*sin(a) + ((s - 2*pow(mh2,2))*(4*hhh*hHH*mW*SW + 6*g*HHH*SW*pow(mH,2)*sin(a)))/2.,2);
	den=8*pow(mW,2)*pow(SW,2)*(pow(mH,2)*pow(wh,2) + pow(s - pow(mH,2),2))*(pow(mh2,2)*pow(wh2,2) + pow(s - pow(mh2,2),2));
	amp=num/den;
	}
	else if (id==2){
	num=-3*pow(g,2)*pow(mqb,2)*(-pow(mh2,2) + (-s + 2*pow(mh2,2))/2. + 2*pow(mqb,2))*pow(hHH*cos(a)*pow(mh2,2) + HHH*(pow(mH,2) - 2*pow(mh2,2))*sin(a) + (s - 2*pow(mh2,2))*(hHH*cos(a) - HHH*sin(a)),2);
	den=pow(mW,2)*(pow(mH,2)*pow(wh,2) + pow(s - pow(mH,2),2))*(pow(mh2,2)*pow(wh2,2) + pow(s - pow(mh2,2),2));
	amp=num/den;
	}
	else if (id==3){
	num=-3*pow(g,2)*pow(mqt,2)*(-pow(mh2,2) + (-s + 2*pow(mh2,2))/2. + 2*pow(mqt,2))*pow(hHH*cos(a)*pow(mh2,2) + HHH*(pow(mH,2) - 2*pow(mh2,2))*sin(a) + (s - 2*pow(mh2,2))*(hHH*cos(a) - HHH*sin(a)),2);
	den=pow(mW,2)*(pow(mH,2)*pow(wh,2) + pow(s - pow(mH,2),2))*(pow(mh2,2)*pow(wh2,2) + pow(s - pow(mh2,2),2));
	amp=num/den;
	}
	else if (id==4){
	num=-3*pow(g,2)*pow(mqs,2)*(-pow(mh2,2) + (-s + 2*pow(mh2,2))/2. + 2*pow(mqs,2))*pow(hHH*cos(a)*pow(mh2,2) + HHH*(pow(mH,2) - 2*pow(mh2,2))*sin(a) + (s - 2*pow(mh2,2))*(hHH*cos(a) - HHH*sin(a)),2);
	den=pow(mW,2)*pow(SW,2)*(pow(mH,2)*pow(wh,2) + pow(s - pow(mH,2),2))*(pow(mh2,2)*pow(wh2,2) + pow(s - pow(mh2,2),2));
	amp=num/den;
	}
	else if (id==5){
	num=-3*pow(g,2)*pow(mqc,2)*(-pow(mh2,2) + (-s + 2*pow(mh2,2))/2. + 2*pow(mqc,2))*pow(hHH*cos(a)*pow(mh2,2) + HHH*(pow(mH,2) - 2*pow(mh2,2))*sin(a) + (s - 2*pow(mh2,2))*(hHH*cos(a) - HHH*sin(a)),2);
	den=pow(mW,2)*(pow(mH,2)*pow(wh,2) + pow(s - pow(mH,2),2))*(pow(mh2,2)*pow(wh2,2) + pow(s - pow(mh2,2),2));
	amp=num/den;
	}
	else if (id==6){
	num=-3*pow(g,2)*pow(mqd,2)*(-pow(mh2,2) + (-s + 2*pow(mh2,2))/2. + 2*pow(mqd,2))*pow(hHH*cos(a)*pow(mh2,2) + HHH*(pow(mH,2) - 2*pow(mh2,2))*sin(a) + (s - 2*pow(mh2,2))*(hHH*cos(a) - HHH*sin(a)),2);
	den=pow(mW,2)*(pow(mH,2)*pow(wh,2) + pow(s - pow(mH,2),2))*(pow(mh2,2)*pow(wh2,2) + pow(s - pow(mh2,2),2));
	amp=num/den;
	}
	else if (id==7){
	num=-3*pow(g,2)*pow(mqu,2)*(-pow(mh2,2) + (-s + 2*pow(mh2,2))/2. + 2*pow(mqu,2))*pow(hHH*cos(a)*pow(mh2,2) + HHH*(pow(mH,2) - 2*pow(mh2,2))*sin(a) + (s - 2*pow(mh2,2))*(hHH*cos(a) - HHH*sin(a)),2);
	den=pow(mW,2)*(pow(mH,2)*pow(wh,2) + pow(s - pow(mH,2),2))*(pow(mh2,2)*pow(wh2,2) + pow(s - pow(mh2,2),2));
	amp=num/den;
	}
	else if (id==8){
	num=-(pow(g,2)*pow(mtau,2)*(-pow(mh2,2) + (-s + 2*pow(mh2,2))/2. + 2*pow(mtau,2))*pow(hHH*cos(a)*pow(mh2,2) + HHH*(pow(mH,2) - 2*pow(mh2,2))*sin(a) + (s - 2*pow(mh2,2))*(hHH*cos(a) - HHH*sin(a)),2));
	den=pow(mW,2)*(pow(mH,2)*pow(wh,2) + pow(s - pow(mH,2),2))*(pow(mh2,2)*pow(wh2,2) + pow(s - pow(mh2,2),2));
	amp=num/den;
	}
	else if (id==9){
	num=-(pow(g,2)*pow(mmu,2)*(-pow(mh2,2) + (-s + 2*pow(mh2,2))/2. + 2*pow(mmu,2))*pow(hHH*cos(a)*pow(mh2,2) + HHH*(pow(mH,2) - 2*pow(mh2,2))*sin(a) + (s - 2*pow(mh2,2))*(hHH*cos(a) - HHH*sin(a)),2));
	den=pow(mW,2)*(pow(mH,2)*pow(wh,2) + pow(s - pow(mH,2),2))*(pow(mh2,2)*pow(wh2,2) + pow(s - pow(mh2,2),2));
	amp=num/den;
	}
	else if (id==10){
	num=-(pow(g,2)*pow(me,2)*(2*pow(me,2) - pow(mh2,2) + (-s + 2*pow(mh2,2))/2.)*pow(hHH*cos(a)*pow(mh2,2) + HHH*(pow(mH,2) - 2*pow(mh2,2))*sin(a) + (s - 2*pow(mh2,2))*(hHH*cos(a) - HHH*sin(a)),2));
	den=pow(mW,2)*(pow(mH,2)*pow(wh,2) + pow(s - pow(mH,2),2))*(pow(mh2,2)*pow(wh2,2) + pow(s - pow(mh2,2),2));
	amp=num/den;
	}
	else if (id==11){
	num=pow(g,2)*(pow(mh2,4) + (s - 2*pow(mh2,2))*(pow(mh2,2) - pow(mW,2)) - 2*pow(mh2,2)*pow(mW,2) + 3*pow(mW,4) + pow(s - 2*pow(mh2,2),2)/4.)*pow(hHH*cos(a)*pow(mh2,2) + HHH*(pow(mH,2) - 2*pow(mh2,2))*sin(a) + (s - 2*pow(mh2,2))*(hHH*cos(a) - HHH*sin(a)),2);
	den=pow(mW,2)*(pow(mH,2)*pow(wh,2) + pow(s - pow(mH,2),2))*(pow(mh2,2)*pow(wh2,2) + pow(s - pow(mh2,2),2));
	amp=num/den;
	}
	else if (id==12){
	num=pow(g,2)*(3*pow(mW,4) + 2*pow(mh2,2)*pow(mW,2)*(-1 + pow(SW,2)) + (s - 2*pow(mh2,2))*(pow(mW,2) + pow(mh2,2)*(-1 + pow(SW,2)))*(-1 + pow(SW,2)) + pow(mh2,4)*pow(-1 + pow(SW,2),2) + (pow(s - 2*pow(mh2,2),2)*pow(-1 + pow(SW,2),2))/4.)*pow(hHH*cos(a)*pow(mh2,2) + HHH*(pow(mH,2) - 2*pow(mh2,2))*sin(a) + (s - 2*pow(mh2,2))*(hHH*cos(a) - HHH*sin(a)),2);
	den=2*pow(mW,2)*(pow(mH,2)*pow(wh,2) + pow(s - pow(mH,2),2))*(pow(mh2,2)*pow(wh2,2) + pow(s - pow(mh2,2),2))*pow(-1 + pow(SW,2),2);
	amp=num/den;
	}
	else if (id==13){
	num=pow(pow(mh2,2)*(-(HHH*cos(a)*(pow(mH,2) - 2*pow(mh2,2))) + hHH*pow(mH,2)*sin(a)) + (s - 2*pow(mh2,2))*(HHH*cos(a)*pow(mh2,2) + hHH*pow(mH,2)*sin(a)),2);
	den=2*pow(vsigma,2)*(pow(mH,2)*pow(wh,2) + pow(s - pow(mH,2),2))*(pow(mh2,2)*pow(wh2,2) + pow(s - pow(mh2,2),2));
	amp=num/den;
	}
    else
    {
        amp = 0.l;
    }
    return amp;
}

										/*RATE OF INTERACTION / HUBBLE CONSTANT*/

// FI rate with MAXWELL-BOLTZMANN approximation
//For 2to2 Process, we need to integrate over s and cc(=cosine) to get the collision rate
int IntegrandForRate(const int *ndim, const long double x[], const int *ncomp, long double ff[], void *_params)
{
    Rate *params = (Rate*)_params;
    long double s, smin, smax, JacS;
    long double cc, ccmin, ccmax, JacCos;
    long double m1, m2, m3, m4;
    long double S12, S34;
    long double func;



    if (params->iproc == 1)// id = 1: theta theta -> h h
    {
        m1 = params->m_h2; m2 = params->m_h2; m3 = params->mH; m4 = params->mH; S12 = 0.5l; S34=0.5l;
    }
    else if (params->iproc == 2)// id = 2: theta theta -> b B
    {
        m1 = params->m_h2; m2 = params->m_h2; m3 = mqb; m4 = mqb; S12 = 0.5l; S34=1.0l;
    }
    else if (params->iproc == 3)// id = 3: theta theta -> t T
    {
        m1 = params->m_h2; m2 = params->m_h2; m3 = mqt; m4 = mqt; S12 = 0.5l; S34=1.0l;
    }
    else if (params->iproc == 4)// id = 4: theta theta -> s S
    {
        m1 = params->m_h2; m2 = params->m_h2; m3 = mqs; m4 = mqs; S12 = 0.5l; S34=1.0l;
    }
    else if (params->iproc == 5)// id = 5: theta theta -> c C
    {
        m1 = params->m_h2; m2 = params->m_h2; m3 = mqc; m4 = mqc; S12 = 0.5l; S34=1.0l;
    }
    else if (params->iproc == 6)// id = 6: theta theta -> d D
    {
        m1 = params->m_h2; m2 = params->m_h2; m3 = mqd; m4 = mqd; S12 = 0.5l; S34=1.0l;
    }
    else if (params->iproc == 7)// id = 7: theta theta -> u U
    {
        m1 = params->m_h2; m2 = params->m_h2; m3 = mqu; m4 = mqu; S12 = 0.5l; S34=1.0l;
    }
    else if (params->iproc == 8)// id = 8: theta theta -> l L
    {
        m1 = params->m_h2; m2 = params->m_h2; m3 = mtau; m4 = mtau; S12 = 0.5l; S34=1.0l;
    }
    else if (params->iproc == 9)// id = 9: theta theta -> m M
    {
        m1 = params->m_h2; m2 = params->m_h2; m3 = mmu; m4 = mmu; S12 = 0.5l; S34=1.0l;
    }
    else if (params->iproc == 10)// id = 10: theta theta -> e E
    {
        m1 = params->m_h2; m2 = params->m_h2; m3 = me; m4 = me; S12 = 0.5l; S34=1.0l;
    }
    else if (params->iproc == 11)// id = 11: theta theta -> W+ W-
    {
        m1 = params->m_h2; m2 = params->m_h2; m3 = mW; m4 = mW; S12 = 0.5l; S34=1.0l;
    }
    else if (params->iproc == 12)// id = 12: theta theta -> Z Z
    {
        m1 = params->m_h2; m2 = params->m_h2; m3 = MZ; m4 = MZ; S12 = 0.5l; S34=0.5l;
    }
    else if (params->iproc == 13)// id = 13: h2 h2 -> theta theta
    {
        m1 = params->m_h2; m2 = params->m_h2; m3 = params->m_theta; m4 = params->m_theta; S12 = 0.5l; S34=0.5l;
    }
    else
    {
        m1 = 0.l; m2 = 0.l; m3 = 0.l; m4 = 0.l; S12 = 0.l; S34=0.l;
    }
    smin = pow(std::max(m1+m2,m3+m4)+1e-4,2.0l);
    smax = params->temp*params->temp*1.0e10;
    s = smin*exp(x[0]*log(smax/smin));
    JacS = smin*log(smax/smin)*exp(x[0]*log(smax/smin));
    if ( s < smin ) return 0.0l;

    cc   = -1.0l + 2.0l*x[1];
    JacCos = 2.0l;

    func = JacS*JacCos;
    func *= params->temp*SQkaller(s,m3*m3,m4*m4)/s;
    func *= SQkaller(s,m1*m1,m2*m2)/sqrt(s)*boost::math::cyl_bessel_k(1,sqrt(s)/params->temp);
    func *= S12*S34 *2.0l*M_PI *params->SqAmp(s,cc,params->iproc);

    ff[0] = func*pow(10.0l,1000);// 10^1000

    return 0;
}
/*GET PROD RATE*/
void Rate::GetProdRate2to2(int id, long double &rate_, double &error_, double &prob_)
{
    int fail, comp, nregions, neval;
    const int ndim = 2, ncomp = 1, seed_CB =0;
    int flags = 4;
    const int nvec = 1, mineval = 10000, maxeval = 100000;
    const int nstart = 2000, nincrease = 500, nbatch = 2000, gridno = 0;
    
    long double integral[ncomp], error[ncomp], prob[ncomp];
    iproc = id;
    Vegas(ndim, ncomp, IntegrandForRate, this, nvec, epsrel_CB, epsabs_CB, flags, seed_CB, mineval, maxeval, nstart, nincrease, nbatch, gridno, NULL, NULL, &neval, &fail, integral, error, prob);
    rate_ = integral[0]/(32.0l*pow(2.0l*M_PI,6.0l))*pow(10.0l,-1000);
    error_ = error[0]/integral[0]*100.0;
    error_ = error_ > 100.0?0:error_;
    prob_ = prob[0];
}
