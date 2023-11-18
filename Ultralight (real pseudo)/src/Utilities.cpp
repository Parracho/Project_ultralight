#include "Utilities.hpp"

long double SQkaller(long double x, long double y, long double z){
    return sqrt( (x-pow(sqrt(y)+sqrt(z),2.l))*(x-pow(sqrt(y)-sqrt(z),2.l)));
}

long double yi(long double mi, long double s){
    return mi*mi/s;
}

long double ki(long double mi, long double s){
    return sqrt(1.l-4.l*mi*mi/s);
}

long double Aij(long double mi, long double mj, long double mN, long double s){
    return 1.l-2.l*(mN*mN/s + mi*mi/s - mj*mj/s);
}

long double rij2(long double mi, long double mj){
    return mi*mi/mj/mj;
}


long double MH(long double lambdaHphi, long double lambdaH, long double vsigma, long double lambdaphi, long double vh)
{
    return pow(2,-0.5)*pow(lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) - 
     pow(pow(lambdaH,2)*pow(vh,4) + 
       2*(-(lambdaH*lambdaphi) + 2*pow(lambdaHphi,2))*pow(vh,2)*
        pow(vsigma,2) + pow(lambdaphi,2)*pow(vsigma,4),0.5),0.5);
}

long double MH2(long double lambdaHphi, long double lambdaH, long double vsigma, long double lambdaphi, long double vh)
{
    return pow(2,-0.5)*pow(lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) + 
     pow(pow(lambdaH,2)*pow(vh,4) + 
       2*(-(lambdaH*lambdaphi) + 2*pow(lambdaHphi,2))*pow(vh,2)*
        pow(vsigma,2) + pow(lambdaphi,2)*pow(vsigma,4),0.5),0.5);
}

long double Alpha(long double vh,long double lambdaHphi,long double mh2,long double mH)
{
    return 1./2.*(asin((2.*lambdaHphi*vsigma*vh)/(pow(mH,2.)-pow(mh2,2.))));
}


// step functio
int StepF(long double x) {

  int step;
    if (x > 0.l) {step = 1;}
    else { step = 0;}

  return step;
}


long double dof(long double _temp, long double mdm, long double mS, long double mN, int i)
// function computing the 2 effective d.o.f >>    i=1: ge, i=2: gs 
{
double TnuoT, res;
long double mH = 125.;

  if (_temp > me) {
    TnuoT = 1.; 
   }
  else {
    if (i == 1) {TnuoT = pow(4./11,4./3);}
    else if (i == 2) {TnuoT = 4./11;}
   } 

  res = StepF(_temp-Lqcd)*( gmasslessboson + 8.*ggluon + 7./8.*gquark*( StepF(_temp-mqt) + StepF(_temp-mqb) + StepF(_temp-mqc) + StepF(_temp-mqs) + StepF(_temp-mqd) + StepF(_temp-mqu) ) + gscalar*(StepF(_temp-mH) + StepF(_temp-mS)) + gmassiveboson*( StepF(_temp-MZ) + 2.*StepF(_temp-mW) ) + 7./8.*( gchargedlepton*( StepF(_temp-mdm) + StepF(_temp-me) + StepF(_temp-mmu) + StepF(_temp-mtau) ) + 3.*gneutrino + 3.*gneutrino*StepF(_temp-mN) ) ) + 
        StepF(Lqcd-_temp)*( gmasslessboson + gscalar*StepF(_temp-mS) + gpion*( 2.*StepF(_temp-mpionC) + StepF(_temp-mpion0)) + 7./8.*(gchargedlepton*( StepF(_temp-mdm) + StepF(_temp-me) + StepF(_temp-mmu)) + 3.*gneutrino*TnuoT + 3.*gneutrino*StepF(_temp-mN) ) ); 

  return res;
}

long double dof_fac(long double temper, long double mdm, long double mS, long double mN, int era)
{
long double ge = dof(temper, mdm, mS, mN, 1), gs = dof(temper, mdm, mS, mN, 2);
long double dof_fac;
long double doffac = 0.;
if (era == 0 || era == 3) {doffac = 1./gs/sqrt(ge);}    
else if (era == 1) {doffac = pow(gs, -1.5);}    
else if (era == 2) {doffac = pow(ge, -2.);}
return doffac;
}


