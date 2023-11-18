#include "Rate.hpp"
#include "gs.hpp"
#include "ge.hpp"
#include <algorithm>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "gsl/gsl_sf_bessel.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_roots.h"
#include <boost/math/special_functions/bessel.hpp>

Rate::Rate(){
    gsstar = 1.0l; gestar = 1.0l;
    alpha = 1.0l/alpha1;
    ee = sqrt(4.0l*M_PI*alpha);
    vev = pow(sqrt(2.0l)*GF,-0.5l);
    long double A = sqrt(M_PI*alpha)*vev;
    thetaW = asin(2.0*A/MZ)/2.0l;
    MW = A/sin(thetaW);
    SW = sin(thetaW);
    g= ee/sin(thetaW);
    wh2 =0.15727878782654084;
    long double vs = 1000;

}
Rate::~Rate(){};

void Rate::SetParam(long double mtheta, long double v_h, long double lam_H, long double v_s, long double lam_phi, long double lam_Hphi, long double mu_s)
{
    m_theta = mtheta; 
    vh = v_h;
    lambdaH = lam_H;
    vsigma = v_s;
    lambdaphi = lam_phi;
    lambdaHphi = lam_Hphi;
    mius = mu_s;
    mh2 = MH2(lam_Hphi, lam_H, v_s, lam_phi, v_h);
    m_h2=mh2;
    mH = MH(lambdaHphi, lam_H, v_s, lam_phi, v_h);
    a = Alpha(lam_Hphi,v_s,v_h,mH,mh2);
    HHH=HHH1(lam_Hphi, lam_H, v_s, lam_phi, v_h, a, mH, mh2);
    hHH=hHH2(lam_Hphi, lam_H, v_s, lam_phi, v_h, a, mH, mh2);
    hhh=hhh3(lam_Hphi, lam_H, v_s, lam_phi, v_h, a, mH, mh2);
    iproc = -1;
    ZpWidth = -1;
}

void Rate::SetTemp(long double _temp)
{
    temp = _temp;
    ZpWidth = -1;
}
void Rate::SetTemp(long double _temp, long double _fac)
{
    temp = _temp;
    powerfac = _fac;
    ZpWidth = -1;
}

long double Rate::Get_dof(long double _temp, int i)
// function computing the 2 effective d.o.f >>    i=1: ge, i=2: gs 
{
    switch (i)
    {
    case 1:
        return ge(_temp);
    case 2: 
        return gs(_temp);
    default:
        return gs(_temp);
    }
}

long double Rate::Get_Hubble(long double _temp)
{
long double ge = Get_dof(_temp,1), gs = Get_dof(_temp,2);
return M_PI*sqrt(ge)/3.0l/sqrt(10.0l)*_temp*_temp/MP;		
}

long double Rate::Get_Nh2(long double _temp)
{
    return 2.0l/2.0l/M_PI/M_PI*m_h2*m_h2*_temp*boost::math::cyl_bessel_k(2,(long double)m_theta/_temp);
}
