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

														/*CONSTANTS*/
Rate::Rate(){
    gsstar = 1.0l; gestar = 1.0l;
    a = 1.0l/alpha1;
    ee = sqrt(4.0l*M_PI*a);
    vev = pow(sqrt(2.0l)*GF,-0.5l);
    long double A = sqrt(M_PI*a)*vev;
    thetaW = asin(2.0*A/MZ)/2.0l;
    MW = A/sin(thetaW);
    g= ee/sin(thetaW);
    lambdaH = 0.25;
    vh = 246.;
    lambdaphi = 0.1;
    lambda12=0.1;
    vsigma = 1000.;
    lambdaHphi = 1.0E-51;
    mius = 0.1;
    wh = 4.0E-3;
    wh2 = 1.0E-5;
    SW = sqrt(0.23);
}
Rate::~Rate(){};

void Rate::SetParam(long double mtheta, long double mh2)
{
    m_theta = mtheta;
    m_h2 = mh2; 
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



long double Rate::Get_Ntheta(long double _temp)
{
    return 1.0l/2.0l/M_PI/M_PI*m_theta*m_theta*_temp*boost::math::cyl_bessel_k(2,(long double)m_theta/_temp);
}

long double Rate::Get_Nh2(long double _temp)
{
   return 1.0l/2.0l/M_PI/M_PI*m_h2*m_h2*_temp*boost::math::cyl_bessel_k(2,(long double)m_h2/_temp);
}




