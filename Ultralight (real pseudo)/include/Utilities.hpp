#ifndef Utilities_H
#define Utilities_H

#include <cmath>

//Kaller Function
long double SQkaller(long double x, long double y, long double z);

long double yi(long double mi, long double s);

long double ki(long double mi, long double s);

long double Aij(long double mi, long double mj, long double m_N, long double s);

long double rij2(long double mi, long double mj);

long double MH(long double lambdaHphi, long double lambdaH, long double vsigma, long double lambdaphi, long double vh);

long double MH2(long double lambdaHphi, long double lambdaH, long double vsigma, long double lambdaphi, long double vh);

long double Alpha(long double lambdaHphi,long double vsigma, long double vh, long double mH, long double mh2);

int StepF(long double x); // step function 

long double dof(long double temper, long double m_dm, long double m_S, long double m_N, int i);

long double dof_fac(long double temper, long double m_dm, long double m_S, long double m_N, int era);

// Constants:
const long double MP = 2.435e18;
//const long double mH = 125.09l;
const long double ml = 1.e-300;
const long double mnu = 1.e-300;


#define alpha1 127.9l
#define MZ 91.1876l
#define GF 0.0000116637l
#define MT 173.5l

// SM masses and dof
#define mqt 173.21
#define mqb 4.18
#define mqc 1.27
#define mqs 96.e-3
#define mqd 4.7e-3
#define mqu 2.2e-3

#define mpionC 139.57e-3
#define mpion0 134.98e-3

#define mtau 1.78
#define mmu 105.66e-3
#define me 0.00051
#define mnmu 8.6e-9
#define mntau 5.02e-8
#define mne 8.6e-12


//#define MZ 91.19
#define mW 80.38
//#define mH 125.09


#define gquark 12.
#define ggluon 2.
#define Lqcd 0.217  // lambdaQCD, GeV
#define gpion 3.
#define gchargedlepton 4.
#define gneutrino 2. // massless neutrino
#define gmassiveboson 3.
#define gmasslessboson 2.
#define gscalar 1.


#endif
