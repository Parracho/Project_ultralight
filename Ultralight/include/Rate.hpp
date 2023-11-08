#ifndef RATE_H
#define RATE_H

#include "Utilities.hpp"
#include "cubal.h"
#include <iostream>
#include <string>

class Rate
{
public:
    Rate();
    ~Rate();
    
    //Settings
    void SetParam(long double mtheta, long double mh2); // Setting Model Parameters.
    void SetTemp(long double temp); // This is the temperature at which we want to calculate the collision term (Rate);
    void SetTemp(long double temp, long double fac); // This is the temperature at which we want to calculate the collision term (Rate);
    void PrintParameters();
    
    //Calculations
    // Rate/Collision Term
    long double SqAmp(long double s, long double x, int id);
   // void GetProdRate1to2(long double &rate);
    void GetProdRate2to2(int processid, long double &Omega, double &error, double &prob);

    // Relic Density
    void GetOmegah2From2to2FI(int processid, long double &Omega, double &error, double &prob);
    void GetOmegah2FromDecayFI(long double &Omega, double &error, double &prob);
    
    void GetOmegah2From2to2FO(int processid, long double &Omega, double &error, double &prob);
    void GetOmegah2FromDecayFO(long double &Omega, double &error, double &prob);

    // Identifier for process
    int iproc;
    
    //long double Width_Zp();
    long double ZpWidth;

// private:
    // SM parameters:
    // long double g_weak;
    // long double gp_hyper;
    // long double yt;
    // long double muh2;
    // long double lamh;
    // long double ch;

    // long double vev;
    // long double alpha;
    // long double ee;
    // long double thetaW;
    // long double MW;
    // long double gdv;
    // long double guv;

// Model Parameters:
    long double m_theta; // DM mass
    long double m_h2; // S2 mass


    long double vev;
    long double alpha;
    long double ee;
    long double thetaW;
    long double MW;
    long double g;
    long double lamH;
    long double vh;
    long double lamphi;
    long double vs;
    long double lamHphi;
    long double mus;
    long double wh2;
    long double wh;
    long double SW;


//------
    long double temp; // temperature at which we calculate the Rate
    long double powerfac; // a power factor to improve the precision.
         
    long double Get_Hubble(long double _temp);
    long double Get_Ntheta(long double _temp); 
    long double Get_Nh2(long double _temp);
    
// The effective dof at different temperature
    long double Get_dof(long double _temp, int i); // dof of energy (i=1) and of entropy (i=2) 

    long double gsstar;// = 1.0l; 
    long double gestar;// = 1.0l; 

// The parameters related to thermal history
    long double Tf;  // final temperature, below which DM yield does not change
    long double TRH; // inflationary reheat temperature

    long double m1, m2, m3, m4; // Temperory variables for external masses of each process
    long double T_min, T_max;

};


#endif
