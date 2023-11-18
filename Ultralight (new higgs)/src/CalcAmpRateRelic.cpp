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
    if (m_Zp > 2.*m_theta)
    {
        // Zp -> X X
        PSdm = 1.;
    }
    rate = m_Zp/(12.l*M_PI)*PSdm*Get_NZp(temp);
}

long double Rate::Width_Zp()
{
    long double tw = 0.0l;
    long double PSdm = 0.0l;
    long double PSsm = 0.0l, PSe = 0.0l;
    if (m_Zp > 2.*m_theta)
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
    
    if (id == 1)
    {
    	num=(81*pow(g,4)*pow(mW,-4)*pow(s - 2*pow(m_h2,2),2)*
    		     pow(lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) -
    		       pow(pow(lambdaH,2)*pow(vh,4) -
    		         2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		         4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		         lambdaphi*pow(vsigma,4),0.5),4)*
    		     pow(-2*(mH + ((1 - cc)*s)/2. + pow(m_h2,2)) + lambdaH*pow(vh,2) +
    		       lambdaphi*pow(vsigma,2) +
    		       pow(pow(lambdaH,2)*pow(vh,4) -
    		         2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		         4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		         lambdaphi*pow(vsigma,4),0.5),-2)*
    		     pow(2*(s - 2*pow(m_h2,2)) - 2*(mH + ((1 - cc)*s)/2. + pow(m_h2,2)) +
    		       lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		       pow(pow(lambdaH,2)*pow(vh,4) -
    		         2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		         4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		         lambdaphi*pow(vsigma,4),0.5),-2)*pow(sin(a),4));
    	den=32;
    	amp=num/den;
    }
    else if (id == 2)
    {
    	num=-3*pow(g,4)*pow(mqb,4)*(8*pow(lambdaH,4)*pow(vh,8) +
    		     32*pow(lambdaH,2)*pow(lambdaHphi,2)*pow(vh,6)*pow(vsigma,2) +
    		     8*lambdaphi*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,4) +
    		     32*lambdaH*lambdaphi*pow(lambdaHphi,2)*pow(vh,4)*pow(vsigma,4) +
    		     16*pow(lambdaHphi,4)*pow(vh,4)*pow(vsigma,4) -
    		     8*pow(lambdaH,2)*pow(lambdaphi,2)*pow(vh,4)*pow(vsigma,4) +
    		     8*lambdaphi*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,6) +
    		     8*lambdaH*pow(lambdaphi,2)*pow(vh,2)*pow(vsigma,6) +
    		     24*pow(lambdaHphi,2)*pow(lambdaphi,2)*pow(vh,2)*pow(vsigma,6) -
    		     8*lambdaH*pow(lambdaphi,3)*pow(vh,2)*pow(vsigma,6) +
    		     pow(lambdaphi,2)*pow(vsigma,8) + 6*pow(lambdaphi,3)*pow(vsigma,8) +
    		     pow(lambdaphi,4)*pow(vsigma,8) +
    		     8*pow(mqb + ((1 - cc)*s)/2. + pow(m_h2,2),4) +
    		     8*pow(lambdaH,3)*pow(vh,6)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 8*lambdaphi*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,2)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 16*lambdaH*pow(lambdaHphi,2)*pow(vh,4)*pow(vsigma,2)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 4*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,4)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 16*lambdaphi*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,4)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 4*lambdaH*pow(lambdaphi,2)*pow(vh,2)*pow(vsigma,4)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 4*pow(lambdaphi,2)*pow(vsigma,6)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 4*pow(lambdaphi,3)*pow(vsigma,6)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) - 16*pow(mqb + ((1 - cc)*s)/2. + pow(m_h2,2),3)*
    		      (lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		        pow(pow(lambdaH,2)*pow(vh,4) -
    		          2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		          4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		          lambdaphi*pow(vsigma,4),0.5)) +
    		     pow(s - 2*pow(m_h2,2),3)*(-2*(mqb + ((1 - cc)*s)/2. + pow(m_h2,2)) +
    		        lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		        pow(pow(lambdaH,2)*pow(vh,4) -
    		          2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		          4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		          lambdaphi*pow(vsigma,4),0.5)) +
    		     14*pow(mqb + ((1 - cc)*s)/2. + pow(m_h2,2),2)*
    		      (2*pow(lambdaH,2)*pow(vh,4) +
    		        2*lambdaH*pow(vh,2)*pow(pow(lambdaH,2)*pow(vh,4) -
    		           2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		           4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		           lambdaphi*pow(vsigma,4),0.5) +
    		        pow(vsigma,2)*(4*pow(lambdaHphi,2)*pow(vh,2) +
    		           lambdaphi*(1 + lambdaphi)*pow(vsigma,2) +
    		           2*lambdaphi*pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5))) -
    		     6*(mqb + ((1 - cc)*s)/2. + pow(m_h2,2))*
    		      (4*pow(lambdaH,3)*pow(vh,6) +
    		        4*pow(lambdaH,2)*pow(vh,4)*
    		         pow(pow(lambdaH,2)*pow(vh,4) -
    		           2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		           4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		           lambdaphi*pow(vsigma,4),0.5) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2)*
    		         (3*lambdaphi*pow(vsigma,2) +
    		           pow(pow(lambdaH,2)*pow(vh,4) -
    		             2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		             4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		             lambdaphi*pow(vsigma,4),0.5)) +
    		        lambdaH*pow(vh,2)*pow(vsigma,2)*
    		         (12*pow(lambdaHphi,2)*pow(vh,2) -
    		           3*(-1 + lambdaphi)*lambdaphi*pow(vsigma,2) +
    		           4*lambdaphi*pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5)) +
    		        lambdaphi*pow(vsigma,4)*
    		         (pow(lambdaphi,2)*pow(vsigma,2) +
    		           pow(pow(lambdaH,2)*pow(vh,4) -
    		             2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		             4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		             lambdaphi*pow(vsigma,4),0.5) +
    		           3*lambdaphi*(pow(vsigma,2) +
    		              pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5)))) +
    		     pow(s - 2*pow(m_h2,2),2)*(32*pow(mqb,4) + 6*pow(lambdaH,2)*pow(vh,4) +
    		        10*pow(mqb + ((1 - cc)*s)/2. + pow(m_h2,2),2) +
    		        6*lambdaH*pow(vh,2)*pow(pow(lambdaH,2)*pow(vh,4) -
    		           2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		           4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		           lambdaphi*pow(vsigma,4),0.5) +
    		        8*pow(mqb,2)*(lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		           pow(pow(lambdaH,2)*pow(vh,4) -
    		             2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		             4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		             lambdaphi*pow(vsigma,4),0.5)) +
    		        3*pow(vsigma,2)*(4*pow(lambdaHphi,2)*pow(vh,2) +
    		           lambdaphi*(1 + lambdaphi)*pow(vsigma,2) +
    		           2*lambdaphi*pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5)) -
    		        2*(mqb + ((1 - cc)*s)/2. + pow(m_h2,2))*
    		         (16*pow(mqb,2) + 5*(lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		              pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5)))) +
    		     (s - 2*pow(m_h2,2))*(12*pow(lambdaH,3)*pow(vh,6) -
    		        16*pow(mqb + ((1 - cc)*s)/2. + pow(m_h2,2),3) +
    		        4*pow(lambdaH,2)*pow(vh,4)*
    		         (4*pow(mqb,2) + 3*pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5)) +
    		        lambdaH*pow(vh,2)*(36*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		           9*lambdaphi*pow(vsigma,4) - 9*pow(lambdaphi,2)*pow(vsigma,4) +
    		           16*pow(mqb,2)*pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5) +
    		           12*lambdaphi*pow(vsigma,2)*
    		            pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5)) +
    		        8*pow(mqb + ((1 - cc)*s)/2. + pow(m_h2,2),2)*
    		         (4*pow(mqb,2) + 3*(lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		              pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5))) -
    		        2*(mqb + ((1 - cc)*s)/2. + pow(m_h2,2))*
    		         (14*pow(lambdaH,2)*pow(vh,4) + 7*pow(lambdaphi,2)*pow(vsigma,4) +
    		           2*lambdaH*pow(vh,2)*
    		            (8*pow(mqb,2) + 7*pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5)) +
    		           lambdaphi*pow(vsigma,2)*
    		            (16*pow(mqb,2) + 7*pow(vsigma,2) +
    		              14*pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5)) +
    		           4*(7*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(mqb,2)*pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5))) +
    		        pow(vsigma,2)*(8*lambdaphi*pow(mqb,2)*
    		            ((1 + lambdaphi)*pow(vsigma,2) +
    		              2*pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5)) +
    		           4*pow(lambdaHphi,2)*pow(vh,2)*
    		            (8*pow(mqb,2) + 9*lambdaphi*pow(vsigma,2) +
    		              3*pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5)) +
    		           3*lambdaphi*pow(vsigma,2)*
    		            (pow(lambdaphi,2)*pow(vsigma,2) +
    		              pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5) +
    		              3*lambdaphi*(pow(vsigma,2) +
    		                 pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5))))))*pow(sin(a),4);
    	den=16*pow(mW,4)*pow(-0.5*((1 - cc)*s) - pow(mqb,2),2)*
    			   pow(-0.5*((1 + cc)*s) - pow(mqb,2),2);
    	amp=num/den;
    }
    else if (id == 3)
    {
    	num=-3*pow(g,4)*pow(mqt,4)*(8*pow(lambdaH,4)*pow(vh,8) +
    		     32*pow(lambdaH,2)*pow(lambdaHphi,2)*pow(vh,6)*pow(vsigma,2) +
    		     8*lambdaphi*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,4) +
    		     32*lambdaH*lambdaphi*pow(lambdaHphi,2)*pow(vh,4)*pow(vsigma,4) +
    		     16*pow(lambdaHphi,4)*pow(vh,4)*pow(vsigma,4) -
    		     8*pow(lambdaH,2)*pow(lambdaphi,2)*pow(vh,4)*pow(vsigma,4) +
    		     8*lambdaphi*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,6) +
    		     8*lambdaH*pow(lambdaphi,2)*pow(vh,2)*pow(vsigma,6) +
    		     24*pow(lambdaHphi,2)*pow(lambdaphi,2)*pow(vh,2)*pow(vsigma,6) -
    		     8*lambdaH*pow(lambdaphi,3)*pow(vh,2)*pow(vsigma,6) +
    		     pow(lambdaphi,2)*pow(vsigma,8) + 6*pow(lambdaphi,3)*pow(vsigma,8) +
    		     pow(lambdaphi,4)*pow(vsigma,8) +
    		     8*pow(mqt + ((1 - cc)*s)/2. + pow(m_h2,2),4) +
    		     8*pow(lambdaH,3)*pow(vh,6)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 8*lambdaphi*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,2)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 16*lambdaH*pow(lambdaHphi,2)*pow(vh,4)*pow(vsigma,2)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 4*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,4)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 16*lambdaphi*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,4)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 4*lambdaH*pow(lambdaphi,2)*pow(vh,2)*pow(vsigma,4)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 4*pow(lambdaphi,2)*pow(vsigma,6)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 4*pow(lambdaphi,3)*pow(vsigma,6)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) - 16*pow(mqt + ((1 - cc)*s)/2. + pow(m_h2,2),3)*
    		      (lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		        pow(pow(lambdaH,2)*pow(vh,4) -
    		          2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		          4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		          lambdaphi*pow(vsigma,4),0.5)) +
    		     pow(s - 2*pow(m_h2,2),3)*(-2*(mqt + ((1 - cc)*s)/2. + pow(m_h2,2)) +
    		        lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		        pow(pow(lambdaH,2)*pow(vh,4) -
    		          2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		          4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		          lambdaphi*pow(vsigma,4),0.5)) +
    		     14*pow(mqt + ((1 - cc)*s)/2. + pow(m_h2,2),2)*
    		      (2*pow(lambdaH,2)*pow(vh,4) +
    		        2*lambdaH*pow(vh,2)*pow(pow(lambdaH,2)*pow(vh,4) -
    		           2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		           4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		           lambdaphi*pow(vsigma,4),0.5) +
    		        pow(vsigma,2)*(4*pow(lambdaHphi,2)*pow(vh,2) +
    		           lambdaphi*(1 + lambdaphi)*pow(vsigma,2) +
    		           2*lambdaphi*pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5))) -
    		     6*(mqt + ((1 - cc)*s)/2. + pow(m_h2,2))*
    		      (4*pow(lambdaH,3)*pow(vh,6) +
    		        4*pow(lambdaH,2)*pow(vh,4)*
    		         pow(pow(lambdaH,2)*pow(vh,4) -
    		           2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		           4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		           lambdaphi*pow(vsigma,4),0.5) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2)*
    		         (3*lambdaphi*pow(vsigma,2) +
    		           pow(pow(lambdaH,2)*pow(vh,4) -
    		             2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		             4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		             lambdaphi*pow(vsigma,4),0.5)) +
    		        lambdaH*pow(vh,2)*pow(vsigma,2)*
    		         (12*pow(lambdaHphi,2)*pow(vh,2) -
    		           3*(-1 + lambdaphi)*lambdaphi*pow(vsigma,2) +
    		           4*lambdaphi*pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5)) +
    		        lambdaphi*pow(vsigma,4)*
    		         (pow(lambdaphi,2)*pow(vsigma,2) +
    		           pow(pow(lambdaH,2)*pow(vh,4) -
    		             2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		             4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		             lambdaphi*pow(vsigma,4),0.5) +
    		           3*lambdaphi*(pow(vsigma,2) +
    		              pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5)))) +
    		     pow(s - 2*pow(m_h2,2),2)*(32*pow(mqt,4) + 6*pow(lambdaH,2)*pow(vh,4) +
    		        10*pow(mqt + ((1 - cc)*s)/2. + pow(m_h2,2),2) +
    		        6*lambdaH*pow(vh,2)*pow(pow(lambdaH,2)*pow(vh,4) -
    		           2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		           4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		           lambdaphi*pow(vsigma,4),0.5) +
    		        8*pow(mqt,2)*(lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		           pow(pow(lambdaH,2)*pow(vh,4) -
    		             2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		             4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		             lambdaphi*pow(vsigma,4),0.5)) +
    		        3*pow(vsigma,2)*(4*pow(lambdaHphi,2)*pow(vh,2) +
    		           lambdaphi*(1 + lambdaphi)*pow(vsigma,2) +
    		           2*lambdaphi*pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5)) -
    		        2*(mqt + ((1 - cc)*s)/2. + pow(m_h2,2))*
    		         (16*pow(mqt,2) + 5*(lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		              pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5)))) +
    		     (s - 2*pow(m_h2,2))*(12*pow(lambdaH,3)*pow(vh,6) -
    		        16*pow(mqt + ((1 - cc)*s)/2. + pow(m_h2,2),3) +
    		        4*pow(lambdaH,2)*pow(vh,4)*
    		         (4*pow(mqt,2) + 3*pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5)) +
    		        lambdaH*pow(vh,2)*(36*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		           9*lambdaphi*pow(vsigma,4) - 9*pow(lambdaphi,2)*pow(vsigma,4) +
    		           16*pow(mqt,2)*pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5) +
    		           12*lambdaphi*pow(vsigma,2)*
    		            pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5)) +
    		        8*pow(mqt + ((1 - cc)*s)/2. + pow(m_h2,2),2)*
    		         (4*pow(mqt,2) + 3*(lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		              pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5))) -
    		        2*(mqt + ((1 - cc)*s)/2. + pow(m_h2,2))*
    		         (14*pow(lambdaH,2)*pow(vh,4) + 7*pow(lambdaphi,2)*pow(vsigma,4) +
    		           2*lambdaH*pow(vh,2)*
    		            (8*pow(mqt,2) + 7*pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5)) +
    		           lambdaphi*pow(vsigma,2)*
    		            (16*pow(mqt,2) + 7*pow(vsigma,2) +
    		              14*pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5)) +
    		           4*(7*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(mqt,2)*pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5))) +
    		        pow(vsigma,2)*(8*lambdaphi*pow(mqt,2)*
    		            ((1 + lambdaphi)*pow(vsigma,2) +
    		              2*pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5)) +
    		           4*pow(lambdaHphi,2)*pow(vh,2)*
    		            (8*pow(mqt,2) + 9*lambdaphi*pow(vsigma,2) +
    		              3*pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5)) +
    		           3*lambdaphi*pow(vsigma,2)*
    		            (pow(lambdaphi,2)*pow(vsigma,2) +
    		              pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5) +
    		              3*lambdaphi*(pow(vsigma,2) +
    		                 pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5))))))*pow(sin(a),4);
    	den=16*pow(mW,4)*pow(-0.5*((1 - cc)*s) - pow(mqt,2),2)*
    			   pow(-0.5*((1 + cc)*s) - pow(mqt,2),2);
    	amp=num/den;
    }
    else if (id == 4)
    {
    	num=-3*pow(g,4)*pow(mqs,4)*(8*pow(lambdaH,4)*pow(vh,8) +
    		     32*pow(lambdaH,2)*pow(lambdaHphi,2)*pow(vh,6)*pow(vsigma,2) +
    		     8*lambdaphi*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,4) +
    		     32*lambdaH*lambdaphi*pow(lambdaHphi,2)*pow(vh,4)*pow(vsigma,4) +
    		     16*pow(lambdaHphi,4)*pow(vh,4)*pow(vsigma,4) -
    		     8*pow(lambdaH,2)*pow(lambdaphi,2)*pow(vh,4)*pow(vsigma,4) +
    		     8*lambdaphi*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,6) +
    		     8*lambdaH*pow(lambdaphi,2)*pow(vh,2)*pow(vsigma,6) +
    		     24*pow(lambdaHphi,2)*pow(lambdaphi,2)*pow(vh,2)*pow(vsigma,6) -
    		     8*lambdaH*pow(lambdaphi,3)*pow(vh,2)*pow(vsigma,6) +
    		     pow(lambdaphi,2)*pow(vsigma,8) + 6*pow(lambdaphi,3)*pow(vsigma,8) +
    		     pow(lambdaphi,4)*pow(vsigma,8) +
    		     8*pow(mqs + ((1 - cc)*s)/2. + pow(m_h2,2),4) +
    		     8*pow(lambdaH,3)*pow(vh,6)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 8*lambdaphi*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,2)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 16*lambdaH*pow(lambdaHphi,2)*pow(vh,4)*pow(vsigma,2)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 4*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,4)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 16*lambdaphi*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,4)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 4*lambdaH*pow(lambdaphi,2)*pow(vh,2)*pow(vsigma,4)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 4*pow(lambdaphi,2)*pow(vsigma,6)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 4*pow(lambdaphi,3)*pow(vsigma,6)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) - 16*pow(mqs + ((1 - cc)*s)/2. + pow(m_h2,2),3)*
    		      (lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		        pow(pow(lambdaH,2)*pow(vh,4) -
    		          2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		          4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		          lambdaphi*pow(vsigma,4),0.5)) +
    		     pow(s - 2*pow(m_h2,2),3)*(-2*(mqs + ((1 - cc)*s)/2. + pow(m_h2,2)) +
    		        lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		        pow(pow(lambdaH,2)*pow(vh,4) -
    		          2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		          4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		          lambdaphi*pow(vsigma,4),0.5)) +
    		     14*pow(mqs + ((1 - cc)*s)/2. + pow(m_h2,2),2)*
    		      (2*pow(lambdaH,2)*pow(vh,4) +
    		        2*lambdaH*pow(vh,2)*pow(pow(lambdaH,2)*pow(vh,4) -
    		           2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		           4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		           lambdaphi*pow(vsigma,4),0.5) +
    		        pow(vsigma,2)*(4*pow(lambdaHphi,2)*pow(vh,2) +
    		           lambdaphi*(1 + lambdaphi)*pow(vsigma,2) +
    		           2*lambdaphi*pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5))) -
    		     6*(mqs + ((1 - cc)*s)/2. + pow(m_h2,2))*
    		      (4*pow(lambdaH,3)*pow(vh,6) +
    		        4*pow(lambdaH,2)*pow(vh,4)*
    		         pow(pow(lambdaH,2)*pow(vh,4) -
    		           2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		           4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		           lambdaphi*pow(vsigma,4),0.5) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2)*
    		         (3*lambdaphi*pow(vsigma,2) +
    		           pow(pow(lambdaH,2)*pow(vh,4) -
    		             2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		             4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		             lambdaphi*pow(vsigma,4),0.5)) +
    		        lambdaH*pow(vh,2)*pow(vsigma,2)*
    		         (12*pow(lambdaHphi,2)*pow(vh,2) -
    		           3*(-1 + lambdaphi)*lambdaphi*pow(vsigma,2) +
    		           4*lambdaphi*pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5)) +
    		        lambdaphi*pow(vsigma,4)*
    		         (pow(lambdaphi,2)*pow(vsigma,2) +
    		           pow(pow(lambdaH,2)*pow(vh,4) -
    		             2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		             4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		             lambdaphi*pow(vsigma,4),0.5) +
    		           3*lambdaphi*(pow(vsigma,2) +
    		              pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5)))) +
    		     pow(s - 2*pow(m_h2,2),2)*(32*pow(mqs,4) + 6*pow(lambdaH,2)*pow(vh,4) +
    		        10*pow(mqs + ((1 - cc)*s)/2. + pow(m_h2,2),2) +
    		        6*lambdaH*pow(vh,2)*pow(pow(lambdaH,2)*pow(vh,4) -
    		           2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		           4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		           lambdaphi*pow(vsigma,4),0.5) +
    		        8*pow(mqs,2)*(lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		           pow(pow(lambdaH,2)*pow(vh,4) -
    		             2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		             4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		             lambdaphi*pow(vsigma,4),0.5)) +
    		        3*pow(vsigma,2)*(4*pow(lambdaHphi,2)*pow(vh,2) +
    		           lambdaphi*(1 + lambdaphi)*pow(vsigma,2) +
    		           2*lambdaphi*pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5)) -
    		        2*(mqs + ((1 - cc)*s)/2. + pow(m_h2,2))*
    		         (16*pow(mqs,2) + 5*(lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		              pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5)))) +
    		     (s - 2*pow(m_h2,2))*(12*pow(lambdaH,3)*pow(vh,6) -
    		        16*pow(mqs + ((1 - cc)*s)/2. + pow(m_h2,2),3) +
    		        4*pow(lambdaH,2)*pow(vh,4)*
    		         (4*pow(mqs,2) + 3*pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5)) +
    		        lambdaH*pow(vh,2)*(36*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		           9*lambdaphi*pow(vsigma,4) - 9*pow(lambdaphi,2)*pow(vsigma,4) +
    		           16*pow(mqs,2)*pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5) +
    		           12*lambdaphi*pow(vsigma,2)*
    		            pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5)) +
    		        8*pow(mqs + ((1 - cc)*s)/2. + pow(m_h2,2),2)*
    		         (4*pow(mqs,2) + 3*(lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		              pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5))) -
    		        2*(mqs + ((1 - cc)*s)/2. + pow(m_h2,2))*
    		         (14*pow(lambdaH,2)*pow(vh,4) + 7*pow(lambdaphi,2)*pow(vsigma,4) +
    		           2*lambdaH*pow(vh,2)*
    		            (8*pow(mqs,2) + 7*pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5)) +
    		           lambdaphi*pow(vsigma,2)*
    		            (16*pow(mqs,2) + 7*pow(vsigma,2) +
    		              14*pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5)) +
    		           4*(7*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(mqs,2)*pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5))) +
    		        pow(vsigma,2)*(8*lambdaphi*pow(mqs,2)*
    		            ((1 + lambdaphi)*pow(vsigma,2) +
    		              2*pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5)) +
    		           4*pow(lambdaHphi,2)*pow(vh,2)*
    		            (8*pow(mqs,2) + 9*lambdaphi*pow(vsigma,2) +
    		              3*pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5)) +
    		           3*lambdaphi*pow(vsigma,2)*
    		            (pow(lambdaphi,2)*pow(vsigma,2) +
    		              pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5) +
    		              3*lambdaphi*(pow(vsigma,2) +
    		                 pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5))))))*pow(sin(a),4);
    	den=16*pow(mW,4)*pow(-0.5*((1 - cc)*s) - pow(mqs,2),2)*
    			   pow(-0.5*((1 + cc)*s) - pow(mqs,2),2);
    	amp=num/den;
    }
    else if (id == 5)
    {
    	num=-3*pow(g,4)*pow(mqc,4)*(8*pow(lambdaH,4)*pow(vh,8) +
    		     32*pow(lambdaH,2)*pow(lambdaHphi,2)*pow(vh,6)*pow(vsigma,2) +
    		     8*lambdaphi*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,4) +
    		     32*lambdaH*lambdaphi*pow(lambdaHphi,2)*pow(vh,4)*pow(vsigma,4) +
    		     16*pow(lambdaHphi,4)*pow(vh,4)*pow(vsigma,4) -
    		     8*pow(lambdaH,2)*pow(lambdaphi,2)*pow(vh,4)*pow(vsigma,4) +
    		     8*lambdaphi*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,6) +
    		     8*lambdaH*pow(lambdaphi,2)*pow(vh,2)*pow(vsigma,6) +
    		     24*pow(lambdaHphi,2)*pow(lambdaphi,2)*pow(vh,2)*pow(vsigma,6) -
    		     8*lambdaH*pow(lambdaphi,3)*pow(vh,2)*pow(vsigma,6) +
    		     pow(lambdaphi,2)*pow(vsigma,8) + 6*pow(lambdaphi,3)*pow(vsigma,8) +
    		     pow(lambdaphi,4)*pow(vsigma,8) +
    		     8*pow(mqc + ((1 - cc)*s)/2. + pow(m_h2,2),4) +
    		     8*pow(lambdaH,3)*pow(vh,6)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 8*lambdaphi*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,2)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 16*lambdaH*pow(lambdaHphi,2)*pow(vh,4)*pow(vsigma,2)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 4*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,4)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 16*lambdaphi*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,4)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 4*lambdaH*pow(lambdaphi,2)*pow(vh,2)*pow(vsigma,4)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 4*pow(lambdaphi,2)*pow(vsigma,6)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 4*pow(lambdaphi,3)*pow(vsigma,6)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) - 16*pow(mqc + ((1 - cc)*s)/2. + pow(m_h2,2),3)*
    		      (lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		        pow(pow(lambdaH,2)*pow(vh,4) -
    		          2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		          4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		          lambdaphi*pow(vsigma,4),0.5)) +
    		     pow(s - 2*pow(m_h2,2),3)*(-2*(mqc + ((1 - cc)*s)/2. + pow(m_h2,2)) +
    		        lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		        pow(pow(lambdaH,2)*pow(vh,4) -
    		          2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		          4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		          lambdaphi*pow(vsigma,4),0.5)) +
    		     14*pow(mqc + ((1 - cc)*s)/2. + pow(m_h2,2),2)*
    		      (2*pow(lambdaH,2)*pow(vh,4) +
    		        2*lambdaH*pow(vh,2)*pow(pow(lambdaH,2)*pow(vh,4) -
    		           2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		           4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		           lambdaphi*pow(vsigma,4),0.5) +
    		        pow(vsigma,2)*(4*pow(lambdaHphi,2)*pow(vh,2) +
    		           lambdaphi*(1 + lambdaphi)*pow(vsigma,2) +
    		           2*lambdaphi*pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5))) -
    		     6*(mqc + ((1 - cc)*s)/2. + pow(m_h2,2))*
    		      (4*pow(lambdaH,3)*pow(vh,6) +
    		        4*pow(lambdaH,2)*pow(vh,4)*
    		         pow(pow(lambdaH,2)*pow(vh,4) -
    		           2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		           4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		           lambdaphi*pow(vsigma,4),0.5) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2)*
    		         (3*lambdaphi*pow(vsigma,2) +
    		           pow(pow(lambdaH,2)*pow(vh,4) -
    		             2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		             4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		             lambdaphi*pow(vsigma,4),0.5)) +
    		        lambdaH*pow(vh,2)*pow(vsigma,2)*
    		         (12*pow(lambdaHphi,2)*pow(vh,2) -
    		           3*(-1 + lambdaphi)*lambdaphi*pow(vsigma,2) +
    		           4*lambdaphi*pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5)) +
    		        lambdaphi*pow(vsigma,4)*
    		         (pow(lambdaphi,2)*pow(vsigma,2) +
    		           pow(pow(lambdaH,2)*pow(vh,4) -
    		             2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		             4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		             lambdaphi*pow(vsigma,4),0.5) +
    		           3*lambdaphi*(pow(vsigma,2) +
    		              pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5)))) +
    		     pow(s - 2*pow(m_h2,2),2)*(32*pow(mqc,4) + 6*pow(lambdaH,2)*pow(vh,4) +
    		        10*pow(mqc + ((1 - cc)*s)/2. + pow(m_h2,2),2) +
    		        6*lambdaH*pow(vh,2)*pow(pow(lambdaH,2)*pow(vh,4) -
    		           2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		           4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		           lambdaphi*pow(vsigma,4),0.5) +
    		        8*pow(mqc,2)*(lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		           pow(pow(lambdaH,2)*pow(vh,4) -
    		             2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		             4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		             lambdaphi*pow(vsigma,4),0.5)) +
    		        3*pow(vsigma,2)*(4*pow(lambdaHphi,2)*pow(vh,2) +
    		           lambdaphi*(1 + lambdaphi)*pow(vsigma,2) +
    		           2*lambdaphi*pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5)) -
    		        2*(mqc + ((1 - cc)*s)/2. + pow(m_h2,2))*
    		         (16*pow(mqc,2) + 5*(lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		              pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5)))) +
    		     (s - 2*pow(m_h2,2))*(12*pow(lambdaH,3)*pow(vh,6) -
    		        16*pow(mqc + ((1 - cc)*s)/2. + pow(m_h2,2),3) +
    		        4*pow(lambdaH,2)*pow(vh,4)*
    		         (4*pow(mqc,2) + 3*pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5)) +
    		        lambdaH*pow(vh,2)*(36*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		           9*lambdaphi*pow(vsigma,4) - 9*pow(lambdaphi,2)*pow(vsigma,4) +
    		           16*pow(mqc,2)*pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5) +
    		           12*lambdaphi*pow(vsigma,2)*
    		            pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5)) +
    		        8*pow(mqc + ((1 - cc)*s)/2. + pow(m_h2,2),2)*
    		         (4*pow(mqc,2) + 3*(lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		              pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5))) -
    		        2*(mqc + ((1 - cc)*s)/2. + pow(m_h2,2))*
    		         (14*pow(lambdaH,2)*pow(vh,4) + 7*pow(lambdaphi,2)*pow(vsigma,4) +
    		           2*lambdaH*pow(vh,2)*
    		            (8*pow(mqc,2) + 7*pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5)) +
    		           lambdaphi*pow(vsigma,2)*
    		            (16*pow(mqc,2) + 7*pow(vsigma,2) +
    		              14*pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5)) +
    		           4*(7*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(mqc,2)*pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5))) +
    		        pow(vsigma,2)*(8*lambdaphi*pow(mqc,2)*
    		            ((1 + lambdaphi)*pow(vsigma,2) +
    		              2*pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5)) +
    		           4*pow(lambdaHphi,2)*pow(vh,2)*
    		            (8*pow(mqc,2) + 9*lambdaphi*pow(vsigma,2) +
    		              3*pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5)) +
    		           3*lambdaphi*pow(vsigma,2)*
    		            (pow(lambdaphi,2)*pow(vsigma,2) +
    		              pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5) +
    		              3*lambdaphi*(pow(vsigma,2) +
    		                 pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5))))))*pow(sin(a),4);
    	den=16*pow(mW,4)*pow(-0.5*((1 - cc)*s) - pow(mqc,2),2)*
    			   pow(-0.5*((1 + cc)*s) - pow(mqc,2),2);
    	amp=num/den;
    }
    else if (id == 6)
    {
    	num=-3*pow(g,4)*pow(mqd,4)*(8*pow(lambdaH,4)*pow(vh,8) +
    		     32*pow(lambdaH,2)*pow(lambdaHphi,2)*pow(vh,6)*pow(vsigma,2) +
    		     8*lambdaphi*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,4) +
    		     32*lambdaH*lambdaphi*pow(lambdaHphi,2)*pow(vh,4)*pow(vsigma,4) +
    		     16*pow(lambdaHphi,4)*pow(vh,4)*pow(vsigma,4) -
    		     8*pow(lambdaH,2)*pow(lambdaphi,2)*pow(vh,4)*pow(vsigma,4) +
    		     8*lambdaphi*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,6) +
    		     8*lambdaH*pow(lambdaphi,2)*pow(vh,2)*pow(vsigma,6) +
    		     24*pow(lambdaHphi,2)*pow(lambdaphi,2)*pow(vh,2)*pow(vsigma,6) -
    		     8*lambdaH*pow(lambdaphi,3)*pow(vh,2)*pow(vsigma,6) +
    		     pow(lambdaphi,2)*pow(vsigma,8) + 6*pow(lambdaphi,3)*pow(vsigma,8) +
    		     pow(lambdaphi,4)*pow(vsigma,8) +
    		     8*pow(mqd + ((1 - cc)*s)/2. + pow(m_h2,2),4) +
    		     8*pow(lambdaH,3)*pow(vh,6)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 8*lambdaphi*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,2)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 16*lambdaH*pow(lambdaHphi,2)*pow(vh,4)*pow(vsigma,2)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 4*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,4)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 16*lambdaphi*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,4)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 4*lambdaH*pow(lambdaphi,2)*pow(vh,2)*pow(vsigma,4)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 4*pow(lambdaphi,2)*pow(vsigma,6)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 4*pow(lambdaphi,3)*pow(vsigma,6)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) - 16*pow(mqd + ((1 - cc)*s)/2. + pow(m_h2,2),3)*
    		      (lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		        pow(pow(lambdaH,2)*pow(vh,4) -
    		          2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		          4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		          lambdaphi*pow(vsigma,4),0.5)) +
    		     pow(s - 2*pow(m_h2,2),3)*(-2*(mqd + ((1 - cc)*s)/2. + pow(m_h2,2)) +
    		        lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		        pow(pow(lambdaH,2)*pow(vh,4) -
    		          2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		          4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		          lambdaphi*pow(vsigma,4),0.5)) +
    		     14*pow(mqd + ((1 - cc)*s)/2. + pow(m_h2,2),2)*
    		      (2*pow(lambdaH,2)*pow(vh,4) +
    		        2*lambdaH*pow(vh,2)*pow(pow(lambdaH,2)*pow(vh,4) -
    		           2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		           4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		           lambdaphi*pow(vsigma,4),0.5) +
    		        pow(vsigma,2)*(4*pow(lambdaHphi,2)*pow(vh,2) +
    		           lambdaphi*(1 + lambdaphi)*pow(vsigma,2) +
    		           2*lambdaphi*pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5))) -
    		     6*(mqd + ((1 - cc)*s)/2. + pow(m_h2,2))*
    		      (4*pow(lambdaH,3)*pow(vh,6) +
    		        4*pow(lambdaH,2)*pow(vh,4)*
    		         pow(pow(lambdaH,2)*pow(vh,4) -
    		           2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		           4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		           lambdaphi*pow(vsigma,4),0.5) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2)*
    		         (3*lambdaphi*pow(vsigma,2) +
    		           pow(pow(lambdaH,2)*pow(vh,4) -
    		             2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		             4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		             lambdaphi*pow(vsigma,4),0.5)) +
    		        lambdaH*pow(vh,2)*pow(vsigma,2)*
    		         (12*pow(lambdaHphi,2)*pow(vh,2) -
    		           3*(-1 + lambdaphi)*lambdaphi*pow(vsigma,2) +
    		           4*lambdaphi*pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5)) +
    		        lambdaphi*pow(vsigma,4)*
    		         (pow(lambdaphi,2)*pow(vsigma,2) +
    		           pow(pow(lambdaH,2)*pow(vh,4) -
    		             2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		             4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		             lambdaphi*pow(vsigma,4),0.5) +
    		           3*lambdaphi*(pow(vsigma,2) +
    		              pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5)))) +
    		     pow(s - 2*pow(m_h2,2),2)*(32*pow(mqd,4) + 6*pow(lambdaH,2)*pow(vh,4) +
    		        10*pow(mqd + ((1 - cc)*s)/2. + pow(m_h2,2),2) +
    		        6*lambdaH*pow(vh,2)*pow(pow(lambdaH,2)*pow(vh,4) -
    		           2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		           4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		           lambdaphi*pow(vsigma,4),0.5) +
    		        8*pow(mqd,2)*(lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		           pow(pow(lambdaH,2)*pow(vh,4) -
    		             2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		             4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		             lambdaphi*pow(vsigma,4),0.5)) +
    		        3*pow(vsigma,2)*(4*pow(lambdaHphi,2)*pow(vh,2) +
    		           lambdaphi*(1 + lambdaphi)*pow(vsigma,2) +
    		           2*lambdaphi*pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5)) -
    		        2*(mqd + ((1 - cc)*s)/2. + pow(m_h2,2))*
    		         (16*pow(mqd,2) + 5*(lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		              pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5)))) +
    		     (s - 2*pow(m_h2,2))*(12*pow(lambdaH,3)*pow(vh,6) -
    		        16*pow(mqd + ((1 - cc)*s)/2. + pow(m_h2,2),3) +
    		        4*pow(lambdaH,2)*pow(vh,4)*
    		         (4*pow(mqd,2) + 3*pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5)) +
    		        lambdaH*pow(vh,2)*(36*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		           9*lambdaphi*pow(vsigma,4) - 9*pow(lambdaphi,2)*pow(vsigma,4) +
    		           16*pow(mqd,2)*pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5) +
    		           12*lambdaphi*pow(vsigma,2)*
    		            pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5)) +
    		        8*pow(mqd + ((1 - cc)*s)/2. + pow(m_h2,2),2)*
    		         (4*pow(mqd,2) + 3*(lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		              pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5))) -
    		        2*(mqd + ((1 - cc)*s)/2. + pow(m_h2,2))*
    		         (14*pow(lambdaH,2)*pow(vh,4) + 7*pow(lambdaphi,2)*pow(vsigma,4) +
    		           2*lambdaH*pow(vh,2)*
    		            (8*pow(mqd,2) + 7*pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5)) +
    		           lambdaphi*pow(vsigma,2)*
    		            (16*pow(mqd,2) + 7*pow(vsigma,2) +
    		              14*pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5)) +
    		           4*(7*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(mqd,2)*pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5))) +
    		        pow(vsigma,2)*(8*lambdaphi*pow(mqd,2)*
    		            ((1 + lambdaphi)*pow(vsigma,2) +
    		              2*pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5)) +
    		           4*pow(lambdaHphi,2)*pow(vh,2)*
    		            (8*pow(mqd,2) + 9*lambdaphi*pow(vsigma,2) +
    		              3*pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5)) +
    		           3*lambdaphi*pow(vsigma,2)*
    		            (pow(lambdaphi,2)*pow(vsigma,2) +
    		              pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5) +
    		              3*lambdaphi*(pow(vsigma,2) +
    		                 pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5))))))*pow(sin(a),4);
    	den=16*pow(mW,4)*pow(-0.5*((1 - cc)*s) - pow(mqd,2),2)*
    			   pow(-0.5*((1 + cc)*s) - pow(mqd,2),2);
    	amp=num/den;
    }
    else if (id == 7)
    {
    	num=-3*pow(g,4)*pow(mqu,4)*(8*pow(lambdaH,4)*pow(vh,8) +
    		     32*pow(lambdaH,2)*pow(lambdaHphi,2)*pow(vh,6)*pow(vsigma,2) +
    		     8*lambdaphi*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,4) +
    		     32*lambdaH*lambdaphi*pow(lambdaHphi,2)*pow(vh,4)*pow(vsigma,4) +
    		     16*pow(lambdaHphi,4)*pow(vh,4)*pow(vsigma,4) -
    		     8*pow(lambdaH,2)*pow(lambdaphi,2)*pow(vh,4)*pow(vsigma,4) +
    		     8*lambdaphi*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,6) +
    		     8*lambdaH*pow(lambdaphi,2)*pow(vh,2)*pow(vsigma,6) +
    		     24*pow(lambdaHphi,2)*pow(lambdaphi,2)*pow(vh,2)*pow(vsigma,6) -
    		     8*lambdaH*pow(lambdaphi,3)*pow(vh,2)*pow(vsigma,6) +
    		     pow(lambdaphi,2)*pow(vsigma,8) + 6*pow(lambdaphi,3)*pow(vsigma,8) +
    		     pow(lambdaphi,4)*pow(vsigma,8) +
    		     8*pow(mqu + ((1 - cc)*s)/2. + pow(m_h2,2),4) +
    		     8*pow(lambdaH,3)*pow(vh,6)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 8*lambdaphi*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,2)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 16*lambdaH*pow(lambdaHphi,2)*pow(vh,4)*pow(vsigma,2)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 4*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,4)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 16*lambdaphi*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,4)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 4*lambdaH*pow(lambdaphi,2)*pow(vh,2)*pow(vsigma,4)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 4*pow(lambdaphi,2)*pow(vsigma,6)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 4*pow(lambdaphi,3)*pow(vsigma,6)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) - 16*pow(mqu + ((1 - cc)*s)/2. + pow(m_h2,2),3)*
    		      (lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		        pow(pow(lambdaH,2)*pow(vh,4) -
    		          2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		          4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		          lambdaphi*pow(vsigma,4),0.5)) +
    		     pow(s - 2*pow(m_h2,2),3)*(-2*(mqu + ((1 - cc)*s)/2. + pow(m_h2,2)) +
    		        lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		        pow(pow(lambdaH,2)*pow(vh,4) -
    		          2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		          4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		          lambdaphi*pow(vsigma,4),0.5)) +
    		     14*pow(mqu + ((1 - cc)*s)/2. + pow(m_h2,2),2)*
    		      (2*pow(lambdaH,2)*pow(vh,4) +
    		        2*lambdaH*pow(vh,2)*pow(pow(lambdaH,2)*pow(vh,4) -
    		           2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		           4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		           lambdaphi*pow(vsigma,4),0.5) +
    		        pow(vsigma,2)*(4*pow(lambdaHphi,2)*pow(vh,2) +
    		           lambdaphi*(1 + lambdaphi)*pow(vsigma,2) +
    		           2*lambdaphi*pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5))) -
    		     6*(mqu + ((1 - cc)*s)/2. + pow(m_h2,2))*
    		      (4*pow(lambdaH,3)*pow(vh,6) +
    		        4*pow(lambdaH,2)*pow(vh,4)*
    		         pow(pow(lambdaH,2)*pow(vh,4) -
    		           2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		           4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		           lambdaphi*pow(vsigma,4),0.5) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2)*
    		         (3*lambdaphi*pow(vsigma,2) +
    		           pow(pow(lambdaH,2)*pow(vh,4) -
    		             2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		             4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		             lambdaphi*pow(vsigma,4),0.5)) +
    		        lambdaH*pow(vh,2)*pow(vsigma,2)*
    		         (12*pow(lambdaHphi,2)*pow(vh,2) -
    		           3*(-1 + lambdaphi)*lambdaphi*pow(vsigma,2) +
    		           4*lambdaphi*pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5)) +
    		        lambdaphi*pow(vsigma,4)*
    		         (pow(lambdaphi,2)*pow(vsigma,2) +
    		           pow(pow(lambdaH,2)*pow(vh,4) -
    		             2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		             4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		             lambdaphi*pow(vsigma,4),0.5) +
    		           3*lambdaphi*(pow(vsigma,2) +
    		              pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5)))) +
    		     pow(s - 2*pow(m_h2,2),2)*(32*pow(mqu,4) + 6*pow(lambdaH,2)*pow(vh,4) +
    		        10*pow(mqu + ((1 - cc)*s)/2. + pow(m_h2,2),2) +
    		        6*lambdaH*pow(vh,2)*pow(pow(lambdaH,2)*pow(vh,4) -
    		           2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		           4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		           lambdaphi*pow(vsigma,4),0.5) +
    		        8*pow(mqu,2)*(lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		           pow(pow(lambdaH,2)*pow(vh,4) -
    		             2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		             4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		             lambdaphi*pow(vsigma,4),0.5)) +
    		        3*pow(vsigma,2)*(4*pow(lambdaHphi,2)*pow(vh,2) +
    		           lambdaphi*(1 + lambdaphi)*pow(vsigma,2) +
    		           2*lambdaphi*pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5)) -
    		        2*(mqu + ((1 - cc)*s)/2. + pow(m_h2,2))*
    		         (16*pow(mqu,2) + 5*(lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		              pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5)))) +
    		     (s - 2*pow(m_h2,2))*(12*pow(lambdaH,3)*pow(vh,6) -
    		        16*pow(mqu + ((1 - cc)*s)/2. + pow(m_h2,2),3) +
    		        4*pow(lambdaH,2)*pow(vh,4)*
    		         (4*pow(mqu,2) + 3*pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5)) +
    		        lambdaH*pow(vh,2)*(36*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		           9*lambdaphi*pow(vsigma,4) - 9*pow(lambdaphi,2)*pow(vsigma,4) +
    		           16*pow(mqu,2)*pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5) +
    		           12*lambdaphi*pow(vsigma,2)*
    		            pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5)) +
    		        8*pow(mqu + ((1 - cc)*s)/2. + pow(m_h2,2),2)*
    		         (4*pow(mqu,2) + 3*(lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		              pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5))) -
    		        2*(mqu + ((1 - cc)*s)/2. + pow(m_h2,2))*
    		         (14*pow(lambdaH,2)*pow(vh,4) + 7*pow(lambdaphi,2)*pow(vsigma,4) +
    		           2*lambdaH*pow(vh,2)*
    		            (8*pow(mqu,2) + 7*pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5)) +
    		           lambdaphi*pow(vsigma,2)*
    		            (16*pow(mqu,2) + 7*pow(vsigma,2) +
    		              14*pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5)) +
    		           4*(7*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(mqu,2)*pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5))) +
    		        pow(vsigma,2)*(8*lambdaphi*pow(mqu,2)*
    		            ((1 + lambdaphi)*pow(vsigma,2) +
    		              2*pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5)) +
    		           4*pow(lambdaHphi,2)*pow(vh,2)*
    		            (8*pow(mqu,2) + 9*lambdaphi*pow(vsigma,2) +
    		              3*pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5)) +
    		           3*lambdaphi*pow(vsigma,2)*
    		            (pow(lambdaphi,2)*pow(vsigma,2) +
    		              pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5) +
    		              3*lambdaphi*(pow(vsigma,2) +
    		                 pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5))))))*pow(sin(a),4);
    	den=16*pow(mW,4)*pow(-0.5*((1 - cc)*s) - pow(mqu,2),2)*
    			   pow(-0.5*((1 + cc)*s) - pow(mqu,2),2);
    	amp=num/den;
       
    }
    else if (id == 8)
    {
    	num=-(pow(g,4)*pow(mtau,4)*(8*pow(lambdaH,4)*pow(vh,8) +
    		       32*pow(lambdaH,2)*pow(lambdaHphi,2)*pow(vh,6)*pow(vsigma,2) +
    		       8*lambdaphi*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,4) +
    		       32*lambdaH*lambdaphi*pow(lambdaHphi,2)*pow(vh,4)*pow(vsigma,4) +
    		       16*pow(lambdaHphi,4)*pow(vh,4)*pow(vsigma,4) -
    		       8*pow(lambdaH,2)*pow(lambdaphi,2)*pow(vh,4)*pow(vsigma,4) +
    		       8*lambdaphi*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,6) +
    		       8*lambdaH*pow(lambdaphi,2)*pow(vh,2)*pow(vsigma,6) +
    		       24*pow(lambdaHphi,2)*pow(lambdaphi,2)*pow(vh,2)*pow(vsigma,6) -
    		       8*lambdaH*pow(lambdaphi,3)*pow(vh,2)*pow(vsigma,6) +
    		       pow(lambdaphi,2)*pow(vsigma,8) + 6*pow(lambdaphi,3)*pow(vsigma,8) +
    		       pow(lambdaphi,4)*pow(vsigma,8) +
    		       8*pow(mtau + ((1 - cc)*s)/2. + pow(m_h2,2),4) +
    		       8*pow(lambdaH,3)*pow(vh,6)*
    		        pow(pow(lambdaH,2)*pow(vh,4) -
    		          2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		          4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		          lambdaphi*pow(vsigma,4),0.5) +
    		       8*lambdaphi*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,2)*
    		        pow(pow(lambdaH,2)*pow(vh,4) -
    		          2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		          4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		          lambdaphi*pow(vsigma,4),0.5) +
    		       16*lambdaH*pow(lambdaHphi,2)*pow(vh,4)*pow(vsigma,2)*
    		        pow(pow(lambdaH,2)*pow(vh,4) -
    		          2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		          4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		          lambdaphi*pow(vsigma,4),0.5) +
    		       4*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,4)*
    		        pow(pow(lambdaH,2)*pow(vh,4) -
    		          2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		          4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		          lambdaphi*pow(vsigma,4),0.5) +
    		       16*lambdaphi*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,4)*
    		        pow(pow(lambdaH,2)*pow(vh,4) -
    		          2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		          4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		          lambdaphi*pow(vsigma,4),0.5) +
    		       4*lambdaH*pow(lambdaphi,2)*pow(vh,2)*pow(vsigma,4)*
    		        pow(pow(lambdaH,2)*pow(vh,4) -
    		          2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		          4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		          lambdaphi*pow(vsigma,4),0.5) +
    		       4*pow(lambdaphi,2)*pow(vsigma,6)*
    		        pow(pow(lambdaH,2)*pow(vh,4) -
    		          2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		          4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		          lambdaphi*pow(vsigma,4),0.5) +
    		       4*pow(lambdaphi,3)*pow(vsigma,6)*
    		        pow(pow(lambdaH,2)*pow(vh,4) -
    		          2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		          4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		          lambdaphi*pow(vsigma,4),0.5) -
    		       16*pow(mtau + ((1 - cc)*s)/2. + pow(m_h2,2),3)*
    		        (lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		          pow(pow(lambdaH,2)*pow(vh,4) -
    		            2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		            4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		            lambdaphi*pow(vsigma,4),0.5)) +
    		       pow(s - 2*pow(m_h2,2),3)*
    		        (-2*(mtau + ((1 - cc)*s)/2. + pow(m_h2,2)) + lambdaH*pow(vh,2) +
    		          lambdaphi*pow(vsigma,2) +
    		          pow(pow(lambdaH,2)*pow(vh,4) -
    		            2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		            4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		            lambdaphi*pow(vsigma,4),0.5)) +
    		       14*pow(mtau + ((1 - cc)*s)/2. + pow(m_h2,2),2)*
    		        (2*pow(lambdaH,2)*pow(vh,4) +
    		          2*lambdaH*pow(vh,2)*pow(pow(lambdaH,2)*pow(vh,4) -
    		             2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		             4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		             lambdaphi*pow(vsigma,4),0.5) +
    		          pow(vsigma,2)*(4*pow(lambdaHphi,2)*pow(vh,2) +
    		             lambdaphi*(1 + lambdaphi)*pow(vsigma,2) +
    		             2*lambdaphi*pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5))) -
    		       6*(mtau + ((1 - cc)*s)/2. + pow(m_h2,2))*
    		        (4*pow(lambdaH,3)*pow(vh,6) +
    		          4*pow(lambdaH,2)*pow(vh,4)*
    		           pow(pow(lambdaH,2)*pow(vh,4) -
    		             2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		             4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		             lambdaphi*pow(vsigma,4),0.5) +
    		          4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2)*
    		           (3*lambdaphi*pow(vsigma,2) +
    		             pow(pow(lambdaH,2)*pow(vh,4) -
    		               2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		               4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		               lambdaphi*pow(vsigma,4),0.5)) +
    		          lambdaH*pow(vh,2)*pow(vsigma,2)*
    		           (12*pow(lambdaHphi,2)*pow(vh,2) -
    		             3*(-1 + lambdaphi)*lambdaphi*pow(vsigma,2) +
    		             4*lambdaphi*pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5)) +
    		          lambdaphi*pow(vsigma,4)*
    		           (pow(lambdaphi,2)*pow(vsigma,2) +
    		             pow(pow(lambdaH,2)*pow(vh,4) -
    		               2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		               4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		               lambdaphi*pow(vsigma,4),0.5) +
    		             3*lambdaphi*(pow(vsigma,2) +
    		                pow(pow(lambdaH,2)*pow(vh,4) -
    		                  2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                  4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                  lambdaphi*pow(vsigma,4),0.5)))) +
    		       pow(s - 2*pow(m_h2,2),2)*
    		        (32*pow(mtau,4) + 6*pow(lambdaH,2)*pow(vh,4) +
    		          10*pow(mtau + ((1 - cc)*s)/2. + pow(m_h2,2),2) +
    		          6*lambdaH*pow(vh,2)*pow(pow(lambdaH,2)*pow(vh,4) -
    		             2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		             4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		             lambdaphi*pow(vsigma,4),0.5) +
    		          8*pow(mtau,2)*(lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		             pow(pow(lambdaH,2)*pow(vh,4) -
    		               2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		               4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		               lambdaphi*pow(vsigma,4),0.5)) +
    		          3*pow(vsigma,2)*(4*pow(lambdaHphi,2)*pow(vh,2) +
    		             lambdaphi*(1 + lambdaphi)*pow(vsigma,2) +
    		             2*lambdaphi*pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5)) -
    		          2*(mtau + ((1 - cc)*s)/2. + pow(m_h2,2))*
    		           (16*pow(mtau,2) + 5*
    		              (lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		                pow(pow(lambdaH,2)*pow(vh,4) -
    		                  2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                  4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                  lambdaphi*pow(vsigma,4),0.5)))) +
    		       (s - 2*pow(m_h2,2))*(12*pow(lambdaH,3)*pow(vh,6) -
    		          16*pow(mtau + ((1 - cc)*s)/2. + pow(m_h2,2),3) +
    		          4*pow(lambdaH,2)*pow(vh,4)*
    		           (4*pow(mtau,2) + 3*pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5)) +
    		          lambdaH*pow(vh,2)*(36*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		             9*lambdaphi*pow(vsigma,4) - 9*pow(lambdaphi,2)*pow(vsigma,4) +
    		             16*pow(mtau,2)*pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5) +
    		             12*lambdaphi*pow(vsigma,2)*
    		              pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5)) +
    		          8*pow(mtau + ((1 - cc)*s)/2. + pow(m_h2,2),2)*
    		           (4*pow(mtau,2) + 3*(lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		                pow(pow(lambdaH,2)*pow(vh,4) -
    		                  2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                  4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                  lambdaphi*pow(vsigma,4),0.5))) -
    		          2*(mtau + ((1 - cc)*s)/2. + pow(m_h2,2))*
    		           (14*pow(lambdaH,2)*pow(vh,4) + 7*pow(lambdaphi,2)*pow(vsigma,4) +
    		             2*lambdaH*pow(vh,2)*
    		              (8*pow(mtau,2) +
    		                7*pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5)) +
    		             lambdaphi*pow(vsigma,2)*
    		              (16*pow(mtau,2) + 7*pow(vsigma,2) +
    		                14*pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5)) +
    		             4*(7*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(mtau,2)*pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5))) +
    		          pow(vsigma,2)*(8*lambdaphi*pow(mtau,2)*
    		              ((1 + lambdaphi)*pow(vsigma,2) +
    		                2*pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5)) +
    		             4*pow(lambdaHphi,2)*pow(vh,2)*
    		              (8*pow(mtau,2) + 9*lambdaphi*pow(vsigma,2) +
    		                3*pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5)) +
    		             3*lambdaphi*pow(vsigma,2)*
    		              (pow(lambdaphi,2)*pow(vsigma,2) +
    		                pow(pow(lambdaH,2)*pow(vh,4) -
    		                  2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                  4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                  lambdaphi*pow(vsigma,4),0.5) +
    		                3*lambdaphi*(pow(vsigma,2) +
    		                   pow(pow(lambdaH,2)*pow(vh,4) -
    		                     2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                     4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                     lambdaphi*pow(vsigma,4),0.5))))))*pow(sin(a),4));
    	den=16*pow(mW,4)*pow(-0.5*((1 - cc)*s) - pow(mtau,2),2)*
    			   pow(-0.5*((1 + cc)*s) - pow(mtau,2),2);
    	amp=num/den;
    }
    else if (id == 9)
    {
    	num=-(pow(g,4)*pow(mmu,4)*(8*pow(lambdaH,4)*pow(vh,8) +
    		       32*pow(lambdaH,2)*pow(lambdaHphi,2)*pow(vh,6)*pow(vsigma,2) +
    		       8*lambdaphi*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,4) +
    		       32*lambdaH*lambdaphi*pow(lambdaHphi,2)*pow(vh,4)*pow(vsigma,4) +
    		       16*pow(lambdaHphi,4)*pow(vh,4)*pow(vsigma,4) -
    		       8*pow(lambdaH,2)*pow(lambdaphi,2)*pow(vh,4)*pow(vsigma,4) +
    		       8*lambdaphi*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,6) +
    		       8*lambdaH*pow(lambdaphi,2)*pow(vh,2)*pow(vsigma,6) +
    		       24*pow(lambdaHphi,2)*pow(lambdaphi,2)*pow(vh,2)*pow(vsigma,6) -
    		       8*lambdaH*pow(lambdaphi,3)*pow(vh,2)*pow(vsigma,6) +
    		       pow(lambdaphi,2)*pow(vsigma,8) + 6*pow(lambdaphi,3)*pow(vsigma,8) +
    		       pow(lambdaphi,4)*pow(vsigma,8) +
    		       8*pow(mmu + ((1 - cc)*s)/2. + pow(m_h2,2),4) +
    		       8*pow(lambdaH,3)*pow(vh,6)*
    		        pow(pow(lambdaH,2)*pow(vh,4) -
    		          2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		          4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		          lambdaphi*pow(vsigma,4),0.5) +
    		       8*lambdaphi*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,2)*
    		        pow(pow(lambdaH,2)*pow(vh,4) -
    		          2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		          4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		          lambdaphi*pow(vsigma,4),0.5) +
    		       16*lambdaH*pow(lambdaHphi,2)*pow(vh,4)*pow(vsigma,2)*
    		        pow(pow(lambdaH,2)*pow(vh,4) -
    		          2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		          4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		          lambdaphi*pow(vsigma,4),0.5) +
    		       4*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,4)*
    		        pow(pow(lambdaH,2)*pow(vh,4) -
    		          2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		          4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		          lambdaphi*pow(vsigma,4),0.5) +
    		       16*lambdaphi*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,4)*
    		        pow(pow(lambdaH,2)*pow(vh,4) -
    		          2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		          4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		          lambdaphi*pow(vsigma,4),0.5) +
    		       4*lambdaH*pow(lambdaphi,2)*pow(vh,2)*pow(vsigma,4)*
    		        pow(pow(lambdaH,2)*pow(vh,4) -
    		          2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		          4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		          lambdaphi*pow(vsigma,4),0.5) +
    		       4*pow(lambdaphi,2)*pow(vsigma,6)*
    		        pow(pow(lambdaH,2)*pow(vh,4) -
    		          2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		          4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		          lambdaphi*pow(vsigma,4),0.5) +
    		       4*pow(lambdaphi,3)*pow(vsigma,6)*
    		        pow(pow(lambdaH,2)*pow(vh,4) -
    		          2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		          4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		          lambdaphi*pow(vsigma,4),0.5) -
    		       16*pow(mmu + ((1 - cc)*s)/2. + pow(m_h2,2),3)*
    		        (lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		          pow(pow(lambdaH,2)*pow(vh,4) -
    		            2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		            4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		            lambdaphi*pow(vsigma,4),0.5)) +
    		       pow(s - 2*pow(m_h2,2),3)*
    		        (-2*(mmu + ((1 - cc)*s)/2. + pow(m_h2,2)) + lambdaH*pow(vh,2) +
    		          lambdaphi*pow(vsigma,2) +
    		          pow(pow(lambdaH,2)*pow(vh,4) -
    		            2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		            4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		            lambdaphi*pow(vsigma,4),0.5)) +
    		       14*pow(mmu + ((1 - cc)*s)/2. + pow(m_h2,2),2)*
    		        (2*pow(lambdaH,2)*pow(vh,4) +
    		          2*lambdaH*pow(vh,2)*pow(pow(lambdaH,2)*pow(vh,4) -
    		             2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		             4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		             lambdaphi*pow(vsigma,4),0.5) +
    		          pow(vsigma,2)*(4*pow(lambdaHphi,2)*pow(vh,2) +
    		             lambdaphi*(1 + lambdaphi)*pow(vsigma,2) +
    		             2*lambdaphi*pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5))) -
    		       6*(mmu + ((1 - cc)*s)/2. + pow(m_h2,2))*
    		        (4*pow(lambdaH,3)*pow(vh,6) +
    		          4*pow(lambdaH,2)*pow(vh,4)*
    		           pow(pow(lambdaH,2)*pow(vh,4) -
    		             2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		             4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		             lambdaphi*pow(vsigma,4),0.5) +
    		          4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2)*
    		           (3*lambdaphi*pow(vsigma,2) +
    		             pow(pow(lambdaH,2)*pow(vh,4) -
    		               2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		               4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		               lambdaphi*pow(vsigma,4),0.5)) +
    		          lambdaH*pow(vh,2)*pow(vsigma,2)*
    		           (12*pow(lambdaHphi,2)*pow(vh,2) -
    		             3*(-1 + lambdaphi)*lambdaphi*pow(vsigma,2) +
    		             4*lambdaphi*pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5)) +
    		          lambdaphi*pow(vsigma,4)*
    		           (pow(lambdaphi,2)*pow(vsigma,2) +
    		             pow(pow(lambdaH,2)*pow(vh,4) -
    		               2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		               4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		               lambdaphi*pow(vsigma,4),0.5) +
    		             3*lambdaphi*(pow(vsigma,2) +
    		                pow(pow(lambdaH,2)*pow(vh,4) -
    		                  2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                  4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                  lambdaphi*pow(vsigma,4),0.5)))) +
    		       pow(s - 2*pow(m_h2,2),2)*
    		        (32*pow(mmu,4) + 6*pow(lambdaH,2)*pow(vh,4) +
    		          10*pow(mmu + ((1 - cc)*s)/2. + pow(m_h2,2),2) +
    		          6*lambdaH*pow(vh,2)*pow(pow(lambdaH,2)*pow(vh,4) -
    		             2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		             4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		             lambdaphi*pow(vsigma,4),0.5) +
    		          8*pow(mmu,2)*(lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		             pow(pow(lambdaH,2)*pow(vh,4) -
    		               2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		               4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		               lambdaphi*pow(vsigma,4),0.5)) +
    		          3*pow(vsigma,2)*(4*pow(lambdaHphi,2)*pow(vh,2) +
    		             lambdaphi*(1 + lambdaphi)*pow(vsigma,2) +
    		             2*lambdaphi*pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5)) -
    		          2*(mmu + ((1 - cc)*s)/2. + pow(m_h2,2))*
    		           (16*pow(mmu,2) + 5*(lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		                pow(pow(lambdaH,2)*pow(vh,4) -
    		                  2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                  4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                  lambdaphi*pow(vsigma,4),0.5)))) +
    		       (s - 2*pow(m_h2,2))*(12*pow(lambdaH,3)*pow(vh,6) -
    		          16*pow(mmu + ((1 - cc)*s)/2. + pow(m_h2,2),3) +
    		          4*pow(lambdaH,2)*pow(vh,4)*
    		           (4*pow(mmu,2) + 3*pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5)) +
    		          lambdaH*pow(vh,2)*(36*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		             9*lambdaphi*pow(vsigma,4) - 9*pow(lambdaphi,2)*pow(vsigma,4) +
    		             16*pow(mmu,2)*pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5) +
    		             12*lambdaphi*pow(vsigma,2)*
    		              pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5)) +
    		          8*pow(mmu + ((1 - cc)*s)/2. + pow(m_h2,2),2)*
    		           (4*pow(mmu,2) + 3*(lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		                pow(pow(lambdaH,2)*pow(vh,4) -
    		                  2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                  4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                  lambdaphi*pow(vsigma,4),0.5))) -
    		          2*(mmu + ((1 - cc)*s)/2. + pow(m_h2,2))*
    		           (14*pow(lambdaH,2)*pow(vh,4) + 7*pow(lambdaphi,2)*pow(vsigma,4) +
    		             2*lambdaH*pow(vh,2)*
    		              (8*pow(mmu,2) + 7*
    		                 pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5)) +
    		             lambdaphi*pow(vsigma,2)*
    		              (16*pow(mmu,2) + 7*pow(vsigma,2) +
    		                14*pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5)) +
    		             4*(7*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(mmu,2)*pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5))) +
    		          pow(vsigma,2)*(8*lambdaphi*pow(mmu,2)*
    		              ((1 + lambdaphi)*pow(vsigma,2) +
    		                2*pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5)) +
    		             4*pow(lambdaHphi,2)*pow(vh,2)*
    		              (8*pow(mmu,2) + 9*lambdaphi*pow(vsigma,2) +
    		                3*pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5)) +
    		             3*lambdaphi*pow(vsigma,2)*
    		              (pow(lambdaphi,2)*pow(vsigma,2) +
    		                pow(pow(lambdaH,2)*pow(vh,4) -
    		                  2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                  4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                  lambdaphi*pow(vsigma,4),0.5) +
    		                3*lambdaphi*(pow(vsigma,2) +
    		                   pow(pow(lambdaH,2)*pow(vh,4) -
    		                     2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                     4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                     lambdaphi*pow(vsigma,4),0.5))))))*pow(sin(a),4));
    	den=16*pow(mW,4)*pow(-0.5*((1 - cc)*s) - pow(mmu,2),2)*
    			   pow(-0.5*((1 + cc)*s) - pow(mmu,2),2);
    	amp=num/den;
    }
    else if (id == 10)
    {
    	num=-(pow(g,4)*pow(me,4)*(8*pow(lambdaH,4)*pow(vh,8) +
    		       32*pow(lambdaH,2)*pow(lambdaHphi,2)*pow(vh,6)*pow(vsigma,2) +
    		       8*lambdaphi*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,4) +
    		       32*lambdaH*lambdaphi*pow(lambdaHphi,2)*pow(vh,4)*pow(vsigma,4) +
    		       16*pow(lambdaHphi,4)*pow(vh,4)*pow(vsigma,4) -
    		       8*pow(lambdaH,2)*pow(lambdaphi,2)*pow(vh,4)*pow(vsigma,4) +
    		       8*lambdaphi*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,6) +
    		       8*lambdaH*pow(lambdaphi,2)*pow(vh,2)*pow(vsigma,6) +
    		       24*pow(lambdaHphi,2)*pow(lambdaphi,2)*pow(vh,2)*pow(vsigma,6) -
    		       8*lambdaH*pow(lambdaphi,3)*pow(vh,2)*pow(vsigma,6) +
    		       pow(lambdaphi,2)*pow(vsigma,8) + 6*pow(lambdaphi,3)*pow(vsigma,8) +
    		       pow(lambdaphi,4)*pow(vsigma,8) +
    		       8*pow(me + ((1 - cc)*s)/2. + pow(m_h2,2),4) +
    		       8*pow(lambdaH,3)*pow(vh,6)*
    		        pow(pow(lambdaH,2)*pow(vh,4) -
    		          2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		          4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		          lambdaphi*pow(vsigma,4),0.5) +
    		       8*lambdaphi*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,2)*
    		        pow(pow(lambdaH,2)*pow(vh,4) -
    		          2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		          4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		          lambdaphi*pow(vsigma,4),0.5) +
    		       16*lambdaH*pow(lambdaHphi,2)*pow(vh,4)*pow(vsigma,2)*
    		        pow(pow(lambdaH,2)*pow(vh,4) -
    		          2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		          4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		          lambdaphi*pow(vsigma,4),0.5) +
    		       4*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,4)*
    		        pow(pow(lambdaH,2)*pow(vh,4) -
    		          2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		          4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		          lambdaphi*pow(vsigma,4),0.5) +
    		       16*lambdaphi*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,4)*
    		        pow(pow(lambdaH,2)*pow(vh,4) -
    		          2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		          4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		          lambdaphi*pow(vsigma,4),0.5) +
    		       4*lambdaH*pow(lambdaphi,2)*pow(vh,2)*pow(vsigma,4)*
    		        pow(pow(lambdaH,2)*pow(vh,4) -
    		          2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		          4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		          lambdaphi*pow(vsigma,4),0.5) +
    		       4*pow(lambdaphi,2)*pow(vsigma,6)*
    		        pow(pow(lambdaH,2)*pow(vh,4) -
    		          2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		          4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		          lambdaphi*pow(vsigma,4),0.5) +
    		       4*pow(lambdaphi,3)*pow(vsigma,6)*
    		        pow(pow(lambdaH,2)*pow(vh,4) -
    		          2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		          4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		          lambdaphi*pow(vsigma,4),0.5) -
    		       16*pow(me + ((1 - cc)*s)/2. + pow(m_h2,2),3)*
    		        (lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		          pow(pow(lambdaH,2)*pow(vh,4) -
    		            2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		            4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		            lambdaphi*pow(vsigma,4),0.5)) +
    		       pow(s - 2*pow(m_h2,2),3)*
    		        (-2*(me + ((1 - cc)*s)/2. + pow(m_h2,2)) + lambdaH*pow(vh,2) +
    		          lambdaphi*pow(vsigma,2) +
    		          pow(pow(lambdaH,2)*pow(vh,4) -
    		            2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		            4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		            lambdaphi*pow(vsigma,4),0.5)) +
    		       14*pow(me + ((1 - cc)*s)/2. + pow(m_h2,2),2)*
    		        (2*pow(lambdaH,2)*pow(vh,4) +
    		          2*lambdaH*pow(vh,2)*pow(pow(lambdaH,2)*pow(vh,4) -
    		             2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		             4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		             lambdaphi*pow(vsigma,4),0.5) +
    		          pow(vsigma,2)*(4*pow(lambdaHphi,2)*pow(vh,2) +
    		             lambdaphi*(1 + lambdaphi)*pow(vsigma,2) +
    		             2*lambdaphi*pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5))) -
    		       6*(me + ((1 - cc)*s)/2. + pow(m_h2,2))*
    		        (4*pow(lambdaH,3)*pow(vh,6) +
    		          4*pow(lambdaH,2)*pow(vh,4)*
    		           pow(pow(lambdaH,2)*pow(vh,4) -
    		             2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		             4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		             lambdaphi*pow(vsigma,4),0.5) +
    		          4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2)*
    		           (3*lambdaphi*pow(vsigma,2) +
    		             pow(pow(lambdaH,2)*pow(vh,4) -
    		               2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		               4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		               lambdaphi*pow(vsigma,4),0.5)) +
    		          lambdaH*pow(vh,2)*pow(vsigma,2)*
    		           (12*pow(lambdaHphi,2)*pow(vh,2) -
    		             3*(-1 + lambdaphi)*lambdaphi*pow(vsigma,2) +
    		             4*lambdaphi*pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5)) +
    		          lambdaphi*pow(vsigma,4)*
    		           (pow(lambdaphi,2)*pow(vsigma,2) +
    		             pow(pow(lambdaH,2)*pow(vh,4) -
    		               2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		               4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		               lambdaphi*pow(vsigma,4),0.5) +
    		             3*lambdaphi*(pow(vsigma,2) +
    		                pow(pow(lambdaH,2)*pow(vh,4) -
    		                  2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                  4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                  lambdaphi*pow(vsigma,4),0.5)))) +
    		       pow(s - 2*pow(m_h2,2),2)*
    		        (32*pow(me,4) + 6*pow(lambdaH,2)*pow(vh,4) +
    		          10*pow(me + ((1 - cc)*s)/2. + pow(m_h2,2),2) +
    		          6*lambdaH*pow(vh,2)*pow(pow(lambdaH,2)*pow(vh,4) -
    		             2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		             4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		             lambdaphi*pow(vsigma,4),0.5) +
    		          8*pow(me,2)*(lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		             pow(pow(lambdaH,2)*pow(vh,4) -
    		               2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		               4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		               lambdaphi*pow(vsigma,4),0.5)) +
    		          3*pow(vsigma,2)*(4*pow(lambdaHphi,2)*pow(vh,2) +
    		             lambdaphi*(1 + lambdaphi)*pow(vsigma,2) +
    		             2*lambdaphi*pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5)) -
    		          2*(me + ((1 - cc)*s)/2. + pow(m_h2,2))*
    		           (16*pow(me,2) + 5*(lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		                pow(pow(lambdaH,2)*pow(vh,4) -
    		                  2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                  4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                  lambdaphi*pow(vsigma,4),0.5)))) +
    		       (s - 2*pow(m_h2,2))*(12*pow(lambdaH,3)*pow(vh,6) -
    		          16*pow(me + ((1 - cc)*s)/2. + pow(m_h2,2),3) +
    		          4*pow(lambdaH,2)*pow(vh,4)*
    		           (4*pow(me,2) + 3*pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5)) +
    		          lambdaH*pow(vh,2)*(36*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		             9*lambdaphi*pow(vsigma,4) - 9*pow(lambdaphi,2)*pow(vsigma,4) +
    		             16*pow(me,2)*pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5) +
    		             12*lambdaphi*pow(vsigma,2)*
    		              pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5)) +
    		          8*pow(me + ((1 - cc)*s)/2. + pow(m_h2,2),2)*
    		           (4*pow(me,2) + 3*(lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		                pow(pow(lambdaH,2)*pow(vh,4) -
    		                  2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                  4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                  lambdaphi*pow(vsigma,4),0.5))) -
    		          2*(me + ((1 - cc)*s)/2. + pow(m_h2,2))*
    		           (14*pow(lambdaH,2)*pow(vh,4) + 7*pow(lambdaphi,2)*pow(vsigma,4) +
    		             2*lambdaH*pow(vh,2)*
    		              (8*pow(me,2) + 7*
    		                 pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5)) +
    		             lambdaphi*pow(vsigma,2)*
    		              (16*pow(me,2) + 7*pow(vsigma,2) +
    		                14*pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5)) +
    		             4*(7*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(me,2)*pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5))) +
    		          pow(vsigma,2)*(8*lambdaphi*pow(me,2)*
    		              ((1 + lambdaphi)*pow(vsigma,2) +
    		                2*pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5)) +
    		             4*pow(lambdaHphi,2)*pow(vh,2)*
    		              (8*pow(me,2) + 9*lambdaphi*pow(vsigma,2) +
    		                3*pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5)) +
    		             3*lambdaphi*pow(vsigma,2)*
    		              (pow(lambdaphi,2)*pow(vsigma,2) +
    		                pow(pow(lambdaH,2)*pow(vh,4) -
    		                  2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                  4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                  lambdaphi*pow(vsigma,4),0.5) +
    		                3*lambdaphi*(pow(vsigma,2) +
    		                   pow(pow(lambdaH,2)*pow(vh,4) -
    		                     2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                     4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                     lambdaphi*pow(vsigma,4),0.5))))))*pow(sin(a),4));
    	den=16*pow(mW,4)*pow(-0.5*((1 - cc)*s) - pow(me,2),2)*
    			   pow(-0.5*((1 + cc)*s) - pow(me,2),2);
    	amp=num/den;
    }
    else if (id == 11)
    {
    	num=pow(g,4)*(64*pow(lambdaH,4)*pow(mW,4)*pow(vh,8) -
    		     64*pow(lambdaH,5)*pow(mW,2)*pow(vh,10) + 32*pow(lambdaH,6)*pow(vh,12) +
    		     256*pow(lambdaH,2)*pow(lambdaHphi,2)*pow(mW,4)*pow(vh,6)*
    		      pow(vsigma,2) - 320*pow(lambdaH,3)*pow(lambdaHphi,2)*pow(mW,2)*
    		      pow(vh,8)*pow(vsigma,2) +
    		     192*pow(lambdaH,4)*pow(lambdaHphi,2)*pow(vh,10)*pow(vsigma,2) +
    		     64*lambdaphi*pow(lambdaH,2)*pow(mW,4)*pow(vh,4)*pow(vsigma,4) +
    		     256*lambdaH*lambdaphi*pow(lambdaHphi,2)*pow(mW,4)*pow(vh,4)*
    		      pow(vsigma,4) + 128*pow(lambdaHphi,4)*pow(mW,4)*pow(vh,4)*
    		      pow(vsigma,4) - 64*pow(lambdaH,2)*pow(lambdaphi,2)*pow(mW,4)*pow(vh,4)*
    		      pow(vsigma,4) - 80*lambdaphi*pow(lambdaH,3)*pow(mW,2)*pow(vh,6)*
    		      pow(vsigma,4) - 320*lambdaphi*pow(lambdaH,2)*pow(lambdaHphi,2)*
    		      pow(mW,2)*pow(vh,6)*pow(vsigma,4) -
    		     320*lambdaH*pow(lambdaHphi,4)*pow(mW,2)*pow(vh,6)*pow(vsigma,4) +
    		     80*pow(lambdaH,3)*pow(lambdaphi,2)*pow(mW,2)*pow(vh,6)*pow(vsigma,4) +
    		     48*lambdaphi*pow(lambdaH,4)*pow(vh,8)*pow(vsigma,4) +
    		     192*lambdaphi*pow(lambdaH,3)*pow(lambdaHphi,2)*pow(vh,8)*
    		      pow(vsigma,4) + 288*pow(lambdaH,2)*pow(lambdaHphi,4)*pow(vh,8)*
    		      pow(vsigma,4) - 48*pow(lambdaH,4)*pow(lambdaphi,2)*pow(vh,8)*
    		      pow(vsigma,4) + 64*lambdaphi*pow(lambdaHphi,2)*pow(mW,4)*pow(vh,2)*
    		      pow(vsigma,6) + 64*lambdaH*pow(lambdaphi,2)*pow(mW,4)*pow(vh,2)*
    		      pow(vsigma,6) + 192*pow(lambdaHphi,2)*pow(lambdaphi,2)*pow(mW,4)*
    		      pow(vh,2)*pow(vsigma,6) -
    		     64*lambdaH*pow(lambdaphi,3)*pow(mW,4)*pow(vh,2)*pow(vsigma,6) -
    		     160*lambdaH*lambdaphi*pow(lambdaHphi,2)*pow(mW,2)*pow(vh,4)*
    		      pow(vsigma,6) - 320*lambdaphi*pow(lambdaHphi,4)*pow(mW,2)*pow(vh,4)*
    		      pow(vsigma,6) - 80*pow(lambdaH,2)*pow(lambdaphi,2)*pow(mW,2)*pow(vh,4)*
    		      pow(vsigma,6) - 160*lambdaH*pow(lambdaHphi,2)*pow(lambdaphi,2)*
    		      pow(mW,2)*pow(vh,4)*pow(vsigma,6) +
    		     80*pow(lambdaH,2)*pow(lambdaphi,3)*pow(mW,2)*pow(vh,4)*pow(vsigma,6) +
    		     144*lambdaphi*pow(lambdaH,2)*pow(lambdaHphi,2)*pow(vh,6)*
    		      pow(vsigma,6) + 384*lambdaH*lambdaphi*pow(lambdaHphi,4)*pow(vh,6)*
    		      pow(vsigma,6) + 64*pow(lambdaHphi,6)*pow(vh,6)*pow(vsigma,6) +
    		     48*pow(lambdaH,3)*pow(lambdaphi,2)*pow(vh,6)*pow(vsigma,6) +
    		     48*pow(lambdaH,2)*pow(lambdaHphi,2)*pow(lambdaphi,2)*pow(vh,6)*
    		      pow(vsigma,6) - 48*pow(lambdaH,3)*pow(lambdaphi,3)*pow(vh,6)*
    		      pow(vsigma,6) + 8*pow(lambdaphi,2)*pow(mW,4)*pow(vsigma,8) +
    		     48*pow(lambdaphi,3)*pow(mW,4)*pow(vsigma,8) +
    		     8*pow(lambdaphi,4)*pow(mW,4)*pow(vsigma,8) -
    		     20*lambdaH*pow(lambdaphi,2)*pow(mW,2)*pow(vh,2)*pow(vsigma,8) -
    		     160*pow(lambdaHphi,2)*pow(lambdaphi,2)*pow(mW,2)*pow(vh,2)*
    		      pow(vsigma,8) - 40*lambdaH*pow(lambdaphi,3)*pow(mW,2)*pow(vh,2)*
    		      pow(vsigma,8) - 160*pow(lambdaHphi,2)*pow(lambdaphi,3)*pow(mW,2)*
    		      pow(vh,2)*pow(vsigma,8) +
    		     60*lambdaH*pow(lambdaphi,4)*pow(mW,2)*pow(vh,2)*pow(vsigma,8) +
    		     48*lambdaphi*pow(lambdaHphi,4)*pow(vh,4)*pow(vsigma,8) +
    		     18*pow(lambdaH,2)*pow(lambdaphi,2)*pow(vh,4)*pow(vsigma,8) +
    		     192*lambdaH*pow(lambdaHphi,2)*pow(lambdaphi,2)*pow(vh,4)*
    		      pow(vsigma,8) + 240*pow(lambdaHphi,4)*pow(lambdaphi,2)*pow(vh,4)*
    		      pow(vsigma,8) + 12*pow(lambdaH,2)*pow(lambdaphi,3)*pow(vh,4)*
    		      pow(vsigma,8) - 30*pow(lambdaH,2)*pow(lambdaphi,4)*pow(vh,4)*
    		      pow(vsigma,8) - 20*pow(lambdaphi,3)*pow(mW,2)*pow(vsigma,10) -
    		     40*pow(lambdaphi,4)*pow(mW,2)*pow(vsigma,10) -
    		     4*pow(lambdaphi,5)*pow(mW,2)*pow(vsigma,10) +
    		     12*pow(lambdaHphi,2)*pow(lambdaphi,2)*pow(vh,2)*pow(vsigma,10) +
    		     24*lambdaH*pow(lambdaphi,3)*pow(vh,2)*pow(vsigma,10) +
    		     120*pow(lambdaHphi,2)*pow(lambdaphi,3)*pow(vh,2)*pow(vsigma,10) +
    		     60*pow(lambdaHphi,2)*pow(lambdaphi,4)*pow(vh,2)*pow(vsigma,10) -
    		     24*lambdaH*pow(lambdaphi,5)*pow(vh,2)*pow(vsigma,10) +
    		     pow(lambdaphi,3)*pow(vsigma,12) + 15*pow(lambdaphi,4)*pow(vsigma,12) +
    		     15*pow(lambdaphi,5)*pow(vsigma,12) + pow(lambdaphi,6)*pow(vsigma,12) +
    		     64*pow(lambdaH,3)*pow(mW,4)*pow(vh,6)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) - 64*pow(lambdaH,4)*pow(mW,2)*pow(vh,8)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 32*pow(lambdaH,5)*pow(vh,10)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 64*lambdaphi*pow(lambdaH,2)*pow(mW,4)*pow(vh,4)*
    		      pow(vsigma,2)*pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 128*lambdaH*pow(lambdaHphi,2)*pow(mW,4)*pow(vh,4)*
    		      pow(vsigma,2)*pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) - 64*lambdaphi*pow(lambdaH,3)*pow(mW,2)*pow(vh,6)*
    		      pow(vsigma,2)*pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) - 192*pow(lambdaH,2)*pow(lambdaHphi,2)*pow(mW,2)*pow(vh,6)*
    		      pow(vsigma,2)*pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 32*lambdaphi*pow(lambdaH,4)*pow(vh,8)*pow(vsigma,2)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 128*pow(lambdaH,3)*pow(lambdaHphi,2)*pow(vh,8)*pow(vsigma,2)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 32*lambdaH*lambdaphi*pow(mW,4)*pow(vh,2)*pow(vsigma,4)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 128*lambdaphi*pow(lambdaHphi,2)*pow(mW,4)*pow(vh,2)*
    		      pow(vsigma,4)*pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 32*lambdaH*pow(lambdaphi,2)*pow(mW,4)*pow(vh,2)*
    		      pow(vsigma,4)*pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) - 48*lambdaphi*pow(lambdaH,2)*pow(mW,2)*pow(vh,4)*
    		      pow(vsigma,4)*pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) - 256*lambdaH*lambdaphi*pow(lambdaHphi,2)*pow(mW,2)*pow(vh,4)*
    		      pow(vsigma,4)*pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) - 64*pow(lambdaHphi,4)*pow(mW,2)*pow(vh,4)*pow(vsigma,4)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) - 16*pow(lambdaH,2)*pow(lambdaphi,2)*pow(mW,2)*pow(vh,4)*
    		      pow(vsigma,4)*pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 32*lambdaphi*pow(lambdaH,3)*pow(vh,6)*pow(vsigma,4)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 192*lambdaphi*pow(lambdaH,2)*pow(lambdaHphi,2)*pow(vh,6)*
    		      pow(vsigma,4)*pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 96*lambdaH*pow(lambdaHphi,4)*pow(vh,6)*pow(vsigma,4)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 32*pow(lambdaphi,2)*pow(mW,4)*pow(vsigma,6)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 32*pow(lambdaphi,3)*pow(mW,4)*pow(vsigma,6)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) - 32*lambdaphi*pow(lambdaHphi,2)*pow(mW,2)*pow(vh,2)*
    		      pow(vsigma,6)*pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) - 64*lambdaH*pow(lambdaphi,2)*pow(mW,2)*pow(vh,2)*
    		      pow(vsigma,6)*pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) - 160*pow(lambdaHphi,2)*pow(lambdaphi,2)*pow(mW,2)*pow(vh,2)*
    		      pow(vsigma,6)*pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 48*lambdaH*lambdaphi*pow(lambdaHphi,2)*pow(vh,4)*
    		      pow(vsigma,6)*pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 96*lambdaphi*pow(lambdaHphi,4)*pow(vh,4)*pow(vsigma,6)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 48*pow(lambdaH,2)*pow(lambdaphi,2)*pow(vh,4)*pow(vsigma,6)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 144*lambdaH*pow(lambdaHphi,2)*pow(lambdaphi,2)*pow(vh,4)*
    		      pow(vsigma,6)*pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) - 16*pow(lambdaH,2)*pow(lambdaphi,3)*pow(vh,4)*pow(vsigma,6)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) - 4*pow(lambdaphi,2)*pow(mW,2)*pow(vsigma,8)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) - 40*pow(lambdaphi,3)*pow(mW,2)*pow(vsigma,8)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) - 20*pow(lambdaphi,4)*pow(mW,2)*pow(vsigma,8)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 6*lambdaH*pow(lambdaphi,2)*pow(vh,2)*pow(vsigma,8)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 48*pow(lambdaHphi,2)*pow(lambdaphi,2)*pow(vh,2)*
    		      pow(vsigma,8)*pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 36*lambdaH*pow(lambdaphi,3)*pow(vh,2)*pow(vsigma,8)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 80*pow(lambdaHphi,2)*pow(lambdaphi,3)*pow(vh,2)*
    		      pow(vsigma,8)*pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) - 10*lambdaH*pow(lambdaphi,4)*pow(vh,2)*pow(vsigma,8)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 6*pow(lambdaphi,3)*pow(vsigma,10)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 20*pow(lambdaphi,4)*pow(vsigma,10)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + 6*pow(lambdaphi,5)*pow(vsigma,10)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) + lambdaphi*pow(vsigma,4)
    		        ,0.5) + pow(s - 2*pow(m_h2,2),4)*
    		      (16*pow(mW,4) + 2*pow(lambdaH,2)*pow(vh,4) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		        lambdaphi*pow(vsigma,4) + pow(lambdaphi,2)*pow(vsigma,4) +
    		        4*pow(mW + ((1 - cc)*s)/2. + pow(m_h2,2),2) +
    		        2*lambdaH*pow(vh,2)*pow(pow(lambdaH,2)*pow(vh,4) -
    		           2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		           4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		           lambdaphi*pow(vsigma,4),0.5) +
    		        2*lambdaphi*pow(vsigma,2)*
    		         pow(pow(lambdaH,2)*pow(vh,4) -
    		           2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		           4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		           lambdaphi*pow(vsigma,4),0.5) -
    		        8*pow(mW,2)*(lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		           pow(pow(lambdaH,2)*pow(vh,4) -
    		             2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		             4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		             lambdaphi*pow(vsigma,4),0.5)) -
    		        4*(mW + ((1 - cc)*s)/2. + pow(m_h2,2))*
    		         (-4*pow(mW,2) + lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		           pow(pow(lambdaH,2)*pow(vh,4) -
    		             2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		             4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		             lambdaphi*pow(vsigma,4),0.5))) +
    		     16*pow(mW + ((1 - cc)*s)/2. + pow(m_h2,2),4)*
    		      (2*pow(lambdaH,2)*pow(vh,4) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		        pow(lambdaphi,2)*pow(vsigma,4) -
    		        2*pow(mW,2)*pow(pow(lambdaH,2)*pow(vh,4) -
    		           2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		           4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		           lambdaphi*pow(vsigma,4),0.5) +
    		        2*lambdaH*pow(vh,2)*(-pow(mW,2) +
    		           pow(pow(lambdaH,2)*pow(vh,4) -
    		             2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		             4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		             lambdaphi*pow(vsigma,4),0.5)) +
    		        lambdaphi*pow(vsigma,2)*
    		         (-2*pow(mW,2) + pow(vsigma,2) +
    		           2*pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5))) -
    		     32*pow(mW + ((1 - cc)*s)/2. + pow(m_h2,2),3)*
    		      (4*pow(lambdaH,3)*pow(vh,6) +
    		        4*pow(lambdaH,2)*pow(vh,4)*
    		         (-pow(mW,2) + pow(pow(lambdaH,2)*pow(vh,4) -
    		             2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		             4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		             lambdaphi*pow(vsigma,4),0.5)) +
    		        lambdaH*pow(vh,2)*(12*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		           3*lambdaphi*pow(vsigma,4) - 3*pow(lambdaphi,2)*pow(vsigma,4) -
    		           4*pow(mW,2)*pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5) +
    		           4*lambdaphi*pow(vsigma,2)*
    		            pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5)) +
    		        pow(vsigma,2)*(4*pow(lambdaHphi,2)*pow(vh,2)*
    		            (-2*pow(mW,2) + 3*lambdaphi*pow(vsigma,2) +
    		              pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5)) -
    		           2*lambdaphi*pow(mW,2)*
    		            ((1 + lambdaphi)*pow(vsigma,2) +
    		              2*pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5)) +
    		           lambdaphi*pow(vsigma,2)*
    		            (pow(lambdaphi,2)*pow(vsigma,2) +
    		              pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5) +
    		              3*lambdaphi*(pow(vsigma,2) +
    		                 pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5))))) +
    		     4*pow(s - 2*pow(m_h2,2),3)*
    		      (4*pow(lambdaH,3)*pow(vh,6) + 8*lambdaphi*pow(mW,4)*pow(vsigma,2) -
    		        24*pow(lambdaHphi,2)*pow(mW,2)*pow(vh,2)*pow(vsigma,2) -
    		        6*lambdaphi*pow(mW,2)*pow(vsigma,4) -
    		        6*pow(lambdaphi,2)*pow(mW,2)*pow(vsigma,4) +
    		        12*lambdaphi*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,4) +
    		        3*pow(lambdaphi,2)*pow(vsigma,6) + pow(lambdaphi,3)*pow(vsigma,6) -
    		        2*pow(mW + ((1 - cc)*s)/2. + pow(m_h2,2),3) +
    		        8*pow(mW,4)*pow(pow(lambdaH,2)*pow(vh,4) -
    		           2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		           4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		           lambdaphi*pow(vsigma,4),0.5) -
    		        12*lambdaphi*pow(mW,2)*pow(vsigma,2)*
    		         pow(pow(lambdaH,2)*pow(vh,4) -
    		           2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		           4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		           lambdaphi*pow(vsigma,4),0.5) +
    		        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2)*
    		         pow(pow(lambdaH,2)*pow(vh,4) -
    		           2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		           4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		           lambdaphi*pow(vsigma,4),0.5) +
    		        lambdaphi*pow(vsigma,4)*
    		         pow(pow(lambdaH,2)*pow(vh,4) -
    		           2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		           4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		           lambdaphi*pow(vsigma,4),0.5) +
    		        3*pow(lambdaphi,2)*pow(vsigma,4)*
    		         pow(pow(lambdaH,2)*pow(vh,4) -
    		           2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		           4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		           lambdaphi*pow(vsigma,4),0.5) +
    		        4*pow(lambdaH,2)*pow(vh,4)*
    		         (-3*pow(mW,2) + pow(pow(lambdaH,2)*pow(vh,4) -
    		             2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		             4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		             lambdaphi*pow(vsigma,4),0.5)) +
    		        pow(mW + ((1 - cc)*s)/2. + pow(m_h2,2),2)*
    		         (-12*pow(mW,2) + 7*(lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		              pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5))) +
    		        lambdaH*pow(vh,2)*(8*pow(mW,4) -
    		           12*pow(mW,2)*pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5) +
    		           pow(vsigma,2)*(12*pow(lambdaHphi,2)*pow(vh,2) -
    		              3*(-1 + lambdaphi)*lambdaphi*pow(vsigma,2) +
    		              4*lambdaphi*pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5))) -
    		        (mW + ((1 - cc)*s)/2. + pow(m_h2,2))*
    		         (16*pow(mW,4) - 20*pow(mW,2)*
    		            (lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		              pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5)) +
    		           5*(2*pow(lambdaH,2)*pow(vh,4) +
    		              2*lambdaH*pow(vh,2)*
    		               pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5) +
    		              pow(vsigma,2)*(4*pow(lambdaHphi,2)*pow(vh,2) +
    		                 lambdaphi*(1 + lambdaphi)*pow(vsigma,2) +
    		                 2*lambdaphi*pow(pow(lambdaH,2)*pow(vh,4) -
    		                    2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                    4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                    lambdaphi*pow(vsigma,4),0.5))))) +
    		     2*pow(s - 2*pow(m_h2,2),2)*
    		      (96*pow(mW,8) + 2*pow(mW + ((1 - cc)*s)/2. + pow(m_h2,2),4) -
    		        32*pow(mW,6)*(lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		           pow(pow(lambdaH,2)*pow(vh,4) -
    		             2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		             4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		             lambdaphi*pow(vsigma,4),0.5)) -
    		        4*pow(mW + ((1 - cc)*s)/2. + pow(m_h2,2),3)*
    		         (-8*pow(mW,2) + 5*(lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		              pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5))) +
    		        28*pow(mW,4)*(2*pow(lambdaH,2)*pow(vh,4) +
    		           2*lambdaH*pow(vh,2)*
    		            pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5) +
    		           pow(vsigma,2)*(4*pow(lambdaHphi,2)*pow(vh,2) +
    		              lambdaphi*(1 + lambdaphi)*pow(vsigma,2) +
    		              2*lambdaphi*pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5))) -
    		        14*pow(mW,2)*(4*pow(lambdaH,3)*pow(vh,6) +
    		           4*pow(lambdaH,2)*pow(vh,4)*
    		            pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5) +
    		           4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2)*
    		            (3*lambdaphi*pow(vsigma,2) +
    		              pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5)) +
    		           lambdaH*pow(vh,2)*pow(vsigma,2)*
    		            (12*pow(lambdaHphi,2)*pow(vh,2) -
    		              3*(-1 + lambdaphi)*lambdaphi*pow(vsigma,2) +
    		              4*lambdaphi*pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5)) +
    		           lambdaphi*pow(vsigma,4)*
    		            (pow(lambdaphi,2)*pow(vsigma,2) +
    		              pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5) +
    		              3*lambdaphi*(pow(vsigma,2) +
    		                 pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5)))) +
    		        2*pow(mW + ((1 - cc)*s)/2. + pow(m_h2,2),2)*
    		         (16*pow(mW,4) - 40*pow(mW,2)*
    		            (lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		              pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5)) +
    		           17*(2*pow(lambdaH,2)*pow(vh,4) +
    		              2*lambdaH*pow(vh,2)*
    		               pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5) +
    		              pow(vsigma,2)*(4*pow(lambdaHphi,2)*pow(vh,2) +
    		                 lambdaphi*(1 + lambdaphi)*pow(vsigma,2) +
    		                 2*lambdaphi*pow(pow(lambdaH,2)*pow(vh,4) -
    		                    2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                    4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                    lambdaphi*pow(vsigma,4),0.5)))) -
    		        2*(mW + ((1 - cc)*s)/2. + pow(m_h2,2))*
    		         (36*pow(lambdaH,3)*pow(vh,6) + 9*pow(lambdaphi,3)*pow(vsigma,6) +
    		           12*pow(lambdaH,2)*pow(vh,4)*
    		            (-5*pow(mW,2) + 3*pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5)) +
    		           3*pow(lambdaphi,2)*pow(vsigma,4)*
    		            (-10*pow(mW,2) + 9*
    		               (pow(vsigma,2) +
    		                 pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5))) +
    		           lambdaphi*pow(vsigma,2)*
    		            (32*pow(mW,4) + 9*pow(vsigma,2)*
    		               (12*pow(lambdaHphi,2)*pow(vh,2) +
    		                 pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5)) -
    		              30*pow(mW,2)*(pow(vsigma,2) +
    		                 2*pow(pow(lambdaH,2)*pow(vh,4) -
    		                    2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                    4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                    lambdaphi*pow(vsigma,4),0.5))) +
    		           4*(8*pow(mW,4)*pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5) +
    		              3*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2)*
    		               (-10*pow(mW,2) +
    		                 3*pow(pow(lambdaH,2)*pow(vh,4) -
    		                    2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                    4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                    lambdaphi*pow(vsigma,4),0.5))) +
    		           lambdaH*pow(vh,2)*(32*pow(mW,4) -
    		              60*pow(mW,2)*pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5) +
    		              9*pow(vsigma,2)*(12*pow(lambdaHphi,2)*pow(vh,2) -
    		                 3*(-1 + lambdaphi)*lambdaphi*pow(vsigma,2) +
    		                 4*lambdaphi*pow(pow(lambdaH,2)*pow(vh,4) -
    		                    2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                    4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                    lambdaphi*pow(vsigma,4),0.5)))) +
    		        3*(8*pow(lambdaH,4)*pow(vh,8) +
    		           8*pow(lambdaH,3)*pow(vh,6)*
    		            pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5) +
    		           8*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,2)*
    		            (4*pow(lambdaHphi,2)*pow(vh,2) +
    		              lambdaphi*(pow(vsigma,2) - lambdaphi*pow(vsigma,2) +
    		                 pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5))) +
    		           4*lambdaH*(4*pow(lambdaHphi,2)*pow(vh,4)*pow(vsigma,2)*
    		               (2*lambdaphi*pow(vsigma,2) +
    		                 pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5)) +
    		              lambdaphi*pow(vh,2)*pow(vsigma,4)*
    		               (2*lambdaphi*pow(vsigma,2) -
    		                 2*pow(lambdaphi,2)*pow(vsigma,2) +
    		                 pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5) +
    		                 lambdaphi*pow(pow(lambdaH,2)*pow(vh,4) -
    		                    2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                    4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                    lambdaphi*pow(vsigma,4),0.5))) +
    		           pow(vsigma,4)*(16*pow(lambdaHphi,4)*pow(vh,4) +
    		              8*lambdaphi*pow(lambdaHphi,2)*pow(vh,2)*
    		               ((1 + 3*lambdaphi)*pow(vsigma,2) +
    		                 2*pow(pow(lambdaH,2)*pow(vh,4) -
    		                    2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                    4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                    lambdaphi*pow(vsigma,4),0.5)) +
    		              pow(lambdaphi,2)*pow(vsigma,2)*
    		               ((1 + 6*lambdaphi + pow(lambdaphi,2))*pow(vsigma,2) +
    		                 4*(1 + lambdaphi)*
    		                  pow(pow(lambdaH,2)*pow(vh,4) -
    		                    2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                    4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                    lambdaphi*pow(vsigma,4),0.5))))) +
    		     8*pow(mW + ((1 - cc)*s)/2. + pow(m_h2,2),2)*
    		      (24*pow(lambdaH,4)*pow(vh,8) +
    		        4*pow(lambdaH,3)*pow(vh,6)*
    		         (-7*pow(mW,2) + 6*pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5)) +
    		        lambdaH*pow(vh,2)*(-24*pow(lambdaphi,3)*pow(vsigma,6) +
    		           8*pow(mW,4)*pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5) +
    		           12*lambdaphi*pow(vsigma,4)*
    		            pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5) +
    		           3*pow(lambdaphi,2)*pow(vsigma,4)*
    		            (7*pow(mW,2) + 8*pow(vsigma,2) +
    		              4*pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5)) +
    		           12*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2)*
    		            (-7*pow(mW,2) + 8*lambdaphi*pow(vsigma,2) +
    		              4*pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5)) -
    		           7*lambdaphi*pow(mW,2)*
    		            (3*pow(vsigma,4) +
    		              4*pow(vsigma,2)*pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5))) +
    		        4*pow(lambdaH,2)*pow(vh,4)*
    		         (2*pow(mW,4) - 7*pow(mW,2)*
    		            pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5) +
    		           6*pow(vsigma,2)*(4*pow(lambdaHphi,2)*pow(vh,2) +
    		              lambdaphi*(pow(vsigma,2) - lambdaphi*pow(vsigma,2) +
    		                 pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5)))) +
    		        pow(vsigma,2)*(48*pow(lambdaHphi,4)*pow(vh,4)*pow(vsigma,2) +
    		           4*pow(lambdaHphi,2)*pow(vh,2)*
    		            (4*pow(mW,4) - 7*pow(mW,2)*
    		               (3*lambdaphi*pow(vsigma,2) +
    		                 pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5)) +
    		              6*lambdaphi*pow(vsigma,2)*
    		               (pow(vsigma,2) + 3*lambdaphi*pow(vsigma,2) +
    		                 2*pow(pow(lambdaH,2)*pow(vh,4) -
    		                    2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                    4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                    lambdaphi*pow(vsigma,4),0.5))) +
    		           lambdaphi*(4*pow(mW,4)*
    		               ((1 + lambdaphi)*pow(vsigma,2) +
    		                 2*pow(pow(lambdaH,2)*pow(vh,4) -
    		                    2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                    4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                    lambdaphi*pow(vsigma,4),0.5)) +
    		              3*lambdaphi*pow(vsigma,4)*
    		               ((1 + 6*lambdaphi + pow(lambdaphi,2))*pow(vsigma,2) +
    		                 4*(1 + lambdaphi)*
    		                  pow(pow(lambdaH,2)*pow(vh,4) -
    		                    2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                    4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                    lambdaphi*pow(vsigma,4),0.5)) -
    		              7*pow(mW,2)*pow(vsigma,2)*
    		               (pow(lambdaphi,2)*pow(vsigma,2) +
    		                 pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5) +
    		                 3*lambdaphi*(pow(vsigma,2) +
    		                    pow(pow(lambdaH,2)*pow(vh,4) -
    		                      2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                      4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                      lambdaphi*pow(vsigma,4),0.5)))))) -
    		     8*(mW + ((1 - cc)*s)/2. + pow(m_h2,2))*
    		      (16*pow(lambdaH,5)*pow(vh,10) +
    		        16*pow(lambdaHphi,4)*pow(vh,4)*pow(vsigma,4)*
    		         (-3*pow(mW,2) + 5*lambdaphi*pow(vsigma,2) +
    		           pow(pow(lambdaH,2)*pow(vh,4) -
    		             2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		             4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		             lambdaphi*pow(vsigma,4),0.5)) +
    		        8*pow(lambdaH,4)*pow(vh,8)*
    		         (-3*pow(mW,2) + 2*pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5)) +
    		        4*pow(lambdaH,2)*pow(vh,4)*
    		         (-5*pow(lambdaphi,3)*pow(vsigma,6) +
    		           4*pow(mW,4)*pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5) +
    		           3*lambdaphi*pow(vsigma,4)*
    		            pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5) -
    		           6*lambdaphi*pow(mW,2)*pow(vsigma,2)*
    		            (pow(vsigma,2) + pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5)) +
    		           pow(lambdaphi,2)*pow(vsigma,4)*
    		            (6*pow(mW,2) + 5*pow(vsigma,2) +
    		              pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5)) +
    		           4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2)*
    		            (-6*pow(mW,2) + 5*lambdaphi*pow(vsigma,2) +
    		              3*pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5))) +
    		        4*pow(lambdaH,3)*pow(vh,6)*
    		         (4*pow(mW,4) - 6*pow(mW,2)*
    		            pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5) +
    		           pow(vsigma,2)*(20*pow(lambdaHphi,2)*pow(vh,2) -
    		              5*(-1 + lambdaphi)*lambdaphi*pow(vsigma,2) +
    		              4*lambdaphi*pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5))) +
    		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2)*
    		         (5*pow(lambdaphi,3)*pow(vsigma,6) +
    		           2*pow(mW,4)*pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5) +
    		           pow(lambdaphi,2)*pow(vsigma,4)*
    		            (-9*pow(mW,2) + 5*(pow(vsigma,2) +
    		                 pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5))) +
    		           lambdaphi*pow(vsigma,2)*
    		            (6*pow(mW,4) + pow(vsigma,2)*
    		               pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5) -
    		              3*pow(mW,2)*(pow(vsigma,2) +
    		                 2*pow(pow(lambdaH,2)*pow(vh,4) -
    		                    2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                    4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                    lambdaphi*pow(vsigma,4),0.5)))) +
    		        lambdaphi*pow(vsigma,4)*
    		         (pow(lambdaphi,4)*pow(vsigma,6) +
    		           4*pow(mW,4)*pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5) +
    		           pow(lambdaphi,3)*pow(vsigma,4)*
    		            (-3*pow(mW,2) + 5*(2*pow(vsigma,2) +
    		                 pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5))) +
    		           pow(lambdaphi,2)*pow(vsigma,2)*
    		            (4*pow(mW,4) - 6*pow(mW,2)*
    		               (3*pow(vsigma,2) +
    		                 2*pow(pow(lambdaH,2)*pow(vh,4) -
    		                    2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                    4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                    lambdaphi*pow(vsigma,4),0.5)) +
    		              5*(pow(vsigma,4) +
    		                 2*pow(vsigma,2)*
    		                  pow(pow(lambdaH,2)*pow(vh,4) -
    		                    2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                    4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                    lambdaphi*pow(vsigma,4),0.5))) +
    		           lambdaphi*(pow(vsigma,4)*
    		               pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5) +
    		              12*pow(mW,4)*(pow(vsigma,2) +
    		                 pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5)) -
    		              3*pow(mW,2)*(pow(vsigma,4) +
    		                 4*pow(vsigma,2)*
    		                  pow(pow(lambdaH,2)*pow(vh,4) -
    		                    2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                    4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                    lambdaphi*pow(vsigma,4),0.5)))) +
    		        lambdaH*pow(vh,2)*pow(vsigma,2)*
    		         (80*pow(lambdaHphi,4)*pow(vh,4)*pow(vsigma,2) +
    		           8*pow(lambdaHphi,2)*pow(vh,2)*
    		            (6*pow(mW,4) + 5*lambdaphi*(1 + lambdaphi)*pow(vsigma,4) +
    		              8*lambdaphi*pow(vsigma,2)*
    		               pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5) -
    		              6*pow(mW,2)*(2*lambdaphi*pow(vsigma,2) +
    		                 pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5))) +
    		           lambdaphi*(4*pow(mW,4)*
    		               (-3*(-1 + lambdaphi)*pow(vsigma,2) +
    		                 4*pow(pow(lambdaH,2)*pow(vh,4) -
    		                    2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                    4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                    lambdaphi*pow(vsigma,4),0.5)) +
    		              lambdaphi*pow(vsigma,4)*
    		               ((5 + 10*lambdaphi - 15*pow(lambdaphi,2))*pow(vsigma,2) +
    		                 16*pow(pow(lambdaH,2)*pow(vh,4) -
    		                    2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                    4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                    lambdaphi*pow(vsigma,4),0.5)) -
    		              12*pow(mW,2)*pow(vsigma,2)*
    		               (-2*pow(lambdaphi,2)*pow(vsigma,2) +
    		                 pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5) +
    		                 lambdaphi*(2*pow(vsigma,2) +
    		                    pow(pow(lambdaH,2)*pow(vh,4) -
    		                      2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                      4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                      lambdaphi*pow(vsigma,4),0.5)))))) +
    		     4*(s - 2*pow(m_h2,2))*(16*pow(lambdaH,5)*pow(vh,10) +
    		        16*pow(lambdaH,4)*pow(vh,8)*
    		         (-2*pow(mW,2) + pow(pow(lambdaH,2)*pow(vh,4) -
    		             2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		             4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		             lambdaphi*pow(vsigma,4),0.5)) +
    		        4*pow(mW + ((1 - cc)*s)/2. + pow(m_h2,2),4)*
    		         (-2*pow(mW,2) + lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		           pow(pow(lambdaH,2)*pow(vh,4) -
    		             2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		             4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		             lambdaphi*pow(vsigma,4),0.5)) +
    		        16*pow(lambdaHphi,4)*pow(vh,4)*pow(vsigma,4)*
    		         (-4*pow(mW,2) + 5*lambdaphi*pow(vsigma,2) +
    		           pow(pow(lambdaH,2)*pow(vh,4) -
    		             2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		             4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		             lambdaphi*pow(vsigma,4),0.5)) -
    		        16*pow(mW + ((1 - cc)*s)/2. + pow(m_h2,2),3)*
    		         (2*pow(lambdaH,2)*pow(vh,4) +
    		           4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		           pow(lambdaphi,2)*pow(vsigma,4) -
    		           2*pow(mW,2)*pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5) +
    		           2*lambdaH*pow(vh,2)*
    		            (-pow(mW,2) + pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5)) +
    		           lambdaphi*pow(vsigma,2)*
    		            (-2*pow(mW,2) + pow(vsigma,2) +
    		              2*pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5))) +
    		        4*pow(lambdaH,2)*pow(vh,4)*
    		         (-5*pow(lambdaphi,3)*pow(vsigma,6) +
    		           8*pow(mW,4)*pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5) +
    		           3*lambdaphi*pow(vsigma,4)*
    		            pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5) -
    		           8*lambdaphi*pow(mW,2)*pow(vsigma,2)*
    		            (pow(vsigma,2) + pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5)) +
    		           pow(lambdaphi,2)*pow(vsigma,4)*
    		            (8*pow(mW,2) + 5*pow(vsigma,2) +
    		              pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5)) +
    		           4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2)*
    		            (-8*pow(mW,2) + 5*lambdaphi*pow(vsigma,2) +
    		              3*pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5))) +
    		        4*pow(lambdaH,3)*pow(vh,6)*
    		         (8*pow(mW,4) - 8*pow(mW,2)*
    		            pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5) +
    		           pow(vsigma,2)*(20*pow(lambdaHphi,2)*pow(vh,2) -
    		              5*(-1 + lambdaphi)*lambdaphi*pow(vsigma,2) +
    		              4*lambdaphi*pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5))) +
    		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2)*
    		         (5*pow(lambdaphi,3)*pow(vsigma,6) +
    		           4*pow(mW,4)*pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5) +
    		           pow(lambdaphi,2)*pow(vsigma,4)*
    		            (-12*pow(mW,2) + 5*
    		               (pow(vsigma,2) +
    		                 pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5))) +
    		           lambdaphi*pow(vsigma,2)*
    		            (12*pow(mW,4) + pow(vsigma,2)*
    		               pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5) -
    		              4*pow(mW,2)*(pow(vsigma,2) +
    		                 2*pow(pow(lambdaH,2)*pow(vh,4) -
    		                    2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                    4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                    lambdaphi*pow(vsigma,4),0.5)))) +
    		        pow(mW + ((1 - cc)*s)/2. + pow(m_h2,2),2)*
    		         (68*pow(lambdaH,3)*pow(vh,6) + 17*pow(lambdaphi,3)*pow(vsigma,6) +
    		           4*pow(lambdaH,2)*pow(vh,4)*
    		            (-19*pow(mW,2) + 17*
    		               pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5)) +
    		           pow(lambdaphi,2)*pow(vsigma,4)*
    		            (-38*pow(mW,2) + 51*
    		               (pow(vsigma,2) +
    		                 pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5))) +
    		           lambdaphi*pow(vsigma,2)*
    		            (16*pow(mW,4) + 17*pow(vsigma,2)*
    		               (12*pow(lambdaHphi,2)*pow(vh,2) +
    		                 pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5)) -
    		              38*pow(mW,2)*(pow(vsigma,2) +
    		                 2*pow(pow(lambdaH,2)*pow(vh,4) -
    		                    2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                    4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                    lambdaphi*pow(vsigma,4),0.5))) +
    		           4*(4*pow(mW,4)*pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5) +
    		              pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2)*
    		               (-38*pow(mW,2) +
    		                 17*pow(pow(lambdaH,2)*pow(vh,4) -
    		                    2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                    4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                    lambdaphi*pow(vsigma,4),0.5))) +
    		           lambdaH*pow(vh,2)*(16*pow(mW,4) -
    		              76*pow(mW,2)*pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5) +
    		              17*pow(vsigma,2)*
    		               (12*pow(lambdaHphi,2)*pow(vh,2) -
    		                 3*(-1 + lambdaphi)*lambdaphi*pow(vsigma,2) +
    		                 4*lambdaphi*pow(pow(lambdaH,2)*pow(vh,4) -
    		                    2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                    4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                    lambdaphi*pow(vsigma,4),0.5)))) +
    		        lambdaphi*pow(vsigma,4)*
    		         (pow(lambdaphi,4)*pow(vsigma,6) +
    		           8*pow(mW,4)*pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5) +
    		           pow(lambdaphi,3)*pow(vsigma,4)*
    		            (-4*pow(mW,2) + 5*(2*pow(vsigma,2) +
    		                 pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5))) +
    		           pow(lambdaphi,2)*pow(vsigma,2)*
    		            (8*pow(mW,4) - 8*pow(mW,2)*
    		               (3*pow(vsigma,2) +
    		                 2*pow(pow(lambdaH,2)*pow(vh,4) -
    		                    2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                    4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                    lambdaphi*pow(vsigma,4),0.5)) +
    		              5*(pow(vsigma,4) +
    		                 2*pow(vsigma,2)*
    		                  pow(pow(lambdaH,2)*pow(vh,4) -
    		                    2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                    4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                    lambdaphi*pow(vsigma,4),0.5))) +
    		           lambdaphi*(pow(vsigma,4)*
    		               pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5) +
    		              24*pow(mW,4)*(pow(vsigma,2) +
    		                 pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5)) -
    		              4*pow(mW,2)*(pow(vsigma,4) +
    		                 4*pow(vsigma,2)*
    		                  pow(pow(lambdaH,2)*pow(vh,4) -
    		                    2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                    4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                    lambdaphi*pow(vsigma,4),0.5)))) +
    		        lambdaH*pow(vh,2)*pow(vsigma,2)*
    		         (80*pow(lambdaHphi,4)*pow(vh,4)*pow(vsigma,2) +
    		           8*pow(lambdaHphi,2)*pow(vh,2)*
    		            (12*pow(mW,4) + 5*lambdaphi*(1 + lambdaphi)*pow(vsigma,4) +
    		              8*lambdaphi*pow(vsigma,2)*
    		               pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5) -
    		              8*pow(mW,2)*(2*lambdaphi*pow(vsigma,2) +
    		                 pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5))) +
    		           lambdaphi*(8*pow(mW,4)*
    		               (-3*(-1 + lambdaphi)*pow(vsigma,2) +
    		                 4*pow(pow(lambdaH,2)*pow(vh,4) -
    		                    2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                    4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                    lambdaphi*pow(vsigma,4),0.5)) +
    		              lambdaphi*pow(vsigma,4)*
    		               ((5 + 10*lambdaphi - 15*pow(lambdaphi,2))*pow(vsigma,2) +
    		                 16*pow(pow(lambdaH,2)*pow(vh,4) -
    		                    2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                    4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                    lambdaphi*pow(vsigma,4),0.5)) -
    		              16*pow(mW,2)*pow(vsigma,2)*
    		               (-2*pow(lambdaphi,2)*pow(vsigma,2) +
    		                 pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5) +
    		                 lambdaphi*(2*pow(vsigma,2) +
    		                    pow(pow(lambdaH,2)*pow(vh,4) -
    		                      2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                      4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                      lambdaphi*pow(vsigma,4),0.5))))) -
    		        (mW + ((1 - cc)*s)/2. + pow(m_h2,2))*
    		         (56*pow(lambdaH,4)*pow(vh,8) +
    		           8*pow(lambdaH,3)*pow(vh,6)*
    		            (-10*pow(mW,2) + 7*
    		               pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5)) +
    		           8*pow(lambdaH,2)*pow(vh,4)*
    		            (6*pow(mW,4) - 10*pow(mW,2)*
    		               pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5) +
    		              7*pow(vsigma,2)*(4*pow(lambdaHphi,2)*pow(vh,2) +
    		                 lambdaphi*(pow(vsigma,2) - lambdaphi*pow(vsigma,2) +
    		                    pow(pow(lambdaH,2)*pow(vh,4) -
    		                      2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                      4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                      lambdaphi*pow(vsigma,4),0.5)))) +
    		           4*lambdaH*pow(vh,2)*
    		            (-14*pow(lambdaphi,3)*pow(vsigma,6) +
    		              12*pow(mW,4)*pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5) +
    		              pow(lambdaphi,2)*pow(vsigma,4)*
    		               (15*pow(mW,2) +
    		                 7*(2*pow(vsigma,2) +
    		                    pow(pow(lambdaH,2)*pow(vh,4) -
    		                      2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                      4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                      lambdaphi*pow(vsigma,4),0.5))) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2)*
    		               (-15*pow(mW,2) +
    		                 7*(2*lambdaphi*pow(vsigma,2) +
    		                    pow(pow(lambdaH,2)*pow(vh,4) -
    		                      2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                      4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                      lambdaphi*pow(vsigma,4),0.5))) +
    		              lambdaphi*(7*pow(vsigma,4)*
    		                  pow(pow(lambdaH,2)*pow(vh,4) -
    		                    2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                    4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                    lambdaphi*pow(vsigma,4),0.5) -
    		                 5*pow(mW,2)*(3*pow(vsigma,4) +
    		                    4*pow(vsigma,2)*
    		                     pow(pow(lambdaH,2)*pow(vh,4) -
    		                       2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                       4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                       lambdaphi*pow(vsigma,4),0.5)))) +
    		           pow(vsigma,2)*(112*pow(lambdaHphi,4)*pow(vh,4)*pow(vsigma,2) +
    		              8*pow(lambdaHphi,2)*pow(vh,2)*
    		               (12*pow(mW,4) -
    		                 10*pow(mW,2)*(3*lambdaphi*pow(vsigma,2) +
    		                    pow(pow(lambdaH,2)*pow(vh,4) -
    		                      2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                      4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                      lambdaphi*pow(vsigma,4),0.5)) +
    		                 7*lambdaphi*pow(vsigma,2)*
    		                  (pow(vsigma,2) + 3*lambdaphi*pow(vsigma,2) +
    		                    2*pow(pow(lambdaH,2)*pow(vh,4) -
    		                       2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                       4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                       lambdaphi*pow(vsigma,4),0.5))) +
    		              lambdaphi*(24*pow(mW,4)*
    		                  ((1 + lambdaphi)*pow(vsigma,2) +
    		                    2*pow(pow(lambdaH,2)*pow(vh,4) -
    		                       2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                       4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                       lambdaphi*pow(vsigma,4),0.5)) +
    		                 7*lambdaphi*pow(vsigma,4)*
    		                  ((1 + 6*lambdaphi + pow(lambdaphi,2))*pow(vsigma,2) +
    		                    4*(1 + lambdaphi)*
    		                     pow(pow(lambdaH,2)*pow(vh,4) -
    		                       2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                       4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                       lambdaphi*pow(vsigma,4),0.5)) -
    		                 20*pow(mW,2)*pow(vsigma,2)*
    		                  (pow(lambdaphi,2)*pow(vsigma,2) +
    		                    pow(pow(lambdaH,2)*pow(vh,4) -
    		                      2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                      4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                      lambdaphi*pow(vsigma,4),0.5) +
    		                    3*lambdaphi*
    		                     (pow(vsigma,2) +
    		                       pow(pow(lambdaH,2)*pow(vh,4) -
    		                         2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                         4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                         lambdaphi*pow(vsigma,4),0.5))))))))*pow(sin(a),4);
    	den=64*pow(mW,4)*pow(-0.5*((1 - cc)*s) - pow(mW,2),2)*
    			   pow(-0.5*((1 + cc)*s) - pow(mW,2),2);
    	amp=num/den;
    }
    else if (id == 12)
    {
    	num=-0.03125*(pow(g,4)*pow(MW,-4)*pow(-1 + pow(SW,2),-6)*
    		      (-2*(MZ + ((1 - cc)*s)/2. + pow(m_h2,2)) + lambdaH*pow(vh,2) +
    		        lambdaphi*pow(vsigma,2) +
    		        pow(pow(lambdaH,2)*pow(vh,4) -
    		          2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		          4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		          lambdaphi*pow(vsigma,4),0.5))*
    		      (2*pow(MZ + ((1 - cc)*s)/2. + pow(m_h2,2),4)*pow(-1 + pow(SW,2),6)*
    		         (lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		           pow(pow(lambdaH,2)*pow(vh,4) -
    		             2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		             4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		             lambdaphi*pow(vsigma,4),0.5)) +
    		        pow(MW,2)*pow(s - 2*pow(m_h2,2),3)*pow(-1 + pow(SW,2),4)*
    		         (-4*pow(MW,2) + 2*(MZ + ((1 - cc)*s)/2. + pow(m_h2,2))*
    		            (-1 + pow(SW,2)) -
    		           (-1 + pow(SW,2))*(lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		              pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5))) -
    		        2*pow(MZ + ((1 - cc)*s)/2. + pow(m_h2,2),3)*pow(-1 + pow(SW,2),5)*
    		         (2*pow(lambdaH,2)*(-1 + pow(SW,2))*pow(vh,4) -
    		           4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		           4*pow(lambdaHphi,2)*pow(SW,2)*pow(vh,2)*pow(vsigma,2) +
    		           pow(lambdaphi,2)*(-1 + pow(SW,2))*pow(vsigma,4) -
    		           2*pow(MW,2)*pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5) +
    		           2*lambdaH*pow(vh,2)*
    		            (-pow(MW,2) + (-1 + pow(SW,2))*
    		               pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5)) +
    		           lambdaphi*pow(vsigma,2)*
    		            (-2*pow(MW,2) + (-1 + pow(SW,2))*
    		               (pow(vsigma,2) +
    		                 2*pow(pow(lambdaH,2)*pow(vh,4) -
    		                    2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                    4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                    lambdaphi*pow(vsigma,4),0.5)))) -
    		        4*(MZ + ((1 - cc)*s)/2. + pow(m_h2,2))*pow(MW,2)*
    		         pow(-1 + pow(SW,2),3)*
    		         (-4*pow(lambdaH,3)*pow(vh,6)*pow(-1 + pow(SW,2),2) -
    		           pow(lambdaphi,3)*pow(vsigma,6)*pow(-1 + pow(SW,2),2) +
    		           2*pow(lambdaH,2)*(-1 + pow(SW,2))*pow(vh,4)*
    		            (-3*pow(MW,2) - 2*(-1 + pow(SW,2))*
    		               pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5)) +
    		           4*pow(lambdaHphi,2)*(-1 + pow(SW,2))*pow(vh,2)*pow(vsigma,2)*
    		            (-3*pow(MW,2) - (-1 + pow(SW,2))*
    		               pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5)) +
    		           pow(lambdaphi,2)*pow(vsigma,4)*
    		            (-3*pow(MW,2)*(-1 + pow(SW,2)) -
    		              3*pow(-1 + pow(SW,2),2)*
    		               (pow(vsigma,2) +
    		                 pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5))) +
    		           lambdaphi*(-1 + pow(SW,2))*pow(vsigma,2)*
    		            (-((-1 + pow(SW,2))*pow(vsigma,2)*
    		                 (12*pow(lambdaHphi,2)*pow(vh,2) +
    		                   pow(pow(lambdaH,2)*pow(vh,4) -
    		                     2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                     4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                     lambdaphi*pow(vsigma,4),0.5))) -
    		              3*pow(MW,2)*(pow(vsigma,2) +
    		                 2*pow(pow(lambdaH,2)*pow(vh,4) -
    		                    2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                    4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                    lambdaphi*pow(vsigma,4),0.5))) +
    		           lambdaH*(-1 + pow(SW,2))*pow(vh,2)*
    		            (-6*pow(MW,2)*pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5) +
    		              (-1 + pow(SW,2))*pow(vsigma,2)*
    		               (-12*pow(lambdaHphi,2)*pow(vh,2) +
    		                 3*(-1 + lambdaphi)*lambdaphi*pow(vsigma,2) -
    		                 4*lambdaphi*pow(pow(lambdaH,2)*pow(vh,4) -
    		                    2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                    4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                    lambdaphi*pow(vsigma,4),0.5)))) +
    		        pow(s - 2*pow(m_h2,2),2)*pow(-1 + pow(SW,2),4)*
    		         (-(pow(MZ + ((1 - cc)*s)/2. + pow(m_h2,2),3)*
    		              pow(-1 + pow(SW,2),2)) -
    		           8*pow(MW,4)*(lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		              pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5)) -
    		           2*(MZ + ((1 - cc)*s)/2. + pow(m_h2,2))*pow(MW,2)*
    		            (-8*pow(MW,2) - 5*(-1 + pow(SW,2))*
    		               (lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		                 pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5))) +
    		           ((-1 + pow(SW,2))*pow(MZ + ((1 - cc)*s)/2. + pow(m_h2,2),2)*
    		              (-8*pow(MW,2) + (-1 + pow(SW,2))*
    		                 (lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		                   pow(pow(lambdaH,2)*pow(vh,4) -
    		                     2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                     4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                     lambdaphi*pow(vsigma,4),0.5))))/2. -
    		           3*pow(MW,2)*(-1 + pow(SW,2))*
    		            (2*pow(lambdaH,2)*pow(vh,4) +
    		              2*lambdaH*pow(vh,2)*
    		               pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5) +
    		              pow(vsigma,2)*(4*pow(lambdaHphi,2)*pow(vh,2) +
    		                 lambdaphi*(1 + lambdaphi)*pow(vsigma,2) +
    		                 2*lambdaphi*pow(pow(lambdaH,2)*pow(vh,4) -
    		                    2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                    4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                    lambdaphi*pow(vsigma,4),0.5)))) -
    		        (pow(MZ + ((1 - cc)*s)/2. + pow(m_h2,2),2)*pow(-1 + pow(SW,2),4)*
    		           (16*pow(MW,4)*(lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		                pow(pow(lambdaH,2)*pow(vh,4) -
    		                  2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                  4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                  lambdaphi*pow(vsigma,4),0.5)) +
    		             12*pow(MW,2)*(-1 + pow(SW,2))*
    		              (2*pow(lambdaH,2)*pow(vh,4) +
    		                2*lambdaH*pow(vh,2)*
    		                 pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5) +
    		                pow(vsigma,2)*(4*pow(lambdaHphi,2)*pow(vh,2) +
    		                   lambdaphi*(1 + lambdaphi)*pow(vsigma,2) +
    		                   2*lambdaphi*
    		                    pow(pow(lambdaH,2)*pow(vh,4) -
    		                      2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                      4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                      lambdaphi*pow(vsigma,4),0.5))) +
    		             (1 - pow(SW,2))*(-1 + pow(SW,2))*
    		              (4*pow(lambdaH,3)*pow(vh,6) +
    		                4*pow(lambdaH,2)*pow(vh,4)*
    		                 pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2)*
    		                 (3*lambdaphi*pow(vsigma,2) +
    		                   pow(pow(lambdaH,2)*pow(vh,4) -
    		                     2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                     4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                     lambdaphi*pow(vsigma,4),0.5)) +
    		                lambdaH*pow(vh,2)*pow(vsigma,2)*
    		                 (12*pow(lambdaHphi,2)*pow(vh,2) -
    		                   3*(-1 + lambdaphi)*lambdaphi*pow(vsigma,2) +
    		                   4*lambdaphi*
    		                    pow(pow(lambdaH,2)*pow(vh,4) -
    		                      2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                      4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                      lambdaphi*pow(vsigma,4),0.5)) +
    		                lambdaphi*pow(vsigma,4)*
    		                 (pow(lambdaphi,2)*pow(vsigma,2) +
    		                   pow(pow(lambdaH,2)*pow(vh,4) -
    		                     2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                     4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                     lambdaphi*pow(vsigma,4),0.5) +
    		                   3*lambdaphi*
    		                    (pow(vsigma,2) +
    		                      pow(pow(lambdaH,2)*pow(vh,4) -
    		                        2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                        4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                        lambdaphi*pow(vsigma,4),0.5))))))/2. +
    		        pow(MW,2)*pow(-1 + pow(SW,2),3)*
    		         (2*pow(MW,2)*(2 - 2*pow(SW,2))*
    		            (4*pow(lambdaH,3)*pow(vh,6) +
    		              4*pow(lambdaH,2)*pow(vh,4)*
    		               pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2)*
    		               (3*lambdaphi*pow(vsigma,2) +
    		                 pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5)) +
    		              lambdaH*pow(vh,2)*pow(vsigma,2)*
    		               (12*pow(lambdaHphi,2)*pow(vh,2) -
    		                 3*(-1 + lambdaphi)*lambdaphi*pow(vsigma,2) +
    		                 4*lambdaphi*pow(pow(lambdaH,2)*pow(vh,4) -
    		                    2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                    4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                    lambdaphi*pow(vsigma,4),0.5)) +
    		              lambdaphi*pow(vsigma,4)*
    		               (pow(lambdaphi,2)*pow(vsigma,2) +
    		                 pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5) +
    		                 3*lambdaphi*(pow(vsigma,2) +
    		                    pow(pow(lambdaH,2)*pow(vh,4) -
    		                      2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                      4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                      lambdaphi*pow(vsigma,4),0.5)))) +
    		           (1 - pow(SW,2))*(-1 + pow(SW,2))*
    		            (8*pow(lambdaH,4)*pow(vh,8) +
    		              8*pow(lambdaH,3)*pow(vh,6)*
    		               pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5) +
    		              8*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,2)*
    		               (4*pow(lambdaHphi,2)*pow(vh,2) +
    		                 lambdaphi*(pow(vsigma,2) - lambdaphi*pow(vsigma,2) +
    		                    pow(pow(lambdaH,2)*pow(vh,4) -
    		                      2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                      4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                      lambdaphi*pow(vsigma,4),0.5))) +
    		              4*lambdaH*(4*pow(lambdaHphi,2)*pow(vh,4)*pow(vsigma,2)*
    		                  (2*lambdaphi*pow(vsigma,2) +
    		                    pow(pow(lambdaH,2)*pow(vh,4) -
    		                      2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                      4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                      lambdaphi*pow(vsigma,4),0.5)) +
    		                 lambdaphi*pow(vh,2)*pow(vsigma,4)*
    		                  (2*lambdaphi*pow(vsigma,2) -
    		                    2*pow(lambdaphi,2)*pow(vsigma,2) +
    		                    pow(pow(lambdaH,2)*pow(vh,4) -
    		                      2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                      4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                      lambdaphi*pow(vsigma,4),0.5) +
    		                    lambdaphi*pow(pow(lambdaH,2)*pow(vh,4) -
    		                       2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                       4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                       lambdaphi*pow(vsigma,4),0.5))) +
    		              pow(vsigma,4)*(16*pow(lambdaHphi,4)*pow(vh,4) +
    		                 8*lambdaphi*pow(lambdaHphi,2)*pow(vh,2)*
    		                  ((1 + 3*lambdaphi)*pow(vsigma,2) +
    		                    2*pow(pow(lambdaH,2)*pow(vh,4) -
    		                       2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                       4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                       lambdaphi*pow(vsigma,4),0.5)) +
    		                 pow(lambdaphi,2)*pow(vsigma,2)*
    		                  ((1 + 6*lambdaphi + pow(lambdaphi,2))*pow(vsigma,2) +
    		                    4*(1 + lambdaphi)*
    		                     pow(pow(lambdaH,2)*pow(vh,4) -
    		                       2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                       4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                       lambdaphi*pow(vsigma,4),0.5))))) -
    		        (s - 2*pow(m_h2,2))*(-1 + pow(SW,2))*
    		         (-(pow(MZ + ((1 - cc)*s)/2. + pow(m_h2,2),4)*
    		              pow(-1 + pow(SW,2),5)) +
    		           pow(MZ + ((1 - cc)*s)/2. + pow(m_h2,2),3)*pow(-1 + pow(SW,2),4)*
    		            (-4*pow(MW,2) + 3*(-1 + pow(SW,2))*
    		               (lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		                 pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5))) +
    		           2*pow(MW,4)*pow(-1 + pow(SW,2),2)*
    		            (-6*pow(lambdaH,2)*(-1 + pow(SW,2))*pow(vh,4) -
    		              3*pow(lambdaphi,2)*(-1 + pow(SW,2))*pow(vsigma,4) +
    		              4*(-3*pow(lambdaHphi,2)*(-1 + pow(SW,2))*pow(vh,2)*
    		                  pow(vsigma,2) +
    		                 pow(MW,2)*pow(pow(lambdaH,2)*pow(vh,4) -
    		                    2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                    4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                    lambdaphi*pow(vsigma,4),0.5)) +
    		              2*lambdaH*pow(vh,2)*
    		               (2*pow(MW,2) - 3*(-1 + pow(SW,2))*
    		                  pow(pow(lambdaH,2)*pow(vh,4) -
    		                    2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                    4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                    lambdaphi*pow(vsigma,4),0.5)) +
    		              lambdaphi*pow(vsigma,2)*
    		               (4*pow(MW,2) - 3*(-1 + pow(SW,2))*
    		                  (pow(vsigma,2) +
    		                    2*pow(pow(lambdaH,2)*pow(vh,4) -
    		                       2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                       4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                       lambdaphi*pow(vsigma,4),0.5)))) +
    		           4*(MZ + ((1 - cc)*s)/2. + pow(m_h2,2))*pow(MW,2)*
    		            pow(-1 + pow(SW,2),3)*
    		            (-7*pow(MW,2)*(lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		                 pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5)) -
    		              3*(-1 + pow(SW,2))*
    		               (2*pow(lambdaH,2)*pow(vh,4) +
    		                 2*lambdaH*pow(vh,2)*
    		                  pow(pow(lambdaH,2)*pow(vh,4) -
    		                    2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                    4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                    lambdaphi*pow(vsigma,4),0.5) +
    		                 pow(vsigma,2)*
    		                  (4*pow(lambdaHphi,2)*pow(vh,2) +
    		                    lambdaphi*(1 + lambdaphi)*pow(vsigma,2) +
    		                    2*lambdaphi*
    		                     pow(pow(lambdaH,2)*pow(vh,4) -
    		                       2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                       4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                       lambdaphi*pow(vsigma,4),0.5)))) +
    		           pow(MZ + ((1 - cc)*s)/2. + pow(m_h2,2),2)*pow(-1 + pow(SW,2),3)*
    		            (16*pow(MW,4) + 14*pow(MW,2)*(-1 + pow(SW,2))*
    		               (lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		                 pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5)) -
    		              pow(-1 + pow(SW,2),2)*
    		               (2*pow(lambdaH,2)*pow(vh,4) +
    		                 2*lambdaH*pow(vh,2)*
    		                  pow(pow(lambdaH,2)*pow(vh,4) -
    		                    2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                    4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                    lambdaphi*pow(vsigma,4),0.5) +
    		                 pow(vsigma,2)*
    		                  (4*pow(lambdaHphi,2)*pow(vh,2) +
    		                    lambdaphi*(1 + lambdaphi)*pow(vsigma,2) +
    		                    2*lambdaphi*
    		                     pow(pow(lambdaH,2)*pow(vh,4) -
    		                       2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                       4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                       lambdaphi*pow(vsigma,4),0.5)))) -
    		           pow(MW,2)*(1 - pow(SW,2))*
    		            (48*pow(MW,6) + 8*pow(MW,4)*(-1 + pow(SW,2))*
    		               (lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		                 pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5)) +
    		              18*pow(MW,2)*pow(-1 + pow(SW,2),2)*
    		               (2*pow(lambdaH,2)*pow(vh,4) +
    		                 2*lambdaH*pow(vh,2)*
    		                  pow(pow(lambdaH,2)*pow(vh,4) -
    		                    2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                    4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                    lambdaphi*pow(vsigma,4),0.5) +
    		                 pow(vsigma,2)*
    		                  (4*pow(lambdaHphi,2)*pow(vh,2) +
    		                    lambdaphi*(1 + lambdaphi)*pow(vsigma,2) +
    		                    2*lambdaphi*
    		                     pow(pow(lambdaH,2)*pow(vh,4) -
    		                       2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                       4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                       lambdaphi*pow(vsigma,4),0.5))) +
    		              3*pow(-1 + pow(SW,2),3)*
    		               (4*pow(lambdaH,3)*pow(vh,6) +
    		                 4*pow(lambdaH,2)*pow(vh,4)*
    		                  pow(pow(lambdaH,2)*pow(vh,4) -
    		                    2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                    4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                    lambdaphi*pow(vsigma,4),0.5) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2)*
    		                  (3*lambdaphi*pow(vsigma,2) +
    		                    pow(pow(lambdaH,2)*pow(vh,4) -
    		                      2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                      4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                      lambdaphi*pow(vsigma,4),0.5)) +
    		                 lambdaH*pow(vh,2)*pow(vsigma,2)*
    		                  (12*pow(lambdaHphi,2)*pow(vh,2) -
    		                    3*(-1 + lambdaphi)*lambdaphi*pow(vsigma,2) +
    		                    4*lambdaphi*
    		                     pow(pow(lambdaH,2)*pow(vh,4) -
    		                       2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                       4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                       lambdaphi*pow(vsigma,4),0.5)) +
    		                 lambdaphi*pow(vsigma,4)*
    		                  (pow(lambdaphi,2)*pow(vsigma,2) +
    		                    pow(pow(lambdaH,2)*pow(vh,4) -
    		                      2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                      4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                      lambdaphi*pow(vsigma,4),0.5) +
    		                    3*lambdaphi*
    		                     (pow(vsigma,2) +
    		                       pow(pow(lambdaH,2)*pow(vh,4) -
    		                         2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                         4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                         lambdaphi*pow(vsigma,4),0.5)))))))*
    		      pow(-MZ + s - ((1 - cc)*s)/2. - 3*pow(m_h2,2) +
    		        (lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		           pow(pow(lambdaH,2)*pow(vh,4) +
    		             2*(-(lambdaH*lambdaphi) + 2*pow(lambdaHphi,2))*pow(vh,2)*
    		              pow(vsigma,2) + lambdaphi*pow(vsigma,4),0.5))/2.,-1)*
    		      pow(-MZ - ((1 - cc)*s)/2. - pow(m_h2,2) +
    		        pow(MW,2)*pow(1 - pow(SW,2),-1) + pow(MW,2)*pow(-1 + pow(SW,2),-1) +
    		        (lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		           pow(pow(lambdaH,2)*pow(vh,4) +
    		             2*(-(lambdaH*lambdaphi) + 2*pow(lambdaHphi,2))*pow(vh,2)*
    		              pow(vsigma,2) + lambdaphi*pow(vsigma,4),0.5))/2.,-3)*
    		      pow(sin(a),4)) + (pow(g,4)*pow(MW,-4)*pow(-1 + pow(SW,2),-6)*
    		      (2*(s - 2*pow(m_h2,2)) - 2*(MZ + ((1 - cc)*s)/2. + pow(m_h2,2)) +
    		        lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		        pow(pow(lambdaH,2)*pow(vh,4) -
    		          2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		          4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		          lambdaphi*pow(vsigma,4),0.5))*
    		      (pow(MW,2)*pow(s - 2*pow(m_h2,2),3)*pow(-1 + pow(SW,2),4)*
    		         (-4*pow(MW,2) + 2*(s - 2*pow(m_h2,2))*(-1 + pow(SW,2)) -
    		           2*(MZ + ((1 - cc)*s)/2. + pow(m_h2,2))*(-1 + pow(SW,2)) -
    		           lambdaH*pow(vh,2) + lambdaH*pow(SW,2)*pow(vh,2) -
    		           lambdaphi*pow(vsigma,2) + lambdaphi*pow(SW,2)*pow(vsigma,2) -
    		           pow(pow(lambdaH,2)*pow(vh,4) -
    		             2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		             4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		             lambdaphi*pow(vsigma,4),0.5) +
    		           pow(SW,2)*pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5)) +
    		        pow(MW,2)*pow(-1 + pow(SW,2),3)*
    		         (2*pow(MW,2)*(2 - 2*pow(SW,2))*
    		            (4*pow(lambdaH,3)*pow(vh,6) +
    		              4*pow(lambdaH,2)*pow(vh,4)*
    		               pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2)*
    		               (3*lambdaphi*pow(vsigma,2) +
    		                 pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5)) +
    		              lambdaH*pow(vh,2)*pow(vsigma,2)*
    		               (12*pow(lambdaHphi,2)*pow(vh,2) -
    		                 3*(-1 + lambdaphi)*lambdaphi*pow(vsigma,2) +
    		                 4*lambdaphi*pow(pow(lambdaH,2)*pow(vh,4) -
    		                    2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                    4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                    lambdaphi*pow(vsigma,4),0.5)) +
    		              lambdaphi*pow(vsigma,4)*
    		               (pow(lambdaphi,2)*pow(vsigma,2) +
    		                 pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5) +
    		                 3*lambdaphi*(pow(vsigma,2) +
    		                    pow(pow(lambdaH,2)*pow(vh,4) -
    		                      2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                      4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                      lambdaphi*pow(vsigma,4),0.5)))) +
    		           (1 - pow(SW,2))*(-1 + pow(SW,2))*
    		            (8*pow(lambdaH,4)*pow(vh,8) +
    		              8*pow(lambdaH,3)*pow(vh,6)*
    		               pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5) +
    		              8*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,2)*
    		               (4*pow(lambdaHphi,2)*pow(vh,2) +
    		                 lambdaphi*(pow(vsigma,2) - lambdaphi*pow(vsigma,2) +
    		                    pow(pow(lambdaH,2)*pow(vh,4) -
    		                      2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                      4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                      lambdaphi*pow(vsigma,4),0.5))) +
    		              4*lambdaH*(4*pow(lambdaHphi,2)*pow(vh,4)*pow(vsigma,2)*
    		                  (2*lambdaphi*pow(vsigma,2) +
    		                    pow(pow(lambdaH,2)*pow(vh,4) -
    		                      2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                      4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                      lambdaphi*pow(vsigma,4),0.5)) +
    		                 lambdaphi*pow(vh,2)*pow(vsigma,4)*
    		                  (2*lambdaphi*pow(vsigma,2) -
    		                    2*pow(lambdaphi,2)*pow(vsigma,2) +
    		                    pow(pow(lambdaH,2)*pow(vh,4) -
    		                      2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                      4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                      lambdaphi*pow(vsigma,4),0.5) +
    		                    lambdaphi*pow(pow(lambdaH,2)*pow(vh,4) -
    		                       2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                       4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                       lambdaphi*pow(vsigma,4),0.5))) +
    		              pow(vsigma,4)*(16*pow(lambdaHphi,4)*pow(vh,4) +
    		                 8*lambdaphi*pow(lambdaHphi,2)*pow(vh,2)*
    		                  ((1 + 3*lambdaphi)*pow(vsigma,2) +
    		                    2*pow(pow(lambdaH,2)*pow(vh,4) -
    		                       2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                       4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                       lambdaphi*pow(vsigma,4),0.5)) +
    		                 pow(lambdaphi,2)*pow(vsigma,2)*
    		                  ((1 + 6*lambdaphi + pow(lambdaphi,2))*pow(vsigma,2) +
    		                    4*(1 + lambdaphi)*
    		                     pow(pow(lambdaH,2)*pow(vh,4) -
    		                       2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                       4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                       lambdaphi*pow(vsigma,4),0.5))))) -
    		        8*pow(MW,2)*pow(-1 + pow(SW,2),3)*
    		         (-4*pow(lambdaH,3)*pow(vh,6)*pow(-1 + pow(SW,2),2) -
    		           pow(lambdaphi,3)*pow(vsigma,6)*pow(-1 + pow(SW,2),2) +
    		           2*pow(lambdaH,2)*(-1 + pow(SW,2))*pow(vh,4)*
    		            (-3*pow(MW,2) - 2*(-1 + pow(SW,2))*
    		               pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5)) +
    		           4*pow(lambdaHphi,2)*(-1 + pow(SW,2))*pow(vh,2)*pow(vsigma,2)*
    		            (-3*pow(MW,2) - (-1 + pow(SW,2))*
    		               pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5)) +
    		           pow(lambdaphi,2)*pow(vsigma,4)*
    		            (-3*pow(MW,2)*(-1 + pow(SW,2)) -
    		              3*pow(-1 + pow(SW,2),2)*
    		               (pow(vsigma,2) +
    		                 pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5))) +
    		           lambdaphi*(-1 + pow(SW,2))*pow(vsigma,2)*
    		            (-((-1 + pow(SW,2))*pow(vsigma,2)*
    		                 (12*pow(lambdaHphi,2)*pow(vh,2) +
    		                   pow(pow(lambdaH,2)*pow(vh,4) -
    		                     2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                     4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                     lambdaphi*pow(vsigma,4),0.5))) -
    		              3*pow(MW,2)*(pow(vsigma,2) +
    		                 2*pow(pow(lambdaH,2)*pow(vh,4) -
    		                    2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                    4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                    lambdaphi*pow(vsigma,4),0.5))) +
    		           lambdaH*(-1 + pow(SW,2))*pow(vh,2)*
    		            (-6*pow(MW,2)*pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5) +
    		              (-1 + pow(SW,2))*pow(vsigma,2)*
    		               (-12*pow(lambdaHphi,2)*pow(vh,2) +
    		                 3*(-1 + lambdaphi)*lambdaphi*pow(vsigma,2) -
    		                 4*lambdaphi*pow(pow(lambdaH,2)*pow(vh,4) -
    		                    2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                    4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                    lambdaphi*pow(vsigma,4),0.5))))*
    		         ((s - 2*pow(m_h2,2))/2. + (-MZ - ((1 - cc)*s)/2. - pow(m_h2,2))/2. +
    		           (lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		              pow(pow(lambdaH,2)*pow(vh,4) +
    		                2*(-(lambdaH*lambdaphi) + 2*pow(lambdaHphi,2))*pow(vh,2)*
    		                 pow(vsigma,2) + lambdaphi*pow(vsigma,4),0.5))/2.) +
    		        (pow(s - 2*pow(m_h2,2),2)*pow(-1 + pow(SW,2),4)*
    		           (-16*pow(MW,4)*(lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		                pow(pow(lambdaH,2)*pow(vh,4) -
    		                  2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                  4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                  lambdaphi*pow(vsigma,4),0.5)) -
    		             8*pow(MW,2)*(-8*pow(MW,2) -
    		                5*(-1 + pow(SW,2))*
    		                 (lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		                   pow(pow(lambdaH,2)*pow(vh,4) -
    		                     2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                     4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                     lambdaphi*pow(vsigma,4),0.5)))*
    		              ((s - 2*pow(m_h2,2))/2. +
    		                (-MZ - ((1 - cc)*s)/2. - pow(m_h2,2) + lambdaH*pow(vh,2) +
    		                   lambdaphi*pow(vsigma,2) +
    		                   pow(pow(lambdaH,2)*pow(vh,4) -
    		                     2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                     4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                     lambdaphi*pow(vsigma,4),0.5))/2.) -
    		             6*pow(MW,2)*(-1 + pow(SW,2))*
    		              (2*pow(lambdaH,2)*pow(vh,4) +
    		                2*lambdaH*pow(vh,2)*
    		                 pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5) +
    		                pow(vsigma,2)*(4*pow(lambdaHphi,2)*pow(vh,2) +
    		                   lambdaphi*(1 + lambdaphi)*pow(vsigma,2) +
    		                   2*lambdaphi*
    		                    pow(pow(lambdaH,2)*pow(vh,4) -
    		                      2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                      4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                      lambdaphi*pow(vsigma,4),0.5))) +
    		             (-1 + pow(SW,2))*(-8*pow(MW,2) +
    		                (-1 + pow(SW,2))*
    		                 (lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		                   pow(pow(lambdaH,2)*pow(vh,4) -
    		                     2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                     4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                     lambdaphi*pow(vsigma,4),0.5)))*
    		              pow(-MZ + s - ((1 - cc)*s)/2. - 3*pow(m_h2,2) +
    		                lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		                pow(pow(lambdaH,2)*pow(vh,4) -
    		                  2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                  4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                  lambdaphi*pow(vsigma,4),0.5),2) -
    		             2*pow(-1 + pow(SW,2),2)*
    		              pow(-MZ + s - ((1 - cc)*s)/2. - 3*pow(m_h2,2) +
    		                lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		                pow(pow(lambdaH,2)*pow(vh,4) -
    		                  2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                  4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                  lambdaphi*pow(vsigma,4),0.5),3)))/2. +
    		        2*pow(-1 + pow(SW,2),6)*
    		         (lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		           pow(pow(lambdaH,2)*pow(vh,4) -
    		             2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		             4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		             lambdaphi*pow(vsigma,4),0.5))*
    		         pow(-MZ + s - ((1 - cc)*s)/2. - 3*pow(m_h2,2) + lambdaH*pow(vh,2) +
    		           lambdaphi*pow(vsigma,2) +
    		           pow(pow(lambdaH,2)*pow(vh,4) -
    		             2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		             4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		             lambdaphi*pow(vsigma,4),0.5),4) -
    		        2*pow(-1 + pow(SW,2),4)*
    		         (16*pow(MW,4)*(lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		              pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5)) +
    		           12*pow(MW,2)*(-1 + pow(SW,2))*
    		            (2*pow(lambdaH,2)*pow(vh,4) +
    		              2*lambdaH*pow(vh,2)*
    		               pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5) +
    		              pow(vsigma,2)*(4*pow(lambdaHphi,2)*pow(vh,2) +
    		                 lambdaphi*(1 + lambdaphi)*pow(vsigma,2) +
    		                 2*lambdaphi*pow(pow(lambdaH,2)*pow(vh,4) -
    		                    2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                    4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                    lambdaphi*pow(vsigma,4),0.5))) +
    		           (1 - pow(SW,2))*(-1 + pow(SW,2))*
    		            (4*pow(lambdaH,3)*pow(vh,6) +
    		              4*pow(lambdaH,2)*pow(vh,4)*
    		               pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2)*
    		               (3*lambdaphi*pow(vsigma,2) +
    		                 pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5)) +
    		              lambdaH*pow(vh,2)*pow(vsigma,2)*
    		               (12*pow(lambdaHphi,2)*pow(vh,2) -
    		                 3*(-1 + lambdaphi)*lambdaphi*pow(vsigma,2) +
    		                 4*lambdaphi*pow(pow(lambdaH,2)*pow(vh,4) -
    		                    2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                    4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                    lambdaphi*pow(vsigma,4),0.5)) +
    		              lambdaphi*pow(vsigma,4)*
    		               (pow(lambdaphi,2)*pow(vsigma,2) +
    		                 pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5) +
    		                 3*lambdaphi*(pow(vsigma,2) +
    		                    pow(pow(lambdaH,2)*pow(vh,4) -
    		                      2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                      4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                      lambdaphi*pow(vsigma,4),0.5)))))*
    		         pow((s - 2*pow(m_h2,2))/2. +
    		           (-MZ - ((1 - cc)*s)/2. - pow(m_h2,2))/2. +
    		           (lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		              pow(pow(lambdaH,2)*pow(vh,4) +
    		                2*(-(lambdaH*lambdaphi) + 2*pow(lambdaHphi,2))*pow(vh,2)*
    		                 pow(vsigma,2) + lambdaphi*pow(vsigma,4),0.5))/2.,2) -
    		        16*pow(-1 + pow(SW,2),5)*
    		         (2*pow(lambdaH,2)*(-1 + pow(SW,2))*pow(vh,4) -
    		           4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		           4*pow(lambdaHphi,2)*pow(SW,2)*pow(vh,2)*pow(vsigma,2) +
    		           pow(lambdaphi,2)*(-1 + pow(SW,2))*pow(vsigma,4) -
    		           2*pow(MW,2)*pow(pow(lambdaH,2)*pow(vh,4) -
    		              2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		              4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		              lambdaphi*pow(vsigma,4),0.5) +
    		           2*lambdaH*pow(vh,2)*
    		            (-pow(MW,2) + (-1 + pow(SW,2))*
    		               pow(pow(lambdaH,2)*pow(vh,4) -
    		                 2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                 lambdaphi*pow(vsigma,4),0.5)) +
    		           lambdaphi*pow(vsigma,2)*
    		            (-2*pow(MW,2) + (-1 + pow(SW,2))*
    		               (pow(vsigma,2) +
    		                 2*pow(pow(lambdaH,2)*pow(vh,4) -
    		                    2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                    4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                    lambdaphi*pow(vsigma,4),0.5))))*
    		         pow((s - 2*pow(m_h2,2))/2. +
    		           (-MZ - ((1 - cc)*s)/2. - pow(m_h2,2))/2. +
    		           (lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		              pow(pow(lambdaH,2)*pow(vh,4) +
    		                2*(-(lambdaH*lambdaphi) + 2*pow(lambdaHphi,2))*pow(vh,2)*
    		                 pow(vsigma,2) + lambdaphi*pow(vsigma,4),0.5))/2.,3) -
    		        (s - 2*pow(m_h2,2))*(-1 + pow(SW,2))*
    		         (2*pow(MW,4)*pow(-1 + pow(SW,2),2)*
    		            (-6*pow(lambdaH,2)*(-1 + pow(SW,2))*pow(vh,4) -
    		              3*pow(lambdaphi,2)*(-1 + pow(SW,2))*pow(vsigma,4) +
    		              4*(-3*pow(lambdaHphi,2)*(-1 + pow(SW,2))*pow(vh,2)*
    		                  pow(vsigma,2) +
    		                 pow(MW,2)*pow(pow(lambdaH,2)*pow(vh,4) -
    		                    2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                    4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                    lambdaphi*pow(vsigma,4),0.5)) +
    		              2*lambdaH*pow(vh,2)*
    		               (2*pow(MW,2) - 3*(-1 + pow(SW,2))*
    		                  pow(pow(lambdaH,2)*pow(vh,4) -
    		                    2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                    4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                    lambdaphi*pow(vsigma,4),0.5)) +
    		              lambdaphi*pow(vsigma,2)*
    		               (4*pow(MW,2) - 3*(-1 + pow(SW,2))*
    		                  (pow(vsigma,2) +
    		                    2*pow(pow(lambdaH,2)*pow(vh,4) -
    		                       2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                       4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                       lambdaphi*pow(vsigma,4),0.5)))) -
    		           pow(MW,2)*(1 - pow(SW,2))*
    		            (48*pow(MW,6) + 8*pow(MW,4)*(-1 + pow(SW,2))*
    		               (lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		                 pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5)) +
    		              18*pow(MW,2)*pow(-1 + pow(SW,2),2)*
    		               (2*pow(lambdaH,2)*pow(vh,4) +
    		                 2*lambdaH*pow(vh,2)*
    		                  pow(pow(lambdaH,2)*pow(vh,4) -
    		                    2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                    4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                    lambdaphi*pow(vsigma,4),0.5) +
    		                 pow(vsigma,2)*
    		                  (4*pow(lambdaHphi,2)*pow(vh,2) +
    		                    lambdaphi*(1 + lambdaphi)*pow(vsigma,2) +
    		                    2*lambdaphi*
    		                     pow(pow(lambdaH,2)*pow(vh,4) -
    		                       2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                       4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                       lambdaphi*pow(vsigma,4),0.5))) +
    		              3*pow(-1 + pow(SW,2),3)*
    		               (4*pow(lambdaH,3)*pow(vh,6) +
    		                 4*pow(lambdaH,2)*pow(vh,4)*
    		                  pow(pow(lambdaH,2)*pow(vh,4) -
    		                    2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                    4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                    lambdaphi*pow(vsigma,4),0.5) +
    		                 4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2)*
    		                  (3*lambdaphi*pow(vsigma,2) +
    		                    pow(pow(lambdaH,2)*pow(vh,4) -
    		                      2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                      4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                      lambdaphi*pow(vsigma,4),0.5)) +
    		                 lambdaH*pow(vh,2)*pow(vsigma,2)*
    		                  (12*pow(lambdaHphi,2)*pow(vh,2) -
    		                    3*(-1 + lambdaphi)*lambdaphi*pow(vsigma,2) +
    		                    4*lambdaphi*
    		                     pow(pow(lambdaH,2)*pow(vh,4) -
    		                       2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                       4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                       lambdaphi*pow(vsigma,4),0.5)) +
    		                 lambdaphi*pow(vsigma,4)*
    		                  (pow(lambdaphi,2)*pow(vsigma,2) +
    		                    pow(pow(lambdaH,2)*pow(vh,4) -
    		                      2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                      4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                      lambdaphi*pow(vsigma,4),0.5) +
    		                    3*lambdaphi*
    		                     (pow(vsigma,2) +
    		                       pow(pow(lambdaH,2)*pow(vh,4) -
    		                         2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                         4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                         lambdaphi*pow(vsigma,4),0.5))))) +
    		           8*pow(MW,2)*pow(-1 + pow(SW,2),3)*
    		            (-7*pow(MW,2)*(lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		                 pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5)) -
    		              3*(-1 + pow(SW,2))*
    		               (2*pow(lambdaH,2)*pow(vh,4) +
    		                 2*lambdaH*pow(vh,2)*
    		                  pow(pow(lambdaH,2)*pow(vh,4) -
    		                    2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                    4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                    lambdaphi*pow(vsigma,4),0.5) +
    		                 pow(vsigma,2)*
    		                  (4*pow(lambdaHphi,2)*pow(vh,2) +
    		                    lambdaphi*(1 + lambdaphi)*pow(vsigma,2) +
    		                    2*lambdaphi*
    		                     pow(pow(lambdaH,2)*pow(vh,4) -
    		                       2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                       4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                       lambdaphi*pow(vsigma,4),0.5))))*
    		            ((s - 2*pow(m_h2,2))/2. +
    		              (-MZ - ((1 - cc)*s)/2. - pow(m_h2,2))/2. +
    		              (lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		                 pow(pow(lambdaH,2)*pow(vh,4) +
    		                   2*(-(lambdaH*lambdaphi) + 2*pow(lambdaHphi,2))*pow(vh,2)*
    		                    pow(vsigma,2) + lambdaphi*pow(vsigma,4),0.5))/2.) -
    		           pow(-1 + pow(SW,2),5)*
    		            pow(-MZ + s - ((1 - cc)*s)/2. - 3*pow(m_h2,2) +
    		              lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		              pow(pow(lambdaH,2)*pow(vh,4) -
    		                2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                lambdaphi*pow(vsigma,4),0.5),4) +
    		           4*pow(-1 + pow(SW,2),3)*
    		            (16*pow(MW,4) + 14*pow(MW,2)*(-1 + pow(SW,2))*
    		               (lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		                 pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5)) -
    		              pow(-1 + pow(SW,2),2)*
    		               (2*pow(lambdaH,2)*pow(vh,4) +
    		                 2*lambdaH*pow(vh,2)*
    		                  pow(pow(lambdaH,2)*pow(vh,4) -
    		                    2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                    4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                    lambdaphi*pow(vsigma,4),0.5) +
    		                 pow(vsigma,2)*
    		                  (4*pow(lambdaHphi,2)*pow(vh,2) +
    		                    lambdaphi*(1 + lambdaphi)*pow(vsigma,2) +
    		                    2*lambdaphi*
    		                     pow(pow(lambdaH,2)*pow(vh,4) -
    		                       2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                       4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                       lambdaphi*pow(vsigma,4),0.5))))*
    		            pow((s - 2*pow(m_h2,2))/2. +
    		              (-MZ - ((1 - cc)*s)/2. - pow(m_h2,2))/2. +
    		              (lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		                 pow(pow(lambdaH,2)*pow(vh,4) +
    		                   2*(-(lambdaH*lambdaphi) + 2*pow(lambdaHphi,2))*pow(vh,2)*
    		                    pow(vsigma,2) + lambdaphi*pow(vsigma,4),0.5))/2.,2) +
    		           8*pow(-1 + pow(SW,2),4)*
    		            (-4*pow(MW,2) + 3*(-1 + pow(SW,2))*
    		               (lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		                 pow(pow(lambdaH,2)*pow(vh,4) -
    		                   2*lambdaH*lambdaphi*pow(vh,2)*pow(vsigma,2) +
    		                   4*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		                   lambdaphi*pow(vsigma,4),0.5)))*
    		            pow((s - 2*pow(m_h2,2))/2. +
    		              (-MZ - ((1 - cc)*s)/2. - pow(m_h2,2))/2. +
    		              (lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		                 pow(pow(lambdaH,2)*pow(vh,4) +
    		                   2*(-(lambdaH*lambdaphi) + 2*pow(lambdaHphi,2))*pow(vh,2)*
    		                    pow(vsigma,2) + lambdaphi*pow(vsigma,4),0.5))/2.,3)))*
    		      pow(-MZ + s - ((1 - cc)*s)/2. - 3*pow(m_h2,2) +
    		        (lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		           pow(pow(lambdaH,2)*pow(vh,4) +
    		             2*(-(lambdaH*lambdaphi) + 2*pow(lambdaHphi,2))*pow(vh,2)*
    		              pow(vsigma,2) + lambdaphi*pow(vsigma,4),0.5))/2.,-3)*
    		      pow(-MZ - ((1 - cc)*s)/2. - pow(m_h2,2) +
    		        pow(MW,2)*pow(1 - pow(SW,2),-1) + pow(MW,2)*pow(-1 + pow(SW,2),-1) +
    		        (lambdaH*pow(vh,2) + lambdaphi*pow(vsigma,2) +
    		           pow(pow(lambdaH,2)*pow(vh,4) +
    		             2*(-(lambdaH*lambdaphi) + 2*pow(lambdaHphi,2))*pow(vh,2)*
    		              pow(vsigma,2) + lambdaphi*pow(vsigma,4),0.5))/2.,-1)*
    		      pow(sin(a),4))/32.;
    	den=2;
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
        m1 = params->m_h2; m2 = params->m_h2; m3 = mH; m4 = mH; S12 = 0.5l; S34=0.5l;
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
    else
    {
        m1 = 0.l; m2 = 0.l; m3 = 0.l; m4 = 0.l; S12 = 0.l; S34=0.l;
    }
    smin = pow(std::max(m1+m2,m3+m4)+1e-4,2.0l);
    smax = params->temp*params->temp*1.0e100;
    s = smin*exp(x[0]*log(smax/smin));
    JacS = smin*log(smax/smin)*exp(x[0]*log(smax/smin));
    if ( s < smin ) return 0.0l;

    cc   = -1.0l + 2.0l*x[1];
    JacCos = 2.0l;

    func = JacS*JacCos;
    func *= params->temp*SQkaller(s,m3*m3,m4*m4)/s;
    func *= SQkaller(s,m1*m1,m2*m2)/sqrt(s)*boost::math::cyl_bessel_k(1,sqrt(s)/params->temp);
    func *= S12*S34 *2.0l*M_PI*params->SqAmp(s,cc,params->iproc);

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
