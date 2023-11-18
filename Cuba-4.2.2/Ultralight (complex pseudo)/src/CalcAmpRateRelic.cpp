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
    	num=-9*pow(g,2)*pow(lambdaH*pow(vh,2) + lambda12*pow(vsigma,2) +
    		     lambdaphi*pow(vsigma,2) + pow(pow(lambdaH,2)*pow(vh,4) -
    		       2*(lambdaH*(lambda12 + lambdaphi) - 4*pow(lambdaHphi,2))*pow(vh,2)*
    		        pow(vsigma,2) + pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5),2)*
    		   (-256*pow(lambdaHphi,2)*pow(mius,4)*pow(vh,2) -
    		     16*pow(lambdaHphi,2)*pow(s,2)*pow(vh,2) +
    		     64*lambdaH*pow(lambdaHphi,2)*pow(mius,2)*pow(vh,4) -
    		     6*pow(lambdaH,2)*pow(lambdaHphi,2)*pow(vh,6) -
    		     2*cos(4*a)*pow(lambdaH,2)*pow(lambdaHphi,2)*pow(vh,6) +
    		     64*lambda12*pow(lambdaHphi,2)*pow(mius,2)*pow(vh,2)*pow(vsigma,2) +
    		     64*lambdaphi*pow(lambdaHphi,2)*pow(mius,2)*pow(vh,2)*pow(vsigma,2) -
    		     2*lambda12*lambdaphi*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,2) +
    		     2*lambda12*lambdaphi*cos(4*a)*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,2) -
    		     pow(lambda12,2)*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,2) +
    		     cos(4*a)*pow(lambda12,2)*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,2) -
    		     4*lambda12*lambdaH*pow(lambdaHphi,2)*pow(vh,4)*pow(vsigma,2) -
    		     4*lambdaH*lambdaphi*pow(lambdaHphi,2)*pow(vh,4)*pow(vsigma,2) +
    		     4*lambda12*lambdaH*cos(4*a)*pow(lambdaHphi,2)*pow(vh,4)*pow(vsigma,2) +
    		     4*lambdaH*lambdaphi*cos(4*a)*pow(lambdaHphi,2)*pow(vh,4)*pow(vsigma,2) -
    		     16*pow(lambdaHphi,4)*pow(vh,4)*pow(vsigma,2) -
    		     16*cos(4*a)*pow(lambdaHphi,4)*pow(vh,4)*pow(vsigma,2) -
    		     pow(lambdaH,2)*pow(lambdaphi,2)*pow(vh,4)*pow(vsigma,2) +
    		     cos(4*a)*pow(lambdaH,2)*pow(lambdaphi,2)*pow(vh,4)*pow(vsigma,2) +
    		     6*lambdaH*lambdaphi*pow(lambda12,2)*pow(vh,2)*pow(vsigma,4) -
    		     6*lambdaH*lambdaphi*cos(4*a)*pow(lambda12,2)*pow(vh,2)*pow(vsigma,4) +
    		     2*lambdaH*pow(lambda12,3)*pow(vh,2)*pow(vsigma,4) -
    		     2*lambdaH*cos(4*a)*pow(lambda12,3)*pow(vh,2)*pow(vsigma,4) -
    		     28*lambda12*lambdaphi*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,4) +
    		     12*lambda12*lambdaphi*cos(4*a)*pow(lambdaHphi,2)*pow(vh,2)*
    		      pow(vsigma,4) - 14*pow(lambda12,2)*pow(lambdaHphi,2)*pow(vh,2)*
    		      pow(vsigma,4) + 6*cos(4*a)*pow(lambda12,2)*pow(lambdaHphi,2)*pow(vh,2)*
    		      pow(vsigma,4) + 6*lambda12*lambdaH*pow(lambdaphi,2)*pow(vh,2)*
    		      pow(vsigma,4) - 6*lambda12*lambdaH*cos(4*a)*pow(lambdaphi,2)*pow(vh,2)*
    		      pow(vsigma,4) - 14*pow(lambdaHphi,2)*pow(lambdaphi,2)*pow(vh,2)*
    		      pow(vsigma,4) + 6*cos(4*a)*pow(lambdaHphi,2)*pow(lambdaphi,2)*pow(vh,2)*
    		      pow(vsigma,4) + 2*lambdaH*pow(lambdaphi,3)*pow(vh,2)*pow(vsigma,4) -
    		     2*lambdaH*cos(4*a)*pow(lambdaphi,3)*pow(vh,2)*pow(vsigma,4) -
    		     4*lambdaphi*pow(lambda12,3)*pow(vsigma,6) +
    		     4*lambdaphi*cos(4*a)*pow(lambda12,3)*pow(vsigma,6) -
    		     pow(lambda12,4)*pow(vsigma,6) + cos(4*a)*pow(lambda12,4)*pow(vsigma,6) -
    		     6*pow(lambda12,2)*pow(lambdaphi,2)*pow(vsigma,6) +
    		     6*cos(4*a)*pow(lambda12,2)*pow(lambdaphi,2)*pow(vsigma,6) -
    		     4*lambda12*pow(lambdaphi,3)*pow(vsigma,6) +
    		     4*lambda12*cos(4*a)*pow(lambdaphi,3)*pow(vsigma,6) -
    		     pow(lambdaphi,4)*pow(vsigma,6) +
    		     cos(4*a)*pow(lambdaphi,4)*pow(vsigma,6) +
    		     8*cos(2*a)*pow(lambdaHphi,2)*pow(vh,2)*
    		      (-8*pow(mius,2) + lambdaH*pow(vh,2) +
    		        (lambda12 + lambdaphi)*pow(vsigma,2))*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
    		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5) -
    		     32*lambda12*lambdaHphi*vh*vsigma*pow(2,0.5)*pow(mius,2)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
    		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) -
    		     32*lambdaHphi*lambdaphi*vh*vsigma*pow(2,0.5)*pow(mius,2)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
    		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) +
    		     4*lambda12*lambdaH*lambdaHphi*vsigma*pow(2,0.5)*pow(vh,3)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
    		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) +
    		     4*lambdaH*lambdaHphi*lambdaphi*vsigma*pow(2,0.5)*pow(vh,3)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
    		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) +
    		     8*lambda12*lambdaHphi*lambdaphi*vh*pow(2,0.5)*pow(vsigma,3)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
    		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) +
    		     4*lambdaHphi*vh*pow(2,0.5)*pow(lambda12,2)*pow(vsigma,3)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
    		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) +
    		     4*lambdaHphi*vh*pow(2,0.5)*pow(lambdaphi,2)*pow(vsigma,3)*
    		      pow(pow(lambdaH,2)*pow(vh,4) -
    		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
    		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) -
    		     8*lambdaHphi*s*vh*(-2*lambdaHphi*vh*
    		         (-8*pow(mius,2) + lambdaH*pow(vh,2) +
    		           (lambda12 + lambdaphi)*pow(vsigma,2)) +
    		        2*lambdaHphi*vh*cos(2*a)*
    		         pow(pow(lambdaH,2)*pow(vh,4) -
    		           2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
    		           8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		           pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5) +
    		        (lambda12 + lambdaphi)*vsigma*pow(2,0.5)*
    		         pow(pow(lambdaH,2)*pow(vh,4) -
    		           2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
    		           8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
    		           pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a)) -
    		     2*lambda12*lambdaHphi*vsigma*pow(2,0.5)*pow(lambdaH,2)*pow(vh,5)*
    		      sin(4*a) - 2*lambdaHphi*lambdaphi*vsigma*pow(2,0.5)*pow(lambdaH,2)*
    		      pow(vh,5)*sin(4*a) + 8*lambda12*lambdaH*lambdaHphi*lambdaphi*pow(2,0.5)*
    		      pow(vh,3)*pow(vsigma,3)*sin(4*a) +
    		     4*lambdaH*lambdaHphi*pow(2,0.5)*pow(lambda12,2)*pow(vh,3)*pow(vsigma,3)*
    		      sin(4*a) - 16*lambda12*pow(2,0.5)*pow(lambdaHphi,3)*pow(vh,3)*
    		      pow(vsigma,3)*sin(4*a) - 16*lambdaphi*pow(2,0.5)*pow(lambdaHphi,3)*
    		      pow(vh,3)*pow(vsigma,3)*sin(4*a) +
    		     4*lambdaH*lambdaHphi*pow(2,0.5)*pow(lambdaphi,2)*pow(vh,3)*pow(vsigma,3)*
    		      sin(4*a) - 6*lambdaHphi*lambdaphi*vh*pow(2,0.5)*pow(lambda12,2)*
    		      pow(vsigma,5)*sin(4*a) - 2*lambdaHphi*vh*pow(2,0.5)*pow(lambda12,3)*
    		      pow(vsigma,5)*sin(4*a) - 6*lambda12*lambdaHphi*vh*pow(2,0.5)*
    		      pow(lambdaphi,2)*pow(vsigma,5)*sin(4*a) -
    		     2*lambdaHphi*vh*pow(2,0.5)*pow(lambdaphi,3)*pow(vsigma,5)*sin(4*a));
    	den=512*pow(mW,2)*((pow(wh,2)*(lambdaH*pow(vh,2) +
    	          (lambda12 + lambdaphi)*pow(vsigma,2) +
    	          pow(pow(lambdaH,2)*pow(vh,4) -
    	            2*(lambdaH*(lambda12 + lambdaphi) - 4*pow(lambdaHphi,2))*pow(vh,2)*
    	             pow(vsigma,2) + pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)))/
    	      2. + pow(s + (-(lambdaH*pow(vh,2)) -
    	          (lambda12 + lambdaphi)*pow(vsigma,2) -
    	          pow(pow(lambdaH,2)*pow(vh,4) -
    	            2*(lambdaH*(lambda12 + lambdaphi) - 4*pow(lambdaHphi,2))*pow(vh,2)*
    	             pow(vsigma,2) + pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5))/2.
    	       ,2))*((pow(wh2,2)*(lambdaH*pow(vh,2) +
    	          (lambda12 + lambdaphi)*pow(vsigma,2) -
    	          pow(pow(lambdaH,2)*pow(vh,4) -
    	            2*(lambdaH*(lambda12 + lambdaphi) - 4*pow(lambdaHphi,2))*pow(vh,2)*
    	             pow(vsigma,2) + pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)))/
    	      2. + pow(s + (-(lambdaH*pow(vh,2)) -
    	          (lambda12 + lambdaphi)*pow(vsigma,2) +
    	          pow(pow(lambdaH,2)*pow(vh,4) -
    	            2*(lambdaH*(lambda12 + lambdaphi) - 4*pow(lambdaHphi,2))*pow(vh,2)*
    	             pow(vsigma,2) + pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5))/2.
    	       ,2));
    	amp=num/den;
    }
    else if (id == 2)
    {
      	num=-3*pow(g,2)*pow(mqc,2)*(-0.5*s + 2*pow(mqc,2) - 2*pow(mius,2))*
      		   (256*pow(lambdaHphi,2)*pow(mius,4)*pow(vh,2) +
      		     16*pow(lambdaHphi,2)*pow(s,2)*pow(vh,2) -
      		     64*lambdaH*pow(lambdaHphi,2)*pow(mius,2)*pow(vh,4) +
      		     6*pow(lambdaH,2)*pow(lambdaHphi,2)*pow(vh,6) +
      		     2*cos(4*a)*pow(lambdaH,2)*pow(lambdaHphi,2)*pow(vh,6) -
      		     64*lambda12*pow(lambdaHphi,2)*pow(mius,2)*pow(vh,2)*pow(vsigma,2) -
      		     64*lambdaphi*pow(lambdaHphi,2)*pow(mius,2)*pow(vh,2)*pow(vsigma,2) +
      		     2*lambda12*lambdaphi*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,2) -
      		     2*lambda12*lambdaphi*cos(4*a)*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,2) +
      		     pow(lambda12,2)*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,2) -
      		     cos(4*a)*pow(lambda12,2)*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,2) +
      		     4*lambda12*lambdaH*pow(lambdaHphi,2)*pow(vh,4)*pow(vsigma,2) +
      		     4*lambdaH*lambdaphi*pow(lambdaHphi,2)*pow(vh,4)*pow(vsigma,2) -
      		     4*lambda12*lambdaH*cos(4*a)*pow(lambdaHphi,2)*pow(vh,4)*pow(vsigma,2) -
      		     4*lambdaH*lambdaphi*cos(4*a)*pow(lambdaHphi,2)*pow(vh,4)*
      		      pow(vsigma,2) + 16*pow(lambdaHphi,4)*pow(vh,4)*pow(vsigma,2) +
      		     16*cos(4*a)*pow(lambdaHphi,4)*pow(vh,4)*pow(vsigma,2) +
      		     pow(lambdaH,2)*pow(lambdaphi,2)*pow(vh,4)*pow(vsigma,2) -
      		     cos(4*a)*pow(lambdaH,2)*pow(lambdaphi,2)*pow(vh,4)*pow(vsigma,2) -
      		     6*lambdaH*lambdaphi*pow(lambda12,2)*pow(vh,2)*pow(vsigma,4) +
      		     6*lambdaH*lambdaphi*cos(4*a)*pow(lambda12,2)*pow(vh,2)*pow(vsigma,4) -
      		     2*lambdaH*pow(lambda12,3)*pow(vh,2)*pow(vsigma,4) +
      		     2*lambdaH*cos(4*a)*pow(lambda12,3)*pow(vh,2)*pow(vsigma,4) +
      		     28*lambda12*lambdaphi*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,4) -
      		     12*lambda12*lambdaphi*cos(4*a)*pow(lambdaHphi,2)*pow(vh,2)*
      		      pow(vsigma,4) + 14*pow(lambda12,2)*pow(lambdaHphi,2)*pow(vh,2)*
      		      pow(vsigma,4) - 6*cos(4*a)*pow(lambda12,2)*pow(lambdaHphi,2)*pow(vh,2)*
      		      pow(vsigma,4) - 6*lambda12*lambdaH*pow(lambdaphi,2)*pow(vh,2)*
      		      pow(vsigma,4) + 6*lambda12*lambdaH*cos(4*a)*pow(lambdaphi,2)*pow(vh,2)*
      		      pow(vsigma,4) + 14*pow(lambdaHphi,2)*pow(lambdaphi,2)*pow(vh,2)*
      		      pow(vsigma,4) - 6*cos(4*a)*pow(lambdaHphi,2)*pow(lambdaphi,2)*
      		      pow(vh,2)*pow(vsigma,4) -
      		     2*lambdaH*pow(lambdaphi,3)*pow(vh,2)*pow(vsigma,4) +
      		     2*lambdaH*cos(4*a)*pow(lambdaphi,3)*pow(vh,2)*pow(vsigma,4) +
      		     4*lambdaphi*pow(lambda12,3)*pow(vsigma,6) -
      		     4*lambdaphi*cos(4*a)*pow(lambda12,3)*pow(vsigma,6) +
      		     pow(lambda12,4)*pow(vsigma,6) -
      		     cos(4*a)*pow(lambda12,4)*pow(vsigma,6) +
      		     6*pow(lambda12,2)*pow(lambdaphi,2)*pow(vsigma,6) -
      		     6*cos(4*a)*pow(lambda12,2)*pow(lambdaphi,2)*pow(vsigma,6) +
      		     4*lambda12*pow(lambdaphi,3)*pow(vsigma,6) -
      		     4*lambda12*cos(4*a)*pow(lambdaphi,3)*pow(vsigma,6) +
      		     pow(lambdaphi,4)*pow(vsigma,6) -
      		     cos(4*a)*pow(lambdaphi,4)*pow(vsigma,6) -
      		     8*cos(2*a)*pow(lambdaHphi,2)*pow(vh,2)*
      		      (-8*pow(mius,2) + lambdaH*pow(vh,2) +
      		        (lambda12 + lambdaphi)*pow(vsigma,2))*
      		      pow(pow(lambdaH,2)*pow(vh,4) -
      		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
      		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
      		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5) +
      		     32*lambda12*lambdaHphi*vh*vsigma*pow(2,0.5)*pow(mius,2)*
      		      pow(pow(lambdaH,2)*pow(vh,4) -
      		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
      		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
      		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) +
      		     32*lambdaHphi*lambdaphi*vh*vsigma*pow(2,0.5)*pow(mius,2)*
      		      pow(pow(lambdaH,2)*pow(vh,4) -
      		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
      		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
      		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) -
      		     4*lambda12*lambdaH*lambdaHphi*vsigma*pow(2,0.5)*pow(vh,3)*
      		      pow(pow(lambdaH,2)*pow(vh,4) -
      		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
      		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
      		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) -
      		     4*lambdaH*lambdaHphi*lambdaphi*vsigma*pow(2,0.5)*pow(vh,3)*
      		      pow(pow(lambdaH,2)*pow(vh,4) -
      		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
      		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
      		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) -
      		     8*lambda12*lambdaHphi*lambdaphi*vh*pow(2,0.5)*pow(vsigma,3)*
      		      pow(pow(lambdaH,2)*pow(vh,4) -
      		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
      		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
      		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) -
      		     4*lambdaHphi*vh*pow(2,0.5)*pow(lambda12,2)*pow(vsigma,3)*
      		      pow(pow(lambdaH,2)*pow(vh,4) -
      		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
      		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
      		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) -
      		     4*lambdaHphi*vh*pow(2,0.5)*pow(lambdaphi,2)*pow(vsigma,3)*
      		      pow(pow(lambdaH,2)*pow(vh,4) -
      		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
      		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
      		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) +
      		     8*lambdaHphi*s*vh*(-2*lambdaHphi*vh*
      		         (-8*pow(mius,2) + lambdaH*pow(vh,2) +
      		           (lambda12 + lambdaphi)*pow(vsigma,2)) +
      		        2*lambdaHphi*vh*cos(2*a)*
      		         pow(pow(lambdaH,2)*pow(vh,4) -
      		           2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
      		           8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
      		           pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5) +
      		        (lambda12 + lambdaphi)*vsigma*pow(2,0.5)*
      		         pow(pow(lambdaH,2)*pow(vh,4) -
      		           2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
      		           8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
      		           pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a)) +
      		     2*lambda12*lambdaHphi*vsigma*pow(2,0.5)*pow(lambdaH,2)*pow(vh,5)*
      		      sin(4*a) + 2*lambdaHphi*lambdaphi*vsigma*pow(2,0.5)*pow(lambdaH,2)*
      		      pow(vh,5)*sin(4*a) - 8*lambda12*lambdaH*lambdaHphi*lambdaphi*
      		      pow(2,0.5)*pow(vh,3)*pow(vsigma,3)*sin(4*a) -
      		     4*lambdaH*lambdaHphi*pow(2,0.5)*pow(lambda12,2)*pow(vh,3)*pow(vsigma,3)*
      		      sin(4*a) + 16*lambda12*pow(2,0.5)*pow(lambdaHphi,3)*pow(vh,3)*
      		      pow(vsigma,3)*sin(4*a) +
      		     16*lambdaphi*pow(2,0.5)*pow(lambdaHphi,3)*pow(vh,3)*pow(vsigma,3)*
      		      sin(4*a) - 4*lambdaH*lambdaHphi*pow(2,0.5)*pow(lambdaphi,2)*pow(vh,3)*
      		      pow(vsigma,3)*sin(4*a) +
      		     6*lambdaHphi*lambdaphi*vh*pow(2,0.5)*pow(lambda12,2)*pow(vsigma,5)*
      		      sin(4*a) + 2*lambdaHphi*vh*pow(2,0.5)*pow(lambda12,3)*pow(vsigma,5)*
      		      sin(4*a) + 6*lambda12*lambdaHphi*vh*pow(2,0.5)*pow(lambdaphi,2)*
      		      pow(vsigma,5)*sin(4*a) +
      		     2*lambdaHphi*vh*pow(2,0.5)*pow(lambdaphi,3)*pow(vsigma,5)*sin(4*a));
       	den=16*pow(mW,2)*((pow(wh,2)*(lambdaH*pow(vh,2) +
                (lambda12 + lambdaphi)*pow(vsigma,2) +
                pow(pow(lambdaH,2)*pow(vh,4) -
                  2*(lambdaH*(lambda12 + lambdaphi) - 4*pow(lambdaHphi,2))*
                   pow(vh,2)*pow(vsigma,2) +
                  pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)))/2. +
           pow(s + (-(lambdaH*pow(vh,2)) - (lambda12 + lambdaphi)*pow(vsigma,2) -
                pow(pow(lambdaH,2)*pow(vh,4) -
                  2*(lambdaH*(lambda12 + lambdaphi) - 4*pow(lambdaHphi,2))*
                   pow(vh,2)*pow(vsigma,2) +
                  pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5))/2.,2))*
         ((pow(wh2,2)*(lambdaH*pow(vh,2) + (lambda12 + lambdaphi)*pow(vsigma,2) -
                pow(pow(lambdaH,2)*pow(vh,4) -
                  2*(lambdaH*(lambda12 + lambdaphi) - 4*pow(lambdaHphi,2))*
                   pow(vh,2)*pow(vsigma,2) +
                  pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)))/2. +
           pow(s + (-(lambdaH*pow(vh,2)) - (lambda12 + lambdaphi)*pow(vsigma,2) +
                pow(pow(lambdaH,2)*pow(vh,4) -
                  2*(lambdaH*(lambda12 + lambdaphi) - 4*pow(lambdaHphi,2))*
                   pow(vh,2)*pow(vsigma,2) +
                  pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5))/2.,2));
       	amp=num/den;
    }
    else if (id == 3)
    {
       	num=-3*pow(g,2)*pow(mqt,2)*(-0.5*s - 2*pow(mius,2) + 2*pow(mqt,2))*
       		   (256*pow(lambdaHphi,2)*pow(mius,4)*pow(vh,2) +
       		     16*pow(lambdaHphi,2)*pow(s,2)*pow(vh,2) -
       		     64*lambdaH*pow(lambdaHphi,2)*pow(mius,2)*pow(vh,4) +
       		     6*pow(lambdaH,2)*pow(lambdaHphi,2)*pow(vh,6) +
       		     2*cos(4*a)*pow(lambdaH,2)*pow(lambdaHphi,2)*pow(vh,6) -
       		     64*lambda12*pow(lambdaHphi,2)*pow(mius,2)*pow(vh,2)*pow(vsigma,2) -
       		     64*lambdaphi*pow(lambdaHphi,2)*pow(mius,2)*pow(vh,2)*pow(vsigma,2) +
       		     2*lambda12*lambdaphi*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,2) -
       		     2*lambda12*lambdaphi*cos(4*a)*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,2) +
       		     pow(lambda12,2)*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,2) -
       		     cos(4*a)*pow(lambda12,2)*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,2) +
       		     4*lambda12*lambdaH*pow(lambdaHphi,2)*pow(vh,4)*pow(vsigma,2) +
       		     4*lambdaH*lambdaphi*pow(lambdaHphi,2)*pow(vh,4)*pow(vsigma,2) -
       		     4*lambda12*lambdaH*cos(4*a)*pow(lambdaHphi,2)*pow(vh,4)*pow(vsigma,2) -
       		     4*lambdaH*lambdaphi*cos(4*a)*pow(lambdaHphi,2)*pow(vh,4)*
       		      pow(vsigma,2) + 16*pow(lambdaHphi,4)*pow(vh,4)*pow(vsigma,2) +
       		     16*cos(4*a)*pow(lambdaHphi,4)*pow(vh,4)*pow(vsigma,2) +
       		     pow(lambdaH,2)*pow(lambdaphi,2)*pow(vh,4)*pow(vsigma,2) -
       		     cos(4*a)*pow(lambdaH,2)*pow(lambdaphi,2)*pow(vh,4)*pow(vsigma,2) -
       		     6*lambdaH*lambdaphi*pow(lambda12,2)*pow(vh,2)*pow(vsigma,4) +
       		     6*lambdaH*lambdaphi*cos(4*a)*pow(lambda12,2)*pow(vh,2)*pow(vsigma,4) -
       		     2*lambdaH*pow(lambda12,3)*pow(vh,2)*pow(vsigma,4) +
       		     2*lambdaH*cos(4*a)*pow(lambda12,3)*pow(vh,2)*pow(vsigma,4) +
       		     28*lambda12*lambdaphi*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,4) -
       		     12*lambda12*lambdaphi*cos(4*a)*pow(lambdaHphi,2)*pow(vh,2)*
       		      pow(vsigma,4) + 14*pow(lambda12,2)*pow(lambdaHphi,2)*pow(vh,2)*
       		      pow(vsigma,4) - 6*cos(4*a)*pow(lambda12,2)*pow(lambdaHphi,2)*pow(vh,2)*
       		      pow(vsigma,4) - 6*lambda12*lambdaH*pow(lambdaphi,2)*pow(vh,2)*
       		      pow(vsigma,4) + 6*lambda12*lambdaH*cos(4*a)*pow(lambdaphi,2)*pow(vh,2)*
       		      pow(vsigma,4) + 14*pow(lambdaHphi,2)*pow(lambdaphi,2)*pow(vh,2)*
       		      pow(vsigma,4) - 6*cos(4*a)*pow(lambdaHphi,2)*pow(lambdaphi,2)*
       		      pow(vh,2)*pow(vsigma,4) -
       		     2*lambdaH*pow(lambdaphi,3)*pow(vh,2)*pow(vsigma,4) +
       		     2*lambdaH*cos(4*a)*pow(lambdaphi,3)*pow(vh,2)*pow(vsigma,4) +
       		     4*lambdaphi*pow(lambda12,3)*pow(vsigma,6) -
       		     4*lambdaphi*cos(4*a)*pow(lambda12,3)*pow(vsigma,6) +
       		     pow(lambda12,4)*pow(vsigma,6) -
       		     cos(4*a)*pow(lambda12,4)*pow(vsigma,6) +
       		     6*pow(lambda12,2)*pow(lambdaphi,2)*pow(vsigma,6) -
       		     6*cos(4*a)*pow(lambda12,2)*pow(lambdaphi,2)*pow(vsigma,6) +
       		     4*lambda12*pow(lambdaphi,3)*pow(vsigma,6) -
       		     4*lambda12*cos(4*a)*pow(lambdaphi,3)*pow(vsigma,6) +
       		     pow(lambdaphi,4)*pow(vsigma,6) -
       		     cos(4*a)*pow(lambdaphi,4)*pow(vsigma,6) -
       		     8*cos(2*a)*pow(lambdaHphi,2)*pow(vh,2)*
       		      (-8*pow(mius,2) + lambdaH*pow(vh,2) +
       		        (lambda12 + lambdaphi)*pow(vsigma,2))*
       		      pow(pow(lambdaH,2)*pow(vh,4) -
       		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
       		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
       		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5) +
       		     32*lambda12*lambdaHphi*vh*vsigma*pow(2,0.5)*pow(mius,2)*
       		      pow(pow(lambdaH,2)*pow(vh,4) -
       		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
       		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
       		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) +
       		     32*lambdaHphi*lambdaphi*vh*vsigma*pow(2,0.5)*pow(mius,2)*
       		      pow(pow(lambdaH,2)*pow(vh,4) -
       		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
       		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
       		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) -
       		     4*lambda12*lambdaH*lambdaHphi*vsigma*pow(2,0.5)*pow(vh,3)*
       		      pow(pow(lambdaH,2)*pow(vh,4) -
       		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
       		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
       		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) -
       		     4*lambdaH*lambdaHphi*lambdaphi*vsigma*pow(2,0.5)*pow(vh,3)*
       		      pow(pow(lambdaH,2)*pow(vh,4) -
       		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
       		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
       		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) -
       		     8*lambda12*lambdaHphi*lambdaphi*vh*pow(2,0.5)*pow(vsigma,3)*
       		      pow(pow(lambdaH,2)*pow(vh,4) -
       		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
       		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
       		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) -
       		     4*lambdaHphi*vh*pow(2,0.5)*pow(lambda12,2)*pow(vsigma,3)*
       		      pow(pow(lambdaH,2)*pow(vh,4) -
       		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
       		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
       		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) -
       		     4*lambdaHphi*vh*pow(2,0.5)*pow(lambdaphi,2)*pow(vsigma,3)*
       		      pow(pow(lambdaH,2)*pow(vh,4) -
       		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
       		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
       		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) +
       		     8*lambdaHphi*s*vh*(-2*lambdaHphi*vh*
       		         (-8*pow(mius,2) + lambdaH*pow(vh,2) +
       		           (lambda12 + lambdaphi)*pow(vsigma,2)) +
       		        2*lambdaHphi*vh*cos(2*a)*
       		         pow(pow(lambdaH,2)*pow(vh,4) -
       		           2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
       		           8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
       		           pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5) +
       		        (lambda12 + lambdaphi)*vsigma*pow(2,0.5)*
       		         pow(pow(lambdaH,2)*pow(vh,4) -
       		           2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
       		           8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
       		           pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a)) +
       		     2*lambda12*lambdaHphi*vsigma*pow(2,0.5)*pow(lambdaH,2)*pow(vh,5)*
       		      sin(4*a) + 2*lambdaHphi*lambdaphi*vsigma*pow(2,0.5)*pow(lambdaH,2)*
       		      pow(vh,5)*sin(4*a) - 8*lambda12*lambdaH*lambdaHphi*lambdaphi*
       		      pow(2,0.5)*pow(vh,3)*pow(vsigma,3)*sin(4*a) -
       		     4*lambdaH*lambdaHphi*pow(2,0.5)*pow(lambda12,2)*pow(vh,3)*pow(vsigma,3)*
       		      sin(4*a) + 16*lambda12*pow(2,0.5)*pow(lambdaHphi,3)*pow(vh,3)*
       		      pow(vsigma,3)*sin(4*a) +
       		     16*lambdaphi*pow(2,0.5)*pow(lambdaHphi,3)*pow(vh,3)*pow(vsigma,3)*
       		      sin(4*a) - 4*lambdaH*lambdaHphi*pow(2,0.5)*pow(lambdaphi,2)*pow(vh,3)*
       		      pow(vsigma,3)*sin(4*a) +
       		     6*lambdaHphi*lambdaphi*vh*pow(2,0.5)*pow(lambda12,2)*pow(vsigma,5)*
       		      sin(4*a) + 2*lambdaHphi*vh*pow(2,0.5)*pow(lambda12,3)*pow(vsigma,5)*
       		      sin(4*a) + 6*lambda12*lambdaHphi*vh*pow(2,0.5)*pow(lambdaphi,2)*
       		      pow(vsigma,5)*sin(4*a) +
       		     2*lambdaHphi*vh*pow(2,0.5)*pow(lambdaphi,3)*pow(vsigma,5)*sin(4*a));
           	den=16*pow(mW,2)*((pow(wh,2)*(lambdaH*pow(vh,2) +
                    (lambda12 + lambdaphi)*pow(vsigma,2) +
                    pow(pow(lambdaH,2)*pow(vh,4) -
                      2*(lambdaH*(lambda12 + lambdaphi) - 4*pow(lambdaHphi,2))*
                       pow(vh,2)*pow(vsigma,2) +
                      pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)))/2. +
               pow(s + (-(lambdaH*pow(vh,2)) - (lambda12 + lambdaphi)*pow(vsigma,2) -
                    pow(pow(lambdaH,2)*pow(vh,4) -
                      2*(lambdaH*(lambda12 + lambdaphi) - 4*pow(lambdaHphi,2))*
                       pow(vh,2)*pow(vsigma,2) +
                      pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5))/2.,2))*
             ((pow(wh2,2)*(lambdaH*pow(vh,2) + (lambda12 + lambdaphi)*pow(vsigma,2) -
                    pow(pow(lambdaH,2)*pow(vh,4) -
                      2*(lambdaH*(lambda12 + lambdaphi) - 4*pow(lambdaHphi,2))*
                       pow(vh,2)*pow(vsigma,2) +
                      pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)))/2. +
               pow(s + (-(lambdaH*pow(vh,2)) - (lambda12 + lambdaphi)*pow(vsigma,2) +
                    pow(pow(lambdaH,2)*pow(vh,4) -
                      2*(lambdaH*(lambda12 + lambdaphi) - 4*pow(lambdaHphi,2))*
                       pow(vh,2)*pow(vsigma,2) +
                      pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5))/2.,2));
           	amp=num/den;
        }
    else if (id == 4)
        {
          	num=-3*pow(g,2)*pow(mqs,2)*(-0.5*s - 2*pow(mius,2) + 2*pow(mqs,2))*
          		   (256*pow(lambdaHphi,2)*pow(mius,4)*pow(vh,2) +
          		     16*pow(lambdaHphi,2)*pow(s,2)*pow(vh,2) -
          		     64*lambdaH*pow(lambdaHphi,2)*pow(mius,2)*pow(vh,4) +
          		     6*pow(lambdaH,2)*pow(lambdaHphi,2)*pow(vh,6) +
          		     2*cos(4*a)*pow(lambdaH,2)*pow(lambdaHphi,2)*pow(vh,6) -
          		     64*lambda12*pow(lambdaHphi,2)*pow(mius,2)*pow(vh,2)*pow(vsigma,2) -
          		     64*lambdaphi*pow(lambdaHphi,2)*pow(mius,2)*pow(vh,2)*pow(vsigma,2) +
          		     2*lambda12*lambdaphi*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,2) -
          		     2*lambda12*lambdaphi*cos(4*a)*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,2) +
          		     pow(lambda12,2)*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,2) -
          		     cos(4*a)*pow(lambda12,2)*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,2) +
          		     4*lambda12*lambdaH*pow(lambdaHphi,2)*pow(vh,4)*pow(vsigma,2) +
          		     4*lambdaH*lambdaphi*pow(lambdaHphi,2)*pow(vh,4)*pow(vsigma,2) -
          		     4*lambda12*lambdaH*cos(4*a)*pow(lambdaHphi,2)*pow(vh,4)*pow(vsigma,2) -
          		     4*lambdaH*lambdaphi*cos(4*a)*pow(lambdaHphi,2)*pow(vh,4)*
          		      pow(vsigma,2) + 16*pow(lambdaHphi,4)*pow(vh,4)*pow(vsigma,2) +
          		     16*cos(4*a)*pow(lambdaHphi,4)*pow(vh,4)*pow(vsigma,2) +
          		     pow(lambdaH,2)*pow(lambdaphi,2)*pow(vh,4)*pow(vsigma,2) -
          		     cos(4*a)*pow(lambdaH,2)*pow(lambdaphi,2)*pow(vh,4)*pow(vsigma,2) -
          		     6*lambdaH*lambdaphi*pow(lambda12,2)*pow(vh,2)*pow(vsigma,4) +
          		     6*lambdaH*lambdaphi*cos(4*a)*pow(lambda12,2)*pow(vh,2)*pow(vsigma,4) -
          		     2*lambdaH*pow(lambda12,3)*pow(vh,2)*pow(vsigma,4) +
          		     2*lambdaH*cos(4*a)*pow(lambda12,3)*pow(vh,2)*pow(vsigma,4) +
          		     28*lambda12*lambdaphi*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,4) -
          		     12*lambda12*lambdaphi*cos(4*a)*pow(lambdaHphi,2)*pow(vh,2)*
          		      pow(vsigma,4) + 14*pow(lambda12,2)*pow(lambdaHphi,2)*pow(vh,2)*
          		      pow(vsigma,4) - 6*cos(4*a)*pow(lambda12,2)*pow(lambdaHphi,2)*pow(vh,2)*
          		      pow(vsigma,4) - 6*lambda12*lambdaH*pow(lambdaphi,2)*pow(vh,2)*
          		      pow(vsigma,4) + 6*lambda12*lambdaH*cos(4*a)*pow(lambdaphi,2)*pow(vh,2)*
          		      pow(vsigma,4) + 14*pow(lambdaHphi,2)*pow(lambdaphi,2)*pow(vh,2)*
          		      pow(vsigma,4) - 6*cos(4*a)*pow(lambdaHphi,2)*pow(lambdaphi,2)*
          		      pow(vh,2)*pow(vsigma,4) -
          		     2*lambdaH*pow(lambdaphi,3)*pow(vh,2)*pow(vsigma,4) +
          		     2*lambdaH*cos(4*a)*pow(lambdaphi,3)*pow(vh,2)*pow(vsigma,4) +
          		     4*lambdaphi*pow(lambda12,3)*pow(vsigma,6) -
          		     4*lambdaphi*cos(4*a)*pow(lambda12,3)*pow(vsigma,6) +
          		     pow(lambda12,4)*pow(vsigma,6) -
          		     cos(4*a)*pow(lambda12,4)*pow(vsigma,6) +
          		     6*pow(lambda12,2)*pow(lambdaphi,2)*pow(vsigma,6) -
          		     6*cos(4*a)*pow(lambda12,2)*pow(lambdaphi,2)*pow(vsigma,6) +
          		     4*lambda12*pow(lambdaphi,3)*pow(vsigma,6) -
          		     4*lambda12*cos(4*a)*pow(lambdaphi,3)*pow(vsigma,6) +
          		     pow(lambdaphi,4)*pow(vsigma,6) -
          		     cos(4*a)*pow(lambdaphi,4)*pow(vsigma,6) -
          		     8*cos(2*a)*pow(lambdaHphi,2)*pow(vh,2)*
          		      (-8*pow(mius,2) + lambdaH*pow(vh,2) +
          		        (lambda12 + lambdaphi)*pow(vsigma,2))*
          		      pow(pow(lambdaH,2)*pow(vh,4) -
          		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5) +
          		     32*lambda12*lambdaHphi*vh*vsigma*pow(2,0.5)*pow(mius,2)*
          		      pow(pow(lambdaH,2)*pow(vh,4) -
          		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) +
          		     32*lambdaHphi*lambdaphi*vh*vsigma*pow(2,0.5)*pow(mius,2)*
          		      pow(pow(lambdaH,2)*pow(vh,4) -
          		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) -
          		     4*lambda12*lambdaH*lambdaHphi*vsigma*pow(2,0.5)*pow(vh,3)*
          		      pow(pow(lambdaH,2)*pow(vh,4) -
          		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) -
          		     4*lambdaH*lambdaHphi*lambdaphi*vsigma*pow(2,0.5)*pow(vh,3)*
          		      pow(pow(lambdaH,2)*pow(vh,4) -
          		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) -
          		     8*lambda12*lambdaHphi*lambdaphi*vh*pow(2,0.5)*pow(vsigma,3)*
          		      pow(pow(lambdaH,2)*pow(vh,4) -
          		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) -
          		     4*lambdaHphi*vh*pow(2,0.5)*pow(lambda12,2)*pow(vsigma,3)*
          		      pow(pow(lambdaH,2)*pow(vh,4) -
          		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) -
          		     4*lambdaHphi*vh*pow(2,0.5)*pow(lambdaphi,2)*pow(vsigma,3)*
          		      pow(pow(lambdaH,2)*pow(vh,4) -
          		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) +
          		     8*lambdaHphi*s*vh*(-2*lambdaHphi*vh*
          		         (-8*pow(mius,2) + lambdaH*pow(vh,2) +
          		           (lambda12 + lambdaphi)*pow(vsigma,2)) +
          		        2*lambdaHphi*vh*cos(2*a)*
          		         pow(pow(lambdaH,2)*pow(vh,4) -
          		           2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		           8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		           pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5) +
          		        (lambda12 + lambdaphi)*vsigma*pow(2,0.5)*
          		         pow(pow(lambdaH,2)*pow(vh,4) -
          		           2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		           8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		           pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a)) +
          		     2*lambda12*lambdaHphi*vsigma*pow(2,0.5)*pow(lambdaH,2)*pow(vh,5)*
          		      sin(4*a) + 2*lambdaHphi*lambdaphi*vsigma*pow(2,0.5)*pow(lambdaH,2)*
          		      pow(vh,5)*sin(4*a) - 8*lambda12*lambdaH*lambdaHphi*lambdaphi*
          		      pow(2,0.5)*pow(vh,3)*pow(vsigma,3)*sin(4*a) -
          		     4*lambdaH*lambdaHphi*pow(2,0.5)*pow(lambda12,2)*pow(vh,3)*pow(vsigma,3)*
          		      sin(4*a) + 16*lambda12*pow(2,0.5)*pow(lambdaHphi,3)*pow(vh,3)*
          		      pow(vsigma,3)*sin(4*a) +
          		     16*lambdaphi*pow(2,0.5)*pow(lambdaHphi,3)*pow(vh,3)*pow(vsigma,3)*
          		      sin(4*a) - 4*lambdaH*lambdaHphi*pow(2,0.5)*pow(lambdaphi,2)*pow(vh,3)*
          		      pow(vsigma,3)*sin(4*a) +
          		     6*lambdaHphi*lambdaphi*vh*pow(2,0.5)*pow(lambda12,2)*pow(vsigma,5)*
          		      sin(4*a) + 2*lambdaHphi*vh*pow(2,0.5)*pow(lambda12,3)*pow(vsigma,5)*
          		      sin(4*a) + 6*lambda12*lambdaHphi*vh*pow(2,0.5)*pow(lambdaphi,2)*
          		      pow(vsigma,5)*sin(4*a) +
          		     2*lambdaHphi*vh*pow(2,0.5)*pow(lambdaphi,3)*pow(vsigma,5)*sin(4*a));
           	den=16*pow(mW,2)*((pow(wh,2)*(lambdaH*pow(vh,2) +
                    (lambda12 + lambdaphi)*pow(vsigma,2) +
                    pow(pow(lambdaH,2)*pow(vh,4) -
                      2*(lambdaH*(lambda12 + lambdaphi) - 4*pow(lambdaHphi,2))*
                       pow(vh,2)*pow(vsigma,2) +
                      pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)))/2. +
               pow(s + (-(lambdaH*pow(vh,2)) - (lambda12 + lambdaphi)*pow(vsigma,2) -
                    pow(pow(lambdaH,2)*pow(vh,4) -
                      2*(lambdaH*(lambda12 + lambdaphi) - 4*pow(lambdaHphi,2))*
                       pow(vh,2)*pow(vsigma,2) +
                      pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5))/2.,2))*
             ((pow(wh2,2)*(lambdaH*pow(vh,2) + (lambda12 + lambdaphi)*pow(vsigma,2) -
                    pow(pow(lambdaH,2)*pow(vh,4) -
                      2*(lambdaH*(lambda12 + lambdaphi) - 4*pow(lambdaHphi,2))*
                       pow(vh,2)*pow(vsigma,2) +
                      pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)))/2. +
               pow(s + (-(lambdaH*pow(vh,2)) - (lambda12 + lambdaphi)*pow(vsigma,2) +
                    pow(pow(lambdaH,2)*pow(vh,4) -
                      2*(lambdaH*(lambda12 + lambdaphi) - 4*pow(lambdaHphi,2))*
                       pow(vh,2)*pow(vsigma,2) +
                      pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5))/2.,2));
           	amp=num/den;
        }
    else if (id == 5)
        {
          	num=-3*pow(g,2)*pow(mqc,2)*(-0.5*s + 2*pow(mqc,2) - 2*pow(mius,2))*
          		   (256*pow(lambdaHphi,2)*pow(mius,4)*pow(vh,2) +
          		     16*pow(lambdaHphi,2)*pow(s,2)*pow(vh,2) -
          		     64*lambdaH*pow(lambdaHphi,2)*pow(mius,2)*pow(vh,4) +
          		     6*pow(lambdaH,2)*pow(lambdaHphi,2)*pow(vh,6) +
          		     2*cos(4*a)*pow(lambdaH,2)*pow(lambdaHphi,2)*pow(vh,6) -
          		     64*lambda12*pow(lambdaHphi,2)*pow(mius,2)*pow(vh,2)*pow(vsigma,2) -
          		     64*lambdaphi*pow(lambdaHphi,2)*pow(mius,2)*pow(vh,2)*pow(vsigma,2) +
          		     2*lambda12*lambdaphi*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,2) -
          		     2*lambda12*lambdaphi*cos(4*a)*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,2) +
          		     pow(lambda12,2)*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,2) -
          		     cos(4*a)*pow(lambda12,2)*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,2) +
          		     4*lambda12*lambdaH*pow(lambdaHphi,2)*pow(vh,4)*pow(vsigma,2) +
          		     4*lambdaH*lambdaphi*pow(lambdaHphi,2)*pow(vh,4)*pow(vsigma,2) -
          		     4*lambda12*lambdaH*cos(4*a)*pow(lambdaHphi,2)*pow(vh,4)*pow(vsigma,2) -
          		     4*lambdaH*lambdaphi*cos(4*a)*pow(lambdaHphi,2)*pow(vh,4)*
          		      pow(vsigma,2) + 16*pow(lambdaHphi,4)*pow(vh,4)*pow(vsigma,2) +
          		     16*cos(4*a)*pow(lambdaHphi,4)*pow(vh,4)*pow(vsigma,2) +
          		     pow(lambdaH,2)*pow(lambdaphi,2)*pow(vh,4)*pow(vsigma,2) -
          		     cos(4*a)*pow(lambdaH,2)*pow(lambdaphi,2)*pow(vh,4)*pow(vsigma,2) -
          		     6*lambdaH*lambdaphi*pow(lambda12,2)*pow(vh,2)*pow(vsigma,4) +
          		     6*lambdaH*lambdaphi*cos(4*a)*pow(lambda12,2)*pow(vh,2)*pow(vsigma,4) -
          		     2*lambdaH*pow(lambda12,3)*pow(vh,2)*pow(vsigma,4) +
          		     2*lambdaH*cos(4*a)*pow(lambda12,3)*pow(vh,2)*pow(vsigma,4) +
          		     28*lambda12*lambdaphi*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,4) -
          		     12*lambda12*lambdaphi*cos(4*a)*pow(lambdaHphi,2)*pow(vh,2)*
          		      pow(vsigma,4) + 14*pow(lambda12,2)*pow(lambdaHphi,2)*pow(vh,2)*
          		      pow(vsigma,4) - 6*cos(4*a)*pow(lambda12,2)*pow(lambdaHphi,2)*pow(vh,2)*
          		      pow(vsigma,4) - 6*lambda12*lambdaH*pow(lambdaphi,2)*pow(vh,2)*
          		      pow(vsigma,4) + 6*lambda12*lambdaH*cos(4*a)*pow(lambdaphi,2)*pow(vh,2)*
          		      pow(vsigma,4) + 14*pow(lambdaHphi,2)*pow(lambdaphi,2)*pow(vh,2)*
          		      pow(vsigma,4) - 6*cos(4*a)*pow(lambdaHphi,2)*pow(lambdaphi,2)*
          		      pow(vh,2)*pow(vsigma,4) -
          		     2*lambdaH*pow(lambdaphi,3)*pow(vh,2)*pow(vsigma,4) +
          		     2*lambdaH*cos(4*a)*pow(lambdaphi,3)*pow(vh,2)*pow(vsigma,4) +
          		     4*lambdaphi*pow(lambda12,3)*pow(vsigma,6) -
          		     4*lambdaphi*cos(4*a)*pow(lambda12,3)*pow(vsigma,6) +
          		     pow(lambda12,4)*pow(vsigma,6) -
          		     cos(4*a)*pow(lambda12,4)*pow(vsigma,6) +
          		     6*pow(lambda12,2)*pow(lambdaphi,2)*pow(vsigma,6) -
          		     6*cos(4*a)*pow(lambda12,2)*pow(lambdaphi,2)*pow(vsigma,6) +
          		     4*lambda12*pow(lambdaphi,3)*pow(vsigma,6) -
          		     4*lambda12*cos(4*a)*pow(lambdaphi,3)*pow(vsigma,6) +
          		     pow(lambdaphi,4)*pow(vsigma,6) -
          		     cos(4*a)*pow(lambdaphi,4)*pow(vsigma,6) -
          		     8*cos(2*a)*pow(lambdaHphi,2)*pow(vh,2)*
          		      (-8*pow(mius,2) + lambdaH*pow(vh,2) +
          		        (lambda12 + lambdaphi)*pow(vsigma,2))*
          		      pow(pow(lambdaH,2)*pow(vh,4) -
          		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5) +
          		     32*lambda12*lambdaHphi*vh*vsigma*pow(2,0.5)*pow(mius,2)*
          		      pow(pow(lambdaH,2)*pow(vh,4) -
          		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) +
          		     32*lambdaHphi*lambdaphi*vh*vsigma*pow(2,0.5)*pow(mius,2)*
          		      pow(pow(lambdaH,2)*pow(vh,4) -
          		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) -
          		     4*lambda12*lambdaH*lambdaHphi*vsigma*pow(2,0.5)*pow(vh,3)*
          		      pow(pow(lambdaH,2)*pow(vh,4) -
          		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) -
          		     4*lambdaH*lambdaHphi*lambdaphi*vsigma*pow(2,0.5)*pow(vh,3)*
          		      pow(pow(lambdaH,2)*pow(vh,4) -
          		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) -
          		     8*lambda12*lambdaHphi*lambdaphi*vh*pow(2,0.5)*pow(vsigma,3)*
          		      pow(pow(lambdaH,2)*pow(vh,4) -
          		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) -
          		     4*lambdaHphi*vh*pow(2,0.5)*pow(lambda12,2)*pow(vsigma,3)*
          		      pow(pow(lambdaH,2)*pow(vh,4) -
          		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) -
          		     4*lambdaHphi*vh*pow(2,0.5)*pow(lambdaphi,2)*pow(vsigma,3)*
          		      pow(pow(lambdaH,2)*pow(vh,4) -
          		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) +
          		     8*lambdaHphi*s*vh*(-2*lambdaHphi*vh*
          		         (-8*pow(mius,2) + lambdaH*pow(vh,2) +
          		           (lambda12 + lambdaphi)*pow(vsigma,2)) +
          		        2*lambdaHphi*vh*cos(2*a)*
          		         pow(pow(lambdaH,2)*pow(vh,4) -
          		           2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		           8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		           pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5) +
          		        (lambda12 + lambdaphi)*vsigma*pow(2,0.5)*
          		         pow(pow(lambdaH,2)*pow(vh,4) -
          		           2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		           8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		           pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a)) +
          		     2*lambda12*lambdaHphi*vsigma*pow(2,0.5)*pow(lambdaH,2)*pow(vh,5)*
          		      sin(4*a) + 2*lambdaHphi*lambdaphi*vsigma*pow(2,0.5)*pow(lambdaH,2)*
          		      pow(vh,5)*sin(4*a) - 8*lambda12*lambdaH*lambdaHphi*lambdaphi*
          		      pow(2,0.5)*pow(vh,3)*pow(vsigma,3)*sin(4*a) -
          		     4*lambdaH*lambdaHphi*pow(2,0.5)*pow(lambda12,2)*pow(vh,3)*pow(vsigma,3)*
          		      sin(4*a) + 16*lambda12*pow(2,0.5)*pow(lambdaHphi,3)*pow(vh,3)*
          		      pow(vsigma,3)*sin(4*a) +
          		     16*lambdaphi*pow(2,0.5)*pow(lambdaHphi,3)*pow(vh,3)*pow(vsigma,3)*
          		      sin(4*a) - 4*lambdaH*lambdaHphi*pow(2,0.5)*pow(lambdaphi,2)*pow(vh,3)*
          		      pow(vsigma,3)*sin(4*a) +
          		     6*lambdaHphi*lambdaphi*vh*pow(2,0.5)*pow(lambda12,2)*pow(vsigma,5)*
          		      sin(4*a) + 2*lambdaHphi*vh*pow(2,0.5)*pow(lambda12,3)*pow(vsigma,5)*
          		      sin(4*a) + 6*lambda12*lambdaHphi*vh*pow(2,0.5)*pow(lambdaphi,2)*
          		      pow(vsigma,5)*sin(4*a) +
          		     2*lambdaHphi*vh*pow(2,0.5)*pow(lambdaphi,3)*pow(vsigma,5)*sin(4*a));
           	den=16*pow(mW,2)*((pow(wh,2)*(lambdaH*pow(vh,2) +
                    (lambda12 + lambdaphi)*pow(vsigma,2) +
                    pow(pow(lambdaH,2)*pow(vh,4) -
                      2*(lambdaH*(lambda12 + lambdaphi) - 4*pow(lambdaHphi,2))*
                       pow(vh,2)*pow(vsigma,2) +
                      pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)))/2. +
               pow(s + (-(lambdaH*pow(vh,2)) - (lambda12 + lambdaphi)*pow(vsigma,2) -
                    pow(pow(lambdaH,2)*pow(vh,4) -
                      2*(lambdaH*(lambda12 + lambdaphi) - 4*pow(lambdaHphi,2))*
                       pow(vh,2)*pow(vsigma,2) +
                      pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5))/2.,2))*
             ((pow(wh2,2)*(lambdaH*pow(vh,2) + (lambda12 + lambdaphi)*pow(vsigma,2) -
                    pow(pow(lambdaH,2)*pow(vh,4) -
                      2*(lambdaH*(lambda12 + lambdaphi) - 4*pow(lambdaHphi,2))*
                       pow(vh,2)*pow(vsigma,2) +
                      pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)))/2. +
               pow(s + (-(lambdaH*pow(vh,2)) - (lambda12 + lambdaphi)*pow(vsigma,2) +
                    pow(pow(lambdaH,2)*pow(vh,4) -
                      2*(lambdaH*(lambda12 + lambdaphi) - 4*pow(lambdaHphi,2))*
                       pow(vh,2)*pow(vsigma,2) +
                      pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5))/2.,2));
           	amp=num/den;
        }
    else if (id == 6)
        {
          	num=-3*pow(g,2)*pow(mqd,2)*(-0.5*s + 2*pow(mqd,2) - 2*pow(mius,2))*
          		   (256*pow(lambdaHphi,2)*pow(mius,4)*pow(vh,2) +
          		     16*pow(lambdaHphi,2)*pow(s,2)*pow(vh,2) -
          		     64*lambdaH*pow(lambdaHphi,2)*pow(mius,2)*pow(vh,4) +
          		     6*pow(lambdaH,2)*pow(lambdaHphi,2)*pow(vh,6) +
          		     2*cos(4*a)*pow(lambdaH,2)*pow(lambdaHphi,2)*pow(vh,6) -
          		     64*lambda12*pow(lambdaHphi,2)*pow(mius,2)*pow(vh,2)*pow(vsigma,2) -
          		     64*lambdaphi*pow(lambdaHphi,2)*pow(mius,2)*pow(vh,2)*pow(vsigma,2) +
          		     2*lambda12*lambdaphi*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,2) -
          		     2*lambda12*lambdaphi*cos(4*a)*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,2) +
          		     pow(lambda12,2)*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,2) -
          		     cos(4*a)*pow(lambda12,2)*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,2) +
          		     4*lambda12*lambdaH*pow(lambdaHphi,2)*pow(vh,4)*pow(vsigma,2) +
          		     4*lambdaH*lambdaphi*pow(lambdaHphi,2)*pow(vh,4)*pow(vsigma,2) -
          		     4*lambda12*lambdaH*cos(4*a)*pow(lambdaHphi,2)*pow(vh,4)*pow(vsigma,2) -
          		     4*lambdaH*lambdaphi*cos(4*a)*pow(lambdaHphi,2)*pow(vh,4)*
          		      pow(vsigma,2) + 16*pow(lambdaHphi,4)*pow(vh,4)*pow(vsigma,2) +
          		     16*cos(4*a)*pow(lambdaHphi,4)*pow(vh,4)*pow(vsigma,2) +
          		     pow(lambdaH,2)*pow(lambdaphi,2)*pow(vh,4)*pow(vsigma,2) -
          		     cos(4*a)*pow(lambdaH,2)*pow(lambdaphi,2)*pow(vh,4)*pow(vsigma,2) -
          		     6*lambdaH*lambdaphi*pow(lambda12,2)*pow(vh,2)*pow(vsigma,4) +
          		     6*lambdaH*lambdaphi*cos(4*a)*pow(lambda12,2)*pow(vh,2)*pow(vsigma,4) -
          		     2*lambdaH*pow(lambda12,3)*pow(vh,2)*pow(vsigma,4) +
          		     2*lambdaH*cos(4*a)*pow(lambda12,3)*pow(vh,2)*pow(vsigma,4) +
          		     28*lambda12*lambdaphi*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,4) -
          		     12*lambda12*lambdaphi*cos(4*a)*pow(lambdaHphi,2)*pow(vh,2)*
          		      pow(vsigma,4) + 14*pow(lambda12,2)*pow(lambdaHphi,2)*pow(vh,2)*
          		      pow(vsigma,4) - 6*cos(4*a)*pow(lambda12,2)*pow(lambdaHphi,2)*pow(vh,2)*
          		      pow(vsigma,4) - 6*lambda12*lambdaH*pow(lambdaphi,2)*pow(vh,2)*
          		      pow(vsigma,4) + 6*lambda12*lambdaH*cos(4*a)*pow(lambdaphi,2)*pow(vh,2)*
          		      pow(vsigma,4) + 14*pow(lambdaHphi,2)*pow(lambdaphi,2)*pow(vh,2)*
          		      pow(vsigma,4) - 6*cos(4*a)*pow(lambdaHphi,2)*pow(lambdaphi,2)*
          		      pow(vh,2)*pow(vsigma,4) -
          		     2*lambdaH*pow(lambdaphi,3)*pow(vh,2)*pow(vsigma,4) +
          		     2*lambdaH*cos(4*a)*pow(lambdaphi,3)*pow(vh,2)*pow(vsigma,4) +
          		     4*lambdaphi*pow(lambda12,3)*pow(vsigma,6) -
          		     4*lambdaphi*cos(4*a)*pow(lambda12,3)*pow(vsigma,6) +
          		     pow(lambda12,4)*pow(vsigma,6) -
          		     cos(4*a)*pow(lambda12,4)*pow(vsigma,6) +
          		     6*pow(lambda12,2)*pow(lambdaphi,2)*pow(vsigma,6) -
          		     6*cos(4*a)*pow(lambda12,2)*pow(lambdaphi,2)*pow(vsigma,6) +
          		     4*lambda12*pow(lambdaphi,3)*pow(vsigma,6) -
          		     4*lambda12*cos(4*a)*pow(lambdaphi,3)*pow(vsigma,6) +
          		     pow(lambdaphi,4)*pow(vsigma,6) -
          		     cos(4*a)*pow(lambdaphi,4)*pow(vsigma,6) -
          		     8*cos(2*a)*pow(lambdaHphi,2)*pow(vh,2)*
          		      (-8*pow(mius,2) + lambdaH*pow(vh,2) +
          		        (lambda12 + lambdaphi)*pow(vsigma,2))*
          		      pow(pow(lambdaH,2)*pow(vh,4) -
          		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5) +
          		     32*lambda12*lambdaHphi*vh*vsigma*pow(2,0.5)*pow(mius,2)*
          		      pow(pow(lambdaH,2)*pow(vh,4) -
          		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) +
          		     32*lambdaHphi*lambdaphi*vh*vsigma*pow(2,0.5)*pow(mius,2)*
          		      pow(pow(lambdaH,2)*pow(vh,4) -
          		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) -
          		     4*lambda12*lambdaH*lambdaHphi*vsigma*pow(2,0.5)*pow(vh,3)*
          		      pow(pow(lambdaH,2)*pow(vh,4) -
          		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) -
          		     4*lambdaH*lambdaHphi*lambdaphi*vsigma*pow(2,0.5)*pow(vh,3)*
          		      pow(pow(lambdaH,2)*pow(vh,4) -
          		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) -
          		     8*lambda12*lambdaHphi*lambdaphi*vh*pow(2,0.5)*pow(vsigma,3)*
          		      pow(pow(lambdaH,2)*pow(vh,4) -
          		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) -
          		     4*lambdaHphi*vh*pow(2,0.5)*pow(lambda12,2)*pow(vsigma,3)*
          		      pow(pow(lambdaH,2)*pow(vh,4) -
          		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) -
          		     4*lambdaHphi*vh*pow(2,0.5)*pow(lambdaphi,2)*pow(vsigma,3)*
          		      pow(pow(lambdaH,2)*pow(vh,4) -
          		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) +
          		     8*lambdaHphi*s*vh*(-2*lambdaHphi*vh*
          		         (-8*pow(mius,2) + lambdaH*pow(vh,2) +
          		           (lambda12 + lambdaphi)*pow(vsigma,2)) +
          		        2*lambdaHphi*vh*cos(2*a)*
          		         pow(pow(lambdaH,2)*pow(vh,4) -
          		           2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		           8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		           pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5) +
          		        (lambda12 + lambdaphi)*vsigma*pow(2,0.5)*
          		         pow(pow(lambdaH,2)*pow(vh,4) -
          		           2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		           8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		           pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a)) +
          		     2*lambda12*lambdaHphi*vsigma*pow(2,0.5)*pow(lambdaH,2)*pow(vh,5)*
          		      sin(4*a) + 2*lambdaHphi*lambdaphi*vsigma*pow(2,0.5)*pow(lambdaH,2)*
          		      pow(vh,5)*sin(4*a) - 8*lambda12*lambdaH*lambdaHphi*lambdaphi*
          		      pow(2,0.5)*pow(vh,3)*pow(vsigma,3)*sin(4*a) -
          		     4*lambdaH*lambdaHphi*pow(2,0.5)*pow(lambda12,2)*pow(vh,3)*pow(vsigma,3)*
          		      sin(4*a) + 16*lambda12*pow(2,0.5)*pow(lambdaHphi,3)*pow(vh,3)*
          		      pow(vsigma,3)*sin(4*a) +
          		     16*lambdaphi*pow(2,0.5)*pow(lambdaHphi,3)*pow(vh,3)*pow(vsigma,3)*
          		      sin(4*a) - 4*lambdaH*lambdaHphi*pow(2,0.5)*pow(lambdaphi,2)*pow(vh,3)*
          		      pow(vsigma,3)*sin(4*a) +
          		     6*lambdaHphi*lambdaphi*vh*pow(2,0.5)*pow(lambda12,2)*pow(vsigma,5)*
          		      sin(4*a) + 2*lambdaHphi*vh*pow(2,0.5)*pow(lambda12,3)*pow(vsigma,5)*
          		      sin(4*a) + 6*lambda12*lambdaHphi*vh*pow(2,0.5)*pow(lambdaphi,2)*
          		      pow(vsigma,5)*sin(4*a) +
          		     2*lambdaHphi*vh*pow(2,0.5)*pow(lambdaphi,3)*pow(vsigma,5)*sin(4*a));
           	den=16*pow(mW,2)*((pow(wh,2)*(lambdaH*pow(vh,2) +
                    (lambda12 + lambdaphi)*pow(vsigma,2) +
                    pow(pow(lambdaH,2)*pow(vh,4) -
                      2*(lambdaH*(lambda12 + lambdaphi) - 4*pow(lambdaHphi,2))*
                       pow(vh,2)*pow(vsigma,2) +
                      pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)))/2. +
               pow(s + (-(lambdaH*pow(vh,2)) - (lambda12 + lambdaphi)*pow(vsigma,2) -
                    pow(pow(lambdaH,2)*pow(vh,4) -
                      2*(lambdaH*(lambda12 + lambdaphi) - 4*pow(lambdaHphi,2))*
                       pow(vh,2)*pow(vsigma,2) +
                      pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5))/2.,2))*
             ((pow(wh2,2)*(lambdaH*pow(vh,2) + (lambda12 + lambdaphi)*pow(vsigma,2) -
                    pow(pow(lambdaH,2)*pow(vh,4) -
                      2*(lambdaH*(lambda12 + lambdaphi) - 4*pow(lambdaHphi,2))*
                       pow(vh,2)*pow(vsigma,2) +
                      pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)))/2. +
               pow(s + (-(lambdaH*pow(vh,2)) - (lambda12 + lambdaphi)*pow(vsigma,2) +
                    pow(pow(lambdaH,2)*pow(vh,4) -
                      2*(lambdaH*(lambda12 + lambdaphi) - 4*pow(lambdaHphi,2))*
                       pow(vh,2)*pow(vsigma,2) +
                      pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5))/2.,2));
           	amp=num/den;
        }
    else if (id == 7)
        {
          	num=-3*pow(g,2)*pow(mqu,2)*(-0.5*s - 2*pow(mius,2) + 2*pow(mqu,2))*
          		   (256*pow(lambdaHphi,2)*pow(mius,4)*pow(vh,2) +
          		     16*pow(lambdaHphi,2)*pow(s,2)*pow(vh,2) -
          		     64*lambdaH*pow(lambdaHphi,2)*pow(mius,2)*pow(vh,4) +
          		     6*pow(lambdaH,2)*pow(lambdaHphi,2)*pow(vh,6) +
          		     2*cos(4*a)*pow(lambdaH,2)*pow(lambdaHphi,2)*pow(vh,6) -
          		     64*lambda12*pow(lambdaHphi,2)*pow(mius,2)*pow(vh,2)*pow(vsigma,2) -
          		     64*lambdaphi*pow(lambdaHphi,2)*pow(mius,2)*pow(vh,2)*pow(vsigma,2) +
          		     2*lambda12*lambdaphi*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,2) -
          		     2*lambda12*lambdaphi*cos(4*a)*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,2) +
          		     pow(lambda12,2)*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,2) -
          		     cos(4*a)*pow(lambda12,2)*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,2) +
          		     4*lambda12*lambdaH*pow(lambdaHphi,2)*pow(vh,4)*pow(vsigma,2) +
          		     4*lambdaH*lambdaphi*pow(lambdaHphi,2)*pow(vh,4)*pow(vsigma,2) -
          		     4*lambda12*lambdaH*cos(4*a)*pow(lambdaHphi,2)*pow(vh,4)*pow(vsigma,2) -
          		     4*lambdaH*lambdaphi*cos(4*a)*pow(lambdaHphi,2)*pow(vh,4)*
          		      pow(vsigma,2) + 16*pow(lambdaHphi,4)*pow(vh,4)*pow(vsigma,2) +
          		     16*cos(4*a)*pow(lambdaHphi,4)*pow(vh,4)*pow(vsigma,2) +
          		     pow(lambdaH,2)*pow(lambdaphi,2)*pow(vh,4)*pow(vsigma,2) -
          		     cos(4*a)*pow(lambdaH,2)*pow(lambdaphi,2)*pow(vh,4)*pow(vsigma,2) -
          		     6*lambdaH*lambdaphi*pow(lambda12,2)*pow(vh,2)*pow(vsigma,4) +
          		     6*lambdaH*lambdaphi*cos(4*a)*pow(lambda12,2)*pow(vh,2)*pow(vsigma,4) -
          		     2*lambdaH*pow(lambda12,3)*pow(vh,2)*pow(vsigma,4) +
          		     2*lambdaH*cos(4*a)*pow(lambda12,3)*pow(vh,2)*pow(vsigma,4) +
          		     28*lambda12*lambdaphi*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,4) -
          		     12*lambda12*lambdaphi*cos(4*a)*pow(lambdaHphi,2)*pow(vh,2)*
          		      pow(vsigma,4) + 14*pow(lambda12,2)*pow(lambdaHphi,2)*pow(vh,2)*
          		      pow(vsigma,4) - 6*cos(4*a)*pow(lambda12,2)*pow(lambdaHphi,2)*pow(vh,2)*
          		      pow(vsigma,4) - 6*lambda12*lambdaH*pow(lambdaphi,2)*pow(vh,2)*
          		      pow(vsigma,4) + 6*lambda12*lambdaH*cos(4*a)*pow(lambdaphi,2)*pow(vh,2)*
          		      pow(vsigma,4) + 14*pow(lambdaHphi,2)*pow(lambdaphi,2)*pow(vh,2)*
          		      pow(vsigma,4) - 6*cos(4*a)*pow(lambdaHphi,2)*pow(lambdaphi,2)*
          		      pow(vh,2)*pow(vsigma,4) -
          		     2*lambdaH*pow(lambdaphi,3)*pow(vh,2)*pow(vsigma,4) +
          		     2*lambdaH*cos(4*a)*pow(lambdaphi,3)*pow(vh,2)*pow(vsigma,4) +
          		     4*lambdaphi*pow(lambda12,3)*pow(vsigma,6) -
          		     4*lambdaphi*cos(4*a)*pow(lambda12,3)*pow(vsigma,6) +
          		     pow(lambda12,4)*pow(vsigma,6) -
          		     cos(4*a)*pow(lambda12,4)*pow(vsigma,6) +
          		     6*pow(lambda12,2)*pow(lambdaphi,2)*pow(vsigma,6) -
          		     6*cos(4*a)*pow(lambda12,2)*pow(lambdaphi,2)*pow(vsigma,6) +
          		     4*lambda12*pow(lambdaphi,3)*pow(vsigma,6) -
          		     4*lambda12*cos(4*a)*pow(lambdaphi,3)*pow(vsigma,6) +
          		     pow(lambdaphi,4)*pow(vsigma,6) -
          		     cos(4*a)*pow(lambdaphi,4)*pow(vsigma,6) -
          		     8*cos(2*a)*pow(lambdaHphi,2)*pow(vh,2)*
          		      (-8*pow(mius,2) + lambdaH*pow(vh,2) +
          		        (lambda12 + lambdaphi)*pow(vsigma,2))*
          		      pow(pow(lambdaH,2)*pow(vh,4) -
          		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5) +
          		     32*lambda12*lambdaHphi*vh*vsigma*pow(2,0.5)*pow(mius,2)*
          		      pow(pow(lambdaH,2)*pow(vh,4) -
          		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) +
          		     32*lambdaHphi*lambdaphi*vh*vsigma*pow(2,0.5)*pow(mius,2)*
          		      pow(pow(lambdaH,2)*pow(vh,4) -
          		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) -
          		     4*lambda12*lambdaH*lambdaHphi*vsigma*pow(2,0.5)*pow(vh,3)*
          		      pow(pow(lambdaH,2)*pow(vh,4) -
          		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) -
          		     4*lambdaH*lambdaHphi*lambdaphi*vsigma*pow(2,0.5)*pow(vh,3)*
          		      pow(pow(lambdaH,2)*pow(vh,4) -
          		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) -
          		     8*lambda12*lambdaHphi*lambdaphi*vh*pow(2,0.5)*pow(vsigma,3)*
          		      pow(pow(lambdaH,2)*pow(vh,4) -
          		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) -
          		     4*lambdaHphi*vh*pow(2,0.5)*pow(lambda12,2)*pow(vsigma,3)*
          		      pow(pow(lambdaH,2)*pow(vh,4) -
          		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) -
          		     4*lambdaHphi*vh*pow(2,0.5)*pow(lambdaphi,2)*pow(vsigma,3)*
          		      pow(pow(lambdaH,2)*pow(vh,4) -
          		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) +
          		     8*lambdaHphi*s*vh*(-2*lambdaHphi*vh*
          		         (-8*pow(mius,2) + lambdaH*pow(vh,2) +
          		           (lambda12 + lambdaphi)*pow(vsigma,2)) +
          		        2*lambdaHphi*vh*cos(2*a)*
          		         pow(pow(lambdaH,2)*pow(vh,4) -
          		           2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		           8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		           pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5) +
          		        (lambda12 + lambdaphi)*vsigma*pow(2,0.5)*
          		         pow(pow(lambdaH,2)*pow(vh,4) -
          		           2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		           8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		           pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a)) +
          		     2*lambda12*lambdaHphi*vsigma*pow(2,0.5)*pow(lambdaH,2)*pow(vh,5)*
          		      sin(4*a) + 2*lambdaHphi*lambdaphi*vsigma*pow(2,0.5)*pow(lambdaH,2)*
          		      pow(vh,5)*sin(4*a) - 8*lambda12*lambdaH*lambdaHphi*lambdaphi*
          		      pow(2,0.5)*pow(vh,3)*pow(vsigma,3)*sin(4*a) -
          		     4*lambdaH*lambdaHphi*pow(2,0.5)*pow(lambda12,2)*pow(vh,3)*pow(vsigma,3)*
          		      sin(4*a) + 16*lambda12*pow(2,0.5)*pow(lambdaHphi,3)*pow(vh,3)*
          		      pow(vsigma,3)*sin(4*a) +
          		     16*lambdaphi*pow(2,0.5)*pow(lambdaHphi,3)*pow(vh,3)*pow(vsigma,3)*
          		      sin(4*a) - 4*lambdaH*lambdaHphi*pow(2,0.5)*pow(lambdaphi,2)*pow(vh,3)*
          		      pow(vsigma,3)*sin(4*a) +
          		     6*lambdaHphi*lambdaphi*vh*pow(2,0.5)*pow(lambda12,2)*pow(vsigma,5)*
          		      sin(4*a) + 2*lambdaHphi*vh*pow(2,0.5)*pow(lambda12,3)*pow(vsigma,5)*
          		      sin(4*a) + 6*lambda12*lambdaHphi*vh*pow(2,0.5)*pow(lambdaphi,2)*
          		      pow(vsigma,5)*sin(4*a) +
          		     2*lambdaHphi*vh*pow(2,0.5)*pow(lambdaphi,3)*pow(vsigma,5)*sin(4*a));
           	den=16*pow(mW,2)*((pow(wh,2)*(lambdaH*pow(vh,2) +
                    (lambda12 + lambdaphi)*pow(vsigma,2) +
                    pow(pow(lambdaH,2)*pow(vh,4) -
                      2*(lambdaH*(lambda12 + lambdaphi) - 4*pow(lambdaHphi,2))*
                       pow(vh,2)*pow(vsigma,2) +
                      pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)))/2. +
               pow(s + (-(lambdaH*pow(vh,2)) - (lambda12 + lambdaphi)*pow(vsigma,2) -
                    pow(pow(lambdaH,2)*pow(vh,4) -
                      2*(lambdaH*(lambda12 + lambdaphi) - 4*pow(lambdaHphi,2))*
                       pow(vh,2)*pow(vsigma,2) +
                      pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5))/2.,2))*
             ((pow(wh2,2)*(lambdaH*pow(vh,2) + (lambda12 + lambdaphi)*pow(vsigma,2) -
                    pow(pow(lambdaH,2)*pow(vh,4) -
                      2*(lambdaH*(lambda12 + lambdaphi) - 4*pow(lambdaHphi,2))*
                       pow(vh,2)*pow(vsigma,2) +
                      pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)))/2. +
               pow(s + (-(lambdaH*pow(vh,2)) - (lambda12 + lambdaphi)*pow(vsigma,2) +
                    pow(pow(lambdaH,2)*pow(vh,4) -
                      2*(lambdaH*(lambda12 + lambdaphi) - 4*pow(lambdaHphi,2))*
                       pow(vh,2)*pow(vsigma,2) +
                      pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5))/2.,2));
           	amp=num/den;
        }
    else if (id == 8)
        {
          	num=-(pow(g,2)*pow(mtau,2)*(-0.5*s - 2*pow(mius,2) + 2*pow(mtau,2))*
          		     (256*pow(lambdaHphi,2)*pow(mius,4)*pow(vh,2) +
          		       16*pow(lambdaHphi,2)*pow(s,2)*pow(vh,2) -
          		       64*lambdaH*pow(lambdaHphi,2)*pow(mius,2)*pow(vh,4) +
          		       6*pow(lambdaH,2)*pow(lambdaHphi,2)*pow(vh,6) +
          		       2*cos(4*a)*pow(lambdaH,2)*pow(lambdaHphi,2)*pow(vh,6) -
          		       64*lambda12*pow(lambdaHphi,2)*pow(mius,2)*pow(vh,2)*pow(vsigma,2) -
          		       64*lambdaphi*pow(lambdaHphi,2)*pow(mius,2)*pow(vh,2)*pow(vsigma,2) +
          		       2*lambda12*lambdaphi*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,2) -
          		       2*lambda12*lambdaphi*cos(4*a)*pow(lambdaH,2)*pow(vh,4)*
          		        pow(vsigma,2) + pow(lambda12,2)*pow(lambdaH,2)*pow(vh,4)*
          		        pow(vsigma,2) - cos(4*a)*pow(lambda12,2)*pow(lambdaH,2)*pow(vh,4)*
          		        pow(vsigma,2) + 4*lambda12*lambdaH*pow(lambdaHphi,2)*pow(vh,4)*
          		        pow(vsigma,2) + 4*lambdaH*lambdaphi*pow(lambdaHphi,2)*pow(vh,4)*
          		        pow(vsigma,2) - 4*lambda12*lambdaH*cos(4*a)*pow(lambdaHphi,2)*
          		        pow(vh,4)*pow(vsigma,2) -
          		       4*lambdaH*lambdaphi*cos(4*a)*pow(lambdaHphi,2)*pow(vh,4)*
          		        pow(vsigma,2) + 16*pow(lambdaHphi,4)*pow(vh,4)*pow(vsigma,2) +
          		       16*cos(4*a)*pow(lambdaHphi,4)*pow(vh,4)*pow(vsigma,2) +
          		       pow(lambdaH,2)*pow(lambdaphi,2)*pow(vh,4)*pow(vsigma,2) -
          		       cos(4*a)*pow(lambdaH,2)*pow(lambdaphi,2)*pow(vh,4)*pow(vsigma,2) -
          		       6*lambdaH*lambdaphi*pow(lambda12,2)*pow(vh,2)*pow(vsigma,4) +
          		       6*lambdaH*lambdaphi*cos(4*a)*pow(lambda12,2)*pow(vh,2)*
          		        pow(vsigma,4) - 2*lambdaH*pow(lambda12,3)*pow(vh,2)*pow(vsigma,4) +
          		       2*lambdaH*cos(4*a)*pow(lambda12,3)*pow(vh,2)*pow(vsigma,4) +
          		       28*lambda12*lambdaphi*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,4) -
          		       12*lambda12*lambdaphi*cos(4*a)*pow(lambdaHphi,2)*pow(vh,2)*
          		        pow(vsigma,4) + 14*pow(lambda12,2)*pow(lambdaHphi,2)*pow(vh,2)*
          		        pow(vsigma,4) - 6*cos(4*a)*pow(lambda12,2)*pow(lambdaHphi,2)*
          		        pow(vh,2)*pow(vsigma,4) -
          		       6*lambda12*lambdaH*pow(lambdaphi,2)*pow(vh,2)*pow(vsigma,4) +
          		       6*lambda12*lambdaH*cos(4*a)*pow(lambdaphi,2)*pow(vh,2)*
          		        pow(vsigma,4) + 14*pow(lambdaHphi,2)*pow(lambdaphi,2)*pow(vh,2)*
          		        pow(vsigma,4) - 6*cos(4*a)*pow(lambdaHphi,2)*pow(lambdaphi,2)*
          		        pow(vh,2)*pow(vsigma,4) -
          		       2*lambdaH*pow(lambdaphi,3)*pow(vh,2)*pow(vsigma,4) +
          		       2*lambdaH*cos(4*a)*pow(lambdaphi,3)*pow(vh,2)*pow(vsigma,4) +
          		       4*lambdaphi*pow(lambda12,3)*pow(vsigma,6) -
          		       4*lambdaphi*cos(4*a)*pow(lambda12,3)*pow(vsigma,6) +
          		       pow(lambda12,4)*pow(vsigma,6) -
          		       cos(4*a)*pow(lambda12,4)*pow(vsigma,6) +
          		       6*pow(lambda12,2)*pow(lambdaphi,2)*pow(vsigma,6) -
          		       6*cos(4*a)*pow(lambda12,2)*pow(lambdaphi,2)*pow(vsigma,6) +
          		       4*lambda12*pow(lambdaphi,3)*pow(vsigma,6) -
          		       4*lambda12*cos(4*a)*pow(lambdaphi,3)*pow(vsigma,6) +
          		       pow(lambdaphi,4)*pow(vsigma,6) -
          		       cos(4*a)*pow(lambdaphi,4)*pow(vsigma,6) -
          		       8*cos(2*a)*pow(lambdaHphi,2)*pow(vh,2)*
          		        (-8*pow(mius,2) + lambdaH*pow(vh,2) +
          		          (lambda12 + lambdaphi)*pow(vsigma,2))*
          		        pow(pow(lambdaH,2)*pow(vh,4) -
          		          2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		          8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		          pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5) +
          		       32*lambda12*lambdaHphi*vh*vsigma*pow(2,0.5)*pow(mius,2)*
          		        pow(pow(lambdaH,2)*pow(vh,4) -
          		          2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		          8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		          pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) +
          		       32*lambdaHphi*lambdaphi*vh*vsigma*pow(2,0.5)*pow(mius,2)*
          		        pow(pow(lambdaH,2)*pow(vh,4) -
          		          2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		          8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		          pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) -
          		       4*lambda12*lambdaH*lambdaHphi*vsigma*pow(2,0.5)*pow(vh,3)*
          		        pow(pow(lambdaH,2)*pow(vh,4) -
          		          2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		          8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		          pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) -
          		       4*lambdaH*lambdaHphi*lambdaphi*vsigma*pow(2,0.5)*pow(vh,3)*
          		        pow(pow(lambdaH,2)*pow(vh,4) -
          		          2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		          8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		          pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) -
          		       8*lambda12*lambdaHphi*lambdaphi*vh*pow(2,0.5)*pow(vsigma,3)*
          		        pow(pow(lambdaH,2)*pow(vh,4) -
          		          2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		          8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		          pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) -
          		       4*lambdaHphi*vh*pow(2,0.5)*pow(lambda12,2)*pow(vsigma,3)*
          		        pow(pow(lambdaH,2)*pow(vh,4) -
          		          2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		          8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		          pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) -
          		       4*lambdaHphi*vh*pow(2,0.5)*pow(lambdaphi,2)*pow(vsigma,3)*
          		        pow(pow(lambdaH,2)*pow(vh,4) -
          		          2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		          8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		          pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) +
          		       8*lambdaHphi*s*vh*(-2*lambdaHphi*vh*
          		           (-8*pow(mius,2) + lambdaH*pow(vh,2) +
          		             (lambda12 + lambdaphi)*pow(vsigma,2)) +
          		          2*lambdaHphi*vh*cos(2*a)*
          		           pow(pow(lambdaH,2)*pow(vh,4) -
          		             2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		             8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		             pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5) +
          		          (lambda12 + lambdaphi)*vsigma*pow(2,0.5)*
          		           pow(pow(lambdaH,2)*pow(vh,4) -
          		             2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		             8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		             pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a)) +
          		       2*lambda12*lambdaHphi*vsigma*pow(2,0.5)*pow(lambdaH,2)*pow(vh,5)*
          		        sin(4*a) + 2*lambdaHphi*lambdaphi*vsigma*pow(2,0.5)*pow(lambdaH,2)*
          		        pow(vh,5)*sin(4*a) - 8*lambda12*lambdaH*lambdaHphi*lambdaphi*
          		        pow(2,0.5)*pow(vh,3)*pow(vsigma,3)*sin(4*a) -
          		       4*lambdaH*lambdaHphi*pow(2,0.5)*pow(lambda12,2)*pow(vh,3)*
          		        pow(vsigma,3)*sin(4*a) +
          		       16*lambda12*pow(2,0.5)*pow(lambdaHphi,3)*pow(vh,3)*pow(vsigma,3)*
          		        sin(4*a) + 16*lambdaphi*pow(2,0.5)*pow(lambdaHphi,3)*pow(vh,3)*
          		        pow(vsigma,3)*sin(4*a) -
          		       4*lambdaH*lambdaHphi*pow(2,0.5)*pow(lambdaphi,2)*pow(vh,3)*
          		        pow(vsigma,3)*sin(4*a) +
          		       6*lambdaHphi*lambdaphi*vh*pow(2,0.5)*pow(lambda12,2)*pow(vsigma,5)*
          		        sin(4*a) + 2*lambdaHphi*vh*pow(2,0.5)*pow(lambda12,3)*pow(vsigma,5)*
          		        sin(4*a) + 6*lambda12*lambdaHphi*vh*pow(2,0.5)*pow(lambdaphi,2)*
          		        pow(vsigma,5)*sin(4*a) +
          		       2*lambdaHphi*vh*pow(2,0.5)*pow(lambdaphi,3)*pow(vsigma,5)*sin(4*a)));
           	den=16*pow(mW,2)*((pow(wh,2)*(lambdaH*pow(vh,2) +
                    (lambda12 + lambdaphi)*pow(vsigma,2) +
                    pow(pow(lambdaH,2)*pow(vh,4) -
                      2*(lambdaH*(lambda12 + lambdaphi) - 4*pow(lambdaHphi,2))*
                       pow(vh,2)*pow(vsigma,2) +
                      pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)))/2. +
               pow(s + (-(lambdaH*pow(vh,2)) - (lambda12 + lambdaphi)*pow(vsigma,2) -
                    pow(pow(lambdaH,2)*pow(vh,4) -
                      2*(lambdaH*(lambda12 + lambdaphi) - 4*pow(lambdaHphi,2))*
                       pow(vh,2)*pow(vsigma,2) +
                      pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5))/2.,2))*
             ((pow(wh2,2)*(lambdaH*pow(vh,2) + (lambda12 + lambdaphi)*pow(vsigma,2) -
                    pow(pow(lambdaH,2)*pow(vh,4) -
                      2*(lambdaH*(lambda12 + lambdaphi) - 4*pow(lambdaHphi,2))*
                       pow(vh,2)*pow(vsigma,2) +
                      pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)))/2. +
               pow(s + (-(lambdaH*pow(vh,2)) - (lambda12 + lambdaphi)*pow(vsigma,2) +
                    pow(pow(lambdaH,2)*pow(vh,4) -
                      2*(lambdaH*(lambda12 + lambdaphi) - 4*pow(lambdaHphi,2))*
                       pow(vh,2)*pow(vsigma,2) +
                      pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5))/2.,2));
           	amp=num/den;
        }
    else if (id == 9)
        {
          	num=-(pow(g,2)*pow(mmu,2)*(-0.5*s - 2*pow(mius,2) + 2*pow(mmu,2))*
          		     (256*pow(lambdaHphi,2)*pow(mius,4)*pow(vh,2) +
          		       16*pow(lambdaHphi,2)*pow(s,2)*pow(vh,2) -
          		       64*lambdaH*pow(lambdaHphi,2)*pow(mius,2)*pow(vh,4) +
          		       6*pow(lambdaH,2)*pow(lambdaHphi,2)*pow(vh,6) +
          		       2*cos(4*a)*pow(lambdaH,2)*pow(lambdaHphi,2)*pow(vh,6) -
          		       64*lambda12*pow(lambdaHphi,2)*pow(mius,2)*pow(vh,2)*pow(vsigma,2) -
          		       64*lambdaphi*pow(lambdaHphi,2)*pow(mius,2)*pow(vh,2)*pow(vsigma,2) +
          		       2*lambda12*lambdaphi*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,2) -
          		       2*lambda12*lambdaphi*cos(4*a)*pow(lambdaH,2)*pow(vh,4)*
          		        pow(vsigma,2) + pow(lambda12,2)*pow(lambdaH,2)*pow(vh,4)*
          		        pow(vsigma,2) - cos(4*a)*pow(lambda12,2)*pow(lambdaH,2)*pow(vh,4)*
          		        pow(vsigma,2) + 4*lambda12*lambdaH*pow(lambdaHphi,2)*pow(vh,4)*
          		        pow(vsigma,2) + 4*lambdaH*lambdaphi*pow(lambdaHphi,2)*pow(vh,4)*
          		        pow(vsigma,2) - 4*lambda12*lambdaH*cos(4*a)*pow(lambdaHphi,2)*
          		        pow(vh,4)*pow(vsigma,2) -
          		       4*lambdaH*lambdaphi*cos(4*a)*pow(lambdaHphi,2)*pow(vh,4)*
          		        pow(vsigma,2) + 16*pow(lambdaHphi,4)*pow(vh,4)*pow(vsigma,2) +
          		       16*cos(4*a)*pow(lambdaHphi,4)*pow(vh,4)*pow(vsigma,2) +
          		       pow(lambdaH,2)*pow(lambdaphi,2)*pow(vh,4)*pow(vsigma,2) -
          		       cos(4*a)*pow(lambdaH,2)*pow(lambdaphi,2)*pow(vh,4)*pow(vsigma,2) -
          		       6*lambdaH*lambdaphi*pow(lambda12,2)*pow(vh,2)*pow(vsigma,4) +
          		       6*lambdaH*lambdaphi*cos(4*a)*pow(lambda12,2)*pow(vh,2)*
          		        pow(vsigma,4) - 2*lambdaH*pow(lambda12,3)*pow(vh,2)*pow(vsigma,4) +
          		       2*lambdaH*cos(4*a)*pow(lambda12,3)*pow(vh,2)*pow(vsigma,4) +
          		       28*lambda12*lambdaphi*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,4) -
          		       12*lambda12*lambdaphi*cos(4*a)*pow(lambdaHphi,2)*pow(vh,2)*
          		        pow(vsigma,4) + 14*pow(lambda12,2)*pow(lambdaHphi,2)*pow(vh,2)*
          		        pow(vsigma,4) - 6*cos(4*a)*pow(lambda12,2)*pow(lambdaHphi,2)*
          		        pow(vh,2)*pow(vsigma,4) -
          		       6*lambda12*lambdaH*pow(lambdaphi,2)*pow(vh,2)*pow(vsigma,4) +
          		       6*lambda12*lambdaH*cos(4*a)*pow(lambdaphi,2)*pow(vh,2)*
          		        pow(vsigma,4) + 14*pow(lambdaHphi,2)*pow(lambdaphi,2)*pow(vh,2)*
          		        pow(vsigma,4) - 6*cos(4*a)*pow(lambdaHphi,2)*pow(lambdaphi,2)*
          		        pow(vh,2)*pow(vsigma,4) -
          		       2*lambdaH*pow(lambdaphi,3)*pow(vh,2)*pow(vsigma,4) +
          		       2*lambdaH*cos(4*a)*pow(lambdaphi,3)*pow(vh,2)*pow(vsigma,4) +
          		       4*lambdaphi*pow(lambda12,3)*pow(vsigma,6) -
          		       4*lambdaphi*cos(4*a)*pow(lambda12,3)*pow(vsigma,6) +
          		       pow(lambda12,4)*pow(vsigma,6) -
          		       cos(4*a)*pow(lambda12,4)*pow(vsigma,6) +
          		       6*pow(lambda12,2)*pow(lambdaphi,2)*pow(vsigma,6) -
          		       6*cos(4*a)*pow(lambda12,2)*pow(lambdaphi,2)*pow(vsigma,6) +
          		       4*lambda12*pow(lambdaphi,3)*pow(vsigma,6) -
          		       4*lambda12*cos(4*a)*pow(lambdaphi,3)*pow(vsigma,6) +
          		       pow(lambdaphi,4)*pow(vsigma,6) -
          		       cos(4*a)*pow(lambdaphi,4)*pow(vsigma,6) -
          		       8*cos(2*a)*pow(lambdaHphi,2)*pow(vh,2)*
          		        (-8*pow(mius,2) + lambdaH*pow(vh,2) +
          		          (lambda12 + lambdaphi)*pow(vsigma,2))*
          		        pow(pow(lambdaH,2)*pow(vh,4) -
          		          2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		          8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		          pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5) +
          		       32*lambda12*lambdaHphi*vh*vsigma*pow(2,0.5)*pow(mius,2)*
          		        pow(pow(lambdaH,2)*pow(vh,4) -
          		          2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		          8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		          pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) +
          		       32*lambdaHphi*lambdaphi*vh*vsigma*pow(2,0.5)*pow(mius,2)*
          		        pow(pow(lambdaH,2)*pow(vh,4) -
          		          2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		          8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		          pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) -
          		       4*lambda12*lambdaH*lambdaHphi*vsigma*pow(2,0.5)*pow(vh,3)*
          		        pow(pow(lambdaH,2)*pow(vh,4) -
          		          2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		          8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		          pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) -
          		       4*lambdaH*lambdaHphi*lambdaphi*vsigma*pow(2,0.5)*pow(vh,3)*
          		        pow(pow(lambdaH,2)*pow(vh,4) -
          		          2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		          8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		          pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) -
          		       8*lambda12*lambdaHphi*lambdaphi*vh*pow(2,0.5)*pow(vsigma,3)*
          		        pow(pow(lambdaH,2)*pow(vh,4) -
          		          2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		          8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		          pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) -
          		       4*lambdaHphi*vh*pow(2,0.5)*pow(lambda12,2)*pow(vsigma,3)*
          		        pow(pow(lambdaH,2)*pow(vh,4) -
          		          2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		          8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		          pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) -
          		       4*lambdaHphi*vh*pow(2,0.5)*pow(lambdaphi,2)*pow(vsigma,3)*
          		        pow(pow(lambdaH,2)*pow(vh,4) -
          		          2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		          8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		          pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) +
          		       8*lambdaHphi*s*vh*(-2*lambdaHphi*vh*
          		           (-8*pow(mius,2) + lambdaH*pow(vh,2) +
          		             (lambda12 + lambdaphi)*pow(vsigma,2)) +
          		          2*lambdaHphi*vh*cos(2*a)*
          		           pow(pow(lambdaH,2)*pow(vh,4) -
          		             2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		             8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		             pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5) +
          		          (lambda12 + lambdaphi)*vsigma*pow(2,0.5)*
          		           pow(pow(lambdaH,2)*pow(vh,4) -
          		             2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		             8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		             pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a)) +
          		       2*lambda12*lambdaHphi*vsigma*pow(2,0.5)*pow(lambdaH,2)*pow(vh,5)*
          		        sin(4*a) + 2*lambdaHphi*lambdaphi*vsigma*pow(2,0.5)*pow(lambdaH,2)*
          		        pow(vh,5)*sin(4*a) - 8*lambda12*lambdaH*lambdaHphi*lambdaphi*
          		        pow(2,0.5)*pow(vh,3)*pow(vsigma,3)*sin(4*a) -
          		       4*lambdaH*lambdaHphi*pow(2,0.5)*pow(lambda12,2)*pow(vh,3)*
          		        pow(vsigma,3)*sin(4*a) +
          		       16*lambda12*pow(2,0.5)*pow(lambdaHphi,3)*pow(vh,3)*pow(vsigma,3)*
          		        sin(4*a) + 16*lambdaphi*pow(2,0.5)*pow(lambdaHphi,3)*pow(vh,3)*
          		        pow(vsigma,3)*sin(4*a) -
          		       4*lambdaH*lambdaHphi*pow(2,0.5)*pow(lambdaphi,2)*pow(vh,3)*
          		        pow(vsigma,3)*sin(4*a) +
          		       6*lambdaHphi*lambdaphi*vh*pow(2,0.5)*pow(lambda12,2)*pow(vsigma,5)*
          		        sin(4*a) + 2*lambdaHphi*vh*pow(2,0.5)*pow(lambda12,3)*pow(vsigma,5)*
          		        sin(4*a) + 6*lambda12*lambdaHphi*vh*pow(2,0.5)*pow(lambdaphi,2)*
          		        pow(vsigma,5)*sin(4*a) +
          		       2*lambdaHphi*vh*pow(2,0.5)*pow(lambdaphi,3)*pow(vsigma,5)*sin(4*a)));
           	den=16*pow(mW,2)*((pow(wh,2)*(lambdaH*pow(vh,2) +
                    (lambda12 + lambdaphi)*pow(vsigma,2) +
                    pow(pow(lambdaH,2)*pow(vh,4) -
                      2*(lambdaH*(lambda12 + lambdaphi) - 4*pow(lambdaHphi,2))*
                       pow(vh,2)*pow(vsigma,2) +
                      pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)))/2. +
               pow(s + (-(lambdaH*pow(vh,2)) - (lambda12 + lambdaphi)*pow(vsigma,2) -
                    pow(pow(lambdaH,2)*pow(vh,4) -
                      2*(lambdaH*(lambda12 + lambdaphi) - 4*pow(lambdaHphi,2))*
                       pow(vh,2)*pow(vsigma,2) +
                      pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5))/2.,2))*
             ((pow(wh2,2)*(lambdaH*pow(vh,2) + (lambda12 + lambdaphi)*pow(vsigma,2) -
                    pow(pow(lambdaH,2)*pow(vh,4) -
                      2*(lambdaH*(lambda12 + lambdaphi) - 4*pow(lambdaHphi,2))*
                       pow(vh,2)*pow(vsigma,2) +
                      pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)))/2. +
               pow(s + (-(lambdaH*pow(vh,2)) - (lambda12 + lambdaphi)*pow(vsigma,2) +
                    pow(pow(lambdaH,2)*pow(vh,4) -
                      2*(lambdaH*(lambda12 + lambdaphi) - 4*pow(lambdaHphi,2))*
                       pow(vh,2)*pow(vsigma,2) +
                      pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5))/2.,2));
           	amp=num/den;
        }
    else if (id == 10)
        {
          	num=-(pow(g,2)*pow(me,2)*(-0.5*s + 2*pow(me,2) - 2*pow(mius,2))*
          		     (256*pow(lambdaHphi,2)*pow(mius,4)*pow(vh,2) +
          		       16*pow(lambdaHphi,2)*pow(s,2)*pow(vh,2) -
          		       64*lambdaH*pow(lambdaHphi,2)*pow(mius,2)*pow(vh,4) +
          		       6*pow(lambdaH,2)*pow(lambdaHphi,2)*pow(vh,6) +
          		       2*cos(4*a)*pow(lambdaH,2)*pow(lambdaHphi,2)*pow(vh,6) -
          		       64*lambda12*pow(lambdaHphi,2)*pow(mius,2)*pow(vh,2)*pow(vsigma,2) -
          		       64*lambdaphi*pow(lambdaHphi,2)*pow(mius,2)*pow(vh,2)*pow(vsigma,2) +
          		       2*lambda12*lambdaphi*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,2) -
          		       2*lambda12*lambdaphi*cos(4*a)*pow(lambdaH,2)*pow(vh,4)*
          		        pow(vsigma,2) + pow(lambda12,2)*pow(lambdaH,2)*pow(vh,4)*
          		        pow(vsigma,2) - cos(4*a)*pow(lambda12,2)*pow(lambdaH,2)*pow(vh,4)*
          		        pow(vsigma,2) + 4*lambda12*lambdaH*pow(lambdaHphi,2)*pow(vh,4)*
          		        pow(vsigma,2) + 4*lambdaH*lambdaphi*pow(lambdaHphi,2)*pow(vh,4)*
          		        pow(vsigma,2) - 4*lambda12*lambdaH*cos(4*a)*pow(lambdaHphi,2)*
          		        pow(vh,4)*pow(vsigma,2) -
          		       4*lambdaH*lambdaphi*cos(4*a)*pow(lambdaHphi,2)*pow(vh,4)*
          		        pow(vsigma,2) + 16*pow(lambdaHphi,4)*pow(vh,4)*pow(vsigma,2) +
          		       16*cos(4*a)*pow(lambdaHphi,4)*pow(vh,4)*pow(vsigma,2) +
          		       pow(lambdaH,2)*pow(lambdaphi,2)*pow(vh,4)*pow(vsigma,2) -
          		       cos(4*a)*pow(lambdaH,2)*pow(lambdaphi,2)*pow(vh,4)*pow(vsigma,2) -
          		       6*lambdaH*lambdaphi*pow(lambda12,2)*pow(vh,2)*pow(vsigma,4) +
          		       6*lambdaH*lambdaphi*cos(4*a)*pow(lambda12,2)*pow(vh,2)*
          		        pow(vsigma,4) - 2*lambdaH*pow(lambda12,3)*pow(vh,2)*pow(vsigma,4) +
          		       2*lambdaH*cos(4*a)*pow(lambda12,3)*pow(vh,2)*pow(vsigma,4) +
          		       28*lambda12*lambdaphi*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,4) -
          		       12*lambda12*lambdaphi*cos(4*a)*pow(lambdaHphi,2)*pow(vh,2)*
          		        pow(vsigma,4) + 14*pow(lambda12,2)*pow(lambdaHphi,2)*pow(vh,2)*
          		        pow(vsigma,4) - 6*cos(4*a)*pow(lambda12,2)*pow(lambdaHphi,2)*
          		        pow(vh,2)*pow(vsigma,4) -
          		       6*lambda12*lambdaH*pow(lambdaphi,2)*pow(vh,2)*pow(vsigma,4) +
          		       6*lambda12*lambdaH*cos(4*a)*pow(lambdaphi,2)*pow(vh,2)*
          		        pow(vsigma,4) + 14*pow(lambdaHphi,2)*pow(lambdaphi,2)*pow(vh,2)*
          		        pow(vsigma,4) - 6*cos(4*a)*pow(lambdaHphi,2)*pow(lambdaphi,2)*
          		        pow(vh,2)*pow(vsigma,4) -
          		       2*lambdaH*pow(lambdaphi,3)*pow(vh,2)*pow(vsigma,4) +
          		       2*lambdaH*cos(4*a)*pow(lambdaphi,3)*pow(vh,2)*pow(vsigma,4) +
          		       4*lambdaphi*pow(lambda12,3)*pow(vsigma,6) -
          		       4*lambdaphi*cos(4*a)*pow(lambda12,3)*pow(vsigma,6) +
          		       pow(lambda12,4)*pow(vsigma,6) -
          		       cos(4*a)*pow(lambda12,4)*pow(vsigma,6) +
          		       6*pow(lambda12,2)*pow(lambdaphi,2)*pow(vsigma,6) -
          		       6*cos(4*a)*pow(lambda12,2)*pow(lambdaphi,2)*pow(vsigma,6) +
          		       4*lambda12*pow(lambdaphi,3)*pow(vsigma,6) -
          		       4*lambda12*cos(4*a)*pow(lambdaphi,3)*pow(vsigma,6) +
          		       pow(lambdaphi,4)*pow(vsigma,6) -
          		       cos(4*a)*pow(lambdaphi,4)*pow(vsigma,6) -
          		       8*cos(2*a)*pow(lambdaHphi,2)*pow(vh,2)*
          		        (-8*pow(mius,2) + lambdaH*pow(vh,2) +
          		          (lambda12 + lambdaphi)*pow(vsigma,2))*
          		        pow(pow(lambdaH,2)*pow(vh,4) -
          		          2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		          8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		          pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5) +
          		       32*lambda12*lambdaHphi*vh*vsigma*pow(2,0.5)*pow(mius,2)*
          		        pow(pow(lambdaH,2)*pow(vh,4) -
          		          2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		          8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		          pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) +
          		       32*lambdaHphi*lambdaphi*vh*vsigma*pow(2,0.5)*pow(mius,2)*
          		        pow(pow(lambdaH,2)*pow(vh,4) -
          		          2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		          8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		          pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) -
          		       4*lambda12*lambdaH*lambdaHphi*vsigma*pow(2,0.5)*pow(vh,3)*
          		        pow(pow(lambdaH,2)*pow(vh,4) -
          		          2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		          8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		          pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) -
          		       4*lambdaH*lambdaHphi*lambdaphi*vsigma*pow(2,0.5)*pow(vh,3)*
          		        pow(pow(lambdaH,2)*pow(vh,4) -
          		          2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		          8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		          pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) -
          		       8*lambda12*lambdaHphi*lambdaphi*vh*pow(2,0.5)*pow(vsigma,3)*
          		        pow(pow(lambdaH,2)*pow(vh,4) -
          		          2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		          8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		          pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) -
          		       4*lambdaHphi*vh*pow(2,0.5)*pow(lambda12,2)*pow(vsigma,3)*
          		        pow(pow(lambdaH,2)*pow(vh,4) -
          		          2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		          8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		          pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) -
          		       4*lambdaHphi*vh*pow(2,0.5)*pow(lambdaphi,2)*pow(vsigma,3)*
          		        pow(pow(lambdaH,2)*pow(vh,4) -
          		          2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		          8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		          pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) +
          		       8*lambdaHphi*s*vh*(-2*lambdaHphi*vh*
          		           (-8*pow(mius,2) + lambdaH*pow(vh,2) +
          		             (lambda12 + lambdaphi)*pow(vsigma,2)) +
          		          2*lambdaHphi*vh*cos(2*a)*
          		           pow(pow(lambdaH,2)*pow(vh,4) -
          		             2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		             8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		             pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5) +
          		          (lambda12 + lambdaphi)*vsigma*pow(2,0.5)*
          		           pow(pow(lambdaH,2)*pow(vh,4) -
          		             2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		             8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		             pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a)) +
          		       2*lambda12*lambdaHphi*vsigma*pow(2,0.5)*pow(lambdaH,2)*pow(vh,5)*
          		        sin(4*a) + 2*lambdaHphi*lambdaphi*vsigma*pow(2,0.5)*pow(lambdaH,2)*
          		        pow(vh,5)*sin(4*a) - 8*lambda12*lambdaH*lambdaHphi*lambdaphi*
          		        pow(2,0.5)*pow(vh,3)*pow(vsigma,3)*sin(4*a) -
          		       4*lambdaH*lambdaHphi*pow(2,0.5)*pow(lambda12,2)*pow(vh,3)*
          		        pow(vsigma,3)*sin(4*a) +
          		       16*lambda12*pow(2,0.5)*pow(lambdaHphi,3)*pow(vh,3)*pow(vsigma,3)*
          		        sin(4*a) + 16*lambdaphi*pow(2,0.5)*pow(lambdaHphi,3)*pow(vh,3)*
          		        pow(vsigma,3)*sin(4*a) -
          		       4*lambdaH*lambdaHphi*pow(2,0.5)*pow(lambdaphi,2)*pow(vh,3)*
          		        pow(vsigma,3)*sin(4*a) +
          		       6*lambdaHphi*lambdaphi*vh*pow(2,0.5)*pow(lambda12,2)*pow(vsigma,5)*
          		        sin(4*a) + 2*lambdaHphi*vh*pow(2,0.5)*pow(lambda12,3)*pow(vsigma,5)*
          		        sin(4*a) + 6*lambda12*lambdaHphi*vh*pow(2,0.5)*pow(lambdaphi,2)*
          		        pow(vsigma,5)*sin(4*a) +
          		       2*lambdaHphi*vh*pow(2,0.5)*pow(lambdaphi,3)*pow(vsigma,5)*sin(4*a)));
           	den=16*pow(mW,2)*((pow(wh,2)*(lambdaH*pow(vh,2) +
                    (lambda12 + lambdaphi)*pow(vsigma,2) +
                    pow(pow(lambdaH,2)*pow(vh,4) -
                      2*(lambdaH*(lambda12 + lambdaphi) - 4*pow(lambdaHphi,2))*
                       pow(vh,2)*pow(vsigma,2) +
                      pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)))/2. +
               pow(s + (-(lambdaH*pow(vh,2)) - (lambda12 + lambdaphi)*pow(vsigma,2) -
                    pow(pow(lambdaH,2)*pow(vh,4) -
                      2*(lambdaH*(lambda12 + lambdaphi) - 4*pow(lambdaHphi,2))*
                       pow(vh,2)*pow(vsigma,2) +
                      pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5))/2.,2))*
             ((pow(wh2,2)*(lambdaH*pow(vh,2) + (lambda12 + lambdaphi)*pow(vsigma,2) -
                    pow(pow(lambdaH,2)*pow(vh,4) -
                      2*(lambdaH*(lambda12 + lambdaphi) - 4*pow(lambdaHphi,2))*
                       pow(vh,2)*pow(vsigma,2) +
                      pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)))/2. +
               pow(s + (-(lambdaH*pow(vh,2)) - (lambda12 + lambdaphi)*pow(vsigma,2) +
                    pow(pow(lambdaH,2)*pow(vh,4) -
                      2*(lambdaH*(lambda12 + lambdaphi) - 4*pow(lambdaHphi,2))*
                       pow(vh,2)*pow(vsigma,2) +
                      pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5))/2.,2));
           	amp=num/den;
        }
    else if (id == 11)
        {
          	num=pow(g,2)*(4*pow(mius,4) + (s*(4*pow(mius,2) - 2*pow(mW,2)))/2. -
          		     4*pow(mius,2)*pow(mW,2) + 3*pow(mW,4) + pow(s,2)/4.)*
          		   (256*pow(lambdaHphi,2)*pow(mius,4)*pow(vh,2) +
          		     16*pow(lambdaHphi,2)*pow(s,2)*pow(vh,2) -
          		     64*lambdaH*pow(lambdaHphi,2)*pow(mius,2)*pow(vh,4) +
          		     6*pow(lambdaH,2)*pow(lambdaHphi,2)*pow(vh,6) +
          		     2*cos(4*a)*pow(lambdaH,2)*pow(lambdaHphi,2)*pow(vh,6) -
          		     64*lambda12*pow(lambdaHphi,2)*pow(mius,2)*pow(vh,2)*pow(vsigma,2) -
          		     64*lambdaphi*pow(lambdaHphi,2)*pow(mius,2)*pow(vh,2)*pow(vsigma,2) +
          		     2*lambda12*lambdaphi*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,2) -
          		     2*lambda12*lambdaphi*cos(4*a)*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,2) +
          		     pow(lambda12,2)*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,2) -
          		     cos(4*a)*pow(lambda12,2)*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,2) +
          		     4*lambda12*lambdaH*pow(lambdaHphi,2)*pow(vh,4)*pow(vsigma,2) +
          		     4*lambdaH*lambdaphi*pow(lambdaHphi,2)*pow(vh,4)*pow(vsigma,2) -
          		     4*lambda12*lambdaH*cos(4*a)*pow(lambdaHphi,2)*pow(vh,4)*pow(vsigma,2) -
          		     4*lambdaH*lambdaphi*cos(4*a)*pow(lambdaHphi,2)*pow(vh,4)*
          		      pow(vsigma,2) + 16*pow(lambdaHphi,4)*pow(vh,4)*pow(vsigma,2) +
          		     16*cos(4*a)*pow(lambdaHphi,4)*pow(vh,4)*pow(vsigma,2) +
          		     pow(lambdaH,2)*pow(lambdaphi,2)*pow(vh,4)*pow(vsigma,2) -
          		     cos(4*a)*pow(lambdaH,2)*pow(lambdaphi,2)*pow(vh,4)*pow(vsigma,2) -
          		     6*lambdaH*lambdaphi*pow(lambda12,2)*pow(vh,2)*pow(vsigma,4) +
          		     6*lambdaH*lambdaphi*cos(4*a)*pow(lambda12,2)*pow(vh,2)*pow(vsigma,4) -
          		     2*lambdaH*pow(lambda12,3)*pow(vh,2)*pow(vsigma,4) +
          		     2*lambdaH*cos(4*a)*pow(lambda12,3)*pow(vh,2)*pow(vsigma,4) +
          		     28*lambda12*lambdaphi*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,4) -
          		     12*lambda12*lambdaphi*cos(4*a)*pow(lambdaHphi,2)*pow(vh,2)*
          		      pow(vsigma,4) + 14*pow(lambda12,2)*pow(lambdaHphi,2)*pow(vh,2)*
          		      pow(vsigma,4) - 6*cos(4*a)*pow(lambda12,2)*pow(lambdaHphi,2)*pow(vh,2)*
          		      pow(vsigma,4) - 6*lambda12*lambdaH*pow(lambdaphi,2)*pow(vh,2)*
          		      pow(vsigma,4) + 6*lambda12*lambdaH*cos(4*a)*pow(lambdaphi,2)*pow(vh,2)*
          		      pow(vsigma,4) + 14*pow(lambdaHphi,2)*pow(lambdaphi,2)*pow(vh,2)*
          		      pow(vsigma,4) - 6*cos(4*a)*pow(lambdaHphi,2)*pow(lambdaphi,2)*
          		      pow(vh,2)*pow(vsigma,4) -
          		     2*lambdaH*pow(lambdaphi,3)*pow(vh,2)*pow(vsigma,4) +
          		     2*lambdaH*cos(4*a)*pow(lambdaphi,3)*pow(vh,2)*pow(vsigma,4) +
          		     4*lambdaphi*pow(lambda12,3)*pow(vsigma,6) -
          		     4*lambdaphi*cos(4*a)*pow(lambda12,3)*pow(vsigma,6) +
          		     pow(lambda12,4)*pow(vsigma,6) -
          		     cos(4*a)*pow(lambda12,4)*pow(vsigma,6) +
          		     6*pow(lambda12,2)*pow(lambdaphi,2)*pow(vsigma,6) -
          		     6*cos(4*a)*pow(lambda12,2)*pow(lambdaphi,2)*pow(vsigma,6) +
          		     4*lambda12*pow(lambdaphi,3)*pow(vsigma,6) -
          		     4*lambda12*cos(4*a)*pow(lambdaphi,3)*pow(vsigma,6) +
          		     pow(lambdaphi,4)*pow(vsigma,6) -
          		     cos(4*a)*pow(lambdaphi,4)*pow(vsigma,6) -
          		     8*cos(2*a)*pow(lambdaHphi,2)*pow(vh,2)*
          		      (-8*pow(mius,2) + lambdaH*pow(vh,2) +
          		        (lambda12 + lambdaphi)*pow(vsigma,2))*
          		      pow(pow(lambdaH,2)*pow(vh,4) -
          		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5) +
          		     32*lambda12*lambdaHphi*vh*vsigma*pow(2,0.5)*pow(mius,2)*
          		      pow(pow(lambdaH,2)*pow(vh,4) -
          		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) +
          		     32*lambdaHphi*lambdaphi*vh*vsigma*pow(2,0.5)*pow(mius,2)*
          		      pow(pow(lambdaH,2)*pow(vh,4) -
          		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) -
          		     4*lambda12*lambdaH*lambdaHphi*vsigma*pow(2,0.5)*pow(vh,3)*
          		      pow(pow(lambdaH,2)*pow(vh,4) -
          		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) -
          		     4*lambdaH*lambdaHphi*lambdaphi*vsigma*pow(2,0.5)*pow(vh,3)*
          		      pow(pow(lambdaH,2)*pow(vh,4) -
          		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) -
          		     8*lambda12*lambdaHphi*lambdaphi*vh*pow(2,0.5)*pow(vsigma,3)*
          		      pow(pow(lambdaH,2)*pow(vh,4) -
          		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) -
          		     4*lambdaHphi*vh*pow(2,0.5)*pow(lambda12,2)*pow(vsigma,3)*
          		      pow(pow(lambdaH,2)*pow(vh,4) -
          		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) -
          		     4*lambdaHphi*vh*pow(2,0.5)*pow(lambdaphi,2)*pow(vsigma,3)*
          		      pow(pow(lambdaH,2)*pow(vh,4) -
          		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) +
          		     8*lambdaHphi*s*vh*(-2*lambdaHphi*vh*
          		         (-8*pow(mius,2) + lambdaH*pow(vh,2) +
          		           (lambda12 + lambdaphi)*pow(vsigma,2)) +
          		        2*lambdaHphi*vh*cos(2*a)*
          		         pow(pow(lambdaH,2)*pow(vh,4) -
          		           2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		           8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		           pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5) +
          		        (lambda12 + lambdaphi)*vsigma*pow(2,0.5)*
          		         pow(pow(lambdaH,2)*pow(vh,4) -
          		           2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		           8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		           pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a)) +
          		     2*lambda12*lambdaHphi*vsigma*pow(2,0.5)*pow(lambdaH,2)*pow(vh,5)*
          		      sin(4*a) + 2*lambdaHphi*lambdaphi*vsigma*pow(2,0.5)*pow(lambdaH,2)*
          		      pow(vh,5)*sin(4*a) - 8*lambda12*lambdaH*lambdaHphi*lambdaphi*
          		      pow(2,0.5)*pow(vh,3)*pow(vsigma,3)*sin(4*a) -
          		     4*lambdaH*lambdaHphi*pow(2,0.5)*pow(lambda12,2)*pow(vh,3)*pow(vsigma,3)*
          		      sin(4*a) + 16*lambda12*pow(2,0.5)*pow(lambdaHphi,3)*pow(vh,3)*
          		      pow(vsigma,3)*sin(4*a) +
          		     16*lambdaphi*pow(2,0.5)*pow(lambdaHphi,3)*pow(vh,3)*pow(vsigma,3)*
          		      sin(4*a) - 4*lambdaH*lambdaHphi*pow(2,0.5)*pow(lambdaphi,2)*pow(vh,3)*
          		      pow(vsigma,3)*sin(4*a) +
          		     6*lambdaHphi*lambdaphi*vh*pow(2,0.5)*pow(lambda12,2)*pow(vsigma,5)*
          		      sin(4*a) + 2*lambdaHphi*vh*pow(2,0.5)*pow(lambda12,3)*pow(vsigma,5)*
          		      sin(4*a) + 6*lambda12*lambdaHphi*vh*pow(2,0.5)*pow(lambdaphi,2)*
          		      pow(vsigma,5)*sin(4*a) +
          		     2*lambdaHphi*vh*pow(2,0.5)*pow(lambdaphi,3)*pow(vsigma,5)*sin(4*a));
           	den=16*pow(mW,2)*((pow(wh,2)*(lambdaH*pow(vh,2) +
                    (lambda12 + lambdaphi)*pow(vsigma,2) +
                    pow(pow(lambdaH,2)*pow(vh,4) -
                      2*(lambdaH*(lambda12 + lambdaphi) - 4*pow(lambdaHphi,2))*
                       pow(vh,2)*pow(vsigma,2) +
                      pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)))/2. +
               pow(s + (-(lambdaH*pow(vh,2)) - (lambda12 + lambdaphi)*pow(vsigma,2) -
                    pow(pow(lambdaH,2)*pow(vh,4) -
                      2*(lambdaH*(lambda12 + lambdaphi) - 4*pow(lambdaHphi,2))*
                       pow(vh,2)*pow(vsigma,2) +
                      pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5))/2.,2))*
             ((pow(wh2,2)*(lambdaH*pow(vh,2) + (lambda12 + lambdaphi)*pow(vsigma,2) -
                    pow(pow(lambdaH,2)*pow(vh,4) -
                      2*(lambdaH*(lambda12 + lambdaphi) - 4*pow(lambdaHphi,2))*
                       pow(vh,2)*pow(vsigma,2) +
                      pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)))/2. +
               pow(s + (-(lambdaH*pow(vh,2)) - (lambda12 + lambdaphi)*pow(vsigma,2) +
                    pow(pow(lambdaH,2)*pow(vh,4) -
                      2*(lambdaH*(lambda12 + lambdaphi) - 4*pow(lambdaHphi,2))*
                       pow(vh,2)*pow(vsigma,2) +
                      pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5))/2.,2));
           	amp=num/den;
        }
    else if (id == 12)
        {
          	num=pow(g,2)*(3*pow(mW,4) + 4*pow(mius,2)*pow(mW,2)*(-1 + pow(SW,2)) +
          		     s*(pow(mW,2) + 2*pow(mius,2)*(-1 + pow(SW,2)))*(-1 + pow(SW,2)) +
          		     4*pow(mius,4)*pow(-1 + pow(SW,2),2) +
          		     (pow(s,2)*pow(-1 + pow(SW,2),2))/4.)*
          		   (256*pow(lambdaHphi,2)*pow(mius,4)*pow(vh,2) +
          		     16*pow(lambdaHphi,2)*pow(s,2)*pow(vh,2) -
          		     64*lambdaH*pow(lambdaHphi,2)*pow(mius,2)*pow(vh,4) +
          		     6*pow(lambdaH,2)*pow(lambdaHphi,2)*pow(vh,6) +
          		     2*cos(4*a)*pow(lambdaH,2)*pow(lambdaHphi,2)*pow(vh,6) -
          		     64*lambda12*pow(lambdaHphi,2)*pow(mius,2)*pow(vh,2)*pow(vsigma,2) -
          		     64*lambdaphi*pow(lambdaHphi,2)*pow(mius,2)*pow(vh,2)*pow(vsigma,2) +
          		     2*lambda12*lambdaphi*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,2) -
          		     2*lambda12*lambdaphi*cos(4*a)*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,2) +
          		     pow(lambda12,2)*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,2) -
          		     cos(4*a)*pow(lambda12,2)*pow(lambdaH,2)*pow(vh,4)*pow(vsigma,2) +
          		     4*lambda12*lambdaH*pow(lambdaHphi,2)*pow(vh,4)*pow(vsigma,2) +
          		     4*lambdaH*lambdaphi*pow(lambdaHphi,2)*pow(vh,4)*pow(vsigma,2) -
          		     4*lambda12*lambdaH*cos(4*a)*pow(lambdaHphi,2)*pow(vh,4)*pow(vsigma,2) -
          		     4*lambdaH*lambdaphi*cos(4*a)*pow(lambdaHphi,2)*pow(vh,4)*
          		      pow(vsigma,2) + 16*pow(lambdaHphi,4)*pow(vh,4)*pow(vsigma,2) +
          		     16*cos(4*a)*pow(lambdaHphi,4)*pow(vh,4)*pow(vsigma,2) +
          		     pow(lambdaH,2)*pow(lambdaphi,2)*pow(vh,4)*pow(vsigma,2) -
          		     cos(4*a)*pow(lambdaH,2)*pow(lambdaphi,2)*pow(vh,4)*pow(vsigma,2) -
          		     6*lambdaH*lambdaphi*pow(lambda12,2)*pow(vh,2)*pow(vsigma,4) +
          		     6*lambdaH*lambdaphi*cos(4*a)*pow(lambda12,2)*pow(vh,2)*pow(vsigma,4) -
          		     2*lambdaH*pow(lambda12,3)*pow(vh,2)*pow(vsigma,4) +
          		     2*lambdaH*cos(4*a)*pow(lambda12,3)*pow(vh,2)*pow(vsigma,4) +
          		     28*lambda12*lambdaphi*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,4) -
          		     12*lambda12*lambdaphi*cos(4*a)*pow(lambdaHphi,2)*pow(vh,2)*
          		      pow(vsigma,4) + 14*pow(lambda12,2)*pow(lambdaHphi,2)*pow(vh,2)*
          		      pow(vsigma,4) - 6*cos(4*a)*pow(lambda12,2)*pow(lambdaHphi,2)*pow(vh,2)*
          		      pow(vsigma,4) - 6*lambda12*lambdaH*pow(lambdaphi,2)*pow(vh,2)*
          		      pow(vsigma,4) + 6*lambda12*lambdaH*cos(4*a)*pow(lambdaphi,2)*pow(vh,2)*
          		      pow(vsigma,4) + 14*pow(lambdaHphi,2)*pow(lambdaphi,2)*pow(vh,2)*
          		      pow(vsigma,4) - 6*cos(4*a)*pow(lambdaHphi,2)*pow(lambdaphi,2)*
          		      pow(vh,2)*pow(vsigma,4) -
          		     2*lambdaH*pow(lambdaphi,3)*pow(vh,2)*pow(vsigma,4) +
          		     2*lambdaH*cos(4*a)*pow(lambdaphi,3)*pow(vh,2)*pow(vsigma,4) +
          		     4*lambdaphi*pow(lambda12,3)*pow(vsigma,6) -
          		     4*lambdaphi*cos(4*a)*pow(lambda12,3)*pow(vsigma,6) +
          		     pow(lambda12,4)*pow(vsigma,6) -
          		     cos(4*a)*pow(lambda12,4)*pow(vsigma,6) +
          		     6*pow(lambda12,2)*pow(lambdaphi,2)*pow(vsigma,6) -
          		     6*cos(4*a)*pow(lambda12,2)*pow(lambdaphi,2)*pow(vsigma,6) +
          		     4*lambda12*pow(lambdaphi,3)*pow(vsigma,6) -
          		     4*lambda12*cos(4*a)*pow(lambdaphi,3)*pow(vsigma,6) +
          		     pow(lambdaphi,4)*pow(vsigma,6) -
          		     cos(4*a)*pow(lambdaphi,4)*pow(vsigma,6) -
          		     8*cos(2*a)*pow(lambdaHphi,2)*pow(vh,2)*
          		      (-8*pow(mius,2) + lambdaH*pow(vh,2) +
          		        (lambda12 + lambdaphi)*pow(vsigma,2))*
          		      pow(pow(lambdaH,2)*pow(vh,4) -
          		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5) +
          		     32*lambda12*lambdaHphi*vh*vsigma*pow(2,0.5)*pow(mius,2)*
          		      pow(pow(lambdaH,2)*pow(vh,4) -
          		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) +
          		     32*lambdaHphi*lambdaphi*vh*vsigma*pow(2,0.5)*pow(mius,2)*
          		      pow(pow(lambdaH,2)*pow(vh,4) -
          		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) -
          		     4*lambda12*lambdaH*lambdaHphi*vsigma*pow(2,0.5)*pow(vh,3)*
          		      pow(pow(lambdaH,2)*pow(vh,4) -
          		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) -
          		     4*lambdaH*lambdaHphi*lambdaphi*vsigma*pow(2,0.5)*pow(vh,3)*
          		      pow(pow(lambdaH,2)*pow(vh,4) -
          		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) -
          		     8*lambda12*lambdaHphi*lambdaphi*vh*pow(2,0.5)*pow(vsigma,3)*
          		      pow(pow(lambdaH,2)*pow(vh,4) -
          		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) -
          		     4*lambdaHphi*vh*pow(2,0.5)*pow(lambda12,2)*pow(vsigma,3)*
          		      pow(pow(lambdaH,2)*pow(vh,4) -
          		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) -
          		     4*lambdaHphi*vh*pow(2,0.5)*pow(lambdaphi,2)*pow(vsigma,3)*
          		      pow(pow(lambdaH,2)*pow(vh,4) -
          		        2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		        8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		        pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a) +
          		     8*lambdaHphi*s*vh*(-2*lambdaHphi*vh*
          		         (-8*pow(mius,2) + lambdaH*pow(vh,2) +
          		           (lambda12 + lambdaphi)*pow(vsigma,2)) +
          		        2*lambdaHphi*vh*cos(2*a)*
          		         pow(pow(lambdaH,2)*pow(vh,4) -
          		           2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		           8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		           pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5) +
          		        (lambda12 + lambdaphi)*vsigma*pow(2,0.5)*
          		         pow(pow(lambdaH,2)*pow(vh,4) -
          		           2*lambdaH*(lambda12 + lambdaphi)*pow(vh,2)*pow(vsigma,2) +
          		           8*pow(lambdaHphi,2)*pow(vh,2)*pow(vsigma,2) +
          		           pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)*sin(2*a)) +
          		     2*lambda12*lambdaHphi*vsigma*pow(2,0.5)*pow(lambdaH,2)*pow(vh,5)*
          		      sin(4*a) + 2*lambdaHphi*lambdaphi*vsigma*pow(2,0.5)*pow(lambdaH,2)*
          		      pow(vh,5)*sin(4*a) - 8*lambda12*lambdaH*lambdaHphi*lambdaphi*
          		      pow(2,0.5)*pow(vh,3)*pow(vsigma,3)*sin(4*a) -
          		     4*lambdaH*lambdaHphi*pow(2,0.5)*pow(lambda12,2)*pow(vh,3)*pow(vsigma,3)*
          		      sin(4*a) + 16*lambda12*pow(2,0.5)*pow(lambdaHphi,3)*pow(vh,3)*
          		      pow(vsigma,3)*sin(4*a) +
          		     16*lambdaphi*pow(2,0.5)*pow(lambdaHphi,3)*pow(vh,3)*pow(vsigma,3)*
          		      sin(4*a) - 4*lambdaH*lambdaHphi*pow(2,0.5)*pow(lambdaphi,2)*pow(vh,3)*
          		      pow(vsigma,3)*sin(4*a) +
          		     6*lambdaHphi*lambdaphi*vh*pow(2,0.5)*pow(lambda12,2)*pow(vsigma,5)*
          		      sin(4*a) + 2*lambdaHphi*vh*pow(2,0.5)*pow(lambda12,3)*pow(vsigma,5)*
          		      sin(4*a) + 6*lambda12*lambdaHphi*vh*pow(2,0.5)*pow(lambdaphi,2)*
          		      pow(vsigma,5)*sin(4*a) +
          		     2*lambdaHphi*vh*pow(2,0.5)*pow(lambdaphi,3)*pow(vsigma,5)*sin(4*a));
           	den=32*pow(mW,2)*pow(-1 + pow(SW,2),2)*
           		   ((pow(wh,2)*(lambdaH*pow(vh,2) + (lambda12 + lambdaphi)*pow(vsigma,2) +
           		          pow(pow(lambdaH,2)*pow(vh,4) -
           		            2*(lambdaH*(lambda12 + lambdaphi) - 4*pow(lambdaHphi,2))*
           		             pow(vh,2)*pow(vsigma,2) +
           		            pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)))/2. +
           		     pow(s + (-(lambdaH*pow(vh,2)) - (lambda12 + lambdaphi)*pow(vsigma,2) -
           		          pow(pow(lambdaH,2)*pow(vh,4) -
           		            2*(lambdaH*(lambda12 + lambdaphi) - 4*pow(lambdaHphi,2))*
           		             pow(vh,2)*pow(vsigma,2) +
           		            pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5))/2.,2))*
           		   ((pow(wh2,2)*(lambdaH*pow(vh,2) + (lambda12 + lambdaphi)*pow(vsigma,2) -
           		          pow(pow(lambdaH,2)*pow(vh,4) -
           		            2*(lambdaH*(lambda12 + lambdaphi) - 4*pow(lambdaHphi,2))*
           		             pow(vh,2)*pow(vsigma,2) +
           		            pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5)))/2. +
           		     pow(s + (-(lambdaH*pow(vh,2)) - (lambda12 + lambdaphi)*pow(vsigma,2) +
           		          pow(pow(lambdaH,2)*pow(vh,4) -
           		            2*(lambdaH*(lambda12 + lambdaphi) - 4*pow(lambdaHphi,2))*
           		             pow(vh,2)*pow(vsigma,2) +
           		            pow(lambda12 + lambdaphi,2)*pow(vsigma,4),0.5))/2.,2));
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
        m1 = params->m_theta; m2 = params->m_theta; m3 = mH; m4 = mH; S12 = 1.0l; S34=0.5l;
    }
    else if (params->iproc == 2)// id = 2: theta theta -> b B
    {
        m1 = params->m_theta; m2 = params->m_theta; m3 = mqb; m4 = mqb; S12 = 1.0l; S34=1.0l;
    }
    else if (params->iproc == 3)// id = 3: theta theta -> t T
    {
        m1 = params->m_theta; m2 = params->m_theta; m3 = mqt; m4 = mqt; S12 =1.0l; S34=1.0l;
    }
    else if (params->iproc == 4)// id = 4: theta theta -> s S
    {
        m1 = params->m_theta; m2 = params->m_theta; m3 = mqs; m4 = mqs; S12 = 1.0l; S34=1.0l;
    }
    else if (params->iproc == 5)// id = 5: theta theta -> c C
    {
        m1 = params->m_theta; m2 = params->m_theta; m3 = mqc; m4 = mqc; S12 = 1.0l; S34=1.0l;
    }
    else if (params->iproc == 6)// id = 6: theta theta -> d D
    {
        m1 = params->m_theta; m2 = params->m_theta; m3 = mqd; m4 = mqd; S12 = 1.0l; S34=1.0l;
    }
    else if (params->iproc == 7)// id = 7: theta theta -> u U
    {
        m1 = params->m_theta; m2 = params->m_theta; m3 = mqu; m4 = mqu; S12 = 1.0l; S34=1.0l;
    }
    else if (params->iproc == 8)// id = 8: theta theta -> l L
    {
        m1 = params->m_theta; m2 = params->m_theta; m3 = mtau; m4 = mtau; S12 = 1.0l; S34=1.0l;
    }
    else if (params->iproc == 9)// id = 9: theta theta -> m M
    {
        m1 = params->m_theta; m2 = params->m_theta; m3 = mmu; m4 = mmu; S12 = 1.0l; S34=1.0l;
    }
    else if (params->iproc == 10)// id = 10: theta theta -> e E
    {
        m1 = params->m_theta; m2 = params->m_theta; m3 = me; m4 = me; S12 = 1.0l; S34=1.0l;
    }
    else if (params->iproc == 11)// id = 11: theta theta -> W+ W-
    {
        m1 = params->m_theta; m2 = params->m_theta; m3 = mW; m4 = mW; S12 = 1.0l; S34=1.0l;
    }
    else if (params->iproc == 12)// id = 12: theta theta -> Z Z
    {
        m1 = params->m_theta; m2 = params->m_theta; m3 = MZ; m4 = MZ; S12 = 1.0l; S34=0.5l;
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
