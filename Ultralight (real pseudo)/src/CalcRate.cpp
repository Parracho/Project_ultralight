#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <list>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <random>
#include <vector>

#include <chrono>  // for high_resolution_clock

#include "Rate.hpp"

using namespace std;

int main(int argc, const char * argv[]) {

// execute in bin: ./Rates idmtheta 

     
    int idmtheta = atoi(argv[1]);
    
    std::ofstream outtime ("./dat/time.dat");           // Record start time
    
    long double Ntemp = 1000., nmin = -3., nmax = 5.;
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<long double> temper;
    long double col;
    for (int i = 0; i < Ntemp; ++i)
    {
        col = pow(10.0l,nmin+i*(nmax-nmin)/(Ntemp-1.));
        temper.push_back(col);
    }
    std::cout << std::endl;
    std::cout << temper.size() << " values of T from " << temper[0] << " to " << temper[temper.size()-1] << std::endl;

    Rate myrate;

    long double mtheta, v_h, lambda_H, v_sigma, lambda_phi, lambda_Hphi, miu_s, a;
    double fac;
    long double rate1 = 0.; double error1 = 0., prob1 = 0.;
    long double rate2 = 0.; double error2 = 0., prob2 = 0.;
    long double rate3 = 0.; double error3 = 0., prob3 = 0.;
    long double rate4 = 0.; double error4 = 0., prob4 = 0.;
    long double rate5 = 0.; double error5 = 0., prob5 = 0.;
    long double rate6 = 0.; double error6 = 0., prob6 = 0.;
    long double rate7 = 0.; double error7 = 0., prob7 = 0.;
    long double rate8 = 0.; double error8 = 0., prob8 = 0.;
    long double rate9 = 0.; double error9 = 0., prob9 = 0.;
    long double rate10 = 0.; double error10 = 0., prob10 = 0.;
    long double rate11 = 0.; double error11 = 0., prob11 = 0.;
    long double rate12 = 0.; double error12 = 0., prob12 = 0.;
    long double rate_total = 0.;
    long double dec1 = 0., dec2 = 0., dec3 = 0., dec4 =0., dec5 = 0.; 
    long double dec6 = 0., dec7 = 0., dec8 = 0., dec9 =0., dec10 = 0.;
    long double dec11 = 0., dec12 = 0.;
    long double freeze_out = 0.;


    long double mthetaSet[4] = {1.e-20,1.e0,1.e4,1.0e6};



    mtheta = mthetaSet[idmtheta];
    v_h = 246.;
    lambda_H = 0.5;
    v_sigma = 1.25E6;
    lambda_phi = 1.0E-8;
    lambda_Hphi = 2.0E-8;
    miu_s = 1.0E-20;
    a=1.0E-2;
    long double   mh2 = MH2(v_h, lambda_H, v_sigma, lambda_phi, lambda_Hphi);
    long double   mH = MH(v_h, lambda_H, v_sigma, lambda_phi, lambda_Hphi);
    //long double   a = Alpha(v_h, lambda_H, v_sigma, lambda_phi, lambda_Hphi);
    
    char temp[100];
    sprintf(temp,"Rates-%d.dat",idmtheta);
    ofstream outrate (temp);

                    myrate.SetParam(mtheta, v_h, lambda_H, v_sigma, lambda_phi, lambda_Hphi, miu_s);
        

        outrate << "T" <<"\t"<< "mtheta" <<"\t"<<  "mH2" << "\t" << "mH" << "\t" << "alpha"<< "\t" << "Rate_Total";
                    outrate << "\t" << "dec1"<< "\t" << "error1" << "\t" << "prob1";
                    outrate << "\t" << "dec2" << "\t" << "error2" << "\t" << "prob2";
                    outrate << "\t" << "dec3" << "\t" << "error3" << "\t" << "prob3";
                    outrate << "\t" << "dec4" << "\t" << "error4" << "\t" << "prob4";
                    outrate << "\t" << "dec5" << "\t" << "error5" << "\t" << "prob5";
                    outrate << "\t" << "dec6" << "\t" << "error6" << "\t" << "prob6";
                    outrate << "\t" << "dec7" << "\t" << "error7" << "\t" << "prob7";
                    outrate << "\t" << "dec8" << "\t" << "error8" << "\t" << "prob8";
                    outrate << "\t" << "dec9" << "\t" << "error9" << "\t" << "prob9";
                    outrate << "\t" << "dec10" << "\t" << "error10" << "\t" << "prob10";
                    outrate << "\t" << "dec11" << "\t" << "error11" << "\t" << "prob11";
                    outrate << "\t" << "dec12" << "\t" << "error12" << "\t" << "prob12";
        outrate << endl;

    for(int i = 0; i < temper.size(); i++)
    {

        if(temper[i]>mH) fac= 100.; else fac= 100.;
        myrate.SetTemp(temper[i],fac);

	
//    void Rate::GetRate(int id, long double &rate_, double &error_, double &prob_)
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

// S and T - Channel 
        myrate.GetProdRate2to2(1,rate1, error1, prob1);
        myrate.GetProdRate2to2(2,rate2, error2, prob2);
        myrate.GetProdRate2to2(3,rate3, error3, prob3);
        myrate.GetProdRate2to2(4,rate4, error4, prob4);
        myrate.GetProdRate2to2(5,rate5, error5, prob5);
        myrate.GetProdRate2to2(6,rate6, error6, prob6);
        myrate.GetProdRate2to2(7,rate7, error7, prob7);
        myrate.GetProdRate2to2(8,rate8, error8, prob8);
        myrate.GetProdRate2to2(9,rate9, error9, prob9);
        myrate.GetProdRate2to2(10,rate10, error10, prob10);
        myrate.GetProdRate2to2(11,rate11, error11, prob11);
        myrate.GetProdRate2to2(12,rate12, error12, prob12);

        dec1 = rate1/myrate.Get_Ntheta(temper[i])/myrate.Get_Hubble(temper[i]);
        dec2 = rate2/myrate.Get_Ntheta(temper[i])/myrate.Get_Hubble(temper[i]);
        dec3 = rate3/myrate.Get_Ntheta(temper[i])/myrate.Get_Hubble(temper[i]);
        dec4 = rate4/myrate.Get_Ntheta(temper[i])/myrate.Get_Hubble(temper[i]);
        dec5 = rate5/myrate.Get_Ntheta(temper[i])/myrate.Get_Hubble(temper[i]);
        dec6 = rate6/myrate.Get_Ntheta(temper[i])/myrate.Get_Hubble(temper[i]);
        dec7 = rate7/myrate.Get_Ntheta(temper[i])/myrate.Get_Hubble(temper[i]);
        dec8 = rate8/myrate.Get_Ntheta(temper[i])/myrate.Get_Hubble(temper[i]);
        dec9 = rate9/myrate.Get_Ntheta(temper[i])/myrate.Get_Hubble(temper[i]);
        dec10 = rate10/myrate.Get_Ntheta(temper[i])/myrate.Get_Hubble(temper[i]);
        dec11 = rate11/myrate.Get_Ntheta(temper[i])/myrate.Get_Hubble(temper[i]);
        dec12 = rate12/myrate.Get_Ntheta(temper[i])/myrate.Get_Hubble(temper[i]);
        
        rate_total = dec1 + dec2 + dec3 + dec4 + dec5 + dec6 + dec7 + dec8 + dec9 + dec10 + dec11 + dec12;


        if (outrate.is_open())
        {   
                    
	outrate << std::scientific << temper[i] <<"\t"<< mtheta <<"\t"<<  mh2 <<"\t"<< mH << "\t" << a <<"\t"<<  rate_total;
                    outrate << "\t" << dec1 << "\t" << error1 << "\t" << prob1;
                    outrate << "\t" << dec2 << "\t" << error2 << "\t" << prob2;
                    outrate << "\t" << dec3 << "\t" << error3 << "\t" << prob3;
                    outrate << "\t" << dec4 << "\t" << error4 << "\t" << prob4;
                    outrate << "\t" << dec5 << "\t" << error5 << "\t" << prob5;
                    outrate << "\t" << dec6 << "\t" << error6 << "\t" << prob6;
                    outrate << "\t" << dec7 << "\t" << error7 << "\t" << prob7;
                    outrate << "\t" << dec8 << "\t" << error8 << "\t" << prob8;
                    outrate << "\t" << dec9 << "\t" << error9 << "\t" << prob9;
                    outrate << "\t" << dec10 << "\t" << error10 << "\t" << prob10;
                    outrate << "\t" << dec11 << "\t" << error11 << "\t" << prob11;
                    outrate << "\t" << dec12 << "\t" << error12 << "\t" << prob12;
                    outrate << std::endl;
        }
    }

    outrate.close(); 


    // Record end time
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Elapsed time: " << elapsed.count()/60. << " min\n";

    if (outtime.is_open()){ 
        outtime << "Elapsed time: " << elapsed.count() << " s, " <<  elapsed.count()/60. << " min, " << elapsed.count()/3600. << " h\n" ;  
    }

    outtime.close();

    return 0;
}
