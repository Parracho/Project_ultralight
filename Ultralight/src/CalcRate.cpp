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

// execute in bin: ./Rates idmtheta idmh2

    int idmtheta = atoi(argv[1]);
    int idmh2 = atoi(argv[2]);

    
    std::ofstream outtime ("./dat/time.dat");           // Record start time
    
    long double Ntemp = 500., nmin = -3., nmax = 5.; // Variando T 500 pontos entre 10^-3 e 10^5 (GeV)
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

    long double mtheta, mh2;
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
    long double dec1 = 0., dec2 = 0., dec3 = 0., dec4 =0., dec5 = 0.; 
    long double dec6 = 0., dec7 = 0., dec8 = 0., dec9 =0., dec10 = 0.;
    long double dec11 = 0., dec12 = 0.;
    long double freeze_out = 0.;

    long double mthetaSet[5] = {4.e-3,4.e-6,4.e3,4.e4,4.e5};
    long double mh2Set[5] = {150.e-3,200.e-3,34,1000.e6,1.0e6};

    mtheta = mthetaSet[idmtheta];
    mh2 = mh2Set[idmh2];
   

    char temp[100];
    sprintf(temp,"Rates-%d-%d.dat",idmtheta,idmh2);
    ofstream outrate (temp);

	myrate.SetParam( mtheta, mh2);

        outrate << "T" <<"\t"<< "mtheta" <<"\t"<< "mh2";
        outrate << "\t" << "dec1"<< "\t" << "error1" << "\t" << "prob1";
        outrate << "\t" << "dec2"<< "\t" << "error2" << "\t" << "prob2";
        outrate << "\t" << "dec3"<< "\t" << "error3" << "\t" << "prob3";
        outrate << "\t" << "dec4"<< "\t" << "error4" << "\t" << "prob4";
        outrate << "\t" << "dec5"<< "\t" << "error5" << "\t" << "prob5";
        outrate << "\t" << "dec6"<< "\t" << "error6" << "\t" << "prob6";
        outrate << "\t" << "dec7"<< "\t" << "error7" << "\t" << "prob7";
        outrate << "\t" << "dec8"<< "\t" << "error8" << "\t" << "prob8";
        outrate << "\t" << "dec9"<< "\t" << "error9" << "\t" << "prob9";
        outrate << "\t" <<  "dec10"<< "\t" << "error10" << "\t" << "prob10";
        outrate << "\t" << "dec11"<< "\t" << "error11" << "\t" << "prob11";
        outrate << "\t" << "dec12"<< "\t" << "error12" << "\t" << "prob12";
        outrate << endl;

    for(int i = 0; i < temper.size(); i++)
    {

        if(temper[i]>mH) fac= 100.; else fac= 100.;
        myrate.SetTemp(temper[i],fac);

	
//    void Rate::GetRate(int id, long double &rate_, double &error_, double &prob_)
//For 2to2 process


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



        if (outrate.is_open())
        {   
	outrate << std::scientific << temper[i] <<"\t"<< mtheta  <<"\t"<< mh2;
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
