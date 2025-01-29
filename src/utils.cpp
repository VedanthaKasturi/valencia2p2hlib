//
//  utils.cpp
//  FSI
//
//  Created by Joanna Sobczyk on 11/03/16.
//  Copyright (c) 2016 Joanna Sobczyk. All rights reserved.
//

#include "utils.h"



double Rho_GFG(double r, nucleus_str *N, bool nucleon){
    double r_max = N->rmax;
    if (r>r_max) return 0;
    else return pow(N->r_p,3)/3./M_PI/M_PI;
    
}

double Rho_F(double r, nucleus_str *N, bool nucleon){ // Proton if true neutron if false
    if (r>8*(N->r_n)) return 0;
    if (nucleon) {
        double r_ch =N->r_p;
        double a= N->a;
        rho_mod_f(r_ch,a);
        return N->rho0_p/(1.+exp((r-r_ch)/a));
    }
    else {
        double r_ch =N->r_n;
        double a= N->a;
        rho_mod_f(r_ch,a);
        return N->rho0_n/(1.+exp((r-r_ch)/a));
    }
}


double Rho_HO(double r, nucleus_str *N, bool nucleon){
    if (r>5*(N->r_n)) return 0;
    if (nucleon) {
        double r_ch =N->r_p;
        double a= N->a;
        rho_mod_ho(r_ch,a);
        return N->rho0_p*(1+a*pow((r/r_ch),2))*exp(-pow((r/r_ch),2));
    }
    else {
        double r_ch =N->r_n;
        double a= N->a;
        rho_mod_ho(r_ch,a);
        return N->rho0_n*(1+a*pow((r/r_ch),2))*exp(-pow((r/r_ch),2));
    }
    
}


nucleus_str C12 = {"C12", Rho_HO, 1.692, 1.692, 0.076022, 0.076022, 1.082, 27.41, 16.827, 13.8778, 6, 6, 4.797};
nucleus_str O16 = {"O16", Rho_HO, 1.833, 1.833, 0.0555037, 0.0555037, 1.544, 22.33, 14.906, 10.931, 8, 8, 5.295};
nucleus_str GFG = {"GFG", Rho_GFG, 225./197.33, 0,0, 0, 0, 27.41, 16.827, 13.8778, 8, 8, pow(9*M_PI*8/4, 1./3.)/(225./197.33)};
nucleus_str GFG_12 = {"GFG_12", Rho_GFG, 225./197.33, 0,0, 0, 0, 27.41, 16.827, 13.8778, 6, 6, pow(9*M_PI*6/4, 1./3.)/(225./197.33)};
nucleus_str GFG_16 = {"GFG_16", Rho_GFG, 225./197.33, 0,0, 0, 0, 22.33, 14.906, 10.931, 8, 8, pow(9*M_PI*8/4, 1./3.)/(225./197.33)};





/////////////////////////////////////////////////////////
///////////DENSITY////////////////////
/////////////////////////////////////////////////////////

void rho_mod_ho(double &r, double &a){
    double rp = sqrt(r*r-2./3.*0.69);
    double x = a*r*r/(1.+3./2.*a)/rp/rp;
    a = 2*x/(2.-3.*x);
    r = rp;
}

void rho_mod_f(double &r, double &a) {
    double rp=r+5*0.69*r/(15*r*r+7*M_PI*M_PI*a*a);
    a = sqrt((pow(r,3)+M_PI*M_PI*a*a*rp-pow(rp,3))/(M_PI*M_PI*rp));
    r = rp;
}



double getRho(double r, nucleus_str *N, bool nucleon){
    return N->Rho(r, N, nucleon);

}





