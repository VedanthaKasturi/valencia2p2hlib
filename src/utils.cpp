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


//                   name      Rho      r_p    r_n   rho0_p     rho0_n     a      Q_2p2h  Qval    Qs      Z  N   rmax
nucleus_str C12 = {"C12", Rho_HO, 1.692, 1.692, 0.076022, 0.076022, 1.082, 27.185, 16.827, 13.8778, 6, 6, 4.797, 31.841};
nucleus_str O16 = {"O16", Rho_HO, 1.833, 1.833, 0.0555037, 0.0555037, 1.544, 22.335, 14.906, 10.931, 8, 8, 5.295, 28.887};

// Q_2p2h (neutrinos) is the two-proton separation energy S_2p and Q_2p2h_anti (antineutrinos) the
// two-neutron one S_2n, in MeV, from AME2020 atomic mass excesses. Neutrinos produce mostly pp
// final states and antineutrinos nn; the np channel is not given its own Q. For reference the
// S_2n/S_2p values are: C12 31.841/27.185, O16 28.887/22.335, O18 12.188/29.055,
// Na23 23.488/24.060, Ar40 16.468/22.757, Ca40 28.931/14.710, Ca44 19.064/21.624,
// As75 18.224/17.913, Cd112 16.370/16.822, Pb208 14.106/15.381.
// Nuclei from O18 on were ported from the legacy table (temp/2p2h/utils.cpp:59-66): their
// densities and radii carry over unchanged, and the legacy Q and Qb columns land in Qval and Qs.
nucleus_str O18 = {"O18", Rho_HO, 1.881, 1.975, 0.0524594, 0.058584, 1.544, 29.055, 1.144, 14.413, 8, 10, 5.749, 12.188};
nucleus_str Na23 = {"Na23", Rho_F, 2.773, 2.81, 0.0889925, 0.0939948, 0.54, 24.06, 3.546, 4.887, 11, 12, 6.174, 23.488};
nucleus_str Ar40 = {"Ar40", Rho_F, 3.47, 3.64, 0.08098841, 0.0874744, 0.569, 22.757, 0.994, 7.991, 18, 22, 7.1557, 16.468};
nucleus_str Ca40 = {"Ca40", Rho_F, 3.51, 3.43, 0.0877563, 0.0931252, 0.563, 14.71, 13.809, 1.822, 20, 20, 6.943, 28.931};
nucleus_str Ca44 = {"Ca44", Rho_F, 3.573, 3.714, 0.0838054, 0.0909102, 0.563, 21.624, 3.142, 6.170, 20, 24, 7.1997, 19.064};
nucleus_str As75 = {"As75", Rho_F, 4.492, 4.64, 0.0745238, 0.0868426, 0.58, 17.913, 0.353, 1.688, 33, 42, 8.21, 18.224};
nucleus_str Cd112 = {"Cd112", Rho_F, 5.38, 5.58, 0.065956, 0.0793988, 0.58, 16.822, 2.075, 4.462, 48, 64, 9.08, 16.37};
nucleus_str Pb208 = {"Pb208", Rho_F, 6.624, 6.890, 0.0630574, 0.0865184, 0.549, 15.381, 2.368, 5.512, 82, 126, 10.188, 14.106};

nucleus_str GFG = {"GFG", Rho_GFG, 225./197.33, 0,0, 0, 0, 27.41, 16.827, 13.8778, 8, 8, pow(9*M_PI*8/4, 1./3.)/(225./197.33), 27.41};
nucleus_str GFG_12 = {"GFG_12", Rho_GFG, 225./197.33, 0,0, 0, 0, 27.41, 16.827, 13.8778, 6, 6, pow(9*M_PI*6/4, 1./3.)/(225./197.33), 27.41};
nucleus_str GFG_16 = {"GFG_16", Rho_GFG, 225./197.33, 0,0, 0, 0, 22.96, 14.906, 10.931, 8, 8, pow(9*M_PI*8/4, 1./3.)/(225./197.33), 22.96};


nucleus_str* nucleus_by_code(int code) {
    switch (code) {
        case 0: return &C12;
        case 1: return &O16;
        case 2: return &O18;
        case 3: return &Na23;
        case 4: return &Ar40;
        case 5: return &Ca40;
        case 6: return &Ca44;
        case 7: return &As75;
        case 8: return &Cd112;
        case 9: return &Pb208;
        default: return nullptr;
    }
}

nucleus_str* gfg_by_code(int code) { // only C12 and O16 have global-Fermi-gas counterparts
    switch (code) {
        case 0: return &GFG_12;
        case 1: return &GFG_16;
        default: return nullptr;
    }
}





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

FermiGas local_fermi_gas(double r, nucleus_str *N){
    FermiGas fg;
    fg.rho_p = getRho(r, N, true)/fermi3;
    fg.rho_n = getRho(r, N, false)/fermi3;
    fg.rho = fg.rho_p + fg.rho_n;
    fg.kF = pow(3*M_PI*M_PI*fg.rho/2., 1./3.);
    fg.kF_p = pow(3*M_PI*M_PI*fg.rho_p, 1./3.);
    fg.kF_n = pow(3*M_PI*M_PI*fg.rho_n, 1./3.);
    fg.kF_max = std::max(fg.kF_p, fg.kF_n);
    return fg;
}





