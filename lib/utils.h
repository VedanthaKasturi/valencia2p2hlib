//
//  utils.h
//  
//
//  Created by Joanna Sobczyk on 11/03/16.
//
//

#ifndef _utils_h
#define _utils_h
#include <stdlib.h>
#include <iostream>
#include <unistd.h>
#include <cmath>
#include <complex>
#include "constants.h"
#include <random>
#include "generatormt.h"

struct nucleus_str {
    std::string name;
    double (*Rho)(double , nucleus_str*, bool);
    double r_p;
    double r_n;
    double rho0_p;
    double rho0_n;
    double a;
    double Q_2p2h;
    double Qval;
    double Qs;
    int A;
    int Z;
    double rmax;
};
double Rho_F(double r, nucleus_str *N, bool nucleon);
double Rho_HO(double r, nucleus_str *N, bool nucleon);

extern nucleus_str C12;
extern nucleus_str O16;
extern nucleus_str GFG;
extern nucleus_str GFG_12;
extern nucleus_str GFG_16;




void rho_mod_ho(double &r, double &a);

void rho_mod_f(double &r, double &a);

double getRho(double r, nucleus_str*, bool nucleon);





#endif
