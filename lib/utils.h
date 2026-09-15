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

// Sentinel for a nucleus whose two-nucleon removal energy has not been supplied yet.
// The legacy nucleus table used its single Q slot inconsistently (the 2p2h value for C12,
// the 1p1h value for O16 and everything else), so Q_2p2h has no source for the nuclei
// ported from it. select_nucleus() refuses to run such a nucleus.
constexpr double Q_2P2H_UNSET = -1.0;

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
    int A;          // MISNOMER: this is the number of PROTONS (Z), not the mass number
    int Z;          // MISNOMER: this is the number of NEUTRONS (N)
    double rmax;
};
double Rho_F(double r, nucleus_str *N, bool nucleon);
double Rho_HO(double r, nucleus_str *N, bool nucleon);

extern nucleus_str C12;
extern nucleus_str O16;
extern nucleus_str O18;
extern nucleus_str Na23;
extern nucleus_str Ar40;
extern nucleus_str Ca40;
extern nucleus_str Ca44;
extern nucleus_str As75;
extern nucleus_str Cd112;
extern nucleus_str Pb208;
extern nucleus_str GFG;
extern nucleus_str GFG_12;
extern nucleus_str GFG_16;

// Target lookup by the integer code used as argv[8]. Return nullptr for an unknown code
// (nucleus_by_code) or for a nucleus with no global-Fermi-gas counterpart (gfg_by_code).
nucleus_str* nucleus_by_code(int code);
nucleus_str* gfg_by_code(int code);




void rho_mod_ho(double &r, double &a);

void rho_mod_f(double &r, double &a);

double getRho(double r, nucleus_str*, bool nucleon);





#endif
