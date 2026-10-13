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
    double Q_2p2h;      // neutrinos: two-proton separation energy S_2p
    double Qval;
    double Qs;
    int Z;          // number of protons
    int N;          // number of neutrons
    double rmax;
    double Q_2p2h_anti; // antineutrinos: two-neutron separation energy S_2n
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

// Local nucleon densities and Fermi momenta of the nucleus at radius r (internal units).
// kF is the isospin-averaged Fermi momentum (3 pi^2 rho/2)^(1/3) that the code currently uses
// for both species; kF_p and kF_n are the per-species values (3 pi^2 rho_p)^(1/3), (3 pi^2 rho_n)^(1/3).
struct FermiGas {
    double rho_p;
    double rho_n;
    double rho;     // rho_p + rho_n
    double kF;
    double kF_p;
    double kF_n;
    double kF_max;  // max(kF_p, kF_n)
};
FermiGas local_fermi_gas(double r, nucleus_str *N);





#endif
