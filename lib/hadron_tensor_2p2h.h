//
//  hadron_tensor_2p2h.h
//  Pions
//
//  Created by Joanna Sobczyk on 20/01/16.
//  Copyright (c) 2016 Joanna Sobczyk. All rights reserved.
//

#ifndef __Pions__hadron_tensor_2p2h__
#define __Pions__hadron_tensor_2p2h__

#include <stdio.h>
#include "1p1h1pion.h"
#include "utils.h"

#include "event.h"
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TVectorD.h"




extern bool FERMI;
extern int NR_PART;
enum CHANNEL {CC, CCanti};
extern CHANNEL CH_2p2h;
extern bool Global;



void ulind(double QZR,double Q, double XKF,std::complex<double> &CUFUN, std::complex<double> &cudel);
double DUR (double q0, double dq, double pf1, double pf2);
std::complex<double> delta_lind (std::complex<double> q_zero,double q_mod, double rho,double k_fermi);


void hadron_2p2h_part1_Epi_thetapi_phipi_kpi_phip_Ep_isospin(double Ep, double phi_p, double kk, double phi_pi, double thetapi, double Epi, double qq, double w, double rho, D4V<DM>axial_nuc, D4V<DM> vector_nuc, D4V<DM> scalar_nuc, std::complex<double>* A00, std::complex<double>* A11, std::complex<double>* A12, std::complex<double>* A03, std::complex<double>* A33);


#endif /* defined(__Pions__hadron_tensor_2p2h__) */
