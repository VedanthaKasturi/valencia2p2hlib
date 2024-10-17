//
//  3p3h.hpp
//  Pions
//
//  Created by Joanna Sobczyk on 12/02/2019.
//  Copyright © 2019 Joanna Sobczyk. All rights reserved.
//

#ifndef _p3h_hpp
#define _p3h_hpp

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include "vect.h"
#include "event_3p3h.h"
#include "hadron_tensor_2p2h.h"

void test_events(double Ein, nucleus_str* nuc, double& max_value, double &total);
void event_3p3h_gaussian(double k_fermi, double w, double qq,  event_3p3h* e);
void saveEvents_3p3h(double Ein, std::string name, int nr, nucleus_str* nuc);
void threebody (double W, double m1, double m2, double m3, vect &p1, vect &p2, vect &p3);
void hadron_delta_p_phi_theta(double theta, double phi, double mom_nuc, double w, double qq, double rho, bool ca3, bool ca2, bool cq, std::complex<double>& A00, std::complex<double>& A11, std::complex<double>& A12, std::complex<double>& A03, std::complex<double>& A33);
double total_xsec_3p3h(double Ein, int nr, nucleus_str* nuc);
#endif /* _p3h_hpp */
