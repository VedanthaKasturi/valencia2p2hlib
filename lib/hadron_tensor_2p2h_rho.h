#ifndef __Pions__hadron_tensor_2p2h_rho__
#define __Pions__hadron_tensor_2p2h_rho__

#include <stdio.h>
#include "1p1h1pion.h"
#include "hadron_tensor_2p2h.h"
#include "hadron_tensor_2p2h_2b.h"



void nucleon_current_medium_rho(D4V<double>&kPi, D4V<double>&p, D4V<double>&q, DM& pslash, DM& qslash,D4V<DM>& vector_nuc,D4V<DM>& axial_nuc, D4V<DM>& scalar_nuc, D4T<DM>& nucleon_vec, D4T<DM>& nucleon_a, D4T<DM>& nucleon_sc, double shift);
void nucleon_cross_current_rho(D4V<double>&kPi, D4V<double>&p, D4V<double>&q,DM& pslash, DM& qslash,D4V<DM>& vector_nuc,D4V<DM>& axial_nuc, D4V<DM>& scalar_nuc, D4T<DM>& nucleon_vec, D4T<DM>& nucleon_a, D4T<DM>& nucleon_sc, double shift);
void delta_current_medium_rho(double rho, D4V<double>& p, D4V<double>& q, D4V<double> &kPi, DM& pslash, DM& qslash, D4T<DM>& delta_vec, D4T<DM>& delta_ax);
void delta_current_crossed_medium_rho(double rho, D4V<double>& p, D4V<double>& q, D4V<double> &kPi, DM& pslash, DM& qslash, D4T<DM>& delta_vec, D4T<DM>& delta_ax);


void hadron_2p2h_rho_Epi_thetapi_phipi_kpi_phip_Ep_isospin(double Ep, double phi_p, double kk, double phi_pi, double thetapi, double Epi, double qq, double w, double rho, D4V<DM>axial_nuc, D4V<DM> vector_nuc, D4V<DM> scalar_nuc, std::complex<double>* A00, std::complex<double>* A11, std::complex<double>* A12, std::complex<double>* A03, std::complex<double>* A33);
#endif
