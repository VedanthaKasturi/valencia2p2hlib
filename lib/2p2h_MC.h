#ifndef __Pions__2p2h_MC__
#define __Pions__2p2h_MC__

#include <stdio.h>
#include "1p1h1pion.h"
#include "hadron_tensor_2p2h.h"
#include "hadron_tensor_2p2h_2b.h"
#include "hadron_tensor_2p2h_rho.h"


void saveEvents_2p2h(double Ein, std::string name, int nr, nucleus_str* nuc);
void event_2p2h_rho_gaussian(double r, double Ein, double theta_lepton,  double Ep2, double Ep, double phi_p, double kk, double phi_pi, double thetapi, double Epi, double w, double qq, double rho, double c, double range, event* e, nucleus_str* nuc, double corr_in_medium);
void total_xsec_2p2h(double Ein,  int nr, nucleus_str* nuc, double& total_xsec);
void MCpoint_2p2h_gaussian(double r, double Ein, double theta_lepton, double Ep2, double Ep, double phi_p, double kk, double phi_pi, double thetapi, double Epi, double w, double qq, double rho, double cos_theta_nuc, double E_range, double& weight, nucleus_str* nuc, double corr);
void saveEvents_2p2h_variedEnergy(double Enu_low, double Enu_high, std::string name, int nr, nucleus_str* nuc);
void saveEvents_2p2h_fluxWeighted(std::string flux_file, std::string flux_hist_name, std::string name, int nr, nucleus_str* nuc);
#endif
