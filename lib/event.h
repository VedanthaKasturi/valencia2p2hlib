#ifndef _event1_h_
#define _event1_h_
#include <iostream>
#include "TObject.h"
#include "units.h"
#include "constants.h"



class event:public TObject
{
public:

	double q0;
    double qz;
	double theta_lepton;
	double r;
    double rho;
	double E_pi;
	double theta_pi;
	double phi_pi;
	double p_pi;
	double E_nuc;
	double E_nuc2;
	double phi_nuc;
	double weight_pp;
	double weight_pn;
	double weight_np;
	double weight_anti_nn;
	double weight_anti_pn;
	double weight_anti_np;
    double energy_neutrino;
    double fermi_momentum;
    double lindhard_Value;
    double transverse_response;
    double mixed_response;
    double parity_violating_response;
    double interference_response;
    double axial_vector_response;
    std::complex<double> rho_transversal_00[3];
    std::complex<double> rho_transversal_11[3];
    std::complex<double> rho_transversal_12[3];
    std::complex<double> rho_transversal_03[3];
    std::complex<double> rho_transversal_33[3];
    std::complex<double> pi_longitudinal_00[3];
    std::complex<double> pi_longitudinal_11[3];
    std::complex<double> pi_longitudinal_12[3];
    std::complex<double> pi_longitudinal_03[3];
    std::complex<double> pi_longitudinal_33[3];
    std::complex<double> both_interference_00[3];
    std::complex<double> both_interference_11[3];
    std::complex<double> both_interference_12[3];
    std::complex<double> both_interference_03[3];
    std::complex<double> both_interference_33[3];
    bool isCC;


 
    
    event ():weight_pp(0),weight_pn(0),weight_np(0),weight_anti_nn(0),weight_anti_pn(0),weight_anti_np(0){}///< default constructor

    
    ClassDef(event,1);
};



#endif
