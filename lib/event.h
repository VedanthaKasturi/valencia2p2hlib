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
    double fermi_momentum; // isospin-averaged kF
    double kF_p;           // proton Fermi momentum at r
    double kF_n;           // neutron Fermi momentum at r
    double lindhard_Value;
    double transverse_response;
    double mixed_response;
    double parity_violating_response;
    double interference_response;
    double axial_vector_response;
    double int_vol_E_pi;
    double int_vol_theta_lepton;
    double int_vol_E_nuc;
    double int_vol_p_pi;
    double int_vol_phi_nuc;
    double int_vol_E_nuc2;
    double int_vol_E_lepton;
    double int_vol_theta_pi;
    double int_vol_phi_pi;
    double int_vol_r;
    bool isCC;


 
    
    event ():weight_pp(0),weight_pn(0),weight_np(0),weight_anti_nn(0),weight_anti_pn(0),weight_anti_np(0),kF_p(0),kF_n(0){}///< default constructor

    
    ClassDef(event,2);
};



#endif
