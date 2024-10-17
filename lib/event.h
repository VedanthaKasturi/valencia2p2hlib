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


 
    
    event ():weight_pp(0),weight_pn(0),weight_np(0),weight_anti_nn(0),weight_anti_pn(0),weight_anti_np(0){}///< default constructor

    
    ClassDef(event,1);
};



#endif
