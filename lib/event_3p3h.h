#ifndef _event2_h_
#define _event2_h_
#include <iostream>
#include "TObject.h"
#include "units.h"
#include "constants.h"
#include "vect.h"


class event_3p3h:public TObject
{
public:

	double q0;
	double theta_lepton;
	double r;
    double rho;

	vect out1;
	vect out2;
	vect out3;
    
    event_3p3h (){}///< default constructor
    
    ClassDef(event_3p3h,1);
};




#endif
