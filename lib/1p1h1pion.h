#ifndef _section_h_
#define _section_h_


#include "units.h"
#include "constants.h"
#include "utils.h"


extern double ReDelta;
extern double width;
extern double width_cross;


extern int channel;
extern double lepton_mass;
extern double mpi;
extern double M;


//formfactors:

extern double deltaffv[3];
extern double deltaffa[3];


void rotate(D4V<double>& pf, D4V<double> qf);
double lambda12_(const double &a,const double &b,const double &c);
void init_deltaFF(double Q2,double W);
double dwidth_(const double &s);
std::complex<double> denominator_shifted(double s, double rho);
void projection_delta(D4V<double> pdelta, D4T<DM>& dvt);
void vertex_delta(D4V<double> p, D4V<double> q, D4T<DM>& dvt_a, D4T<DM>& dvt_v);
double F_CT(D4V<double> q);
double F_rho(double t);
void Wnn(D4V<double> q, D4V<DM>& vector_nuc, D4V<DM>& axial_nuc, D4V<DM>& scalar_nuc);
void contact_term(D4V<double>& q, D4V<double>& kPi, D4V<DM>& vector_contact, D4V<DM>& axial_contact);
void pion_pole(D4V<double>& q,D4V<double>& kPi, DM& qslash, D4V<DM>& pion_pole_cur);
void pion_in_flight(D4V<double>& q,D4V<double>& kPi, D4V<DM>& pion_in_flight_cur);

void delta_current_medium(double rho, D4V<double>& p, D4V<double>& q, D4V<double> &kPi, DM& pslash, DM& qslash, D4V<DM>& delta_vec, D4V<DM>& delta_ax);
void delta_current_crossed_medium(double rho, D4V<double>& p, D4V<double>& q, D4V<double> &kPi, DM& pslash, DM& qslash, D4V<DM>& delta_vec, D4V<DM>& delta_ax);
void nucleon_current_medium(D4V<double>&kPi, D4V<double>&p, D4V<double>&q, DM& pslash, DM& qslash,D4V<DM>& vector_nuc,D4V<DM>& axial_nuc, D4V<DM>& scalar_nuc, D4V<DM>& nucleon_vec, D4V<DM>& nucleon_a, D4V<DM>& nucleon_sc, double shift);
void nucleon_cross_current_medium(D4V<double>&kPi, D4V<double>&p, D4V<double>&q,DM& pslash, DM& qslash,D4V<DM>& vector_nuc,D4V<DM>& axial_nuc, D4V<DM>& scalar_nuc, D4V<DM>& nucleon_vec, D4V<DM>& nucleon_a, D4V<DM>& nucleon_sc, double shift);


void boost (D4V<double>& v, D4V<double>& vect);
double im_delta_denominator(double s, double rho, bool ca3, bool ca2, bool cq);
double gprime_l(double kPi2, double k2);
double gprime_t(double kPi2, double k2);
void delta_32(D4V<double> pdelta, D4T<DM>& dvt);

#endif
