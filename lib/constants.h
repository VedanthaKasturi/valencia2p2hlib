#ifndef _constants_h_
#define _constants_h_
#include "units.h"
#include "MatrixDir.h"
#include <cmath>
//typical PDG staticants in HEP units, put all the particle masses, couplings etc. here
//Masses:

const double Mp=938.27203*MeV; //Proton mass in MeV
const double Mn=939.56536*MeV; //Neutron mass in MeV (optionally precise, different masses)
const double M_N=940*MeV;//0.5*(Mp+Mn);//=938.9186795
const double M2=M_N*M_N;
const double M4=M2*M2;
const double Md=1232.0*MeV; // Delta mass in Mev
const double Md13=1520*MeV; //D13 mass in MeV
//const double mpip=139.57*MeV;  //pi+- mass
//const double mpi0=134.98*MeV;  //pi0 mass
//const double mpi=0.138*GeV;
//const double mpi=mpip;//(2*mpip+mpi0)*onethird;//=138.04
//const double mpi2=mpip*mpip;
//const double mpi4=mpi2*mpi2;
extern double mpip;  //pi+- mass
extern double mpi0;  //pi0 mass
//const double mpi=0.138*GeV;
//const double mpi=mpip;//(2*mpip+mpi0)*onethird;//=138.04
extern double mpi2;
extern double mpi4;

const double minus_=-1;
//const double zero=0;
const std::complex<double> delusive(0,1);



const double crho=2.;
const double mrho= 770*MeV;//775.8*MeV; //rho meson mass
const double mrho2= mrho*mrho;
const double mrho4= mrho2*mrho2;
const double Md2=Md*Md;
const double Md13_2=Md13*Md13;
const double Md4=Md2*Md2;
const double mmu=105.658357*MeV;
const double me=0.511*MeV;
const double mtau=1776.82*MeV;
const bool PROTON=true;
const bool NEUTRON=false;
const int CARBON=1;
const int OXYGEN=2;


const double MA=1.05*GeV;//nucleon axial mass
const double MA2=MA*MA;
const double Mv=0.84*GeV;
const double Mn24=4*M2;
const double Mv2=Mv*Mv;
const double MaD=948.2*MeV;//0.985*GeV;// //BEFORE: 1050*MeV
const double MaD2=MaD*MaD;



//couplings
const double alpha= 1.0/137.0;
//const double e_charge= 2*sqrt(alpha*Pi);// electron charge
const double GF=1.1664*1e-11/MeV/MeV;


//pi N Delta coupling
const double fstar= 2.14;//<- this one like in 2p2h 2.151;//2.151; //BEFORE: 2.13
const double fpi=93*MeV;//93.2
const double gA=1.26;//axial nucleon coupling
const double gS=-0.15;
const double gD=20/GeV;
const double lambda_pi=1.2*GeV;
const double fpNN2 = 0.08*4*M_PI;
const double gprime=0.63;
const double lambda_rho=2.5*GeV;

const double cabbibo=0.974;
const double cabbibo2=cabbibo*cabbibo;
const double mup= 2.792847351;
const double mun= -1.91304273;
const double sin2wein=0.231;
const double v_us=0.2257;







const double CGDP[11]={0,1,-onethird*sqrt2,onethird,onethird*sqrt2,sqrt2/3.,-onethird, sqrt3inv};
const double CGDPC[11]={0,1,sqrt2,3,-sqrt2, sqrt2,1.};
const double CGNP[11]={0,0,sqrt2inv,1,-sqrt2inv, sqrt2inv,1.};
const double CGNPC[11]={0,1,-sqrt2inv,0,sqrt2inv, sqrt2inv,1.};
const double CGOTH[11]={0,1,-sqrt2,-1,sqrt2,0.,-1.};
const double CG2ND[11]={0,onethird,2*onethird,2*onethird,onethird};
//  n --> p + pi0  channel = 3
//  n --> n + pi+  channel = 2
//  p --> p + pi+  channel = 1

//  p --> p + pi-  channel = 4 ==  2
//  n --> n + pi-  channel = 5 ==  1
//  p --> n + pi0  channel = 6 == -3
// For antineutrinos the code uses the p<->n mirror of the neutrino channels: slot 1 (final nn) takes
// channel 5, slot 2 (final pn) channel 4 and slot 3 (final np) channel 6 (see chans[] in the hadron_tensor files).

const double CGDP_v[16]={cabbibo,   cabbibo,  cabbibo*onethird,   -cabbibo*onethird*sqrt2,    cabbibo*onethird,   cabbibo,  cabbibo*onethird*sqrt2,
    2./3.*sqrt2*(1-2*sin2wein),     2./3.*sqrt2*(1-2*sin2wein),     -2./3.*(1-2*sin2wein), 2./3.*(1-2*sin2wein),
    sqrt2*onethird, -onethird, sqrt2*onethird, onethird,
    1./sqrt(3.)};

const double CGDPC_v[16]={cabbibo,  cabbibo,  3*cabbibo,          cabbibo*sqrt2,              cabbibo*3.,         cabbibo, -cabbibo*sqrt2,
    4.*sqrt2inv*(1-2*sin2wein),     4.*sqrt2inv*(1-2*sin2wein),     2.*(1-2*sin2wein), -2.*(1-2*sin2wein),
    sqrt2, 1., sqrt2, -1.,
    1./sqrt(3.)};

const double CGNP_v[16]={cabbibo,   0,  cabbibo,          cabbibo*sqrt2inv,           cabbibo,          0,  -cabbibo*sqrt2inv,
    sqrt2inv*(1-2*sin2wein),     sqrt2inv*(1-2*sin2wein),   (1-2*sin2wein), -(1-2*sin2wein),
    sqrt2inv/2., 1./2., sqrt2inv/2., -1./2.,
    0};

const double CGNP_s[16]={cabbibo, 0, 0, 0, 0, 0, 0,
    -2*sqrt2*sin2wein, 2*sqrt2*sin2wein, -4*sin2wein, -4*sin2wein,
    sqrt2inv, 1., -sqrt2inv, 1.,
    0};


const double CGNPC_v[16]={cabbibo,  cabbibo,  0,          -cabbibo*sqrt2inv,          0,          cabbibo,  cabbibo*sqrt2inv,
    sqrt2inv*(1-2*sin2wein),  sqrt2inv*(1-2*sin2wein),   -(1-2*sin2wein), (1-2*sin2wein),
    sqrt2inv/2., -1./2., sqrt2inv/2., 1./2.,
    1./sqrt(3.)};
const double CGNPC_s[16]={cabbibo, 0,  0, 0,  0,  0, 0,
    -2*sqrt2*sin2wein, 2*sqrt2*sin2wein, -4*sin2wein, -4*sin2wein,
    sqrt2inv, 1., -sqrt2inv, 1.,
    0};


const double CGCT_v[16]={cabbibo,  cabbibo,  -cabbibo,         -cabbibo*sqrt2,             -cabbibo,         cabbibo,  cabbibo*sqrt2,
    0,     0,     -2.*(1-2*sin2wein), 2.*(1-2*sin2wein),
    0, -1., 0, 1.,
    1./sqrt(3.)}; //contact term
const double CGPF_v[16]={cabbibo,  cabbibo,  -cabbibo,         -cabbibo*sqrt2,             -cabbibo,         cabbibo,  cabbibo*sqrt2,
    0,     0,     -2.*(1-2*sin2wein), 2.*(1-2*sin2wein),
    0, -1., 0, 1.,
    1./sqrt(3.) }; //pion-in-flight term


const double CGDP_a[15]={cabbibo,   cabbibo,  cabbibo*onethird,   -cabbibo*onethird*sqrt2,    cabbibo*onethird,   cabbibo,  cabbibo*onethird*sqrt2,
    2./3.*sqrt2,     2./3.*sqrt2,     -2./3.,    2./3.,
    0, 0, 0, 0};
const double CGDPC_a[15]={cabbibo,  cabbibo,  cabbibo*3,          cabbibo*sqrt2,              cabbibo*3.,         cabbibo, -cabbibo*sqrt2,
    4.*sqrt2inv,     4.*sqrt2inv,     2., -2.,
    0, 0, 0, 0};
const double CGNP_a[15]={cabbibo,   0,  cabbibo,          cabbibo*sqrt2inv,           cabbibo,          0,  -cabbibo*sqrt2inv,
    sqrt2inv,     sqrt2inv,    1.,      -1.,
    0, 0, 0, 0};
const double CGNPC_a[15]={cabbibo,  cabbibo,  0,          -cabbibo*sqrt2inv,          0,          cabbibo,  cabbibo*sqrt2inv,
    sqrt2inv,     sqrt2inv,    -1.,      1.,
    0, 0, 0, 0};

const double CGCT_a[15]={cabbibo,  cabbibo,  -cabbibo,         -cabbibo*sqrt2,             -cabbibo,         cabbibo,  cabbibo*sqrt2,
    0,     0,     -2., 2.,
    0, 0, 0, 0}; //contact term

const double CGPP_a[15]={cabbibo,  cabbibo,  -cabbibo,         -cabbibo*sqrt2,             -cabbibo,         cabbibo,  cabbibo*sqrt2,
    0,     0,     -2., 2.,
    0, 0, 0, 0}; //pion pole






extern std::complex<double> DP_const;
extern std::complex<double> CDP_const;
const std::complex<double> NP_const=-delusive*gA/fpi*sqrt2inv;
const std::complex<double> CNP_const=NP_const;
const std::complex<double> CT_const=-delusive/fpi*sqrt2inv;
const std::complex<double> PP_const=CT_const;
const std::complex<double> PF_const=-delusive*gA/fpi*sqrt2inv;
//
//Dirac algebra
//these will be constant elements of algebra:

//Dirac matrices and 1_{4x4}:
const DM gamma0(0);
const DM gamma1(1);
const DM gamma2(2);
const DM gamma3(3);
const DM identi(4);
const DM gamma5(5);
const DM nic(6);

const D4V<std::complex<double > > zeros(0,0,0,0);
//default 4-vector of Dirac matrices
const D4V<DM> dirac(gamma0,gamma1,gamma2,gamma3);
//default 4-vector of identities
const D4V<DM> identity(identi,identi,identi,identi);
//default contravariant vector of identities
const D4V<DM> contridentity(identi,identi*minus_,identi*minus_,identi*minus_);
//default tensor \gamma^\mu\gamma\nu
const D4T<DM> Dirac2(dirac,dirac);
//default metric tensor in 4x4 matrix space
const D4T<DM> gmunu(identi,nic,nic,nic,nic,(-1*identi),nic,nic,nic,nic,(-1*identi),nic,nic,nic,nic,(-1*identi));
//default metric tensor in 4x4 matrix space
const D4T<double> gmunud(1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,-1);
const D4T<std::complex<double>> gmunudim(1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,-1);
const D4T<std::complex<double>> oneone(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1);
const D4T<std::complex<double>> gmunuc(1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,-1);
// \sigma^{\mu\nu}=\frac{i}{2}\left\{\gamma^\mu,\gamma^\nu\right\}
const D4T<DM> sigma(nic, 0.5*(gamma0*gamma1-gamma1*gamma0)*delusive, 0.5*(gamma0*gamma2-gamma2*gamma0)*delusive, 0.5*(gamma0*gamma3-gamma3*gamma0)*delusive, -0.5*(gamma0*gamma1-gamma1*gamma0)*delusive, nic, 0.5*(gamma1*gamma2-gamma2*gamma1)*delusive, 0.5*(gamma1*gamma3-gamma3*gamma1)*delusive, -0.5*(gamma0*gamma2-gamma2*gamma0)*delusive, -0.5*(gamma1*gamma2-gamma2*gamma1)*delusive, nic, 0.5*(gamma2*gamma3-gamma3*gamma2)*delusive, -0.5*(gamma0*gamma3-gamma3*gamma0)*delusive, -0.5*(gamma1*gamma3-gamma3*gamma1)*delusive, -0.5*(gamma2*gamma3-gamma3*gamma2)*delusive,nic);


 
#endif
