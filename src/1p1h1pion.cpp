#include "1p1h1pion.h"

double ReDelta=0;
int channel=1;
double lepton_mass=mmu;
double mpi=138.04*MeV;
double M=M_N;
double mpip=139.57018*MeV;  //pi+- mass
double mpi0=134.9766*MeV;  //pi0 mass
double mpi2=mpi*mpi;
double mpi4=mpi2*mpi2;

std::complex<double> DP_const=delusive*fstar/mpip*sqrt3; //we define it here because it depends on mpip
std::complex<double> CDP_const=delusive*fstar/mpip/sqrt3;

double width=20*MeV; //an artificial nucleon width, introduced not to hit the pole
double width_cross=20.*MeV;



double deltaffv[3]={0,0,0};
double deltaffa[3]={0,0,0};


double w_max=1e10*GeV;



//pf is defined with respect to qf direction. We rotate it so that it's defined with respect to Z direction
void rotate(D4V<double>& pf, D4V<double> qf) {
	D4V<double> q = qf;
	D4V<double> p = pf;
	D4V<double> zaxis(0,0,0,1);
	q(0)=0;
	//	p(0)=0;
	double length_q = sqrt(-(q*q));
	//	double length_p = sqrt(-p*p);
	//angle between q and z axis
	q(1)/=length_q;
	q(2)/=length_q;
	q(3)/=length_q;
	
	double costheta=(-((zaxis+q)*(zaxis+q))-2.)/2.;
	double sintheta = sqrt(1.-costheta*costheta);
	double qx=-q(2)/sqrt(q(1)*q(1)+q(2)*q(2));
	double qy=q(1)/sqrt(q(1)*q(1)+q(2)*q(2));
	double qx2= qx*qx;
	double qy2= qy*qy;
	double qxy= qx*qy;
	
	pf(1) = (costheta+qx2*(1-costheta))*p(1) + qxy*(1-costheta)*p(2) + qy*sintheta*p(3);
	pf(2) = qxy*(1-costheta)*p(1) + (costheta+qy2*(1-costheta))*p(2) -qx*sintheta*p(3);
	pf(3) = -qy*sintheta*p(1) + qx*sintheta*p(2) + costheta*p(3);
}


double gprime_l(double kPi2, double kk2) {
	

	
	double qc2 = 780*MeV*780*MeV;
	double fpi2_dpi=pow(lambda_pi*lambda_pi - mpi2,2)/pow(lambda_pi*lambda_pi - kPi2+qc2 ,2) / (mpi2-kPi2+qc2 );
	return 1./3.* (
				   fpi2_dpi * (3*kk2+qc2)+
				   
				   pow(lambda_rho*lambda_rho - mrho*mrho,2)/pow(lambda_rho*lambda_rho - kPi2+qc2 ,2)* 2.*crho* qc2 / (mrho*mrho-kPi2+qc2)
				  
				   
				   );
}

double gprime_t(double kPi2, double kk2) {

	
	double qc2 = 780*MeV*780*MeV;
	double frho2_drho = pow(lambda_rho*lambda_rho - mrho*mrho,2)/pow(lambda_rho*lambda_rho - kPi2+qc2 ,2)*crho/ (mrho*mrho-kPi2+qc2);
	return 1./3.* (
				   pow(lambda_pi*lambda_pi - mpi2,2)/pow(lambda_pi*lambda_pi - kPi2+qc2 ,2) * qc2 / (mpi2-kPi2+qc2 )+
				   
				   (2.*qc2+3.*kk2)*frho2_drho
				   
				   );
}


void boost (D4V<double>& v, D4V<double>& vect)  // orig boost
{
	
	D4V<double> velocity = v*(1./v(0));
	velocity(0) = 0;
	double vv = std::abs(velocity*velocity);
	if(vv==0) return;
	if(vv>=1)
	std::cerr<<"vv>1: "<<v(0)<<" "<<v(3)<<"\n";
	
	D4V<double> zmienne =  velocity*(-(vect * velocity)/vv); //0 component = 0

	D4V<double> stale = vect - zmienne;
	stale(0) = 0;
	
	vect -= stale;
	double gamma = 1 / sqrt (1 - vv);
	D4V<double> vect_temp = vect;
	vect(0) = vect_temp(0) - vect_temp(1) * velocity(1) - vect_temp(2) * velocity(2) - vect_temp(3) * velocity(3);
	vect(1) = vect_temp(1) - vect_temp(0) *  velocity(1);
	vect(2)  = vect_temp(2) - vect_temp(0) *  velocity(2);
	vect(3) = vect_temp(3) - vect_temp(0) *  velocity(3);
	vect*= gamma;
	vect+=stale;
	return;
}




double lambda (double a, double b, double c){
    return a*a+b*b+c*c-2*(a*b+b*c+a*c);
}

double lambda12_(const double &a,const double &b,const double &c)
{
    
    double result= sqrt(a*a+b*b+c*c-2*(a*b+b*c+c*a));
    return result;
}

void init_deltaFF(double Q2,double W) {
	

	deltaffv[0]=2.13/pow(1.+Q2/Mv2,2)/(1.+0.25*Q2/Mv2)/M_N; //C3
	deltaffv[1]=-1.51/pow(1.+Q2/Mv2,2)/(1.+0.25*Q2/Mv2)/M2; //C4
	deltaffv[2]=0.48/pow(1.+Q2/Mv2,2)/(1.+Q2/Mv2/0.776)/M2; //C5
 
    
    if (channel<11) {
		deltaffa[0]= 1.183/pow(1+Q2/MaD2,2);
		deltaffa[1]=-deltaffa[0]/4.; //C4A
        deltaffa[2]=deltaffa[0]*M2/(mpi2+Q2); //C6A
        deltaffa[1]/=M2;
        deltaffa[2]/=M2;
        
    }
}



/**
 *
 *
 * DELTA
 *
 *
 **/
//Delta(1232) free width.
double dwidth_(const double &s)
{
    double result=0;
    if(s>pow(M_N+mpip,2))
    {
        double lamb=lambda12_(s,mpi2,M2);
        double W=sqrt(s);
        result=fstar*fstar*pow(lamb/2./W,3)*M/(6*mpi2*W*M_PI);//nominator/denominator; OLD


        double qcm=0.5*lamb/W;
        double qcm2=qcm*qcm;
        double qcm3=qcm2*qcm;
        result=0.25*fstar*fstar*qcm3*(0.5*(s+M2-mpi2)+W*M_N)/(3*mpi2*s*M_PI);//nominator/denominator;

    }
    return result;
}



std::complex<double> denominator_shifted(double s, double rho) {
    

    double DW=dwidth_(s);
    std::complex<double> denom(0,0);
    double shifted_mass = Md+ReDelta*rho/(0.16/fermi3);
    denom=s- pow(shifted_mass,2) + delusive*shifted_mass*DW;
    denom=1.0/denom;
    return denom;
}

//delta projection
void projection_delta(D4V<double> pdelta, D4T<DM>& dvt) {
	
    DM pdslash(-pdelta(0),-pdelta(1),-pdelta(2),-pdelta(3),-Md);
    dvt=gmunu;
    dvt.Add((pdelta*(-onethird*2./Md2)),pdelta);
    dvt.Add((dirac*(-onethird)),dirac);
    dvt.Add((pdelta*(onethird/Md)),dirac);
    dvt.Add((dirac*(-onethird/Md)),pdelta);
    dvt&=pdslash;
	
}

void delta_32(D4V<double> pdelta, D4T<DM>& dvt) {

	DM pdslash(-pdelta(0),-pdelta(1),-pdelta(2),-pdelta(3),-Md);
	DM pdslash_real(pdelta(0),pdelta(1),pdelta(2),pdelta(3));
	double pdelta2 = pdelta*pdelta;
	dvt=gmunu; //ma dodatkowo strukturę spinorowa!
	dvt.Add((dirac*(-onethird)),dirac);
	dvt*=pdelta2/Md2;
	dvt.Add((pdelta*(-onethird/Md2)),dirac*pdslash_real);
	dvt.Add(pdslash_real*(dirac*(-onethird/Md2)),pdelta);
	dvt&=pdslash;

	
}



//vertex N-Delta-W
void vertex_delta(D4V<double> p, D4V<double> q, D4T<DM>& dvt_a, D4T<DM>& dvt_v) {
    DM qslash=q*dirac;
    D4V<double> pdelta=p+q;
    double qp=p*q;
    double qpd=q*pdelta;
    dvt_v=gmunu;
    dvt_v=dvt_v*(qslash*deltaffv[0] + qpd*deltaffv[1] + qp*deltaffv[2]);
    
    dvt_v.Add((q*(-deltaffv[0])),dirac);
    dvt_v.Add((q*(-deltaffv[1])),pdelta);
    dvt_v.Add((q*(-deltaffv[2])),p);
    dvt_v=dvt_v*gamma5;
    
    dvt_a=(gmunu*(qpd*deltaffa[1] + deltaffa[0]));
    dvt_a.Add((q*(-deltaffa[1])),pdelta); //1 - C4A, 2 - C6A, 0 - C5A
    dvt_a.Add((q*(deltaffa[2])),q);
   
}

double F_CT(D4V<double> q){
    
    double Q_2 = q*q*(-1.);
    double dmagp=2.792847;
    double dmagn=-1.913043;
    double FFAC=843*MeV;
    double xln=5.6;
    
    double tau=Q_2/Mn24;
    double gq=1./pow(1.+Q_2/(FFAC*FFAC),2);
    
    double F1P=gq*(1.+tau*dmagp)/(1.+tau);
    double F1N=gq*dmagn*xln*pow(tau,2)/(1.+xln*tau)/(1.+tau);
    return F1P-F1N;
    
}

double F_rho(double t) {
    return 1./(1-t/mrho2);
}


void Wnn(D4V<double> q, D4V<DM>& vector_nuc, D4V<DM>& axial_nuc, D4V<DM>& scalar_nuc) {
    double dmagp=2.792847;
    double dmagn=-1.913043;
    double FFAC=843*MeV;
    double xln=5.6;
    
    double Q_2 = -(q*q);
    
    double tau=Q_2/Mn24;
    double gq=1./pow(1.+Q_2/(FFAC*FFAC),2);
    
    double F1P=gq*(1.+tau*dmagp)/(1.+tau);
    double F1N=gq*dmagn*xln*pow(tau,2)/(1.+xln*tau)/(1.+tau);
    double F2P=gq*(-1.+dmagp)/(1.+tau);
    double F2N=gq*dmagn*(1.+xln*tau+tau)/(1.+xln*tau)/(1.+tau);
    //CC
    double F1V=F1P-F1N;
    double F2V=F2P-F2N;
    
    double F1S = (F1P+F1N)/2.;
    double F2S = (F2P+F2N)/2.;

    double GaGa=-gA/((1+Q_2/MA2)*(1+Q_2/MA2));
    double GaS = -gS/((1+Q_2/MA2)*(1+Q_2/MA2));
    

    vector_nuc = (q*dirac)*dirac-dirac*(q*dirac);
    vector_nuc*=(0.25*F2V/M_N);
    vector_nuc+=(dirac*F1V);
    
    axial_nuc=(((dirac+((q*dirac)*q)*(1.0/(Q_2+mpi2)))*gamma5)*(GaGa));

    scalar_nuc=(q*dirac)*dirac-dirac*(q*dirac);
    scalar_nuc*=(0.25*F2S/M_N);
    scalar_nuc += dirac*F1S;
    if (channel <11) scalar_nuc +=dirac*gamma5*(GaS/(4*sin2wein));
}


void nucleon_current_medium(D4V<double>&kPi, D4V<double>&p, D4V<double>&q, DM& pslash, DM& qslash,D4V<DM>& vector_nuc,D4V<DM>& axial_nuc, D4V<DM>& scalar_nuc, D4V<DM>& nucleon_vec, D4V<DM>& nucleon_a, D4V<DM>& nucleon_sc,double shift){
    
    DM kslash = kPi*dirac;
    double s = (p+q)*(p+q);

    //shift = 0;
	std::complex<double> propagat = 1./(s-pow(M-shift,2)+delusive*M*width); //(s-M2)/( pow(s-M2,2)+ pow(width*M,2) );
    nucleon_vec = NP_const*CGNP_v[channel]*kslash*gamma5*(pslash+qslash+M_N)*vector_nuc*propagat;
    nucleon_a =   NP_const*CGNP_a[channel]*kslash*gamma5*(pslash+qslash+M_N)*axial_nuc*propagat;
    nucleon_sc =  NP_const*CGNP_s[channel]*kslash*gamma5*(pslash+qslash+M_N)*scalar_nuc*propagat;

    
}

void nucleon_cross_current_medium(D4V<double>&kPi, D4V<double>&p, D4V<double>&q,DM& pslash, DM& qslash,D4V<DM>& vector_nuc,D4V<DM>& axial_nuc, D4V<DM>& scalar_nuc, D4V<DM>& nucleon_vec, D4V<DM>& nucleon_a, D4V<DM>& nucleon_sc,double shift){
	D4V<double> pprime = p+q-kPi;
	DM pprimeslash = pprime*dirac;
	DM kslash = kPi*dirac;
	double s = (p-kPi)*(p-kPi);

	std::complex<double> propagat =1./(s-pow(M-shift,2)+delusive*M*width_cross);
	nucleon_vec = vector_nuc*(pprimeslash-qslash+M_N)*kslash*gamma5*propagat;
	nucleon_vec *= CNP_const*CGNPC_v[channel];
	nucleon_a = axial_nuc*(pprimeslash-qslash+M_N)*kslash*gamma5*propagat;
	nucleon_a *= CNP_const*CGNPC_a[channel];
	nucleon_sc = scalar_nuc*(pprimeslash-qslash+M_N)*kslash*gamma5*propagat;
	nucleon_sc *= CNP_const*CGNPC_s[channel];
	
}





double im_delta_denominator(double s, double rho, bool ca3,  bool ca2, bool cq) {
	
	// ca3 = true: calculate ca3
	// cq = true: calculate cq and ca2
	
	double kf = pow(rho*3*M_PI*M_PI/2.,1./3.);
	double lamb=lambda12_(s,mpi2,M2);
	double W=sqrt(s);
	double qcm=0.5*lamb/W;
	
	double qtil = qcm/kf;
	double qtil2=qtil*qtil;
	double I1, I2;
	if (qtil>1) {
		I1 = 1-0.2/qtil2;
		I2 = 1-0.6/qtil2-4./21./pow(qtil2,3)+18./35./pow(qtil2,2);
	}
	else {
		I1 = qtil-0.2*pow(qtil,3);
		I2 = 33./35.*qtil-23./105.*pow(qtil,3);
	}
	double pauliDelta= 0;
	if (s>0 && lamb>0/*W>mpi+M && qcm>0*/) pauliDelta= 1./6./M_PI*pow(fstar/mpi,2)*pow(qcm,3)*(M_N/W)*(I1+I2)/2.;
	
	//MEDIUM CORRECTION
	//function of pion kinetic energy that would excite a delta of W invariant mass
	//(Epion+M)^2-kk^2 = s
	// Epion_kin = (s-M^2-mpi^2)/2M-mpi
	
	double C_Q, C_A2, C_A3;
	
	double kin_pion = (s-M*M-mpi*mpi)/2./M-mpi;
	double Tpi = kin_pion/mpi;
	double ImSigma, ImSigma_tot, alpha, alpha3, beta;
	if (kin_pion<85*MeV && kin_pion>0) {
		
		C_Q= Tpi*(20.2-8.58*Tpi+0.702*Tpi*Tpi)*MeV;
		C_A2=  (1.06*Tpi*Tpi-6.64*Tpi+22.66)*MeV;
		C_A3= kin_pion/85.*3.7;
		alpha = 1.+Tpi*(-0.309-0.315*Tpi+0.151*Tpi*Tpi);
		alpha3=1.+Tpi*(0.984-0.512*Tpi+0.1*Tpi*Tpi);
		beta =-0.038*Tpi*Tpi+0.204*Tpi+0.613;
		ImSigma = cq*C_Q*pow(rho/(0.17/fermi3),alpha)+ca2*C_A2*pow(rho/(0.17/fermi3),beta)+ca3*C_A3*pow(rho/(0.17/fermi3),alpha3);
		ImSigma_tot = C_Q*pow(rho/(0.17/fermi3),alpha)+C_A2*pow(rho/(0.17/fermi3),beta)+C_A3*pow(rho/(0.17/fermi3),alpha3);
	}
	else if (kin_pion>0){
		if (kin_pion>315*MeV) Tpi = 315*MeV/mpi;
		
		C_Q= (-5.19*Tpi*Tpi+15.35*Tpi+2.06)*MeV;
		C_A2= (1.06*Tpi*Tpi-6.64*Tpi+22.66)*MeV;
		C_A3= (-13.46*Tpi*Tpi+46.17*Tpi-20.34)*MeV;
		alpha =0.382*Tpi*Tpi-1.322*Tpi+1.466;
		beta =-0.038*Tpi*Tpi+0.204*Tpi+0.613;
		ImSigma = cq*C_Q*pow(rho/(0.17/fermi3),alpha)+ca2*C_A2*pow(rho/(0.17/fermi3),beta)+ca3*C_A3*pow(rho/(0.17/fermi3),2*beta);
		ImSigma_tot = C_Q*pow(rho/(0.17/fermi3),alpha)+C_A2*pow(rho/(0.17/fermi3),beta)+C_A3*pow(rho/(0.17/fermi3),2*beta);
	}
	else {
		if(s>M*M){
			C_A2= 22.66*MeV;
			beta =0.613;
			ImSigma=ca2*C_A2*pow(rho/(0.17/fermi3),beta);
			ImSigma_tot=C_A2*pow(rho/(0.17/fermi3),beta);
		}
		else    { ImSigma=0.; ImSigma_tot=0;}
		
	}
    double shifted_mass = Md + ReDelta*rho/(0.16/fermi3);
	return -(0*pauliDelta/2.+ImSigma)/(W+shifted_mass)/(pow(W-shifted_mass,2)+ pow(pauliDelta/2.+ImSigma_tot,2));

}



void delta_current_medium(double rho, D4V<double>& p, D4V<double>& q, D4V<double> &kPi, DM& pslash, DM& qslash, D4V<DM>& delta_vec, D4V<DM>& delta_ax){
    
    D4V<double> pdelta=p+q;
    
    
    D4T<DM> projD, del_corr;
    projection_delta(pdelta, projD); //OK
    
    
    
    D4T<DM> verD_a, verD_v;
    vertex_delta(p, q, verD_a, verD_v); //OK
    
    std::complex<double> denom = denominator_shifted(pdelta*pdelta,rho);
	if (denom!=denom) std::cerr<<"pdelta^2: "<<pdelta*pdelta<<" rho="<<rho*fermi3<<"\n";
	

    delta_32(pdelta, projD);
    verD_v&=projD*denom;
    verD_a&=projD*denom;

    
    delta_vec = (kPi*verD_v)*DP_const*CGDP_v[channel];
    delta_ax =  (kPi*verD_a)*DP_const*CGDP_a[channel];

    
}

void delta_current_crossed_medium(double rho, D4V<double>& p, D4V<double>& q, D4V<double> &kPi, DM& pslash, DM& qslash, D4V<DM>& delta_vec, D4V<DM>& delta_ax){
    
    D4V<double> pdelta=p-kPi;
    D4V<double> pprime=pdelta+q;
    
    D4T<DM> projD, del_corr;
    projection_delta(pdelta, projD); //OK
    
    D4T<DM> verD_a, verD_v;
    vertex_delta(pprime, q*(-1.), verD_a, verD_v); //OK
    
    verD_v.hermit();
    verD_v&=gamma0;
    verD_v*=gamma0;
    
    verD_a.hermit();
    verD_a&=gamma0;
    verD_a*=gamma0;
    
    
    std::complex<double> denom = denominator_shifted(pdelta*pdelta,rho);
    

	delta_32(pdelta, projD);
    verD_v*=projD*denom;
    verD_a*=projD*denom;

    
    delta_vec = (verD_v*kPi)*CDP_const*CGDPC_v[channel];
    delta_ax =  (verD_a*kPi)*CDP_const*CGDPC_a[channel];
    

    
}





void contact_term(D4V<double>& q, D4V<double>& kPi, D4V<DM>& vector_contact, D4V<DM>& axial_contact){
    
    double qk = (q-kPi)*(q-kPi);
    
    vector_contact=  dirac*(gA*F_CT(q)*gamma5);
    vector_contact*=CT_const*CGCT_v[channel];

    axial_contact = dirac*F_rho(qk);
    axial_contact*=(-1.)*CT_const*CGCT_a[channel];
    
}

void pion_pole(D4V<double>& q,D4V<double>& kPi, DM& qslash, D4V<DM>& pion_pole_cur){
    double qk = (q-kPi)*(q-kPi);
    pion_pole_cur= (q*qslash)*(1./(q*q-mpi2))*PP_const*CGPP_a[channel]*F_rho(qk);
    
}


void pion_in_flight(D4V<double>& q,D4V<double>& kPi, D4V<DM>& pion_in_flight_cur){

    double qk = (q-kPi)*(q-kPi);
     if (channel>10) pion_in_flight_cur = (kPi*2.-q)*gamma5*(2.*M_N*PF_const*CGPF_v[channel]*F_CT(q)/(qk-mpi*mpi));
    else pion_in_flight_cur = (kPi*2.-q)*gamma5*(2.*M_N*PF_const*CGPF_v[channel]*F_CT(q)/(qk-mpi0*mpi0));
    
}
